#!/usr/bin/env python3
"""
Extract interface residues from PDB files based on inter-chain distance.

Usage:
  # Extract everything within 5A of chain A (from all other chains)
  python extract_interface.py -i complex.pdb -c A -cut 5

  # Extract interface between chains A and B (residues of A near B + residues of B near A)
  python extract_interface.py -i complex.pdb -c A,B -cut 5

  # Batch mode: process multiple PDBs with different chain mappings
  python extract_interface.py -b batch.txt -cut 5 -od output_dir/

  batch.txt format (tab-separated):
    complex1.pdb    A,B
    complex2.pdb    A
    complex3.pdb    B,C
"""

import argparse
import multiprocessing
import os
import sys
import time
from collections import Counter

from Bio.PDB import PDBParser, PDBIO, Select, NeighborSearch


class InterfaceSelect(Select):
    """Select only residues in the interface set."""

    def __init__(self, interface_residues):
        self._residues = interface_residues

    def accept_residue(self, residue):
        return residue in self._residues


def get_interface_residues(structure, chains, cutoff):
    """
    Find interface residues using KDTree neighbor search.

    One chain (e.g. [A]): extract chain A in full, plus residues from
    all other chains that are within cutoff of chain A.

    Two chains (e.g. [A, B]): extract residues of A near B AND
    residues of B near A. Both chains are included.
    """
    model = structure[0]
    chain_ids = {c.id for c in model.get_chains()}

    for c in chains:
        if c not in chain_ids:
            print(
                f"ERROR: chain '{c}' not found in structure. "
                f"Available chains: {', '.join(sorted(chain_ids))}",
                file=sys.stderr,
            )
            return None

    if len(chains) == 1:
        query_chain = chains[0]
        partner_chains = chain_ids - {query_chain}
        if not partner_chains:
            print("ERROR: only one chain in structure, nothing to compare", file=sys.stderr)
            return None

        query_atoms = list(model[query_chain].get_atoms())
        partner_atoms = [a for cid in partner_chains for a in model[cid].get_atoms()]

        # Include all residues of the query chain
        interface = set(model[query_chain].get_residues())

        # Add partner residues near the query chain
        ns_query = NeighborSearch(query_atoms)
        for atom in partner_atoms:
            neighbors = ns_query.search(atom.get_vector().get_array(), cutoff, level='R')
            if neighbors:
                interface.add(atom.get_parent())

    else:
        chain_a, chain_b = chains
        atoms_a = list(model[chain_a].get_atoms())
        atoms_b = list(model[chain_b].get_atoms())

        # Residues of A near B
        ns_b = NeighborSearch(atoms_b)
        interface = set()
        for atom in atoms_a:
            neighbors = ns_b.search(atom.get_vector().get_array(), cutoff, level='R')
            if neighbors:
                interface.add(atom.get_parent())

        # Residues of B near A
        ns_a = NeighborSearch(atoms_a)
        for atom in atoms_b:
            neighbors = ns_a.search(atom.get_vector().get_array(), cutoff, level='R')
            if neighbors:
                interface.add(atom.get_parent())

    return interface


def process_one(parser, io, pdb_path, chains, cutoff, output):
    """Process a single PDB file. Returns True on success."""
    structure = parser.get_structure('s', pdb_path)
    interface = get_interface_residues(structure, chains, cutoff)

    if interface is None:
        return False

    if not interface:
        print(f"WARNING: no interface residues found in {pdb_path}", file=sys.stderr)
        return False

    io.set_structure(structure)
    io.save(output, InterfaceSelect(interface))

    chain_counts = Counter(r.get_parent().id for r in interface)
    summary = ", ".join(f"chain {c}: {n}" for c, n in sorted(chain_counts.items()))
    print(f"  {os.path.basename(pdb_path)}: {len(interface)} residues ({summary}) → {output}", file=sys.stderr)
    return True


def _process_one_worker(args):
    """Worker for multiprocessing — creates its own parser/io per call."""
    pdb_path, chains, cutoff, output = args
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    return process_one(parser, io, pdb_path, chains, cutoff, output)


def main():
    p = argparse.ArgumentParser(
        description="Extract interface residues from PDB files"
    )

    # Single-file or directory mode
    p.add_argument("-i", "--input", default=None, help="Input PDB file or directory")
    p.add_argument(
        "-c", "--chain", default="A,B",
        help="Chain(s) to extract interface for. "
             "One chain (e.g. A): extract other chains' residues within cutoff of A. "
             "Two chains (e.g. A,B): extract interface between A and B."
    )
    p.add_argument(
        "-o", "--output", default=None,
        help="Output: file path (single-file mode) or directory (dir/batch mode). "
             "Defaults: <input>_interface.pdb for files, <input_dir>_interface/ for dirs."
    )

    # Batch mode
    p.add_argument(
        "-b", "--batch", default=None,
        help="Batch file: each line is <pdb_path> [<chains>]"
    )

    # Shared
    p.add_argument(
        "-cut", "--cutoff", type=float, default=5.0,
        help="Distance cutoff in Angstroms (default: 5.0)"
    )
    p.add_argument(
        "-t", "--threads", type=int, default=1,
        help="Number of parallel workers for directory/batch mode (default: 1)"
    )
    args = p.parse_args()

    if args.input is None and args.batch is None:
        p.error("provide either -i (single file) or -b (batch file)")
    if args.input is not None and args.batch is not None:
        p.error("-i and -b are mutually exclusive")

    parser = PDBParser(QUIET=True)
    io = PDBIO()

    chains = [c.strip() for c in args.chain.split(',')]
    if len(chains) > 2:
        p.error("--chain accepts at most 2 chains (e.g. A or A,B)")

    if args.input is not None and os.path.isdir(args.input):
        # Directory mode: process all .pdb files
        import glob
        pdb_files = sorted(glob.glob(os.path.join(args.input, '*.pdb')))
        if not pdb_files:
            print(f"ERROR: no .pdb files found in {args.input}", file=sys.stderr)
            sys.exit(1)

        out_dir = args.output or (args.input.rstrip('/') + '_interface')
        os.makedirs(out_dir, exist_ok=True)

        work_items = []
        for pdb_path in pdb_files:
            base = os.path.basename(pdb_path).rsplit('.pdb', 1)[0]
            output = os.path.join(out_dir, f"{base}_interface.pdb")
            work_items.append((pdb_path, chains, args.cutoff, output))

        t0 = time.time()
        if args.threads > 1:
            with multiprocessing.Pool(args.threads) as pool:
                results = pool.map(_process_one_worker, work_items)
            n_ok = sum(results)
        else:
            n_ok = 0
            for item in work_items:
                if _process_one_worker(item):
                    n_ok += 1

        elapsed = time.time() - t0
        print(f"Done: {n_ok}/{len(pdb_files)} processed in {elapsed:.1f}s → {out_dir}", file=sys.stderr)

    elif args.input is not None:
        # Single file mode
        output = args.output
        if output is None:
            base = args.input.rsplit('.pdb', 1)[0]
            output = f"{base}_interface.pdb"

        process_one(parser, io, args.input, chains, args.cutoff, output)

    else:
        # Batch mode
        t0 = time.time()
        work_items = []

        if args.output:
            os.makedirs(args.output, exist_ok=True)

        with open(args.batch) as fh:
            for lineno, line in enumerate(fh, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split()
                if len(parts) == 1:
                    pdb_path = parts[0]
                    chain_str = args.chain
                elif len(parts) == 2:
                    pdb_path, chain_str = parts
                else:
                    print(f"WARNING: skipping line {lineno}: expected <pdb> [<chains>]", file=sys.stderr)
                    continue

                batch_chains = [c.strip() for c in chain_str.split(',')]
                if len(batch_chains) > 2:
                    print(f"WARNING: skipping line {lineno}: max 2 chains", file=sys.stderr)
                    continue

                if not os.path.isfile(pdb_path):
                    print(f"WARNING: file not found: {pdb_path}", file=sys.stderr)
                    continue

                if args.output:
                    base = os.path.basename(pdb_path).rsplit('.pdb', 1)[0]
                    output = os.path.join(args.output, f"{base}_interface.pdb")
                else:
                    base = pdb_path.rsplit('.pdb', 1)[0]
                    output = f"{base}_interface.pdb"

                work_items.append((pdb_path, batch_chains, args.cutoff, output))

        if args.threads > 1 and len(work_items) > 1:
            with multiprocessing.Pool(args.threads) as pool:
                results = pool.map(_process_one_worker, work_items)
            n_ok = sum(results)
        else:
            n_ok = 0
            for item in work_items:
                if _process_one_worker(item):
                    n_ok += 1

        elapsed = time.time() - t0
        print(f"Batch done: {n_ok}/{len(work_items)} processed in {elapsed:.1f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
