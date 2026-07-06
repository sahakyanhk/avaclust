#!/usr/bin/env python3
"""
Extract interface residues from PDB files based on inter-chain distance.

One chain (e.g. "A"):  keep chain A in full + residues of other chains within
                       cutoff of A.
Two chains (e.g. "A,B"): keep residues of A near B and residues of B near A.

Importable API: extract_one(), extract_dir().
Run standalone for single-file or batch mode (see --help).
"""

import glob
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
    """Find interface residues with a KDTree neighbor search. Returns a set or None."""
    model = structure[0]
    chain_ids = {c.id for c in model.get_chains()}

    for c in chains:
        if c not in chain_ids:
            print(f"ERROR: chain '{c}' not found. Available: "
                  f"{', '.join(sorted(chain_ids))}", file=sys.stderr)
            return None

    if len(chains) == 1:
        query_chain = chains[0]
        partner_chains = chain_ids - {query_chain}
        if not partner_chains:
            print("ERROR: only one chain in structure, nothing to compare", file=sys.stderr)
            return None

        query_atoms = list(model[query_chain].get_atoms())
        partner_atoms = [a for cid in partner_chains for a in model[cid].get_atoms()]
        interface = set(model[query_chain].get_residues())  # full query chain
        ns_query = NeighborSearch(query_atoms)
        for atom in partner_atoms:
            if ns_query.search(atom.get_vector().get_array(), cutoff, level='R'):
                interface.add(atom.get_parent())
    else:
        chain_a, chain_b = chains
        atoms_a = list(model[chain_a].get_atoms())
        atoms_b = list(model[chain_b].get_atoms())
        interface = set()
        ns_b = NeighborSearch(atoms_b)
        for atom in atoms_a:
            if ns_b.search(atom.get_vector().get_array(), cutoff, level='R'):
                interface.add(atom.get_parent())
        ns_a = NeighborSearch(atoms_a)
        for atom in atoms_b:
            if ns_a.search(atom.get_vector().get_array(), cutoff, level='R'):
                interface.add(atom.get_parent())

    return interface


def extract_one(pdb_path, chains, cutoff, output):
    """Extract the interface of one PDB to `output`. Returns True on success."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('s', pdb_path)
    interface = get_interface_residues(structure, chains, cutoff)

    if not interface:
        if interface is not None:
            print(f"WARNING: no interface residues in {pdb_path}", file=sys.stderr)
        return False

    io = PDBIO()
    io.set_structure(structure)
    io.save(output, InterfaceSelect(interface))

    counts = Counter(r.get_parent().id for r in interface)
    summary = ", ".join(f"chain {c}: {n}" for c, n in sorted(counts.items()))
    print(f"  {os.path.basename(pdb_path)}: {len(interface)} residues ({summary}) "
          f"→ {output}", file=sys.stderr)
    return True


def _worker(item):
    return extract_one(*item)


def extract_dir(input_dir, chains, cutoff, output=None, threads=1):
    """Extract interfaces for every .pdb in `input_dir`. Returns the output dir."""
    pdb_files = sorted(glob.glob(os.path.join(input_dir, '*.pdb')))
    if not pdb_files:
        raise SystemExit(f"ERROR: no .pdb files found in {input_dir}")

    out_dir = output or (input_dir.rstrip('/') + '_interface')
    os.makedirs(out_dir, exist_ok=True)

    work = []
    for pdb_path in pdb_files:
        base = os.path.basename(pdb_path).rsplit('.pdb', 1)[0]
        work.append((pdb_path, chains, cutoff, os.path.join(out_dir, f"{base}.pdb")))

    t0 = time.time()
    if threads > 1:
        with multiprocessing.Pool(threads) as pool:
            n_ok = sum(pool.map(_worker, work))
    else:
        n_ok = sum(_worker(item) for item in work)

    print(f"Done: {n_ok}/{len(pdb_files)} processed in {time.time() - t0:.1f}s "
          f"→ {out_dir}", file=sys.stderr)
    return out_dir


def main(argv=None):
    import argparse

    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="Input PDB file or directory")
    p.add_argument("-c", "--chain", default="A,B",
                   help="Chain(s): 'A' or 'A,B' (default: A,B)")
    p.add_argument("-o", "--output", default=None,
                   help="Output file (single) or directory (dir mode)")
    p.add_argument("--cut", "--cutoff", dest="cutoff", type=float, default=5.0,
                   help="Distance cutoff in Angstroms (default: 5.0)")
    p.add_argument("-t", "--threads", type=int, default=1,
                   help="Parallel workers for directory mode (default: 1)")
    args = p.parse_args(argv)

    chains = [c.strip() for c in args.chain.split(',')]
    if len(chains) > 2:
        p.error("--chain accepts at most 2 chains (e.g. A or A,B)")

    if os.path.isdir(args.input):
        extract_dir(args.input, chains, args.cutoff, args.output, args.threads)
    else:
        output = args.output or (args.input.rsplit('.pdb', 1)[0] + '_interface.pdb')
        extract_one(args.input, chains, args.cutoff, output)


if __name__ == "__main__":
    main()
