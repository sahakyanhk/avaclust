#!/usr/bin/env python3
"""
avaclust — orchestrate the all-vs-all → UPGMA → cluster → plot workflow.

Steps (each can be skipped):
  0. (optional) extract inter-chain interface residues          [--interface-cutoff]
  1. all-vs-all USalign                       → <out>/all_vs_all.tsv
  2. hierarchical (UPGMA) clustering          → <out>/clusters.dat (+ tree, matrix)
  3. (optional) align members within clusters → <out>/aligned_clusters/
  4. (optional) dendrogram + heatmap          → <out>/clustermap.png
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

from . import cluster, plot, interface

SCRIPTS = Path(__file__).parent / "scripts"
USALIGN_ALL = SCRIPTS / "usalign_all_vs_all.bash"
ALIGN_CLUSTERS = SCRIPTS / "align_clusters.bash"

ARCHIVE_SUFFIXES = (".tar.gz", ".tgz", ".tar", ".tar.bz2", ".tar.xz")
STRUCT_EXTS = (".pdb", ".ent")


def log(msg):
    print(msg, file=sys.stderr)


def prepare_structures(input_path, out_dir):
    """Return a directory of structures as .pdb files.

    Accepts either a directory or a tar archive (.tar.gz/.tgz/.tar/...), and
    picks up both .pdb and .ent files. A plain flat directory of .pdb files is
    used as-is. Otherwise structures are collected into <out_dir>/structures as
    <stem>.pdb — symlinks for a directory input, extracted copies for an archive
    — so the downstream bash scripts can assume one flat dir with a .pdb suffix.
    """
    input_path = input_path.rstrip('/')
    norm = os.path.join(out_dir, "structures")

    # Archive: extract structure members straight into structures/ as <stem>.pdb.
    # Using only the basename also neutralises any path-traversal in the archive.
    if input_path.endswith(ARCHIVE_SUFFIXES):
        os.makedirs(norm, exist_ok=True)
        log(f"Extracting {input_path} → {norm}")
        n = 0
        with tarfile.open(input_path) as tar:
            for m in tar.getmembers():
                base = os.path.basename(m.name)
                # skip non-structures and hidden/AppleDouble sidecars (._foo.pdb)
                if not (m.isfile() and base.endswith(STRUCT_EXTS) and not base.startswith(".")):
                    continue
                src = tar.extractfile(m)
                if src is None:
                    continue
                dest = os.path.join(norm, os.path.splitext(os.path.basename(m.name))[0] + ".pdb")
                with src, open(dest, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                n += 1
        if n == 0:
            sys.exit(f"ERROR: no {'/'.join(STRUCT_EXTS)} files found in {input_path}")
        log(f"Prepared {n} structures → {norm}")
        return norm

    if not os.path.isdir(input_path):
        sys.exit(f"ERROR: {input_path} is not a directory or a supported archive "
                 f"({', '.join(ARCHIVE_SUFFIXES)})")

    files = []
    for ext in STRUCT_EXTS:
        files += glob.glob(os.path.join(input_path, "**", f"*{ext}"), recursive=True)
    files = sorted(set(files))
    if not files:
        sys.exit(f"ERROR: no {'/'.join(STRUCT_EXTS)} files found in {input_path}")

    # Fast path: a flat directory of only .pdb files needs no normalisation.
    if all(f.endswith(".pdb") and os.path.dirname(f) == input_path for f in files):
        return input_path

    os.makedirs(norm, exist_ok=True)
    for f in files:
        link = os.path.join(norm, os.path.splitext(os.path.basename(f))[0] + ".pdb")
        if os.path.lexists(link):
            os.remove(link)
        os.symlink(os.path.abspath(f), link)
    log(f"Prepared {len(files)} structures → {norm}")
    return norm


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="avaclust", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True,
                   help="Directory of .pdb/.ent files, or a tar archive (.tar.gz/.tgz/...)")
    p.add_argument("-o", "--output", default="avaclust_out", help="Output directory")

    g = p.add_mutually_exclusive_group()
    g.add_argument("-c", "--cutoff", type=float, default=None,
                   help="TM-score cutoff; cut the tree at height 1 - cutoff "
                        "(default: auto-detect from the largest dendrogram gap)")
    g.add_argument("--k", type=int, default=None, help="Force a specific cluster count")

    p.add_argument("--linkage", default="average",
                   choices=["average", "complete", "single"],
                   help="Linkage method (default: average = UPGMA)")
    p.add_argument("--chains", default=None,
                   help="Chains to align, e.g. 'A,B' or 'A' "
                        "(default: USalign uses the first chain of each structure)")
    p.add_argument("--threads", default=None, help="Threads for USalign (default: all)")
    p.add_argument("--interface-cutoff", type=float, default=None,
                   help="Extract interface residues within this many Angstroms first")
    p.add_argument("--min-cluster-size", type=int, default=3,
                   help="Min members for a cluster to be structurally aligned (default: 3)")
    p.add_argument("--no-optimal-leaf-order", action="store_true",
                   help="Skip optimal leaf ordering (faster for large N)")
    p.add_argument("--no-align", action="store_true",
                   help="Skip per-cluster structural alignment")
    p.add_argument("--no-plot", action="store_true",
                   help="Skip the dendrogram + heatmap figure")
    p.add_argument("--figsize", nargs=2, type=float, default=(6.0, 6.0),
                   metavar=("W", "H"), help="Figure size in inches (default: 6 6)")
    p.add_argument("--labels", action="store_true",
                   help="Label the dendrogram leaves with structure names (default: off)")
    p.add_argument("--leaf-fontsize", type=float, default=4,
                   help="Font size for leaf labels when --labels is set (default: 4)")
    return p.parse_args(argv)


def _run_signature(pdb_dir, chains):
    """Signature of an all-vs-all run: its input structures and chain selection.

    Used to decide whether an existing TSV can be reused. We can't compare line
    counts, because unalignable pairs are legitimately dropped from the output
    (an incomplete-but-final matrix would otherwise always look stale).
    """
    names = sorted(f for f in os.listdir(pdb_dir) if f.endswith(".pdb"))
    return f"chains={chains or 'first'}\n" + "\n".join(names) + "\n"


def run_all_vs_all(pdb_dir, out_tsv, chains, threads):
    """Run USalign all-vs-all, reusing results from an identical previous run.

    A sidecar '<tsv>.done' holding the run signature is written only after the
    run succeeds, so an interrupted run is never mistaken for a complete one.
    """
    sig = _run_signature(pdb_dir, chains)
    done = out_tsv + ".done"
    if os.path.isfile(out_tsv) and os.path.isfile(done):
        with open(done) as fh:
            if fh.read() == sig:
                log(f"=== Reusing existing {out_tsv} ===")
                return
    log("=== Running all-vs-all USalign ===")
    cmd = ["bash", str(USALIGN_ALL), pdb_dir, out_tsv,
           chains or "", str(threads) if threads else "100%"]
    subprocess.run(cmd, check=True)
    with open(done, "w") as fh:
        fh.write(sig)


def run_align(clusters_dat, out_dir, pdb_dir, min_size, chains):
    log("=== Aligning structures within clusters ===")
    subprocess.run(
        ["bash", str(ALIGN_CLUSTERS), clusters_dat, out_dir, pdb_dir,
         str(min_size), chains or ""],
        check=True,
    )


def main(argv=None):
    args = parse_args(argv)

    out_dir = args.output.rstrip('/')
    os.makedirs(out_dir, exist_ok=True)
    pdb_dir = prepare_structures(args.input, out_dir)

    log(f"PDB dir:    {pdb_dir}")
    log(f"Output dir: {out_dir}")
    log(f"Cutoff:     {args.cutoff if args.cutoff is not None else 'auto'}")
    log(f"Chains:     {args.chains or 'first chain of each structure'}")

    # 0. Optional interface extraction.
    if args.interface_cutoff is not None:
        if not args.chains:
            sys.exit("ERROR: --interface-cutoff requires --chains (which chain[s] "
                     "define the interface)")
        log(f"=== Extracting interface residues within {args.interface_cutoff} Å ===")
        chains = [c.strip() for c in args.chains.split(',')]
        threads = int(args.threads) if args.threads else 1
        pdb_dir = interface.extract_dir(
            pdb_dir, chains, args.interface_cutoff,
            output=os.path.join(out_dir, f"interface_i{args.interface_cutoff}"),
            threads=threads,
        )

    # 1. All-vs-all alignments.
    all_vs_all = os.path.join(out_dir, "all_vs_all.tsv")
    run_all_vs_all(pdb_dir, all_vs_all, args.chains, args.threads)

    # 2. Hierarchical clustering.
    log("=== Hierarchical (UPGMA) clustering ===")
    paths = cluster.run(
        all_vs_all, os.path.join(out_dir, "clusters.dat"),
        cutoff=args.cutoff, k=args.k, linkage_method=args.linkage,
        optimal_leaf_order=not args.no_optimal_leaf_order,
    )

    # 3. Optional per-cluster structural alignment.
    if not args.no_align:
        run_align(paths["clusters"], os.path.join(out_dir, "aligned_clusters"),
                  pdb_dir, args.min_cluster_size, args.chains)

    # 4. Optional dendrogram + heatmap.
    if not args.no_plot:
        log("=== Plotting dendrogram + heatmap ===")
        plot.plot_clustermap(
            paths["newick"], paths["matrix"],
            output=os.path.join(out_dir, "clustermap.png"),
            cutoff=args.cutoff, figsize=tuple(args.figsize),
            add_labels=args.labels, leaf_fontsize=args.leaf_fontsize,
        )

    log("=== Done ===")


if __name__ == "__main__":
    main()
