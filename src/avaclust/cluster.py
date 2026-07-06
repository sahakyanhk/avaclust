#!/usr/bin/env python3
"""
Hierarchical (UPGMA-style) clustering of pairwise structural similarity.

Input: TSV from usalign_all_vs_all.bash (header lines start with '#').
Each PDB pair appears once with the best TM-score across chain combos.

Pipeline:
  1. Parse pairwise TM-scores into a symmetric distance matrix (1 - TM).
  2. Build a linkage tree (default 'average' = UPGMA).
  3. Cut it: --cutoff (TM threshold), --k (fixed count), or auto (largest gap).
  4. Pick a centroid per cluster (member with highest mean TM to its mates).
  5. Write clusters.dat plus *_summary.tsv, *_tree.nwk, *_matrix.tsv.

Importable API: load_distance_matrix(), build_clusters(), run().
Run standalone for debugging a single step (see --help).
"""

import sys
import time
from collections import defaultdict

import numpy as np
from scipy.cluster.hierarchy import (
    linkage, fcluster, to_tree, leaves_list, optimal_leaf_ordering,
)
from scipy.spatial.distance import squareform


def load_distance_matrix(input_tsv, qcol=0, tcol=1, col=2):
    """Read pairwise TM-scores into a symmetric distance matrix (1 - TM)."""
    name_to_id = {}
    names = []
    scores = {}
    t0 = time.time()
    n_lines = 0

    with open(input_tsv, buffering=1 << 22) as fh:
        for line in fh:
            if not line or line[0] == '#':
                continue
            parts = line.rstrip('\n').split('\t')
            try:
                val = float(parts[col])
            except (IndexError, ValueError):
                continue

            n_lines += 1
            query, target = parts[qcol], parts[tcol]
            for name in (query, target):
                if name not in name_to_id:
                    name_to_id[name] = len(names)
                    names.append(name)
            qid, tid = name_to_id[query], name_to_id[target]
            scores[(qid, tid)] = val
            scores[(tid, qid)] = val

    n = len(names)
    dist = np.ones((n, n), dtype=np.float64)
    np.fill_diagonal(dist, 0.0)
    for (i, j), tm in scores.items():
        dist[i, j] = 1.0 - tm

    print(
        f"Loaded {n_lines:,} pairs, {n:,} structures in {time.time() - t0:.1f}s",
        file=sys.stderr,
    )
    return dist, names


def find_auto_cut(Z, n):
    """Pick k by the largest relative gap in the sorted merge heights."""
    heights = Z[:, 2]
    if len(heights) < 2:
        return 2

    gaps = np.diff(heights)
    rel_gaps = gaps / np.maximum(heights[:-1], 1e-10)
    best_idx = int(np.argmax(rel_gaps))
    k = max(n - best_idx - 1, 2)

    print(
        f"Auto-cut: largest gap at height {heights[best_idx]:.4f} "
        f"(gap={gaps[best_idx]:.4f}, relative={rel_gaps[best_idx]:.4f}), k={k}",
        file=sys.stderr,
    )
    return k


def build_clusters(dist, names, cutoff=None, k=None,
                   linkage_method="average", optimal_leaf_order=False):
    """Cluster the distance matrix and return (clusters, Z).

    clusters: list of member-index lists, centroid first, largest cluster first.
    Z: scipy linkage matrix.
    """
    n = len(names)
    condensed = squareform(dist)
    Z = linkage(condensed, method=linkage_method)

    if optimal_leaf_order:
        print("Applying optimal leaf ordering...", file=sys.stderr)
        Z = optimal_leaf_ordering(Z, condensed)

    if k is not None:
        print(f"Using user-specified k={k}", file=sys.stderr)
        labels = fcluster(Z, t=k, criterion='maxclust')
    elif cutoff is not None:
        cut_height = 1.0 - cutoff
        print(f"Cutting at height={cut_height:.4f} (TM cutoff={cutoff})", file=sys.stderr)
        labels = fcluster(Z, t=cut_height, criterion='distance')
    else:
        labels = fcluster(Z, t=find_auto_cut(Z, n), criterion='maxclust')

    members_by_label = defaultdict(list)
    for idx, label in enumerate(labels):
        members_by_label[label].append(idx)

    tm = 1.0 - dist
    clusters = []
    for label in sorted(members_by_label):
        members = members_by_label[label]
        if len(members) == 1:
            centroid = members[0]
        else:
            centroid = max(
                members,
                key=lambda m: np.mean([tm[m, o] for o in members if o != m]),
            )
        ordered = [centroid] + sorted(
            (m for m in members if m != centroid), key=lambda x: names[x]
        )
        clusters.append(ordered)

    clusters.sort(key=len, reverse=True)
    return clusters, Z


def _linkage_to_newick(Z, names):
    """Convert a scipy linkage matrix to a Newick string with branch lengths."""
    tree = to_tree(Z, rd=False)

    def _build(node, parent_height):
        branch = parent_height - node.dist
        if node.is_leaf():
            return f"{names[node.id]}:{branch:.6f}"
        left = _build(node.get_left(), node.dist)
        right = _build(node.get_right(), node.dist)
        return f"({left},{right}):{branch:.6f}"

    if tree.is_leaf():
        return f"{names[tree.id]};"
    left = _build(tree.get_left(), tree.dist)
    right = _build(tree.get_right(), tree.dist)
    return f"({left},{right});"


def write_clusters(clusters, names, path):
    """Wide format: one line per cluster, members space-separated, centroid first."""
    with open(path, "w") as fh:
        for members in clusters:
            fh.write(" ".join(names[m] for m in members) + "\n")
    print(f"Wrote {path}", file=sys.stderr)


def write_newick(Z, names, path):
    with open(path, "w") as fh:
        fh.write(_linkage_to_newick(Z, names) + "\n")
    print(f"Wrote {path}", file=sys.stderr)


def write_ordered_matrix(Z, dist, names, path):
    """TM-score matrix reordered to match the dendrogram leaf order."""
    order = leaves_list(Z)
    ordered_names = [names[i] for i in order]
    ordered = (1.0 - dist)[np.ix_(order, order)]
    with open(path, "w") as fh:
        fh.write("\t" + "\t".join(ordered_names) + "\n")
        for i, name in enumerate(ordered_names):
            fh.write(name + "\t" + "\t".join(f"{v:.4f}" for v in ordered[i]) + "\n")
    print(f"Wrote {path}", file=sys.stderr)


def write_summary(clusters, names, dist, path):
    """Per-cluster statistics."""
    tm = 1.0 - dist
    with open(path, "w") as fh:
        fh.write("cluster\tsize\tcentroid\tmean_TM\tmin_TM\tmax_TM\n")
        for ci, members in enumerate(clusters, 1):
            centroid = names[members[0]]
            size = len(members)
            if size == 1:
                fh.write(f"{ci}\t{size}\t{centroid}\tNA\tNA\tNA\n")
                continue
            pair_tms = [
                tm[members[i], members[j]]
                for i in range(size) for j in range(i + 1, size)
            ]
            fh.write(
                f"{ci}\t{size}\t{centroid}\t{np.mean(pair_tms):.4f}\t"
                f"{np.min(pair_tms):.4f}\t{np.max(pair_tms):.4f}\n"
            )
    print(f"Wrote {path}", file=sys.stderr)


def run(input_tsv, output, cutoff=None, k=None, linkage_method="average",
        optimal_leaf_order=True, qcol=0, tcol=1, col=2):
    """Full clustering step. Writes output + sibling files and returns a dict
    of the paths produced (keys: clusters, summary, newick, matrix)."""
    dist, names = load_distance_matrix(input_tsv, qcol, tcol, col)
    if len(names) < 2:
        raise ValueError("need at least 2 structures to cluster")

    clusters, Z = build_clusters(
        dist, names, cutoff=cutoff, k=k,
        linkage_method=linkage_method, optimal_leaf_order=optimal_leaf_order,
    )

    sizes = [len(c) for c in clusters]
    print(
        f"Clusters: {len(clusters)}, largest: {max(sizes)}, "
        f"singletons: {sizes.count(1)}, total: {sum(sizes)}",
        file=sys.stderr,
    )

    base = output.rsplit('.', 1)[0]
    paths = {
        "clusters": output,
        "summary": base + "_summary.tsv",
        "newick": base + "_tree.nwk",
        "matrix": base + "_matrix.tsv",
    }
    write_clusters(clusters, names, paths["clusters"])
    write_summary(clusters, names, dist, paths["summary"])
    write_newick(Z, names, paths["newick"])
    write_ordered_matrix(Z, dist, names, paths["matrix"])
    return paths


def main(argv=None):
    import argparse

    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="Pairwise results TSV (from usalign_all_vs_all.bash)")
    p.add_argument("-o", "--output", default="clusters.dat", help="Output file")
    p.add_argument("--cutoff", type=float, default=None,
                   help="TM-score cutoff; cuts the tree at height = 1 - cutoff")
    p.add_argument("--k", type=int, default=None, help="Force a specific cluster count")
    p.add_argument("--linkage", default="average",
                   choices=["average", "complete", "single"], help="Linkage method")
    p.add_argument("--no-optimal-leaf-order", action="store_true",
                   help="Skip optimal leaf ordering (faster for large N)")
    p.add_argument("--qcol", type=int, default=0, help="Query-name column (0-indexed)")
    p.add_argument("--tcol", type=int, default=1, help="Target-name column (0-indexed)")
    p.add_argument("--col", type=int, default=2, help="TM-score column (0-indexed)")
    args = p.parse_args(argv)

    if args.cutoff is not None and args.k is not None:
        p.error("--cutoff and --k are mutually exclusive")

    run(args.input, args.output, cutoff=args.cutoff, k=args.k,
        linkage_method=args.linkage,
        optimal_leaf_order=not args.no_optimal_leaf_order,
        qcol=args.qcol, tcol=args.tcol, col=args.col)


if __name__ == "__main__":
    main()
