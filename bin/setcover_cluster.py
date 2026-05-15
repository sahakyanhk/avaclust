#!/usr/bin/env python3
"""
Greedy set-cover clustering of pairwise structural similarity results.

Input: TSV file (no header, tab-separated)

Algorithm:
  1. Build a graph: edge between A-B if score passes all filters
  2. Pick the node with the most unclustered neighbors → representative
  3. Assign all its unclustered neighbors to its cluster
  4. Repeat until every node is assigned
  5. (Optional) Compute centroid per cluster: member with highest avg score

Optimizations for large files (100M+ lines):
  - Integer node IDs instead of string comparisons
  - Buffered I/O with large read buffer
  - Max-heap for O(1) best-node lookup instead of O(n) scan
  - Lazy deletion from heap (skip stale entries)

Usage:
  # Basic: cluster on column 4 (TM-score >= 0.5)
  python setcover_cluster.py input.tsv --col 4 --cutoff 0.5 -o clusters.tsv

  # With extra filters: TM-score >= 0.5 AND last two columns <= 1.0 and <= 0.95
  python setcover_cluster.py input.tsv --col 4 --cutoff 0.5 \\
      --filter1 5 --cutoff1 1.0 --filter2 6 --cutoff2 0.95 -o clusters.tsv

  # Use centroid (highest avg within-cluster score) instead of greedy rep
  python setcover_cluster.py input.tsv --col 4 --cutoff 0.5 --centroid
"""

import argparse
import heapq
import os
import subprocess
import sys
import time
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def parse_args():
    p = argparse.ArgumentParser(
        description="Greedy set-cover clustering on pairwise similarity results"
    )
    p.add_argument("input", help="Pairwise results TSV")
    p.add_argument(
        "--qcol", type=int, default=0,
        help="0-indexed column for query name (default: 0)"
    )
    p.add_argument(
        "--tcol", type=int, default=1,
        help="0-indexed column for target name (default: 1)"
    )

    # Primary filter
    p.add_argument(
        "--col", type=int, required=True,
        help="0-indexed column for primary score"
    )
    p.add_argument(
        "--cutoff", type=float, required=True,
        help="Primary threshold"
    )
    p.add_argument(
        "--gt", action="store_true", default=True,
        help="Primary filter: value >= cutoff (default)"
    )
    p.add_argument(
        "--lt", action="store_true", default=False,
        help="Primary filter: value <= cutoff (e.g. for E-values)"
    )

    # Optional secondary filters (both use <=)
    p.add_argument(
        "--filter1", type=int, default=None,
        help="0-indexed column for extra filter 1"
    )
    p.add_argument(
        "--cutoff1", type=float, default=None,
        help="Threshold for extra filter 1 (value <= cutoff1)"
    )
    p.add_argument(
        "--filter2", type=int, default=None,
        help="0-indexed column for extra filter 2"
    )
    p.add_argument(
        "--cutoff2", type=float, default=None,
        help="Threshold for extra filter 2 (value <= cutoff2)"
    )

    p.add_argument(
        "-o", "--output", default=None,
        help="Output file (default: stdout)"
    )
    p.add_argument(
        "--fmt", choices=["tsv", "wide"], default="wide",
        help="Output format: 'tsv' = one row per member (centroid<tab>member), "
             "'wide' = one row per cluster (centroid<tab>member1<tab>member2...)"
    )
    p.add_argument(
        "--centroid", action="store_true", default=False,
        help="Compute centroid (highest avg within-cluster score) as cluster rep. "
             "Requires a second pass over the input file."
    )
    p.add_argument(
        "--centroid-col", type=int, default=None,
        help="0-indexed column for centroid score (default: same as --col). "
             "Useful when clustering on one metric but picking centroids by another."
    )
    p.add_argument(
        "-sc", "--sort-clusters", default=None,
        help="Directory for aligned cluster PDBs"
    )
    p.add_argument(
        "--pdb-dir", default=None,
        help="PDB directory (required with --sort-clusters to resolve basenames)"
    )
    args = p.parse_args()

    if (args.filter1 is None) != (args.cutoff1 is None):
        p.error("--filter1 and --cutoff1 must be used together")
    if (args.filter2 is None) != (args.cutoff2 is None):
        p.error("--filter2 and --cutoff2 must be used together")
    if args.centroid_col is not None and not args.centroid:
        p.error("--centroid-col requires --centroid")
    if args.centroid and args.centroid_col is None:
        args.centroid_col = args.col
    if args.sort_clusters and args.fmt == "tsv":
        p.error("--sort-clusters requires wide format (incompatible with --fmt tsv)")
    if args.sort_clusters and not args.pdb_dir:
        p.error("--sort-clusters requires --pdb-dir")

    return args


def load_graph(args):
    """
    Read pairwise file and build adjacency lists using integer node IDs.

    Returns:
        neighbors: list of sets, neighbors[i] = set of neighbor IDs for node i
        id_to_name: list mapping int ID -> string name
        name_to_id: dict mapping string name -> int ID
        n_nodes: total number of unique nodes
    """
    name_to_id = {}
    neighbors = []
    n_lines = 0
    n_edges = 0

    col = args.col
    qcol = args.qcol
    tcol = args.tcol
    cutoff = args.cutoff
    use_lt = args.lt
    has_f1 = args.filter1 is not None
    has_f2 = args.filter2 is not None
    filter1 = args.filter1
    cutoff1 = args.cutoff1
    filter2 = args.filter2
    cutoff2 = args.cutoff2

    t0 = time.time()
    n_skipped = 0

    with open(args.input, buffering=1 << 22) as fh:
        for line in fh:
            if not line or line[0] == '#':
                continue

            parts = line.split('\t')
            try:
                val = float(parts[col])
            except (IndexError, ValueError):
                n_skipped += 1
                if n_skipped == 1:
                    ncols = len(parts)
                    print(
                        f"WARNING: first skipped line has {ncols} columns "
                        f"(indices 0-{ncols-1}), but --col={col} requested. "
                        f"Check your column index!",
                        file=sys.stderr,
                    )
                continue

            n_lines += 1

            if use_lt:
                if val > cutoff:
                    continue
            else:
                if val < cutoff:
                    continue

            if has_f1:
                try:
                    if float(parts[filter1]) > cutoff1:
                        continue
                except (IndexError, ValueError):
                    continue

            if has_f2:
                try:
                    if float(parts[filter2]) > cutoff2:
                        continue
                except (IndexError, ValueError):
                    continue

            query = parts[qcol]
            target = parts[tcol]
            if query == target:
                continue

            qid = name_to_id.get(query)
            if qid is None:
                qid = len(name_to_id)
                name_to_id[query] = qid
                neighbors.append(set())

            tid = name_to_id.get(target)
            if tid is None:
                tid = len(name_to_id)
                name_to_id[target] = tid
                neighbors.append(set())

            neighbors[qid].add(tid)
            neighbors[tid].add(qid)
            n_edges += 1

    t1 = time.time()

    id_to_name = [''] * len(name_to_id)
    for name, nid in name_to_id.items():
        id_to_name[nid] = name

    n_nodes = len(name_to_id)

    desc = f"col{col} {'<=' if use_lt else '>='} {cutoff}"
    if has_f1:
        desc += f", col{filter1} <= {cutoff1}"
    if has_f2:
        desc += f", col{filter2} <= {cutoff2}"

    print(
        f"Loaded {n_lines:,} lines in {t1-t0:.1f}s, "
        f"{n_nodes:,} proteins, {n_edges:,} edges [{desc}]",
        file=sys.stderr,
    )
    return neighbors, id_to_name, name_to_id, n_nodes


def load_all_nodes(args):
    """Quick pass to find all unique node names (for singletons)."""
    all_names = set()
    qcol = args.qcol
    tcol = args.tcol

    with open(args.input, buffering=1 << 22) as fh:
        for line in fh:
            if not line or line[0] == '#':
                continue
            parts = line.split('\t')
            try:
                all_names.add(parts[qcol])
                all_names.add(parts[tcol])
            except IndexError:
                continue

    return all_names


def greedy_setcover(neighbors, n_nodes):
    """
    Greedy set-cover using a max-heap for O(1) best-node selection.
    Returns list of (rep_id, [member_ids]).
    """
    clustered = bytearray(n_nodes)
    clusters = []

    heap = []
    for i in range(n_nodes):
        deg = len(neighbors[i]) + 1
        heap.append((-deg, i))
    heapq.heapify(heap)

    n_remaining = n_nodes
    t0 = time.time()

    while n_remaining > 0:
        while heap:
            neg_deg, node = heapq.heappop(heap)
            if clustered[node]:
                continue

            actual_nbrs = [nb for nb in neighbors[node] if not clustered[nb]]
            actual_deg = len(actual_nbrs) + 1

            if -neg_deg != actual_deg:
                heapq.heappush(heap, (-actual_deg, node))
                continue

            members = [node] + actual_nbrs
            for m in members:
                clustered[m] = 1
            n_remaining -= len(members)
            clusters.append((node, members))
            break

    elapsed = time.time() - t0
    print(f"Clustering took {elapsed:.1f}s", file=sys.stderr)

    clusters.sort(key=lambda c: len(c[1]), reverse=True)
    return clusters


def compute_centroids(clusters, id_to_name, name_to_id, args):
    """
    Second pass over input file to find the centroid of each cluster.
    Centroid = member with highest average score to other cluster members.

    Returns dict: greedy_rep_id -> centroid_id
    """
    # Build member -> cluster index mapping
    member_to_cidx = {}
    for cidx, (_, members) in enumerate(clusters):
        for m in members:
            member_to_cidx[m] = cidx

    # Accumulate scores: (node_id) -> [scores to cluster-mates]
    score_sums = defaultdict(float)
    score_counts = defaultdict(int)

    qcol = args.qcol
    tcol = args.tcol
    scol = args.centroid_col

    t0 = time.time()
    with open(args.input, buffering=1 << 22) as fh:
        for line in fh:
            if not line or line[0] == '#':
                continue
            parts = line.split('\t')
            try:
                a_name = parts[qcol]
                b_name = parts[tcol]
            except IndexError:
                continue

            a_id = name_to_id.get(a_name)
            b_id = name_to_id.get(b_name)
            if a_id is None or b_id is None:
                continue

            ca = member_to_cidx.get(a_id)
            cb = member_to_cidx.get(b_id)
            if ca is None or ca != cb or a_id == b_id:
                continue

            try:
                val = float(parts[scol])
            except (IndexError, ValueError):
                continue

            score_sums[a_id] += val
            score_counts[a_id] += 1
            score_sums[b_id] += val
            score_counts[b_id] += 1

    # Pick centroid per cluster
    centroids = {}
    for cidx, (rep, members) in enumerate(clusters):
        if len(members) == 1:
            centroids[rep] = rep
            continue
        best_id, best_avg = rep, 0.0
        for m in members:
            cnt = score_counts[m]
            if cnt > 0:
                avg = score_sums[m] / cnt
                if avg > best_avg:
                    best_id, best_avg = m, avg
        centroids[rep] = best_id

    elapsed = time.time() - t0
    print(f"Centroid computation took {elapsed:.1f}s", file=sys.stderr)
    return centroids


def write_summary(clusters, id_to_name, name_to_id, singleton_names,
                   summary_path, args, centroids=None):
    """Write cluster_summary.tsv with per-cluster statistics."""
    member_to_cidx = {}
    for cidx, (_, members) in enumerate(clusters):
        for m in members:
            member_to_cidx[m] = cidx

    # Collect intra-cluster pairwise scores
    pair_scores = defaultdict(list)  # cidx -> [scores]
    scol = args.centroid_col if args.centroid_col is not None else args.col

    with open(args.input, buffering=1 << 22) as fh:
        for line in fh:
            if not line or line[0] == '#':
                continue
            parts = line.split('\t')
            try:
                a_id = name_to_id.get(parts[args.qcol])
                b_id = name_to_id.get(parts[args.tcol])
            except IndexError:
                continue
            if a_id is None or b_id is None:
                continue
            ca = member_to_cidx.get(a_id)
            cb = member_to_cidx.get(b_id)
            if ca is None or ca != cb or a_id == b_id:
                continue
            try:
                val = float(parts[scol])
            except (IndexError, ValueError):
                continue
            pair_scores[ca].append(val)

    with open(summary_path, "w") as fh:
        fh.write("cluster\tsize\tcentroid\tmean_TM\tmin_TM\tmax_TM\n")
        ci = 0
        for cidx, (rep, members) in enumerate(clusters):
            ci += 1
            c = centroids[rep] if centroids else rep
            centroid_name = id_to_name[c]
            size = len(members)
            scores = pair_scores.get(cidx, [])
            if size == 1 or not scores:
                fh.write(f"{ci}\t{size}\t{centroid_name}\tNA\tNA\tNA\n")
            else:
                fh.write(
                    f"{ci}\t{size}\t{centroid_name}\t"
                    f"{sum(scores)/len(scores):.4f}\t{min(scores):.4f}\t{max(scores):.4f}\n"
                )
        for i, name in enumerate(sorted(singleton_names), ci + 1):
            fh.write(f"{i}\t1\t{name}\tNA\tNA\tNA\n")

    print(f"Wrote {summary_path}", file=sys.stderr)


def write_output(clusters, id_to_name, singleton_names, out, fmt, centroids=None):
    """Write cluster assignments. If centroids provided, use them as reps."""
    if fmt == "tsv":
        for rep, members in clusters:
            c = centroids[rep] if centroids else rep
            c_name = id_to_name[c]
            # Sort members but put centroid first
            others = sorted(
                [m for m in members if m != c], key=lambda x: id_to_name[x]
            )
            ordered = [c] + others
            for m in ordered:
                out.write(f"{c_name}\t{id_to_name[m]}\n")
        for name in sorted(singleton_names):
            out.write(f"{name}\t{name}\n")
    else:  # wide
        for rep, members in clusters:
            c = centroids[rep] if centroids else rep
            c_name = id_to_name[c]
            others = sorted(
                (id_to_name[m] for m in members if m != c)
            )
            out.write(" ".join([c_name] + others) + "\n")
        for name in sorted(singleton_names):
            out.write(f"{name}\n")


def main():
    args = parse_args()
    t_start = time.time()

    neighbors, id_to_name, name_to_id, n_nodes = load_graph(args)

    print("Scanning for singletons...", file=sys.stderr)
    all_names = load_all_nodes(args)
    graph_names = set(id_to_name)
    singleton_names = all_names - graph_names
    n_total = len(all_names)
    n_isolated = len(singleton_names)
    print(
        f"{n_total:,} total proteins, {n_isolated:,} isolated "
        f"(no edges passing filters)",
        file=sys.stderr,
    )

    clusters = greedy_setcover(neighbors, n_nodes)

    centroids = None
    if args.centroid:
        print("Computing centroids...", file=sys.stderr)
        centroids = compute_centroids(clusters, id_to_name, name_to_id, args)
        n_changed = sum(
            1 for rep, _ in clusters
            if centroids[rep] != rep and len(_) > 1
        )
        print(
            f"Centroid differs from greedy rep in {n_changed}/{len(clusters)} "
            f"multi-member clusters",
            file=sys.stderr,
        )

    sizes = [len(m) for _, m in clusters]
    n_clusters = len(clusters) + n_isolated
    largest = max(sizes) if sizes else 1
    n_graph_singletons = sizes.count(1)
    total_members = sum(sizes) + n_isolated
    print(
        f"Clusters: {n_clusters:,}, largest: {largest}, "
        f"singletons: {n_graph_singletons + n_isolated:,} "
        f"({n_isolated} isolated + {n_graph_singletons} left over after greedy), "
        f"mean size: {total_members/n_clusters:.1f}",
        file=sys.stderr,
    )

    fh = open(args.output, "w") if args.output else sys.stdout
    write_output(clusters, id_to_name, singleton_names, fh, args.fmt, centroids)
    if args.output:
        fh.close()
        print(f"Wrote {args.output}", file=sys.stderr)
        summary_path = args.output.rsplit('.', 1)[0] + '_summary.tsv'
        write_summary(clusters, id_to_name, name_to_id, singleton_names,
                      summary_path, args, centroids)

    print(f"Total time: {time.time()-t_start:.1f}s", file=sys.stderr)

    if args.sort_clusters:
        if not args.output:
            print("ERROR: --sort-clusters requires --output", file=sys.stderr)
            sys.exit(1)
        align_script = os.path.join(SCRIPT_DIR, "align_clusters.bash")
        ret = subprocess.run(
            ["bash", align_script, args.output, args.sort_clusters, args.pdb_dir],
        )
        if ret.returncode != 0:
            print(f"ERROR: align_clusters.bash exited with code {ret.returncode}", file=sys.stderr)
            sys.exit(ret.returncode)

if __name__ == "__main__":
    main()
