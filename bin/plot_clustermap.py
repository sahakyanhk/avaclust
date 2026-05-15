#!/usr/bin/env python3
"""
Plot dendrogram + heatmap from a Newick tree and a matching similarity matrix
using seaborn.clustermap.

The matrix is expected to be a TSV with row and column labels matching the
leaf names in the Newick tree (the format produced by hierarchical_cluster.py).

Usage:
  plot_clustermap.py tree.nwk matrix.tsv -o clustermap.png
"""

import argparse
import math
import sys

import matplotlib
matplotlib.use("Agg")
from matplotlib.transforms import blended_transform_factory
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Phylo


def newick_to_linkage(newick_path):
    """
    Parse a Newick file into a SciPy-style linkage matrix.
    Assumes a binary, rooted tree (as produced by hierarchical_cluster.py).

    Returns:
        Z: linkage matrix, shape (n-1, 4)
        leaf_names: leaf names ordered so that index i is the i-th leaf in Z
    """
    tree = Phylo.read(newick_path, "newick")
    # Ladderize so the larger child is the first clade at every internal node.
    # This carries through into the linkage matrix and ultimately into
    # leaves_list(Z), so big clusters end up at the top/left of the heatmap.
    tree.ladderize(reverse=False)
    leaves = list(tree.get_terminals())
    leaf_names = [leaf.name for leaf in leaves]
    n = len(leaves)
    clade_id = {leaf: i for i, leaf in enumerate(leaves)}

    def height(clade):
        if clade.is_terminal():
            return 0.0
        return max((c.branch_length or 0.0) + height(c) for c in clade.clades)

    internals = sorted(tree.get_nonterminals(), key=height)
    Z = []
    for i, clade in enumerate(internals):
        if len(clade.clades) != 2:
            raise ValueError(
                f"Non-binary node with {len(clade.clades)} children — "
                "only strictly binary trees are supported"
            )
        c1, c2 = clade.clades
        Z.append([clade_id[c1], clade_id[c2], height(clade), len(clade.get_terminals())])
        clade_id[clade] = n + i

    return np.array(Z, dtype=float), leaf_names


def pick_scale(max_height):
    """Pick a 1-2-5 round scale length around 1/5 of the max dendrogram height."""
    target = max_height / 5
    if target <= 0:
        return 1.0
    magnitude = 10 ** math.floor(math.log10(target))
    for mult in (1, 2, 5):
        if mult * magnitude >= target:
            return mult * magnitude
    return 10 * magnitude


def add_tree_scale(ax, scale_length):
    """
    Draw a tree-scale bar just above the row-dendrogram axes.
    x is in data coordinates (so the bar length matches a real tree distance);
    y is in axes coordinates (so we always end up visually above the dendrogram,
    regardless of the dendrogram's inverted y-axis).
    """
    trans = blended_transform_factory(ax.transData, ax.transAxes)

    x0, x1 = sorted(ax.get_xlim())
    bar_x_start = x0 + 0.05 * (x1 - x0)
    bar_x_end = bar_x_start + scale_length

    bar_y = 1.02      # just above the top of the visible axes
    label_y = 1.03    # a clear gap above the bar

    ax.plot(
        [bar_x_start, bar_x_end], [bar_y, bar_y],
        color="black", linewidth=1.2, solid_capstyle="butt",
        transform=trans, clip_on=False,
    )
    ax.text(
        (bar_x_start + bar_x_end) / 2, label_y,
        f"Tree scale: {scale_length:g}",
        ha="center", va="bottom", fontsize=8,
        transform=trans, clip_on=False,
    )


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("newick", help="Newick tree file")
    p.add_argument("matrix", help="TSV similarity matrix with row/col labels")
    p.add_argument(
        "-o", "--output", default="clustermap.png",
        help="Output figure path (default: clustermap.png)"
    )
    p.add_argument("--cmap", default="Reds", help="Colormap (default: Reds)")
    p.add_argument("--vmin", type=float, default=0.0)
    p.add_argument("--vmax", type=float, default=1.0)
    p.add_argument(
        "--figsize", nargs=2, type=float, default=(10, 10),
        metavar=("W", "H"), help="Figure size in inches (default: 10 10)"
    )
    p.add_argument(
        "--add-labels", action="store_true",
        help="Show leaf names at the dendrogram tips (left side of heatmap)"
    )
    p.add_argument(
        "--leaf-fontsize", type=float, default=4,
        help="Font size for leaf-tip labels (default: 4)"
    )
    p.add_argument("--label", default="TM-score", help="Colorbar label")
    p.add_argument(
        "--scale", type=float, default=None,
        help="Tree-scale bar length (default: auto-pick a 1/2/5 round number)"
    )
    p.add_argument(
        "--cutoff", type=float, default=None,
        help="TM-score cutoff to mark on the dendrogram as a dashed line "
             "(drawn at distance = 1 - cutoff)"
    )
    p.add_argument("--dpi", type=int, default=150)
    args = p.parse_args()

    print(f"Parsing {args.newick}...", file=sys.stderr)
    Z, leaf_names = newick_to_linkage(args.newick)
    print(f"  {len(leaf_names)} leaves, {len(Z)} merges", file=sys.stderr)

    print(f"Loading {args.matrix}...", file=sys.stderr)
    df = pd.read_csv(args.matrix, sep="\t", index_col=0)

    missing = set(leaf_names) - set(df.index)
    if missing:
        sys.exit(
            f"ERROR: {len(missing)} tree leaves missing from matrix "
            f"(e.g. {sorted(missing)[0]})"
        )

    df = df.reindex(index=leaf_names, columns=leaf_names)

    print("Plotting...", file=sys.stderr)
    g = sns.clustermap(
        df,
        row_linkage=Z,
        col_cluster=False,  # matrix is already in the right column order
        cmap=args.cmap,
        vmin=args.vmin,
        vmax=args.vmax,
        figsize=tuple(args.figsize),
        xticklabels=False,
        yticklabels=args.add_labels,
        cbar_kws={"label": args.label},
        cbar_pos=(1, 0.01, 0.02, 0.45),  # x kept, y/h overridden below
    )

    # Drop the (now-empty) column dendrogram so nothing draws above the heatmap
    g.ax_col_dendrogram.set_visible(False)

    # Move y labels to the left (heatmap edge nearest the tree tips), and put a
    # white background behind each label so the dendrogram leaf branches don't
    # show through the text.
    if args.add_labels:
        g.ax_heatmap.yaxis.tick_left()
        g.ax_heatmap.tick_params(
            axis="y", labelsize=args.leaf_fontsize, length=0, pad=1, rotation=0,
        )
        for lbl in g.ax_heatmap.get_yticklabels():
            lbl.set_bbox(dict(facecolor="white", edgecolor="none", pad=0.5))

    # Stretch the colorbar to match the heatmap's vertical extent
    hm_pos = g.ax_heatmap.get_position()
    cb_pos = g.ax_cbar.get_position()
    g.ax_cbar.set_position([cb_pos.x0, hm_pos.y0, cb_pos.width, hm_pos.height])

    scale = args.scale if args.scale is not None else pick_scale(Z[:, 2].max())
    add_tree_scale(g.ax_row_dendrogram, scale)

    if args.cutoff is not None:
        x_dist = 1.0 - args.cutoff
        g.ax_row_dendrogram.axvline(
            x=x_dist,
            color="red", linestyle="--", linewidth=1.0,
        )
        cutoff_trans = blended_transform_factory(
            g.ax_row_dendrogram.transData,
            g.ax_row_dendrogram.transAxes,
        )
        g.ax_row_dendrogram.text(
            x_dist, 1.02,
            f"TM={args.cutoff:g}",
            ha="center", va="bottom", fontsize=8, color="red",
            transform=cutoff_trans, clip_on=False,
        )

    g.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    print(f"Wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
