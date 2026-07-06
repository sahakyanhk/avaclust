#!/usr/bin/env python3
"""
Plot a UPGMA dendrogram + TM-score heatmap from a Newick tree and a matching
matrix (both produced by cluster.py) using seaborn.clustermap.

Importable API: plot_clustermap(newick, matrix, output, ...).
"""

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
    """Parse a binary, rooted Newick file into a SciPy-style linkage matrix.

    Returns (Z, leaf_names) where leaf_names[i] is the i-th leaf in Z.
    """
    tree = Phylo.read(newick_path, "newick")
    # Ladderize so big clusters end up at the top/left of the heatmap.
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


def _pick_scale(max_height):
    """Pick a 1-2-5 round scale-bar length near 1/5 of the max tree height."""
    target = max_height / 5
    if target <= 0:
        return 1.0
    magnitude = 10 ** math.floor(math.log10(target))
    for mult in (1, 2, 5):
        if mult * magnitude >= target:
            return mult * magnitude
    return 10 * magnitude


def _add_tree_scale(ax, scale_length, fontsize=8):
    """Draw a tree-scale bar just above the row-dendrogram axes."""
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    x0, x1 = sorted(ax.get_xlim())
    bar_x_start = x0 + 0.05 * (x1 - x0)
    bar_x_end = bar_x_start + scale_length

    ax.plot([bar_x_start, bar_x_end], [1.02, 1.02], color="black", linewidth=fontsize / 6.5,
            solid_capstyle="butt", transform=trans, clip_on=False)
    ax.text((bar_x_start + bar_x_end) / 2, 1.03, f"Tree scale: {scale_length:g}",
            ha="center", va="bottom", fontsize=fontsize, transform=trans, clip_on=False)


def plot_clustermap(newick, matrix, output="clustermap.png", cutoff=None,
                    cmap="Blues", vmin=0.0, vmax=1.0, figsize=(10, 10),
                    add_labels=False, leaf_fontsize=4, label="TM-score",
                    scale=None, linewidth=None, dpi=300):
    """Render the dendrogram + heatmap to `output`.

    figsize is used as given (it controls the output size and the heatmap aspect;
    use a square figsize for square cells).

    linewidth: dendrogram line width. None scales it with the figure so branches
    stay visible on large canvases; the save dpi is raised so it renders crisply.

    add_labels: draw leaf names at the dendrogram tips, centered on each leaf. They
    auto-shrink to the row height (--leaf-fontsize caps it), so enlarge `figsize`
    to make them bigger. Auxiliary text scales with the figure.
    """
    print(f"Parsing {newick}...", file=sys.stderr)
    Z, leaf_names = newick_to_linkage(newick)
    print(f"  {len(leaf_names)} leaves, {len(Z)} merges", file=sys.stderr)

    print(f"Loading {matrix}...", file=sys.stderr)
    df = pd.read_csv(matrix, sep="\t", index_col=0)
    missing = set(leaf_names) - set(df.index)
    if missing:
        raise SystemExit(
            f"ERROR: {len(missing)} tree leaves missing from matrix "
            f"(e.g. {sorted(missing)[0]})"
        )
    df = df.reindex(index=leaf_names, columns=leaf_names)

    n = len(leaf_names)
    figsize = tuple(figsize)               # honoured as given; controls the size
    HEATMAP_FRAC = 0.78                    # measured heatmap height fraction of the figure

    # Scale the small text (scale bar, cutoff, colorbar) with the figure so it
    # stays proportional instead of shrinking on a large --figsize.
    fontscale = min(2.5, max(1.0, max(figsize) / 8.0))

    # Dendrogram line width: scale with the figure so branches stay visible when
    # the whole figure is viewed at once (a fixed thin line vanishes on a big
    # canvas), but never thicker than the row spacing so dense trees don't smear.
    row_pts = min(figsize) * HEATMAP_FRAC * 72.0 / n
    lw = (linewidth if linewidth is not None
          else min(min(1.2, max(0.5, max(figsize) / 26.0)), max(0.15, 0.7 * row_pts)))

    # Ensure the line is at least ~1 px at the save dpi so it doesn't wash to grey.
    render_dpi = min(300, max(dpi, math.ceil(72.0 / lw)))

    print("Plotting...", file=sys.stderr)
    g = sns.clustermap(
        df, row_linkage=Z, col_cluster=False, cmap=cmap, vmin=vmin, vmax=vmax,
        figsize=figsize, xticklabels=False, yticklabels=False,
        cbar_kws={"label": label}, cbar_pos=(1, 0.01, 0.02, 0.45),
        tree_kws={"linewidths": lw},
    )
    g.ax_col_dendrogram.set_visible(False)

    if add_labels:
        # Label the dendrogram tips ourselves. We keep yticklabels off in the
        # clustermap call on purpose: with them on, seaborn reserves a strip on
        # the RIGHT of the heatmap for the labels, which leaves the heatmap
        # non-square and a wide gap before the colorbar. Instead we place the
        # labels on the LEFT (right-aligned at the heatmap edge, centered on each
        # leaf line, on a white halo) and shorten the tree by the *measured*
        # label width so no branch runs under the text. Font auto-shrinks to the
        # per-row height; --leaf-fontsize caps it.
        hm = g.ax_heatmap.get_position()
        fit = hm.height * figsize[1] * 72.0 / n   # points of vertical space per row
        fs = min(leaf_fontsize, fit)
        if fs < 2.0:
            print(f"  note: {n} labels render at {fs:.1f}pt (tiny) — increase "
                  f"--figsize height to enlarge them", file=sys.stderr)
        ordered = [df.index[i] for i in g.dendrogram_row.reordered_ind]
        g.ax_heatmap.set_yticks(np.arange(n) + 0.5)
        g.ax_heatmap.set_yticklabels(ordered, fontsize=fs)
        g.ax_heatmap.yaxis.tick_left()
        g.ax_heatmap.tick_params(axis="y", length=0, pad=1)
        labels = g.ax_heatmap.get_yticklabels()
        for lbl in labels:
            lbl.set_ha("right")     # text ends at the heatmap edge
            lbl.set_va("center")    # centered on the leaf line (default sits low)
            lbl.set_bbox(dict(facecolor="white", edgecolor="none", pad=0.2))

        # Shorten the tree so its tips stop a few points left of the widest label.
        # Measure the real rendered width (needs a draw) rather than estimating,
        # so long names can't overrun the branches.
        fig = g.figure
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        fig_w_px = fig.get_figwidth() * fig.dpi
        label_w = max(t.get_window_extent(renderer).width for t in labels) / fig_w_px
        gap = 3.0 / 72.0 / figsize[0]         # ~3 pt gap between tips and labels
        dend = g.ax_row_dendrogram.get_position()
        tree_w = max(0.02, (hm.x0 - label_w - gap) - dend.x0)
        g.ax_row_dendrogram.set_position([dend.x0, dend.y0, tree_w, dend.height])

    # Stretch the colorbar to match the heatmap's vertical extent, and scale its
    # text with the figure so it doesn't look tiny on a large --figsize.
    hm_pos = g.ax_heatmap.get_position()
    cb_pos = g.ax_cbar.get_position()
    g.ax_cbar.set_position([cb_pos.x0, hm_pos.y0, cb_pos.width, hm_pos.height])
    g.ax_cbar.tick_params(labelsize=9 * fontscale)
    g.ax_cbar.set_ylabel(label, fontsize=11 * fontscale)

    _add_tree_scale(g.ax_row_dendrogram,
                    scale if scale is not None else _pick_scale(Z[:, 2].max()),
                    fontsize=8 * fontscale)

    if cutoff is not None:
        x_dist = 1.0 - cutoff
        g.ax_row_dendrogram.axvline(x=x_dist, color="DarkSlateBlue", linestyle="--",
                                    linewidth=1.0 * fontscale)
        cutoff_trans = blended_transform_factory(
            g.ax_row_dendrogram.transData, g.ax_row_dendrogram.transAxes)
        g.ax_row_dendrogram.text(x_dist, -0.03, f"TM={cutoff:g}", ha="center",
                                 va="bottom", fontsize=7 * fontscale, color="DarkSlateBlue",
                                 transform=cutoff_trans, clip_on=False)

    g.savefig(output, dpi=render_dpi, bbox_inches="tight")
    print(f"Wrote {output}", file=sys.stderr)


def main(argv=None):
    import argparse

    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("newick", help="Newick tree file")
    p.add_argument("matrix", help="TSV similarity matrix with row/col labels")
    p.add_argument("-o", "--output", default="clustermap.png", help="Output figure path")
    p.add_argument("--cmap", default="Blues", help="Colormap")
    p.add_argument("--vmin", type=float, default=0.0)
    p.add_argument("--vmax", type=float, default=1.0)
    p.add_argument("--figsize", nargs=2, type=float, default=(10, 10),
                   metavar=("W", "H"), help="Figure size in inches")
    p.add_argument("--add-labels", action="store_true",
                   help="Show leaf names at the heatmap edge")
    p.add_argument("--leaf-fontsize", type=float, default=4)
    p.add_argument("--label", default="TM-score", help="Colorbar label")
    p.add_argument("--scale", type=float, default=None, help="Tree-scale bar length")
    p.add_argument("--linewidth", type=float, default=None,
                   help="Dendrogram line width (default: auto-thin for large trees)")
    p.add_argument("--cutoff", type=float, default=None,
                   help="Mark a TM cutoff on the dendrogram (at distance 1 - cutoff)")
    p.add_argument("--dpi", type=int, default=300)
    args = p.parse_args(argv)

    plot_clustermap(args.newick, args.matrix, output=args.output, cutoff=args.cutoff,
                    cmap=args.cmap, vmin=args.vmin, vmax=args.vmax,
                    figsize=args.figsize, add_labels=args.add_labels,
                    leaf_fontsize=args.leaf_fontsize, label=args.label,
                    scale=args.scale, linewidth=args.linewidth, dpi=args.dpi)


if __name__ == "__main__":
    main()
