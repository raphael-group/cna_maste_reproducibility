#!/usr/bin/env python
"""Draw the cell-typing figures of a cell-typing output folder to the lab figure rules.

Author (2026-11-04)

Dependencies:
    python 3.11, scanpy 1.11, anndata 0.12, pandas, numpy, matplotlib 3.11 (configs/envs/celltyping.yaml)

Usage:
    python plot_celltypes.py <out_dir> --scalefactors FILE [--figures-dir DIR]
        out_dir         # cell-typing output folder holding the inputs below
        --scalefactors  # Space Ranger spatial/scalefactors_json.json of the sample, any bin size
        --figures-dir   # output folder, created if missing; default <out_dir>/figures-refined

Inputs:
    <out_dir>/*.h5ad: the one annotated AnnData; obs cluster, cell_type, cluster_label; obsm X_umap,
        score_ulm, spatial (full-resolution pixels); uns cluster_colors (Space Ranger web summary colours), if present.
    <out_dir>/cluster_annotation.tsv: one row per cluster: panglaodb_cell_type and support_genes
        (comma-separated, by descending log2 fold change).
    <out_dir>/cluster_celltype_ranking.tsv: decoupler rankby_group per cluster (group, name, stat),
        in decoupler order.
    <out_dir>/cluster_top_features.tsv: Space Ranger top features per cluster (cluster, rank, gene).
    scalefactors: microns_per_pixel, the full-resolution image's pixel side, the same for every bin size.
Parameters:
    N_MARKERS: supporting markers drawn per cell type.
    N_TOP_FEATURES: top features drawn per cluster.
    N_TOP_CELL_TYPES: ranked cell types drawn per cluster.
    FONT_SIZE, LABEL_SIZE: points of every text and of the axis labels, at PAGE_WIDTH.
    PAGE_WIDTH: inches of the paper's text column (cnamaste.figure_style.PAPER_WIDTH).
    PAGE_HEIGHT: inches of a full-page figure, its caption's room left (cnamaste.figure_style.page_size).
    DPI: resolution of the rasterized bin layers.
    KEY_WIDTH: inches of the matrix plot's key column.
    DOT_AREA: points squared of a dot plot's largest dot, under the gene spacing.
    PAD: inches of margin round a cropped page.
    SPATIAL_WIDTH: inches of a spatial map, one scale for every map, its page under PAGE_HEIGHT.
    ROW_PITCH: inches between dot plot rows, room for labels at 45 degrees.
    KEY_COLUMNS: columns of a spatial map's key.
Outputs (figures dir):
    umap_cluster.pdf, umap_cell_type.pdf: Space Ranger UMAP of the bins, coloured by cluster and by cell type.
    spatial_cluster.pdf, spatial_cell_type.pdf: bins at their physical position, coloured by cluster and by cell type.
    dotplot__celltype_markers.pdf: supporting PanglaoDB markers of each cluster's cell type.
    dotplot__10x_top_features.pdf: Space Ranger top features of each cluster.
    matrixplot__ulm_top_celltypes.pdf: mean ULM score of each cluster's top-ranked cell types.
Notes/References:
    Rules: the lab's figure rules, by ID. Text is Computer Modern (cmr10, MP-1), saved as TrueType (pdf.fonttype 42,
    MP-3), FONT_SIZE throughout and LABEL_SIZE for axis labels (FA-2, FL-2, MP-5); keys sit below the dot plots and
    spatial maps, right of the UMAPs and matrix plot (FL-1); bins, > 10K points, are rasterized at DPI and the text
    stays vector (FX-3); pages are cropped, at most PAGE_WIDTH by PAGE_HEIGHT (FP-5).
    A cell type takes the colour of its lowest-numbered cluster, so cluster and cell-type figures share colours
    (FC-1). Dot and matrix plots scale each gene or cell type to [0, 1] (scanpy standard_scale "var",
    https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.dotplot.html) in Reds, zero white (FH-3); dot plot
    values are scanpy's DotPlot tables, drawn here so the keys go below.
    The matrix plot follows the decoupler tutorial
    (https://decoupler.scverse.org/en/latest/notebooks/scell/rna_sc.html).
    A cluster without supporting markers (low confidence) has no marker group.
    Spatial maps place bins as cnamaste's spatial panels do (cnamaste.platforms.visium.pixel_physical): x the pixel
    column, y the pixel row negated, times microns_per_pixel; framed from the slide origin on 1,000 um steps, every
    other tick labelled in thousands (cnamaste.plot.layout.frame_extent). A bin's array row and column, read from its
    Space Ranger HD barcode s_<side>um_<row>_<col>-1, lay out the bin image, placed by an affine map fitted on every
    bin's centre.
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, to_rgba_array
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.transforms import Affine2D
from scanpy.plotting import palettes

N_MARKERS = 5
N_TOP_FEATURES = 3
N_TOP_CELL_TYPES = 3
FONT_SIZE = 10.0
LABEL_SIZE = 12.5
PAGE_WIDTH = 8.5 - 2 * 1.01
PAGE_HEIGHT = 11.0 - 2 * 1.01 - 58.0 / 72.27 - 1.0
DPI = 600
CMAP = "Reds"
KEY_WIDTH = 1.5
DOT_AREA = 80
SIZE_EXPONENT = 1.5
DOT_EDGE = 0.2
ROW_PITCH = 0.25
BRACKET_HEIGHT = 0.12
KEY_GAP = 0.15
KEY_TITLE = 0.25
KEY_HEIGHT = 0.3
KEY_COLUMNS = 2
TICK_STEP_UM = 1000.0
FRAME_PAD = 0.02
SPATIAL_WIDTH = 4.9
CANVAS = 14.0
AXES_ORIGIN = (3.0, 5.0)
CLUSTER_LABEL = "Space Ranger cluster"
BARCODE = r"^s_(\d+)um_(\d+)_(\d+)-\d+$"
PAD = 0.02


def set_figure_style():
    """Set every figure's face, text sizes and PDF font type (MP-1, MP-3, MP-5)."""
    # NB: matplotlib's cmr10.ttf has head timestamps fontTools warns about on every save
    logging.getLogger("fontTools").setLevel(logging.ERROR)
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["cmr10"],
            "mathtext.fontset": "cm",
            "axes.formatter.use_mathtext": True,
            "font.size": FONT_SIZE,
            "axes.labelsize": LABEL_SIZE,
            "axes.titlesize": LABEL_SIZE,
            "xtick.labelsize": FONT_SIZE,
            "ytick.labelsize": FONT_SIZE,
            "legend.fontsize": FONT_SIZE,
            "legend.title_fontsize": FONT_SIZE,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def set_text_sizes(fig):
    """Set every text of a figure to FONT_SIZE and its axis labels to LABEL_SIZE, over scanpy's relative sizes.

    Args:
        fig: the figure.
    """
    for text in fig.findobj(Text):
        text.set_fontsize(FONT_SIZE)
    for ax in fig.axes:
        ax.xaxis.label.set_fontsize(LABEL_SIZE)
        ax.yaxis.label.set_fontsize(LABEL_SIZE)


def save_figure(fig, path):
    """Save a figure cropped to its content, on a transparent page, then close it (FP-5, MP-3).

    Args:
        fig: the figure.
        path: output PDF.
    """
    fig.savefig(path, dpi=DPI, bbox_inches="tight", pad_inches=PAD, transparent=True)
    plt.close(fig)
    print(f"[plot] {path}")


def get_cluster_colours(adata):
    """Get each cluster's colour: uns cluster_colors (Space Ranger's web summary), else scanpy's default palette for
    that many categories.

    Args:
        adata: AnnData with categorical obs cluster.

    Returns:
        Hex colour per cluster, in cluster order.
    """
    clusters = list(adata.obs["cluster"].cat.categories)
    k = len(clusters)
    default = palettes.default_20 if k <= 20 else palettes.default_28 if k <= 28 else palettes.default_102
    colours = adata.uns.get("cluster_colors", default[:k])
    return dict(zip(clusters, colours, strict=True))


def get_cell_type_colours(adata, cluster_colours):
    """Get each cell type's colour: its lowest-numbered cluster's (FC-1).

    Args:
        adata: AnnData with categorical obs cluster and cell_type.
        cluster_colours: output of get_cluster_colours.

    Returns:
        Hex colour per cell type, in the order of its first cluster.
    """
    cell_types = adata.obs.groupby("cluster", observed=True)["cell_type"].first()
    colours = {}
    for cluster, cell_type in cell_types.items():
        colours.setdefault(cell_type, cluster_colours[cluster])
    return colours


def get_marker_groups(markers, n_markers):
    """Get the first n_markers supporting markers of each cell type, from its lowest-numbered cluster.

    Args:
        markers: one row per supporting marker: cluster, panglaodb_cell_type, gene; by cluster, then log2FC.
        n_markers: markers per cell type.

    Returns:
        Marker genes per PanglaoDB cell type, in cluster order.
    """
    first = markers.drop_duplicates("panglaodb_cell_type")["cluster"]
    kept = markers[markers["cluster"].isin(first)]
    return {t: list(g)[:n_markers] for t, g in kept.groupby("panglaodb_cell_type", sort=False)["gene"]}


def get_key_handles(colours):
    """Get a key's square handles, one per category, in the bins' colours.

    Args:
        colours: hex colour per category.

    Returns:
        One handle per category.
    """
    return [
        Line2D([], [], linestyle="none", marker="s", markersize=FONT_SIZE * 0.8, color=c, label=k)
        for k, c in colours.items()
    ]


def add_key(ax, colours, title):
    """Put a categorical key right of the axes, centred on them (FL-1).

    Args:
        ax: the axes.
        colours: hex colour per category.
        title: the key's title.
    """
    ax.legend(
        handles=get_key_handles(colours), title=title, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False
    )


def plot_umap(adata, key, colours, key_title, path):
    """Draw the bins' Space Ranger UMAP coloured by an obs column, its key right of the axes.

    Args:
        adata: AnnData with obsm X_umap and obs key.
        key: obs column, e.g. cluster.
        colours: hex colour per category of key.
        key_title: the key's title.
        path: output PDF.
    """
    fig, ax = plt.subplots(figsize=(PAGE_WIDTH, 4.0), layout="constrained")
    umap = adata.obsm["X_umap"]
    ax.scatter(
        umap[:, 0],
        umap[:, 1],
        c=adata.obs[key].map(colours).astype(str),
        s=120_000 / adata.n_obs,
        linewidths=0,
        rasterized=True,
    )
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2", xticks=[], yticks=[], aspect="equal")
    ax.spines[["top", "right"]].set_visible(False)
    add_key(ax, colours, key_title)
    save_figure(fig, path)


def get_bin_grid(obs_names):
    """Get each bin's side and array position from its Space Ranger HD barcode.

    Args:
        obs_names: barcodes, s_<side>um_<row>_<col>-1.

    Returns:
        DataFrame with integer columns side_um, row, col, one row per barcode.
    """
    grid = obs_names.str.extract(BARCODE)
    if grid.isna().any().any():
        sys.exit("a barcode is not s_<side>um_<row>_<col>-1")
    grid = grid.astype(int)
    grid.columns = ["side_um", "row", "col"]
    if grid["side_um"].nunique() != 1:
        sys.exit("bins of more than one size")
    return grid


def add_axes_inches(fig, left, bottom, width, height):
    """Add axes to a CANVAS-square figure at a position in inches.

    Args:
        fig: the figure, CANVAS inches a side.
        left, bottom, width, height: inches.

    Returns:
        The axes.
    """
    return fig.add_axes((left / CANVAS, bottom / CANVAS, width / CANVAS, height / CANVAS))


def get_bottom_inches(fig, ax):
    """Get the bottom of an axes' x axis, its tick labels and label included, in inches.

    Args:
        fig: the figure.
        ax: the axes.

    Returns:
        Inches from the figure's bottom.
    """
    box = ax.xaxis.get_tightbbox(fig.canvas.get_renderer())
    return float(box.transformed(fig.dpi_scale_trans.inverted()).y0)


def get_physical_um(adata, microns_per_pixel):
    """Get each bin's centre in um as cnamaste places it: x the full-resolution pixel column, y the pixel row
    negated, so y runs up (cnamaste.platforms.visium.pixel_physical).

    Args:
        adata: AnnData with obsm spatial, (pxl_col_in_fullres, pxl_row_in_fullres) per bin.
        microns_per_pixel: Space Ranger's scalefactors_json.json value.

    Returns:
        (n_bins, 2) array of (x, y) in um.
    """
    pixels = np.asarray(adata.obsm["spatial"], dtype=float)
    return np.column_stack([pixels[:, 0], -pixels[:, 1]]) * microns_per_pixel


def get_grid_affine(grid, xy):
    """Get the affine map from a bin's (col, row) on the array to its centre in um, fitted on every bin.

    Args:
        grid: output of get_bin_grid.
        xy: output of get_physical_um.

    Returns:
        The map, col and row to x and y.
    """
    design = np.column_stack([grid["col"], grid["row"], np.ones(len(grid))])
    fit = np.linalg.lstsq(design, xy, rcond=None)[0]
    return Affine2D(np.vstack([fit.T, [0.0, 0.0, 1.0]]))


def get_spatial_ticks(low, high):
    """Get ticks every TICK_STEP_UM from the step at or below low to the one at or above high
    (cnamaste.plot.layout.framed).

    Args:
        low, high: the least and greatest coordinate, um.

    Returns:
        The ticks, um.
    """
    return np.arange(np.floor(low / TICK_STEP_UM), np.ceil(high / TICK_STEP_UM) + 1) * TICK_STEP_UM


def get_spatial_tick_labels(ticks):
    """Get every other tick's label from 0, in thousands, the outermost labelled one ending in K
    (cnamaste.plot.layout.coordinates and get_alternate_tick_labels).

    Args:
        ticks: ticks holding 0, um.

    Returns:
        One label per tick, every other one empty.
    """
    zero = int(np.flatnonzero(ticks == 0)[0])
    shown = [i for i in range(len(ticks)) if (i - zero) % 2 == 0]
    outer = max(shown, key=lambda i: abs(ticks[i]))
    return [f"{ticks[i] / 1000:g}" + ("K" if i == outer else "") if i in shown else "" for i in range(len(ticks))]


def plot_spatial(adata, key, colours, key_title, microns_per_pixel, path):
    """Draw each bin as a square at its physical position coloured by an obs column, framed as cnamaste's spatial
    panels from the slide origin, SPATIAL_WIDTH wide, its key below in KEY_COLUMNS columns as wide as the axes.

    Args:
        adata: AnnData with Space Ranger HD barcodes as obs_names, obsm spatial and obs key.
        key: obs column, e.g. cell_type.
        colours: hex colour per category of key.
        key_title: the key's title.
        microns_per_pixel: Space Ranger's scalefactors_json.json value.
        path: output PDF.
    """
    grid = get_bin_grid(adata.obs_names)
    xy = get_physical_um(adata, microns_per_pixel)
    image = np.zeros((grid["row"].max() + 1, grid["col"].max() + 1, 4))
    image[grid["row"], grid["col"]] = to_rgba_array(adata.obs[key].map(colours).astype(str))
    xs = get_spatial_ticks(min(xy[:, 0].min(), 0.0), max(xy[:, 0].max(), 0.0))
    ys = get_spatial_ticks(min(xy[:, 1].min(), 0.0), max(xy[:, 1].max(), 0.0))
    xlim = (xs[0] - FRAME_PAD * (xs[-1] - xs[0]), xs[-1] + FRAME_PAD * (xs[-1] - xs[0]))
    ylim = (ys[0] - FRAME_PAD * (ys[-1] - ys[0]), ys[-1] + FRAME_PAD * (ys[-1] - ys[0]))
    left, bottom = AXES_ORIGIN
    fig = plt.figure(figsize=(CANVAS, CANVAS))
    ax = add_axes_inches(fig, left, bottom, SPATIAL_WIDTH, SPATIAL_WIDTH * np.diff(ylim)[0] / np.diff(xlim)[0])
    im = ax.imshow(image, origin="lower", interpolation="nearest", rasterized=True)
    im.set_transform(get_grid_affine(grid, xy) + ax.transData)
    ax.set(xlim=xlim, ylim=ylim, aspect="equal", xlabel=r"$x$ ($\mu$m)", ylabel=r"$y$ ($\mu$m)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines["bottom"].set_bounds(xs[0], xs[-1])
    ax.spines["left"].set_bounds(ys[0], ys[-1])
    ax.set_xticks(xs, get_spatial_tick_labels(xs))
    ax.set_yticks(ys, get_spatial_tick_labels(ys))
    set_text_sizes(fig)
    add_key_below(fig, ax, colours, key_title, left, SPATIAL_WIDTH)
    save_figure(fig, path)


def add_key_below(fig, ax, colours, title, left, width):
    """Put a categorical key below an axes' x axis in KEY_COLUMNS columns, each starting at its share of the axes'
    width, filled top to bottom, the title over the first.

    Args:
        fig: the figure.
        ax: the axes.
        colours: hex colour per category.
        title: the key's title.
        left, width: inches of the axes' left edge and width.
    """
    handles = get_key_handles(colours)
    per_column = -(-len(handles) // KEY_COLUMNS)
    top = (get_bottom_inches(fig, ax) - KEY_GAP) / CANVAS
    for i in range(KEY_COLUMNS):
        fig.legend(
            handles=handles[i * per_column : (i + 1) * per_column],
            title=title if i == 0 else " ",
            loc="upper left",
            bbox_to_anchor=((left + i * width / KEY_COLUMNS) / CANVAS, top),
            alignment="left",
            borderpad=0,
            borderaxespad=0,
            handlelength=0.8,
            frameon=False,
        )


def draw_page_fit(draw):
    """Draw a figure whose labels fall outside its axes, then redraw it, its width adjusted by the secant method,
    until the saved page is PAGE_WIDTH wide, or PAGE_HEIGHT tall where that binds first.

    Args:
        draw: function of a width in inches (the figure's or its main axes'), returning the drawn figure.

    Returns:
        The redrawn figure.
    """
    target_width, target_height = PAGE_WIDTH - 2 * PAD, PAGE_HEIGHT - 2 * PAD

    def get_excess(fig):
        box = fig.get_tightbbox()
        return max(box.width - target_width, box.height - target_height)

    widths = [target_width]
    fig = draw(widths[-1])
    excesses = [get_excess(fig)]
    for _ in range(6):
        if abs(excesses[-1]) < 0.01:
            break
        slope = 1.0 if len(widths) < 2 else (excesses[-1] - excesses[-2]) / (widths[-1] - widths[-2])
        plt.close(fig)
        widths.append(widths[-1] - excesses[-1] / (slope if slope > 0.1 else 1.0))
        fig = draw(widths[-1])
        excesses.append(get_excess(fig))
    return fig


def finish_scanpy_plot(plot, xlabel):
    """Render a scanpy dot or matrix plot, label its axes and set its text sizes.

    Args:
        plot: scanpy DotPlot or MatrixPlot, its legend set.
        xlabel: the main axes' x label.

    Returns:
        The figure.
    """
    plot.make_figure()
    plot.get_axes()["mainplot_ax"].set(xlabel=xlabel, ylabel=CLUSTER_LABEL)
    set_text_sizes(plot.fig)
    return plot.fig


def get_dot_areas(fractions, dot_max):
    """Get dot areas as scanpy's dot plot sizes them: the fraction clipped to [0, dot_max], over dot_max, to the
    SIZE_EXPONENT, times DOT_AREA.

    Args:
        fractions: fractions of bins with expression.
        dot_max: the fraction drawn at DOT_AREA.

    Returns:
        Areas in points squared.
    """
    return (np.clip(fractions, 0.0, dot_max) / dot_max) ** SIZE_EXPONENT * DOT_AREA


def get_size_key_fractions(dot_max):
    """Get the size key's fractions as scanpy's: down from dot_max in steps of 0.2, 0.1 or 0.05 as dot_max exceeds
    0.6, 0.3 or not.

    Args:
        dot_max: the fraction drawn at DOT_AREA.

    Returns:
        The fractions, ascending, above 0.
    """
    step = 0.2 if dot_max > 0.6 else 0.1 if dot_max > 0.3 else 0.05
    fractions = np.arange(dot_max, 0.0, -step)[::-1]
    return fractions[fractions > step / 2]


def add_dot_keys(fig, ax, dot_max, left, width):
    """Put a dot plot's size key and colour bar side by side below its gene labels, across the axes' width.

    Args:
        fig: the figure.
        ax: the dot plot's axes.
        dot_max: the fraction drawn at DOT_AREA.
        left, width: inches of the axes' left edge and width.
    """
    top = get_bottom_inches(fig, ax) - KEY_GAP - KEY_TITLE
    size_ax = add_axes_inches(fig, left, top - KEY_HEIGHT, 0.45 * width, KEY_HEIGHT)
    fractions = get_size_key_fractions(dot_max)
    size_ax.scatter(
        range(len(fractions)),
        np.zeros(len(fractions)),
        s=get_dot_areas(fractions, dot_max),
        color="grey",
        edgecolors="black",
        linewidths=DOT_EDGE,
    )
    size_ax.set_xticks(range(len(fractions)), [f"{round(f * 100)}" for f in fractions])
    size_ax.set(xlim=(-0.5, len(fractions) - 0.5), yticks=[], title="Fraction of bins with expression (%)")
    size_ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    colour_ax = add_axes_inches(fig, left + 0.55 * width, top - 0.6 * KEY_HEIGHT, 0.45 * width, 0.6 * KEY_HEIGHT)
    fig.colorbar(ScalarMappable(Normalize(0.0, 1.0), CMAP), cax=colour_ax, orientation="horizontal", ticks=[0, 0.5, 1])
    colour_ax.set_title("Mean expression, scaled per gene")


def plot_dotplot(adata, groups, groupby, xlabel, rotation, path):
    """Draw a dot plot of genes per group of bins from scanpy's dot plot tables, each gene scaled to [0, 1], group
    names over brackets above the genes, keys below.

    Dots are sized and coloured as scanpy's dot plot draws them (DotPlot.dot_size_df, dot_color_df,
    https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.DotPlot.html), its largest dot DOT_AREA.

    Args:
        adata: log1p-normalized AnnData with obs groupby.
        groups: genes per bracket label.
        groupby: obs column of the rows.
        xlabel: the gene axis' label.
        rotation: degrees of the row labels and the bracket labels.
        path: output PDF.
    """
    dp = sc.pl.DotPlot(adata, groups, groupby=groupby, standard_scale="var")
    fractions, means = dp.dot_size_df, dp.dot_color_df
    dot_max = float(np.ceil(fractions.to_numpy().max() * 10) / 10)
    n_rows, n_genes = fractions.shape
    ends = np.cumsum([len(g) for g in groups.values()])
    x, y = np.meshgrid(np.arange(n_genes), np.arange(n_rows))
    left, bottom = AXES_ORIGIN
    height = n_rows * ROW_PITCH

    def draw(width):
        fig = plt.figure(figsize=(CANVAS, CANVAS))
        ax = add_axes_inches(fig, left, bottom, width, height)
        ax.scatter(
            x.ravel(),
            y.ravel(),
            s=get_dot_areas(fractions.to_numpy().ravel(), dot_max),
            c=means.to_numpy().ravel(),
            cmap=CMAP,
            norm=Normalize(0.0, 1.0),
            edgecolors="black",
            linewidths=DOT_EDGE,
        )
        ax.set(xlim=(-0.6, n_genes - 0.4), ylim=(n_rows - 0.4, -0.6), xlabel=xlabel, ylabel=CLUSTER_LABEL)
        ax.set_xticks(range(n_genes), fractions.columns, rotation=90)
        ax.set_yticks(
            range(n_rows), fractions.index, rotation=rotation, ha="right", va="center", rotation_mode="anchor"
        )
        top = add_axes_inches(fig, left, bottom + height + KEY_GAP / 2, width, BRACKET_HEIGHT)
        top.set(xlim=ax.get_xlim(), ylim=(0, 1))
        top.axis("off")
        for label, end, size in zip(groups, ends, [len(g) for g in groups.values()], strict=True):
            start = end - size
            top.plot([start - 0.4, start - 0.4, end - 0.6, end - 0.6], [0, 1, 1, 0], color="black", linewidth=1)
            top.text(
                (start + end - 1) / 2,
                1.5,
                label,
                rotation=rotation,
                ha="left" if rotation else "center",
                va="bottom",
                rotation_mode="anchor",
            )
        set_text_sizes(fig)
        add_dot_keys(fig, ax, dot_max, left, width)
        set_text_sizes(fig)
        return fig

    save_figure(draw_page_fit(draw), path)


def plot_score_matrix(adata, ranking, n_cell_types, path):
    """Draw each cluster's mean ULM score of every cluster's top-ranked cell types, rows by dendrogram.

    Args:
        adata: AnnData with obs cluster and obsm score_ulm (bins x cell types).
        ranking: decoupler rankby_group table (group, name, stat), in decoupler order.
        n_cell_types: top cell types per cluster with stat > 0.
        path: output PDF.
    """
    scores = ad.AnnData(X=adata.obsm["score_ulm"].to_numpy(), obs=adata.obs[["cluster"]])
    scores.var_names = adata.obsm["score_ulm"].columns
    sc.tl.dendrogram(scores, groupby="cluster")
    top = ranking[ranking["stat"] > 0].groupby("group", sort=False).head(n_cell_types)
    groups = {g: list(top.loc[top["group"] == g, "name"]) for g in scores.obs["cluster"].cat.categories}

    def draw(width):
        mp = sc.pl.matrixplot(
            scores,
            groups,
            groupby="cluster",
            dendrogram=True,
            standard_scale="var",
            cmap=CMAP,
            figsize=(width, 2.8),
            return_fig=True,
            show=False,
        )
        mp.legend(width=KEY_WIDTH, title="Mean ULM score,\nscaled per cell type")
        return finish_scanpy_plot(mp, "PanglaoDB cell type")

    save_figure(draw_page_fit(draw), path)


def main():
    """Read a cell-typing output folder and draw its figures."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--scalefactors", type=Path, required=True)
    ap.add_argument("--figures-dir", type=Path, default=None)
    a = ap.parse_args()
    microns_per_pixel = json.loads(a.scalefactors.read_text())["microns_per_pixel"]
    figures_dir = a.figures_dir or a.out_dir / "figures-refined"
    h5ads = sorted(a.out_dir.glob("*.h5ad"))
    if len(h5ads) != 1:
        sys.exit(f"expected one .h5ad in {a.out_dir}, found {len(h5ads)}")
    adata = ad.read_h5ad(h5ads[0])
    print(f"[load] {h5ads[0]}: {adata.n_obs} bins x {adata.n_vars} genes")
    annotation = pd.read_csv(a.out_dir / "cluster_annotation.tsv", sep="\t", dtype={"cluster": str})
    ranking = pd.read_csv(a.out_dir / "cluster_celltype_ranking.tsv", sep="\t", dtype={"group": str})
    features = pd.read_csv(a.out_dir / "cluster_top_features.tsv", sep="\t", dtype={"cluster": str})
    figures_dir.mkdir(parents=True, exist_ok=True)
    set_figure_style()

    cluster_colours = get_cluster_colours(adata)
    cell_type_colours = get_cell_type_colours(adata, cluster_colours)
    for key, colours, title in [
        ("cluster", cluster_colours, CLUSTER_LABEL),
        ("cell_type", cell_type_colours, "Cell type"),
    ]:
        plot_umap(adata, key, colours, title, figures_dir / f"umap_{key}.pdf")
        plot_spatial(adata, key, colours, title, microns_per_pixel, figures_dir / f"spatial_{key}.pdf")

    markers = annotation.assign(gene=annotation["support_genes"].str.split(",")).explode("gene").dropna(subset="gene")
    plot_dotplot(
        adata,
        get_marker_groups(markers, N_MARKERS),
        "cluster_label",
        "PanglaoDB marker gene",
        45,
        figures_dir / "dotplot__celltype_markers.pdf",
    )
    top = features[features["rank"] <= N_TOP_FEATURES]
    plot_dotplot(
        adata,
        {c: list(g) for c, g in top.groupby("cluster", sort=False)["gene"]},
        "cluster",
        "Space Ranger top feature",
        0,
        figures_dir / "dotplot__10x_top_features.pdf",
    )
    plot_score_matrix(adata, ranking, N_TOP_CELL_TYPES, figures_dir / "matrixplot__ulm_top_celltypes.pdf")
    print(f"[done] {figures_dir}")


if __name__ == "__main__":
    main()
