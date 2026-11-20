#!/usr/bin/env python
"""Annotate Space Ranger clusters of a Visium HD sample with PanglaoDB cell types via decoupler.

Author (2026-10-28)

Dependencies:
    python 3.11, scanpy 1.11, decoupler 2.1, anndata 0.12, pandas, numpy, matplotlib, pyarrow

Usage:
    python annotate_celltypes.py <binned_dir> <markers_tsv> <out_dir> [--analysis-dir DIR] [--sample-id ID]
        binned_dir      # Space Ranger bin folder, e.g. binned_outputs/square_016um
        markers_tsv     # OmniPath PanglaoDB annotation table (download command below)
        out_dir         # output folder, created; must not already hold the output h5ad
        --analysis-dir  # Space Ranger secondary analysis folder, default <binned_dir>/analysis
        --sample-id     # prefix of the output h5ad and figure titles, default "sample"

    curl -o omnipath_panglaodb_annotations.tsv \
        "https://omnipathdb.org/annotations?resources=PanglaoDB&format=tsv"

Inputs:
    <binned_dir>/filtered_feature_bc_matrix.h5: filtered bin x gene UMI counts.
    <binned_dir>/spatial/tissue_positions.parquet: full-resolution pixel coordinates per bin.
    <binned_dir>/spatial/tissue_hires_image.png, scalefactors_json.json: H&E image and its scale.
    <analysis_dir>/clustering/gene_expression_graphclust/clusters.csv: graph-based cluster per bin.
    <analysis_dir>/umap/gene_expression_2_components/projection.csv: UMAP per bin.
    <analysis_dir>/diffexp/gene_expression_graphclust/differential_expression.csv: per-cluster
        mean counts, log2 fold change and adjusted p-value (cluster vs rest).
    markers_tsv: OmniPath long table (genesymbol, label, value, record_id) for PanglaoDB.
Parameters:
    MIN_SENSITIVITY: keep PanglaoDB human canonical markers with human_sensitivity above it.
    ULM_TMIN: minimum markers in the data for a cell type to be scored.
    SUPPORT_*: Space Ranger DE cutoffs a marker must pass to support its cluster's label.
    MIN_SUPPORT: minimum supporting markers to keep a label; otherwise "Low confidence (<type>)".
    BROAD: lineage of each of the 60 PanglaoDB cell types, used in the label text.
    N_DOTPLOT: supporting markers per cell type in the dot plot.
Outputs (out_dir/):
    cluster_marker_genes.tsv: one row per supporting marker gene: cluster, cell_type,
        panglaodb_cell_type, gene, mean_counts, log2fc, padj, human_sensitivity; sorted by
        cluster and log2fc.
    cluster_annotation.tsv: one row per cluster: cell_type, cell_type_broad, candidate
        (top-ranked PanglaoDB type), rank_stat, n_support, marker_genes, n_spots.
        cell_type is "<broad> (<PanglaoDB type>)", e.g. "Epithelial (Cholangiocytes)".
    <sample_id>.annotated.h5ad: log1p-normalized expression (counts in layers["counts"]),
        obs cluster, cell_type, cell_type_broad and cluster_label, obsm X_umap, spatial, score_ulm and padj_ulm,
        uns spatial (H&E image).
    spot_annotation.tsv: barcode, cluster, cell_type, cell_type_broad.
    dotplot__celltype_markers.pdf, umap_cluster.pdf, umap_cell_type.pdf, spatial_cluster.pdf,
        spatial_cell_type.pdf, drawn by plot_celltypes.py.
Notes/References:
    Clusters, UMAP and DE are Space Ranger outputs reused as is (no reclustering). Space Ranger
    clusters every filtered bin with >= 1 UMI; bins without a cluster are dropped.
    Markers: PanglaoDB (Franzen et al. 2019, doi:10.1093/database/baz046) via OmniPath
    (Turei et al. 2021, doi:10.15252/msb.20209923), filtered to human canonical markers with
    sensitivity > 0.5 as in the decoupler tutorial
    (https://decoupler.scverse.org/en/latest/notebooks/scell/rna_sc.html).
    Annotation follows that tutorial (decoupler, doi:10.1093/bioadv/vbac016): dc.mt.ulm per bin,
    dc.tl.rankby_group per cluster vs rest (t-test_overestim_var), keep stat > 0, candidate = first
    row in decoupler order.
    A candidate is kept only if >= MIN_SUPPORT of its markers have Space Ranger DE log2FC > 1,
    adjusted p < 0.1 and mean counts >= 0.1 in that cluster (the last two are the Space Ranger web
    summary cutoffs). A cluster failing this is labeled "Low confidence (<candidate>)" with
    cell_type_broad "Low confidence"; a cluster with no positive stat is "Unassigned".
    PanglaoDB has no lineage field, so BROAD assigns one lineage (Epithelial, Mesothelial,
    Stromal, Mural, Endothelial, Immune, Blood, Neural, Germ) to every PanglaoDB cell type,
    fixed before annotation.
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import decoupler as dc
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from plot_celltypes import (
    get_cell_type_colours,
    get_cluster_colours,
    get_marker_groups,
    plot_dotplot,
    plot_spatial,
    plot_umap,
    set_figure_style,
)

MIN_SENSITIVITY = 0.5
ULM_TMIN = 3
SUPPORT_L2FC_GT = 1.0
SUPPORT_PADJ_LT = 0.1
SUPPORT_MIN_MEAN_COUNTS = 0.1
MIN_SUPPORT = 3
N_DOTPLOT = 5
UNASSIGNED = "Unassigned"
LOW_CONFIDENCE = "Low confidence"
CLUSTERING = "gene_expression_graphclust"
BROAD = {
    "Acinar cells": "Epithelial",
    "Adipocytes": "Stromal",
    "Alpha cells": "Epithelial",
    "Astrocytes": "Neural",
    "B cells": "Immune",
    "B cells naive": "Immune",
    "Basal cells": "Epithelial",
    "Beta cells": "Epithelial",
    "Cholangiocytes": "Epithelial",
    "Clara cells": "Epithelial",
    "Delta cells": "Epithelial",
    "Dendritic cells": "Immune",
    "Ductal cells": "Epithelial",
    "Endothelial cells": "Endothelial",
    "Enteric glia cells": "Neural",
    "Enterocytes": "Epithelial",
    "Enteroendocrine cells": "Epithelial",
    "Ependymal cells": "Neural",
    "Epithelial cells": "Epithelial",
    "Erythroid-like and erythroid precursor cells": "Blood",
    "Fibroblasts": "Stromal",
    "Gamma (PP) cells": "Epithelial",
    "Gamma delta T cells": "Immune",
    "Germ cells": "Germ",
    "Goblet cells": "Epithelial",
    "Granulosa cells": "Stromal",
    "Hepatic stellate cells": "Stromal",
    "Hepatocytes": "Epithelial",
    "Keratinocytes": "Epithelial",
    "Leydig cells": "Stromal",
    "Loop of Henle cells": "Epithelial",
    "Luminal epithelial cells": "Epithelial",
    "Macrophages": "Immune",
    "Mammary epithelial cells": "Epithelial",
    "Mast cells": "Immune",
    "Mesangial cells": "Mural",
    "Mesothelial cells": "Mesothelial",
    "Monocytes": "Immune",
    "Myoepithelial cells": "Epithelial",
    "NK cells": "Immune",
    "Neurons": "Neural",
    "Neutrophils": "Immune",
    "Oligodendrocyte progenitor cells": "Neural",
    "Oligodendrocytes": "Neural",
    "Pancreatic stellate cells": "Stromal",
    "Peri-islet Schwann cells": "Neural",
    "Pericytes": "Mural",
    "Peritubular myoid cells": "Stromal",
    "Plasma cells": "Immune",
    "Plasmacytoid dendritic cells": "Immune",
    "Platelets": "Blood",
    "Podocytes": "Epithelial",
    "Proximal tubule cells": "Epithelial",
    "Pulmonary alveolar type II cells": "Epithelial",
    "Retinal ganglion cells": "Neural",
    "Serotonergic neurons": "Neural",
    "Smooth muscle cells": "Mural",
    "T cells": "Immune",
    "T cytotoxic cells": "Immune",
    "T memory cells": "Immune",
}


def load_markers(path):
    """Read the OmniPath PanglaoDB table and keep human canonical markers.

    Args:
        path: OmniPath annotations TSV for resource PanglaoDB.

    Returns:
        DataFrame with columns cell_type, gene, human_sensitivity; one row per (cell_type, gene).
    """
    raw = pd.read_csv(path, sep="\t", dtype=str)
    wide = raw.pivot_table(
        index=["genesymbol", "record_id"],
        columns="label",
        values="value",
        aggfunc="first",
    ).reset_index()
    wide["human_sensitivity"] = wide["human_sensitivity"].astype(float)
    keep = (
        (wide["human"] == "True") & (wide["canonical_marker"] == "True") & (wide["human_sensitivity"] > MIN_SENSITIVITY)
    )
    return (
        wide.loc[keep, ["cell_type", "genesymbol", "human_sensitivity"]]
        .rename(columns={"genesymbol": "gene"})
        .drop_duplicates(["cell_type", "gene"])
        .sort_values(["cell_type", "gene"])
        .reset_index(drop=True)
    )


def load_spaceranger(binned_dir, analysis_dir, sample_id):
    """Read Space Ranger counts, clusters, UMAP, H&E image and DE table.

    Args:
        binned_dir: Space Ranger bin folder (square_XXXum).
        analysis_dir: Space Ranger analysis folder.
        sample_id: Library key for adata.uns["spatial"].

    Returns:
        (adata, de): AnnData of clustered bins with obs cluster, obsm X_umap and spatial,
        uns spatial; DE table indexed by gene name.
    """
    adata = sc.read_10x_h5(binned_dir / "filtered_feature_bc_matrix.h5")
    adata.var_names_make_unique()
    clus = pd.read_csv(analysis_dir / "clustering" / CLUSTERING / "clusters.csv", index_col=0)["Cluster"]
    umap = pd.read_csv(
        analysis_dir / "umap" / "gene_expression_2_components" / "projection.csv",
        index_col=0,
    )
    pos = pd.read_parquet(binned_dir / "spatial" / "tissue_positions.parquet").set_index("barcode")
    if not clus.index.isin(adata.obs_names).all():
        sys.exit("clustered barcodes missing from the count matrix")
    adata = adata[adata.obs_names.isin(clus.index)].copy()
    ids = adata.obs_names
    k = clus.nunique()
    adata.obs["cluster"] = pd.Categorical(clus.loc[ids].astype(str), categories=[str(c) for c in range(1, k + 1)])
    adata.obsm["X_umap"] = umap.loc[ids].to_numpy()
    adata.obsm["spatial"] = pos.loc[ids, ["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
    sf = json.loads((binned_dir / "spatial" / "scalefactors_json.json").read_text())
    img = plt.imread(binned_dir / "spatial" / "tissue_hires_image.png")
    adata.uns["spatial"] = {sample_id: {"images": {"hires": img}, "scalefactors": sf}}
    de = pd.read_csv(analysis_dir / "diffexp" / CLUSTERING / "differential_expression.csv")
    de = de.drop_duplicates("Feature Name").set_index("Feature Name")
    print(f"[load] {adata.n_obs} clustered bins x {adata.n_vars} genes, {k} clusters")
    return adata, de


def annotate(adata, markers, de):
    """Score cell types per bin, rank them per cluster, and keep labels with DE support.

    Args:
        adata: log1p-normalized AnnData with obs cluster.
        markers: Output of load_markers.
        de: Space Ranger DE table indexed by gene name.

    Returns:
        (summary, support): per-cluster DataFrame (cell_type, cell_type_broad, candidate,
        rank_stat, n_support, marker_genes) and long DataFrame of supporting markers.
    """
    net = markers.rename(columns={"cell_type": "source", "gene": "target"})[["source", "target"]]
    dc.mt.ulm(data=adata, net=net, tmin=ULM_TMIN)
    score = dc.pp.get_obsm(adata, key="score_ulm")
    rank = dc.tl.rankby_group(adata=score, groupby="cluster", reference="rest", method="t-test_overestim_var")
    rank = rank[rank["stat"] > 0]
    sens = markers.set_index(["cell_type", "gene"])["human_sensitivity"]
    rows, support = [], []
    for c in adata.obs["cluster"].cat.categories:
        g = rank[rank["group"].astype(str) == c]
        if g.empty:
            rows.append(
                {
                    "cluster": c,
                    "cell_type": UNASSIGNED,
                    "cell_type_broad": UNASSIGNED,
                    "candidate": "",
                    "n_support": 0,
                }
            )
            continue
        cand = g.iloc[0]["name"]
        genes = [x for x in markers.loc[markers["cell_type"] == cand, "gene"] if x in adata.var_names]
        d = de.reindex([x for x in genes if x in de.index])
        d = pd.DataFrame(
            {
                "mean_counts": d[f"Cluster {c} Mean Counts"],
                "log2fc": d[f"Cluster {c} Log2 fold change"],
                "padj": d[f"Cluster {c} Adjusted p value"],
            }
        )
        ok = (
            (d["log2fc"] > SUPPORT_L2FC_GT)
            & (d["padj"] < SUPPORT_PADJ_LT)
            & (d["mean_counts"] >= SUPPORT_MIN_MEAN_COUNTS)
        )
        d = d[ok].sort_values("log2fc", ascending=False)
        passed = len(d) >= MIN_SUPPORT
        lineage = BROAD[cand] if passed else LOW_CONFIDENCE
        label = f"{lineage} ({cand})"
        rows.append(
            {
                "cluster": c,
                "cell_type": label,
                "cell_type_broad": lineage,
                "candidate": cand,
                "rank_stat": g.iloc[0]["stat"],
                "n_support": len(d),
                "marker_genes": ",".join(d.index),
            }
        )
        if passed:
            d = d.rename_axis("gene").reset_index()
            d.insert(0, "panglaodb_cell_type", cand)
            d.insert(0, "cell_type", label)
            d.insert(0, "cluster", c)
            d["human_sensitivity"] = [sens[(cand, x)] for x in d["gene"]]
            support.append(d)
    if not support:
        sys.exit("no cluster passed the marker support cutoffs; nothing to annotate")
    return pd.DataFrame(rows).set_index("cluster"), pd.concat(support, ignore_index=True)


def main():
    """Run the annotation and write tables, h5ad and figures."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("binned_dir", type=Path)
    ap.add_argument("markers_tsv", type=Path)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--analysis-dir", type=Path, default=None)
    ap.add_argument("--sample-id", default="sample")
    a = ap.parse_args()
    analysis_dir = a.analysis_dir or a.binned_dir / "analysis"
    h5ad = a.out_dir / f"{a.sample_id}.annotated.h5ad"
    if h5ad.exists():
        sys.exit(f"output exists, refusing to overwrite: {h5ad}")
    a.out_dir.mkdir(parents=True, exist_ok=True)

    adata, de = load_spaceranger(a.binned_dir, analysis_dir, a.sample_id)
    markers = load_markers(a.markers_tsv)
    print(f"[markers] {markers['cell_type'].nunique()} cell types, {len(markers)} marker rows")
    nobroad = sorted(set(markers["cell_type"]) - set(BROAD))
    if nobroad:
        sys.exit(f"{len(nobroad)} PanglaoDB cell types without a BROAD lineage, e.g. {nobroad[0]}")

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    summary, support = annotate(adata, markers, de)
    summary["n_spots"] = adata.obs["cluster"].value_counts().reindex(summary.index)
    print(summary.to_string())
    summary.to_csv(a.out_dir / "cluster_annotation.tsv", sep="\t")
    support.to_csv(a.out_dir / "cluster_marker_genes.tsv", sep="\t", index=False)

    for col in ["cell_type", "cell_type_broad"]:
        adata.obs[col] = adata.obs["cluster"].map(summary[col]).astype("category")
    adata.obs[["cluster", "cell_type", "cell_type_broad"]].to_csv(
        a.out_dir / "spot_annotation.tsv", sep="\t", index_label="barcode"
    )
    labels = [f"{c}: {t}" for c, t in summary["cell_type"].items()]
    adata.obs["cluster_label"] = pd.Categorical(
        adata.obs["cluster"].map(dict(zip(summary.index, labels))), categories=labels
    )

    set_figure_style()
    plot_dotplot(
        adata,
        get_marker_groups(support, N_DOTPLOT),
        "cluster_label",
        "PanglaoDB marker gene",
        45,
        a.out_dir / "dotplot__celltype_markers.pdf",
    )
    cluster_colours = get_cluster_colours(adata)
    cell_type_colours = get_cell_type_colours(adata, cluster_colours)
    microns_per_pixel = adata.uns["spatial"][a.sample_id]["scalefactors"]["microns_per_pixel"]
    for key, colours, title in [
        ("cluster", cluster_colours, "Space Ranger cluster"),
        ("cell_type", cell_type_colours, "Cell type"),
    ]:
        plot_umap(adata, key, colours, title, a.out_dir / f"umap_{key}.pdf")
        plot_spatial(adata, key, colours, title, microns_per_pixel, a.out_dir / f"spatial_{key}.pdf")

    adata.write_h5ad(h5ad, compression="gzip")
    print(f"[done] {a.out_dir}")


if __name__ == "__main__":
    main()
