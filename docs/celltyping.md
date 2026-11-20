# Cell typing: PanglaoDB cell types of a Visium HD sample's Space Ranger clusters

`scripts/celltyping/annotate_celltypes.py` labels each Space Ranger graph-based cluster of one Visium HD bin size with a
PanglaoDB cell type, scored per bin with decoupler's univariate linear model and kept only where the cluster's own Space
Ranger differential-expression table supports the label's marker genes. The pipeline runs it on the ovarian sample's
16 um bins (rule `celltyping_ovarian`, in `CELLTYPING_ENV`, `configs/envs/celltyping.yaml`), with PanglaoDB's
markers through OmniPath (rule `celltyping_markers`); the analysis reads its `spot_annotation.tsv`
([reference.md](reference.md#method-inputs)).

## Method

1. **Markers**: PanglaoDB entries that are human, canonical, and have human sensitivity > 0.5, as in the decoupler
   scRNA tutorial.
2. **Scoring**: counts are library-size normalized and log1p transformed (`scanpy.pp.normalize_total`,
   `scanpy.pp.log1p`); each bin is scored per cell type with `decoupler.mt.ulm` (cell types with at least 3 markers
   present).
3. **Ranking**: `decoupler.tl.rankby_group` ranks cell types per cluster against the rest (`t-test_overestim_var`);
   among types with a positive statistic, the first in decoupler's order is the cluster's candidate.
4. **Support**: the candidate is kept where at least 3 of its markers have, in that cluster's Space Ranger DE table,
   log2 fold change > 1, adjusted p < 0.1 and mean counts >= 0.1 (the last two Space Ranger's web summary cutoffs);
   otherwise the cluster is `Unassigned`.

The cutoffs are constants at the top of the script (`MIN_SENSITIVITY`, `ULM_TMIN`, `SUPPORT_L2FC_GT`,
`SUPPORT_PADJ_LT`, `SUPPORT_MIN_MEAN_COUNTS`, `MIN_SUPPORT`, `N_DOTPLOT`). Clusters, UMAP and differential expression
are Space Ranger's own, reused without reclustering; bins Space Ranger did not cluster are dropped.

## Inputs

| File | Used for |
| --- | --- |
| `<binned_dir>/filtered_feature_bc_matrix.h5` | bin x gene UMI counts |
| `<binned_dir>/spatial/tissue_positions.parquet` | full-resolution pixel coordinates per bin |
| `<binned_dir>/spatial/tissue_hires_image.png`, `scalefactors_json.json` | H&E image and its scale for the spatial plot |
| `<binned_dir>/analysis/clustering/gene_expression_graphclust/clusters.csv` | cluster per bin |
| `<binned_dir>/analysis/umap/gene_expression_2_components/projection.csv` | UMAP per bin |
| `<binned_dir>/analysis/diffexp/gene_expression_graphclust/differential_expression.csv` | per-cluster mean counts, log2 fold change, adjusted p-value |

`<binned_dir>` is Space Ranger's `outs/binned_outputs/square_016um`.

## Outputs

Under `<OUTDIR>/ovarian/celltyping/<sample>/square_016um/`:

| File | Content |
| --- | --- |
| `spot_annotation.tsv` | `barcode`, `cluster`, `cell_type`, `cell_type_broad`; the analysis' `celltype` |
| `cluster_annotation.tsv` | one row per cluster: `cell_type` (`<broad> (<PanglaoDB type>)` or `Unassigned`), `cell_type_broad`, `candidate`, `rank_stat`, `n_support`, `marker_genes`, `n_spots` |
| `cluster_marker_genes.tsv` | one row per supporting marker: `cluster`, `cell_type`, `panglaodb_cell_type`, `gene`, `mean_counts`, `log2fc`, `padj`, `human_sensitivity` |
| `<sample>.annotated.h5ad` | log1p-normalized expression (counts in `layers["counts"]`); `obs` `cluster`, `cell_type`, `cell_type_broad`, `cluster_label`; `obsm` `X_umap`, `spatial`, `score_ulm`, `padj_ulm`; `uns["spatial"]` |
| `dotplot__celltype_markers.pdf`, `umap_{cluster,cell_type}.pdf`, `spatial_{cluster,cell_type}.pdf` | supporting markers per cluster; the Space Ranger UMAP and the bins in place, by cluster and by cell type |

The figures are drawn by `plot_celltypes.py`. It also redraws them from an existing folder alone,
`python plot_celltypes.py <out_dir> --scalefactors <binned_dir>/spatial/scalefactors_json.json`, into
`<out_dir>/figures-refined/`, adding Space Ranger's top features per cluster and each cluster's mean ULM score of its
top cell types.

## References

- PanglaoDB: Franzen et al. 2019, Database, doi:[10.1093/database/baz046](https://doi.org/10.1093/database/baz046)
- OmniPath: Turei et al. 2021, Mol Syst Biol, doi:[10.15252/msb.20209923](https://doi.org/10.15252/msb.20209923)
- decoupler: Badia-i-Mompel et al. 2022, Bioinform Adv,
  doi:[10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016); tutorial
  <https://decoupler.scverse.org/en/latest/notebooks/scell/rna_sc.html>
- scanpy: Wolf et al. 2018, Genome Biol, doi:[10.1186/s13059-017-1382-x](https://doi.org/10.1186/s13059-017-1382-x)
