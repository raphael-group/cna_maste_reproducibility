# annotate_celltypes: PanglaoDB cell-type labels for Space Ranger clusters of a Visium HD sample

`annotate_celltypes.py` labels each Space Ranger graph-based cluster of one Visium HD bin size with a PanglaoDB cell type, scored per bin with decoupler's univariate linear model and kept only when the cluster's own Space Ranger differential-expression table supports the label's marker genes.

## Files
```
annotate_celltypes.py   the annotation
plot_celltypes.py       the figures, also redrawn from an existing output folder
environment.yaml        pinned conda environment
README.md               this file
```

## 1. Install the environment
```sh
mamba env create -f environment.yaml      # or: conda env create -f environment.yaml
mamba activate celltyping
```
Pinned versions: Python 3.11.17, scanpy 1.11.5, decoupler-py 2.1.6, anndata 0.12.19, pandas 2.3.3, numpy 2.4.6, matplotlib 3.11.2, pyarrow 25.0.0, h5py 3.16.0 (conda-forge).

## 2. Download the marker table
PanglaoDB markers through OmniPath:
```sh
curl -o omnipath_panglaodb_annotations.tsv \
    "https://omnipathdb.org/annotations?resources=PanglaoDB&format=tsv"
```

## 3. Run
```sh
python annotate_celltypes.py <binned_dir> omnipath_panglaodb_annotations.tsv <out_dir> \
    [--analysis-dir <analysis_dir>] [--sample-id <sample_id>]
```

| argument | meaning |
|---|---|
| `binned_dir` | Space Ranger bin folder, e.g. `outs/binned_outputs/square_016um` |
| `markers_tsv` | OmniPath PanglaoDB table from step 2 |
| `out_dir` | output folder, created if missing |
| `--analysis-dir` | Space Ranger `analysis/` folder for that bin size; default `<binned_dir>/analysis` |
| `--sample-id` | prefix of the output `.h5ad` and figure titles; default `sample` |

Example (Space Ranger outputs of the 10x public Visium HD 3' human ovarian cancer sample, 16 um bins):
```sh
python annotate_celltypes.py \
    Visium_HD_3prime_Human_Ovarian_Cancer/outs/binned_outputs/square_016um \
    omnipath_panglaodb_annotations.tsv \
    results/square_016um \
    --sample-id H1-N3KXRGQ
```

> [!IMPORTANT]
> The script refuses to run when `<out_dir>/<sample_id>.annotated.h5ad` already exists. Use a new `out_dir` or remove the old file to rerun.

### Required inputs
| file | used for |
|---|---|
| `<binned_dir>/filtered_feature_bc_matrix.h5` | bin x gene UMI counts |
| `<binned_dir>/spatial/tissue_positions.parquet` | full-resolution pixel coordinates per bin |
| `<binned_dir>/spatial/tissue_hires_image.png`, `scalefactors_json.json` | H&E image and its scale for the spatial plot |
| `<analysis_dir>/clustering/gene_expression_graphclust/clusters.csv` | cluster per bin |
| `<analysis_dir>/umap/gene_expression_2_components/projection.csv` | UMAP per bin |
| `<analysis_dir>/diffexp/gene_expression_graphclust/differential_expression.csv` | per-cluster mean counts, log2 fold change, adjusted p-value |

> [!NOTE]
> Clusters, UMAP and differential expression are Space Ranger's own outputs, reused without reclustering. Bins that Space Ranger did not cluster are dropped.

## 4. Outputs (`out_dir/`)
| file | content |
|---|---|
| `cluster_annotation.tsv` | one row per cluster: `cell_type` (final label or `Unassigned`), `candidate` (top-ranked PanglaoDB type), `rank_stat`, `n_support` (supporting markers), `marker_genes`, `n_spots` |
| `cluster_marker_genes.tsv` | one row per supporting marker: `cluster`, `cell_type`, `gene`, `mean_counts`, `log2fc`, `padj`, `human_sensitivity`; sorted by cluster and log2fc |
| `<sample_id>.annotated.h5ad` | log1p-normalized expression (raw counts in `layers["counts"]`); `obs`: `cluster`, `cell_type`, `cluster_label`; `obsm`: `X_umap`, `spatial`, `score_ulm`, `padj_ulm`; `uns["spatial"]`: H&E image and scale factors |
| `dotplot__celltype_markers.pdf` | supporting markers (up to 5 per cell type) per cluster |
| `umap_cluster.pdf`, `umap_cell_type.pdf` | Space Ranger UMAP colored by cluster and by cell type |
| `spatial_cluster.pdf`, `spatial_cell_type.pdf` | bins at their physical position in um, as cnamaste's spatial panels, colored by cluster and by cell type; key below in two columns |

Figures are drawn by `plot_celltypes.py` to the lab figure rules (`../../docs/figure-rule.md`): at most one text column wide (6.48 in) and one full-page figure tall (7.18 in), Computer Modern text at 10 pt with 12.5 pt axis labels, saved as editable TrueType, bins rasterized at 600 dpi. A cell type takes the color of its lowest-numbered cluster. Dot plots put their keys below, the marker dot plot's row and cell-type labels at 45 degrees.

## 5. Redraw the figures
`plot_celltypes.py` redraws the figures of an existing cell-typing folder from its outputs alone, without re-running the annotation:
```sh
python plot_celltypes.py <out_dir> --scalefactors <spaceranger_dir>/spatial/scalefactors_json.json [--figures-dir <dir>]
```
It reads the folder's one `.h5ad`, `cluster_annotation.tsv` (`panglaodb_cell_type`, `support_genes`), `cluster_celltype_ranking.tsv` and `cluster_top_features.tsv`, and `microns_per_pixel` from the sample's Space Ranger scale factors (any bin size: it is the full-resolution image's pixel), and writes to `<out_dir>/figures-refined/` by default:

| file | content |
|---|---|
| `umap_cluster.pdf`, `umap_cell_type.pdf`, `spatial_cluster.pdf`, `spatial_cell_type.pdf` | as above |
| `dotplot__celltype_markers.pdf` | supporting markers (up to 5 per cell type) per cluster; a cluster without supporting markers has no group |
| `dotplot__10x_top_features.pdf` | Space Ranger top features (3 per cluster) per cluster |
| `matrixplot__ulm_top_celltypes.pdf` | mean ULM score of each cluster's top 3 cell types, scaled per cell type, rows ordered by dendrogram |

## Method
1. **Markers**: PanglaoDB entries that are human, canonical, and have human sensitivity > 0.5, as in the decoupler scRNA tutorial.
2. **Scoring**: counts are library-size normalized and log1p transformed (`scanpy.pp.normalize_total`, `scanpy.pp.log1p`); each bin is scored per cell type with `decoupler.mt.ulm` (cell types with at least 3 markers present).
3. **Ranking**: `decoupler.tl.rankby_group` ranks cell types per cluster against the rest (`t-test_overestim_var`); among types with a positive statistic, the first one in decoupler's order is the cluster's candidate.
4. **Support**: the candidate is kept when at least 3 of its markers have, in that cluster's Space Ranger DE table, log2 fold change > 1, adjusted p < 0.1 and mean counts >= 0.1 (the last two are Space Ranger's web summary cutoffs). Otherwise the cluster is `Unassigned`.

The cutoffs are constants at the top of the script (`MIN_SENSITIVITY`, `ULM_TMIN`, `SUPPORT_L2FC_GT`, `SUPPORT_PADJ_LT`, `SUPPORT_MIN_MEAN_COUNTS`, `MIN_SUPPORT`, `N_DOTPLOT`).

## References
- PanglaoDB: Franzen et al. 2019, Database, doi:[10.1093/database/baz046](https://doi.org/10.1093/database/baz046)
- OmniPath: Turei et al. 2021, Mol Syst Biol, doi:[10.15252/msb.20209923](https://doi.org/10.15252/msb.20209923)
- decoupler: Badia-i-Mompel et al. 2022, Bioinform Adv, doi:[10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016); tutorial <https://decoupler.scverse.org/en/latest/notebooks/scell/rna_sc.html>
- scanpy: Wolf et al. 2018, Genome Biol, doi:[10.1186/s13059-017-1382-x](https://doi.org/10.1186/s13059-017-1382-x)
