# cna_maste_reproducibility: benchmark of cnamaste against CalicoST, Numbat and Clonalscope

![Status](https://img.shields.io/badge/status-in%20development-yellow)

Scores cnamaste, CalicoST, Numbat and Clonalscope on the CalicoST simulations against their truth, compares their runs on a Visium HD 3prime ovarian sample (H1-N3KXRGQ) with one another and with the sample's cell types, and draws the paper's tables and figures from those scores.

## Overview

`scripts/run_benchmark.sh` runs `scripts/run_benchmark.py` on H1-N3KXRGQ. The script runs four steps in order:

- **Inputs, made beforehand.** Each method's runs: Numbat, Clonalscope and CalicoST by [baselines/](baselines/README.md),
  cnamaste by its own `run_cnamaste`. The sample's UGP bins, its Space Ranger outputs, and its per-spot cell types
  ([scripts/celltyping/](scripts/celltyping/README.md)). All are read from `input_data/` ([Configuration](#configuration)).
- **Step 1, simulations** (`write_simulation_metrics`). Each method's run of each CalicoST simulation is scored against
  the truth, on the run's own segments and at the genes: CNA sensitivity, normal specificity, N-CNA recall (exact
  allele copies), CNA recall (CNA state), clone ARI and spot coverage, with each run's runtime and memory.
- **Step 2, real data** (`write_real_data_metrics`). The sample's runs are compared on its in-tissue 16 um spots:
  - each run's spot coverage, and its false-positive tumor spots on normal cell types;
  - its share of the genome in each CNA state;
  - its clone ARI and CNA ARI against the reference run, cnamaste;
  - the clone ARI and CNP similarity of every two runs;
  - the copy-number error: how well each run's copies predict its clones' pseudobulk RDR and BAF.
- **Step 3, LaTeX** (`write_latex_table`, `write_metric_text`). A LaTeX table of the runs' coverage, false-positive
  rate, runtime and memory, plus each metric's definition as LaTeX text, each typeset to a PDF by `pdflatex`.
- **Step 4, figures** (`write_figures`). `combined.pdf` (main figures), `diagnosis.pdf` (diagnostic pages) and
  `appendix.pdf` (supplementary figures), with each page also written as its own PDF.

The code is in `scripts/metrics/` (the metrics), `scripts/plots/` (the figures) and `scripts/platforms/` (the Visium
HD bins). Reference docs:

- [docs/pipeline.md](docs/pipeline.md): each step, input reader and metric.
- [docs/metrics.md](docs/metrics.md): the simulation metrics and the CNP similarity.
- [docs/figures.md](docs/figures.md): every figure page.

## Environment

The benchmark runs in cnamaste's uv environment, with cnamaste at commit `42b96ff` (Python >= 3.12.2, its
[pyproject.toml](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/pyproject.toml)). It needs cnamaste's `luo`
dependency group, which adds cnplot for the profiles and tracks. `pdflatex` (TeX Live) typesets the LaTeX tables.

The setup below takes uv from conda-forge, clones both repositories, and builds cnamaste's environment in
`cna_maste/.venv`. uv downloads a matching Python if none is found. The last line checks that `pdflatex` and the LaTeX
packages the tables use are installed.

> [!NOTE]
> cnamaste's repository is private: cloning it needs read access and a git that holds your GitHub credentials. On
> macOS the system git (`/usr/bin/git`) reads them from the keychain; conda-forge's git has no credential helper set.

```sh
mamba create -y -n cna-benchmark -c conda-forge uv
conda activate cna-benchmark

git clone https://github.com/raphael-group/cna_maste_reproducibility.git
git clone https://github.com/michaelJwilson/cna_maste.git
git -C cna_maste checkout 42b96ff
uv sync --project cna_maste --locked --all-extras --group luo
export CNAMASTE_ROOT=$PWD/cna_maste

command -v pdflatex && kpsewhich geometry.sty amsmath.sty amssymb.sty booktabs.sty
```

Re-running the methods themselves needs their own conda environments, described in [baselines/README.md](baselines/README.md).
Re-making the cell types needs `scripts/celltyping/environment.yaml` ([scripts/celltyping/README.md](scripts/celltyping/README.md)).

## Configuration

A relative path in a yaml is resolved under this directory; an absolute one is read as given
(`metrics.read_method_input.get_benchmark_path`).

- `configs/benchmark.local.yaml` (`--benchmark_yaml`): the runs.
  - `vishd: <sample>:` holds, per run key, a `dir` (what the method's reader takes), `labels` and `profile`. It also
    holds `ugp` (the in-tissue spots, the Space Ranger directory, a `bb` directory per bin side), `cell_type.labels`,
    and `runtime` (the runs' runtime table).
  - `simulation:` holds `groundtruth_dir`, `runtime`, and per method a `dir` template of `{sample_id}`. Each extra
    cnamaste solver is a key `cnamaste-<solver>` with its `solver`.
  - A run key's method is the key before its first `_` or `-`, e.g. `calicost_016um_msu10` is CalicoST.
- `configs/benchmark.remote.yaml`: the cluster paths the local copies were taken from, in the same schema.
- `configs/config.local.yaml` (`--config_yaml`): the reference files, `reference.bed_file` (the region bed) and
  `reference.refdata.gene_list.gtf_file` (the genes).
- `data/hg38.regions.bed`, `data/hg38.chrom.sizes`: the chromosome arms and lengths, copied from
  Universal-Genotyping-Pipeline's `resources/data/`.
- `scripts/run_benchmark.sh`: the sample's paths for the figure step (`$B`, `$H`, `$CNAMASTE`).

> [!IMPORTANT]
> The metric step reads each run from `benchmark.local.yaml`, by its key's `dir`. The figure step reads each run from
> `run_benchmark.sh`'s options (`--numbat_dir`, `--calicost_dir`, `--cnamaste_labels`, ...). To swap a run, change
> both places.

`input_data/` is gitignored and holds the local copies the yamls name:

```text
input_data/
  reference/genes.gtf.gz              the genes, Cell Ranger's refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz
  real_data/
    output_10092026/<sample>/         UGP bins (bb/, bb-032um/, ...) and Space Ranger outputs
    celltyping/<sample>/              spot_annotation.tsv per bin side
    numbat/ clonalscope/ calicost/ cnamaste/    each method's runs
    runtime.<sample>.tsv              the runs' runtime and memory
  simulation_data/
    groundtruth/                      per simulation, truth_clone_labels.tsv and truth_acn_profile.tsv
    simulation.runtime.tsv            each run's runtime and memory
    numbat/ clonalscope/ method_results/calicost/ nomixing_cnamaste_related/
```

The genes come from 10x Genomics' GRCh38-2024-A reference, an 11 GB download; only the GTF is kept:

```sh
curl -LO https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
mkdir -p input_data/reference
tar -xzOf refdata-gex-GRCh38-2024-A.tar.gz refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz > input_data/reference/genes.gtf.gz
```

### Inputs

The command takes Numbat's, Clonalscope's and CalicoST's run directories and finds these files in them:

| Method | Labels | Profile |
| --- | --- | --- |
| Numbat | `clone_post_<i>.tsv`, the last iteration, 32 um bins given to the 16 um bins inside each | `bulk_clones_final.tsv.gz` |
| Clonalscope | `<sample>.<assay>.cells.tsv` of a Visium tutorial export's `postprocess/`, coarser bins given to the 16 um bins inside each | `<sample>.<assay>.profile.tsv` |
| CalicoST | `clone_labels.tsv` of the best rectangle by `total_llf` | `cnv_seglevel.tsv`, unconstrained ploidy |
| cnamaste | `clone_labels_int.tsv.gz` | `cnv_seglevel.tsv` |
| cell type | `spot_annotation.tsv` | |

- The H&E image and the bins' places come from a Space Ranger directory of the sample, at any bin side.
- `scripts/metrics/read_method_input.py` reads each method's run, one reader per method. Every method's normal clone is
  clone 0.
- `read_clonalscope` reads a Clonalscope Visium tutorial export, its normal clusters first.
- A CalicoST run whose unconstrained table leaves out a labelled clone is read at a fallback ploidy
  ([docs/pipeline.md](docs/pipeline.md#calicost-clones-missing-from-the-unconstrained-table)).

## Run

From any directory, with `CNAMASTE_ROOT` set as in [Environment](#environment):

```sh
export CNAMASTE_ROOT=/path/to/cna_maste
bash cna_maste_reproducibility/scripts/run_benchmark.sh
```

Any option given to `run_benchmark.sh` is appended to `run_benchmark.py`'s command line ([Options](#options)):

```sh
bash cna_maste_reproducibility/scripts/run_benchmark.sh --cnamaste_normal_fraction <clone>=<fraction>
```

`run_benchmark.sh` does the following:

1. Stops if `CNAMASTE_ROOT` is unset.
2. Puts `$CNAMASTE_ROOT/python` and this directory on `PYTHONPATH`, and logs cnamaste's commit
   (`git describe --always --dirty`).
3. Deletes the previous run's `output/benchmark_results/`, then makes it again with `tsvs/gene_level_cnp/`, `plots/`,
   `pdf-pages/` and `text/`.
4. Runs `run_benchmark.py` on H1-N3KXRGQ in cnamaste's environment, `uv run --project "$CNAMASTE_ROOT" --group luo`.

`run_benchmark.py` stops if `cnamaste` is imported from anywhere other than `$CNAMASTE_ROOT/python/cnamaste`
(`get_cnamaste_root`). To use a cnamaste run other than `cnamaste`, pass `--run <key> --reference <key>` with that
run's own `--cnamaste_labels` and `--cnamaste_seglevel`.

> [!NOTE]
> `run_benchmark.sh` passes no sweep run (no `--sweep_run`, `--calicost_sweep` or `--false_positive_sweep`). So its
> diagnosis.pdf has no sweep pages, its pages 1 and 2 hold Numbat and cnamaste, and it writes no `cn_error.sweep`. The
> sweep's runs stay in `configs/benchmark.local.yaml`: `calicost_016um` (`H1-N3KXRGQ_016um_msu20`),
> `calicost_016um_msu10`, `_msu15`, `_msu25`, `_msu50` and `calicost_032um`.

### Options

| Option | Default | Meaning |
| --- | --- | --- |
| `--sample` | required | The sample, a key under the benchmark yaml's platform; in the figures' titles and tables. |
| `--outdir` | required | The run's directory, `output/benchmark_results/` in `run_benchmark.sh`, gitignored, its four subdirectories made beforehand. |
| `--spaceranger_dir` | required | A Space Ranger directory of the sample, any bin side: `spatial/tissue_positions.parquet`, `scalefactors_json.json` and the H&E image. |
| `--spots` | required | Every in-tissue spot at `--side_um`, a UGP `barcodes.tsv.gz`. |
| `--celltype` | required | Per-spot cell types at `--side_um`, `spot_annotation.tsv`. |
| `--bb_dir` | none | A UGP `bb` directory and its bins' side, `<side_um>=<dir>`, once per side: the appendix's UMI page, and the diagnosis' UMI page at `--side_um`. |
| `--cnamaste_labels`, `--cnamaste_seglevel` | required | The cnamaste run's `clone_labels_int.tsv[.gz]` and `cnv_<ploidy>_seglevel.tsv`. With `binned_data.npz` and `rdrbaf_final_nstates<n>_smp.npz` beside them, combined.pdf holds cnamaste's own page. |
| `--calicost_dir` | required | A CalicoST sample directory, `[<label>=]<dir>`, read at its best rectangle; its label in its title. |
| `--numbat_dir` | required | A Numbat output directory. |
| `--clonalscope_dir` | none | A Clonalscope export's `postprocess/`; without it, Clonalscope is a placeholder. |
| `--calicost_sweep` | none | A CalicoST run of a parameter sweep, `<label>=<dir>`, once per run: its clones and profile in diagnosis.pdf. |
| `--false_positive_sweep` | none | The label of a `--calicost_sweep` run, once per run: its clones on diagnosis.pdf's false-positive page. |
| `--benchmark_yaml` | `configs/benchmark.local.yaml` | The runs, per sample and for the simulations ([Configuration](#configuration)). |
| `--config_yaml` | `configs/config.local.yaml` | The region bed and the genes. |
| `--method` | every method of `simulation` | A simulation method scored, once per method. |
| `--simulation_dataset` | `calicost-sim` | The simulations' wide table's name, `tsvs/<name>.tsv.gz`. |
| `--run` | `numbat`, `clonalscope`, `calicost`, `cnamaste` | A run key of the sample, once per run, its method the key before its first `_` or `-`. |
| `--reference` | `cnamaste` | The run every other run is compared with. |
| `--sweep_run` | none | A run key scored with the reference on the copy-number error, on the bins they share. |
| `--msr` | 300 | UGP's `min_snp_reads` of the bins scored, their directory `MSR<msr>`. |
| `--real_data_dataset` | `ovarian` | The sample's wide table's name, `tsvs/<name>.tsv.gz`. |
| `--table` | `table` | The LaTeX table's name under `text/`. |
| `--metric_text` | `metric` | The metric text's name under `text/`. |
| `--cnamaste_normal_fraction` | the run's log | A tumour clone of the cnamaste run's fit and its decoder's normal fraction, `<clone>=<fraction>`, once per tumour clone. By default each integer clone's from the run's `cnamaste_<id>.log`, run_cnamaste's `normal fractions per clone` (`plots.cnamaste_page.get_normal_fractions`); the normal clone is held at 1. Used where the run wrote `binned_data.npz`: the cnamaste page's expected lines. |
| `--combined` | `combined.pdf` | The main figures' name under `plots/`, `.pdf` or `.png`; a `.png`'s page k after the first is `<stem>_page<k>.png`. |
| `--diagnosis` | `diagnosis.pdf` | The diagnostic pages' name under `plots/`, `.pdf` or `.png` as `--combined`; its UMI table is `tsvs/<stem>.tsv.gz`. |
| `--appendix` | `appendix.pdf` | The supplementary figures' name under `plots/`. |
| `--side_um` | 16 | The side of the spots compared, drawn and scored; coarser bins are given to each such spot inside them. |
| `--calicost_ploidy` | `unconstrained` | CalicoST's table: `cnv<_ploidy>_seglevel.tsv`. |
| `--calicost_fallback_ploidy` | `diploid tetraploid` | The tables tried, in order, where `--calicost_ploidy`'s leaves out a labelled tumour clone: the first holding them all is read, else `--calicost_ploidy`'s; `none` for none. |
| `--region_bed` | `data/hg38.regions.bed` | The regions the diagnosis' sweep profiles' genome axis draws, a chromosome's arms. |
| `--chrom_sizes` | `data/hg38.chrom.sizes` | The chromosome lengths of the genome axis. |
| `--assay` | `VISIUM3prime` | UGP's assay type, the bins' subdirectory and the stem of a `--bb_dir`'s AnnData and allele matrices' directory. |
| `--bins` | 100 | Each appendix UMI histogram's bin count. |

## Output

`<outdir>` is `--outdir`, `output/benchmark_results/` for `run_benchmark.sh` (gitignored). Names below are at their
defaults: `calicost-sim` (`--simulation_dataset`), `ovarian` (`--real_data_dataset`), `<sample>` (`--sample`), `table`
and `metric` (`--table`, `--metric_text`), and the three PDFs. Tables are tab-separated and gzipped, with the gzip
header's time set to 0, so equal tables give equal bytes (`scripts/metrics/outputs.py`).

| File | Description |
| --- | --- |
| `<outdir>/tsvs/calicost-sim.tsv.gz` | Simulations, a row per simulation and method: `#spots`, `#called`, `spot-coverage-rate`, each metric per level, `<metric> (segment)` and `<metric> (gene)`, then `Runtime [mins]` and `Memory [GB]`. |
| `<outdir>/tsvs/calicost-sim.median.tsv.gz` | A row per simulated class and method: `sim_id`, `method`, `#replicates`, then each metric's median over the class' replicates, NaN skipped. |
| `<outdir>/tsvs/gene_level_cnp/<sample_id>.gene_level_cnp.tsv.gz` | Per simulation, a row per gene: the truth's and each run's fields per clone. |
| `<outdir>/tsvs/ovarian.tsv.gz` | The sample's wide table, a row per run: spot coverage and false-positive counts, genome composition, clone ARI and CNA ARI against the reference, copy-number error, runtime and memory. |
| `<outdir>/tsvs/ovarian.json` | The wide table's summary: per run, its summary metrics and their provenance counts. |
| `<outdir>/tsvs/<sample>.gene_level_cnp.tsv.gz` | A row per gene: each run's fields per clone. |
| `<outdir>/tsvs/<sample>.segment_level_cnp.tsv.gz` | A row per joint segment of the runs giving copies (CalicoST, cnamaste, Clonalscope): each run's fields per clone. |
| `<outdir>/tsvs/<sample>.genome_composition.tsv.gz` | Per run and clone, and its `all` row, the share of the genome in each CNA state and `uncalled`. |
| `<outdir>/tsvs/<sample>.cna_confusion.tsv.gz` | Per run, the cells by `(reference state, run state)` at the run's `level`. |
| `<outdir>/tsvs/<sample>.clone_ari.tsv.gz` | Every ordered pair of runs: `run_a`, `run_b`, `#shared-spots`, `clone_ari`. |
| `<outdir>/tsvs/<sample>.cnp_similarity.tsv.gz` | Every ordered pair of runs: `run_a`, `run_b`, `clone_a`, `clone_b`, `#altered-genes`, `cnp_similarity`. |
| `<outdir>/tsvs/<sample>.cn_error.sweep.tsv.gz` | With `--sweep_run` only: the copy-number error of the reference and sweep runs, scored on the bins they share. |
| `<outdir>/tsvs/diagnosis.tsv.gz` | The numbers behind diagnosis.pdf's UMI page, a row per run and UMI measure. |
| `<outdir>/text/table.tex`, `table.pdf`, `table.log` | The runs' spot coverage rate, false-positive rate, runtime and memory as a booktabs `table` float, typeset on letter pages, and pdflatex's log. |
| `<outdir>/text/metric.tex`, `metric.pdf`, `metric.log` | Each metric's definition as computed, a paragraph per metric with the pages that draw it, typeset, and pdflatex's log. |
| `<outdir>/plots/combined.pdf` | Main figures: the simulations' metrics at the genes; cnamaste's own page; each method's clone profiles beside every two runs' CNP similarity and clone ARI and each run's skill scores; last, the panel names. |
| `<outdir>/plots/diagnosis.pdf` | Diagnostic pages: each run's false-positive tumor spots and their UMI; each run's copy-number error; with `--calicost_sweep`, the sweep runs; the simulation pages for every class; last, the panel names. |
| `<outdir>/plots/appendix.pdf` | Supplementary figures: per-spot UMI, H&E with cell types and each method's clones, cells per altered state, state Jaccard index, the simulations' runtime and memory, cnamaste's solvers; last, the panel names. |
| `<outdir>/pdf-pages/<stem>-P<k>.pdf` | Page k of each of the three PDFs alone, cropped to its content, unless the page is a paper figure, named below. |
| `<outdir>/pdf-pages/main-fig-zenodo-sim-result.pdf` | combined.pdf's simulation metrics page. |
| `<outdir>/pdf-pages/main-fig-cnamaste-result.pdf` | combined.pdf's cnamaste page. |
| `<outdir>/pdf-pages/main-fig-cnp-compare.pdf` | combined.pdf's clone profiles and run comparison page. |
| `<outdir>/pdf-pages/suppl-fig-umi.pdf` | appendix.pdf's per-spot UMI page. |
| `<outdir>/pdf-pages/suppl-fig-ovarian-label.pdf` | appendix.pdf's spatial labels page. |
| `<outdir>/pdf-pages/suppl-fig-cna-composition.pdf` | appendix.pdf's cells per altered state page. |
| `<outdir>/pdf-pages/suppl-fig-cna-jaccard.pdf` | appendix.pdf's state Jaccard index page. |
| `<outdir>/pdf-pages/suppl-fig-zenodo-sim-runtime.pdf` | appendix.pdf's simulations' runtime and memory page. |
| `<outdir>/pdf-pages/suppl-fig-zenodo-sim-ablation.pdf` | appendix.pdf's cnamaste solvers page. |

Every table's columns are defined in [docs/pipeline.md](docs/pipeline.md), and every page's panels in
[docs/figures.md](docs/figures.md).

## Test

`tests/` holds the benchmark's tests. They need no `input_data/`. `pytest.ini` puts this directory on the import path.
Run from this directory, in cnamaste's environment:

```sh
cd cna_maste_reproducibility
uv run --project "$CNAMASTE_ROOT" --group luo pytest
```
