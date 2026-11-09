# Reference

`run_cna_maste_reproducibility.smk`, one Snakemake workflow, every step a rule of it. The ovarian sample's cell types are
in [celltyping.md](celltyping.md).

## Table of Contents
- [Environments](#environments)
- [Configuration](#configuration)
  - [Paths and resources](#paths-and-resources)
  - [Fixed settings](#fixed-settings)
  - [Repositories](#repositories)
  - [Method parameters](#method-parameters)
  - [Analysis](#analysis)
  - [Analysis options](#analysis-options)
- [Steps](#steps)
  - [Baselines](#baselines)
  - [cnamaste](#cnamaste)
  - [Analysis step](#analysis-step)
- [Metrics](#metrics)
  - [Method inputs](#method-inputs)
  - [CNA states and fields](#cna-states-and-fields)
  - [Genes and joint segments](#genes-and-joint-segments)
  - [Against the truth (simulations)](#against-the-truth-simulations)
  - [Between runs (real data)](#between-runs-real-data)
  - [CNP similarity](#cnp-similarity)
  - [Genome composition](#genome-composition)
  - [Copy-number error](#copy-number-error)
  - [CalicoST clones missing from the unconstrained table](#calicost-clones-missing-from-the-unconstrained-table)
- [Outputs](#outputs)
  - [Run layout](#run-layout)
  - [Tables](#tables)
  - [LaTeX tables](#latex-tables)
  - [Figures](#figures)
- [Code](#code)
- [Test](#test)

---

## Environments

`CNA_MASTE_REPRODUCIBILITY_ENV` runs the workflow; the setup step builds each other one that is missing from its
`ENV_RECIPES` recipe: a conda yaml, or `uv:<repository>[:<group>]`, that repository's `uv.lock`
(`uv sync --locked --all-extras` into the prefix).

| Environment | Recipe | Runs |
|---|---|---|
| `CNA_MASTE_REPRODUCIBILITY_ENV` | `configs/envs/cna_maste_reproducibility.yaml` | Snakemake 9, pyyaml, pandas, pysam, uv, git: the workflow, its config writers, `convert_bam_visium3prime.py`, the runtime tables. |
| `CALICOST_ENV` | `configs/envs/calicost.yaml` | CalicoST, installed editable from `CALICOST_REPO`; the simulations' 10x matrices. |
| `NUMBAT_ENV` | `configs/envs/numbat.yaml` | Numbat, bioconda `r-numbat` 1.5.2 (CRAN). |
| `CLONALSCOPE_ENV` | `configs/envs/clonalscope.yaml` | Clonalscope, installed from `CLONALSCOPE_REPO`. |
| `CNAMASTE_ENV` | `uv:CNAMASTE_REPO` | cnamaste at `cleanup3`: `run_cnamaste`, `run_cnamaste_plot`, `run_sim_qa`. |
| `CNAMASTE_ANALYSIS_ENV` | `uv:CNAMASTE_ANALYSIS_REPO:luo` | The analysis step and the tests: cnamaste at `42b96ff`, its `luo` group adding cnplot. |
| `UNIVERSAL_GENOTYPING_PIPELINE_ENV` | `{UGP_REPO}/environment.yaml` | Universal-Genotyping-Pipeline (UGP) on the ovarian sample. |
| `CELLTYPING_ENV` | `configs/envs/celltyping.yaml` | The cell types. |

> [!IMPORTANT]
> - The analysis imports cnamaste at `42b96ff` and the runs use `cleanup3`: `scripts/metrics/` and `scripts/plots/`
>   import modules that `cleanup3` no longer has.
> - cnamaste's repository is private: cloning it needs read access and a git holding your GitHub credentials. On macOS
>   the system git (`/usr/bin/git`) reads them from the keychain; conda-forge's git has no credential helper set.
> - The method runs need Linux x86_64: bioconda builds `eagle2` (UGP's phaser) and `r-numbat` for `linux-64` and
>   `osx-64` alone, and CalicoST's runs need GNU time and systemd user scopes with the memory controller delegated,
>   for their `MemoryMax`. The LaTeX tables need `pdflatex` (TeX Live). The README's first commands check each.

---

## Configuration

`configs/cna_maste_reproducibility.yaml` is the one config. Its "Set" sections are what a run sets: the `/path/to/`
values, `CORES` and `MEM_GB`, in a copy of those sections passed with `--configfile`, whose keys override the
default ([README](../README.md#configuration--run)). Its "Fixed" sections pin the benchmark: changing them changes
what is run or scored. A relative path is under the repository, and
`--config KEY=VALUE` overrides one key.

### Paths and resources

| Key | Description |
|---|---|
| `CORES`, `MEM_GB` | The run's CPUs and memory; pass the same to `snakemake --cores <CORES> --resources mem_mb=<MEM_GB * 1024>`. |
| `CONDA_ROOT` | Its `etc/profile.d/conda.sh` activates the method environments; its mamba builds them. |
| `CALICOST_ZENODO_SIM_ROOT` | `CalicoST_simulation_deposit.zip` of [doi:10.5281/zenodo.14175627](https://doi.org/10.5281/zenodo.14175627), unzipped: `simulated_data_related/`, `nomixing_calicost_related/`, `numbat_related/`. |
| `OVARIAN_SPACERANGER_ROOT` | The Visium HD 3' ovarian sample's Space Ranger `outs/`, with `possorted_genome_bam.bam` (`spaceranger count --create-bam=true`) and `binned_outputs/square_016um/`. |
| `PATH_TO_REFDATA` | `refdata-gex-GRCh38-2024-A`, 10x's [Space Ranger downloads](https://www.10xgenomics.com/support/software/space-ranger/downloads). |
| `PATH_TO_SNP_PANEL`, `PATH_TO_PHASE_PANEL` | UGP's SNP panel (bgzipped, tabix-indexed VCF) and phasing panel, UGP `resources/README.md` ("SNP Panels", "Population-haplotype Panels"). |
| `PATH_TO_EAGLE_MAP` | Eagle2's hg38 genetic map, `Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz` ([Eagle](https://alkesgroup.broadinstitute.org/Eagle/)). |
| `PATH_TO_ORDERED_GENES` | Cell Ranger GRCh38-2020-A `genes/ordered_genes.tsv`: the simulations' gene symbols to Ensembl IDs. |
| `FORMAT_BB_SCRIPT`, `FORMAT_BB_PYTHON` | The script writing UGP's bins in CalicoST's input layout, and its python; null: the workflow's. |
| `OUTDIR` | Every output ([Run layout](#run-layout)). |
| `TMPDIR` | Temporary files; null: `<OUTDIR>/tmp`. |
| `<METHOD>_REPO` | Each clone's root. |
| `<METHOD>_ENV` | Each environment's prefix. |

### Fixed settings

| Key | Description |
|---|---|
| `OVARIAN_SAMPLE_ID` | The ovarian sample's id, `H1-N3KXRGQ`; `scripts/metrics/runtime.py` names its runs by it. |
| `REPOS` | Each clone's `url` and `rev` ([Repositories](#repositories)). |
| `ENV_RECIPES` | Each environment's recipe ([Environments](#environments)). |
| `CALICOST_CONFIG`, `NUMBAT_CONFIG`, `CLONALSCOPE_CONFIG` | The baselines' parameters per case, `configs/baselines/` ([Method parameters](#method-parameters)). |
| `CNAMASTE` | cnamaste's parameters per case ([Method parameters](#method-parameters)). |
| `ANALYSIS` | The CalicoST runs the analysis scores ([Analysis](#analysis)). |
| `UGP_SAMPLE_TEMPLATE`, `UGP_CONFIG_TEMPLATE` | UGP's sample sheet and config, `configs/templates/`, filled with the paths above. |
| `CELLTYPING_MARKERS_URL` | PanglaoDB's markers through OmniPath ([celltyping.md](celltyping.md)). |

### Repositories

A missing clone is cloned at its `rev`; an existing one is checked against it and left as it is.

| Key | Repository | Revision |
|---|---|---|
| `CALICOST_REPO` | [raphael-group/CalicoST](https://github.com/raphael-group/CalicoST) | `c1abcae3e3657e01e547ee4529e3b9d039221453` (main) |
| `NUMBAT_REPO` | [kharchenkolab/numbat](https://github.com/kharchenkolab/numbat) | `54d8268f3d14bd6afaa6784eeb839db7479a8947`, the commit setting `Version: 1.5.2`; GitHub has no `v1.5.2` tag |
| `CLONALSCOPE_REPO` | [seasoncloud/Clonalscope](https://github.com/seasoncloud/Clonalscope) | `07dd0c576b6754d18dda86ec8905d8a9be5e77b0` (main; package version 1.1.0) |
| `CNAMASTE_REPO` | michaelJwilson/cna_maste (private) | `cleanup3` |
| `CNAMASTE_ANALYSIS_REPO` | michaelJwilson/cna_maste (private) | `42b96ff` |
| `UGP_REPO` | [raphael-group/Universal-Genotyping-Pipeline](https://github.com/raphael-group/Universal-Genotyping-Pipeline) | `c5c0b9b9ccb84268f4d29f9840dd91b71d793d7d` |

### Method parameters

A key left out keeps the method's default; each config's header cites the source of every non-default value.

| Config | Case | Keys |
|---|---|---|
| `configs/baselines/numbat.yaml` | `calicost-sim`, `ovarian` | `run_label`, `ncores`, `max_entropy` (0.8, the [spatial vignette](https://kharchenkolab.github.io/numbat/articles/spatial-rna.html)'s; default 0.5); ovarian `bin_um: 32`. |
| `configs/baselines/clonalscope.yaml` | `calicost-sim`, `ovarian` | `run_label`, `ploidy`; ovarian `bins_um`, the bin sides run. The method's arguments are fixed in `scripts/baselines/clonalscope/tutorial.R`. |
| `configs/baselines/calicost.yaml` | `calicost-sim` | `sims` (`all`, every `numcnas*` folder of the deposit, or a list), `rects` (the initializations run at once per simulation, one core each). |
| | `ovarian` | `mem_max` (cgroup `MemoryMax` per run; empty: `MEM_GB`), `inputs` (bin side), `runs` (tag: its input and the `configuration_cna` keys it sets; `msu<N>` sets `min_snpumi_perspot = N`), `configuration_cna` (CalicoST's default). |
| `CNAMASTE` | `zenodo_sim` | `solvers` (`<label_solver>-<initializer>-<decoder>`), `random_seeds`, `qa_solver`, `calicost_ploidy`, `threads`, `config` (cnamaste's `zenodo_sim_config.yaml` at `cleanup3`, less its paths and references). |
| | `ovarian` | `threads`, `config` (`hmrf.lattice_type: square`, which has no default). |

### Analysis

The analysis runs on this pipeline's runs. Rule `benchmark_yaml` writes `<OUTDIR>/configs/benchmark.yaml`, the runs,
and `config.yaml`, the references; rule `analysis_options` writes `<OUTDIR>/configs/analysis_options.yaml`, the options
each step reads. `ANALYSIS` picks the CalicoST runs:

| Key | Description |
|---|---|
| `calicost_run` | The ovarian CalicoST run, a tag of `CALICOST_CONFIG`'s `runs`, the analysis scores and draws as CalicoST's; its tag's last word, e.g. `msu20`, names it in the figures. |
| `calicost_sim` | The simulations' CalicoST runs scored: `deposit`, the deposit's `nomixing_calicost_related/`; `rerun`, this pipeline's. |

The benchmark yaml names the runs:

- `vishd: <sample>:` per run key a `dir` (what the method's reader takes); `ugp` (`spots`, the in-tissue spots;
  `spaceranger_dir`; `bb_dir`, a directory per bin side); `cell_type.labels`; and `runtime`, the runs' runtime table.
  The ovarian CalicoST runs are `calicost`, `calicost_run`, and `calicost_<tag less the sample>`. The figure step's
  options are read from it.
- `simulation:` `groundtruth_dir`, `runtime`, and per method a `dir` template of `{sample_id}`. Each cnamaste solver
  but `cnamaste`'s is a key `cnamaste-<solver>` with its `solver`.
- A run key's method is the key before its first `_` or `-` (`read_method_input.get_run_method`):
  `calicost_016um_msu10` is CalicoST, `cnamaste-opengm_icm-lattice-lattice` cnamaste.

The config yaml (`config_yaml`) names `reference.bed_file`, the region bed, `data/hg38.regions.bed`, and
`reference.refdata.gene_list.gtf_file`, the genes, `PATH_TO_REFDATA`'s `genes/genes.gtf.gz`. `data/hg38.regions.bed` and
`data/hg38.chrom.sizes`, the chromosome arms and lengths, are copied from UGP's `resources/data/`.

### Analysis options

The fields of `scripts/analysis/options.py`'s `AnalysisOptions`, a list option a yaml list. The pipeline sets `sample`,
`outdir`, `benchmark_yaml`, `config_yaml` and the figure step's inputs from the benchmark yaml; the others keep their
defaults.

| Option | Default | Description |
|---|---|---|
| `sample` | required | The sample, a key under the benchmark yaml's `vishd`. |
| `outdir` | required | `<OUTDIR>/benchmark_results/` in the pipeline. |
| `spaceranger_dir` | required | A Space Ranger directory of the sample, any bin side: `spatial/tissue_positions.parquet`, `scalefactors_json.json` and the H&E image. |
| `spots` | required | Every in-tissue spot at `side_um`, a UGP `barcodes.tsv.gz`. |
| `celltype` | required | Per-spot cell types at `side_um`, `spot_annotation.tsv`. |
| `bb_dir` | none | A UGP `bb` directory and its bins' side, `<side_um>=<dir>`, once per side: the UMI pages. |
| `cnamaste_labels`, `cnamaste_seglevel` | required | The cnamaste run's `clone_labels_int.tsv[.gz]` and `cnv_<ploidy>_seglevel.tsv`. With `binned_data.npz` and `rdrbaf_final_nstates<n>_smp.npz` beside them, combined.pdf holds cnamaste's own page. |
| `calicost_dir` | required | A CalicoST sample directory, `[<label>=]<dir>`, read at its best rectangle. |
| `numbat_dir` | required | A Numbat output directory. |
| `clonalscope_dir` | none | A Clonalscope export's `postprocess/`; without it, Clonalscope is a placeholder. |
| `calicost_sweep` | none | A CalicoST run of a parameter sweep, `<label>=<dir>`, once per run: its clones and profile in diagnosis.pdf. |
| `false_positive_sweep` | none | The label of a `calicost_sweep` run, once per run: its clones on diagnosis.pdf's false-positive page. |
| `sweep_runs` | none | Run keys scored with the reference on the copy-number error, on the bins they share. |
| `benchmark_yaml` | required | The runs ([Analysis](#analysis)). |
| `config_yaml` | required | The region bed and the genes. |
| `methods` | every method of `simulation` | The simulation methods scored. |
| `runs` | `numbat`, `clonalscope`, `calicost`, `cnamaste` | The sample's run keys. |
| `reference` | `cnamaste` | The run every other run is compared with. |
| `simulation_dataset` | `calicost-sim` | The simulations' wide table, `tsvs/<name>.tsv.gz`. |
| `real_data_dataset` | `ovarian` | The sample's wide table, `tsvs/<name>.tsv.gz`. |
| `table`, `metric_text` | `table`, `metric` | The LaTeX table's and metric text's names under `text/`. |
| `combined`, `diagnosis`, `appendix` | `combined.pdf`, `diagnosis.pdf`, `appendix.pdf` | The three PDFs' names under `plots/`; `combined` and `diagnosis` may be `.png`, page k after the first `<stem>_page<k>.png`. |
| `msr` | 300 | UGP's `min_snp_reads` of the bins scored, `MSR<msr>`. |
| `side_um` | 16 | The side of the spots compared, drawn and scored; a coarser bin is given to each such spot inside it. |
| `cnamaste_normal_fraction` | the run's log | `<clone>=<fraction>` per tumour clone of the cnamaste fit, its decoder's normal fraction: the cnamaste page's expected lines. By default each integer clone's `normal fractions per clone` from `cnamaste_<id>.log` (`plots.cnamaste_page.get_normal_fractions`); the normal clone held at 1. |
| `calicost_ploidy` | `unconstrained` | CalicoST's table, `cnv<_ploidy>_seglevel.tsv`. |
| `calicost_fallback_ploidy` | `diploid tetraploid` | The tables tried in order where `calicost_ploidy`'s leaves out a labelled tumour clone; `none` for none ([CalicoST clones missing from the unconstrained table](#calicost-clones-missing-from-the-unconstrained-table)). |
| `region_bed`, `chrom_sizes` | `data/hg38.regions.bed`, `data/hg38.chrom.sizes` | The arms and lengths of the sweep profiles' genome axis. |
| `assay` | `VISIUM3prime` | UGP's assay type, the bins' subdirectory. |
| `bins` | 100 | Each appendix UMI histogram's bin count. |

---

## Steps

A target runs the steps up to it: `snakemake -s run_cna_maste_reproducibility.smk ... <target>`.

| Target | Rules |
|---|---|
| `setup` | `check_inputs` (every input path); `clone_repo` (a missing repository cloned, an existing one checked); `build_env` (a missing environment built). |
| `prepare` | Per simulation: its 10x matrix, Numbat count matrix and allele counts, and cnamaste sample sheet. On the ovarian sample: a 16 um BAM where Space Ranger has none (UGP's `convert_bam_visium3prime.py`, `ovarian_bam_016um`); UGP's genotyping, phasing and 16 um bins (`ugp_ovarian`); then Numbat's 32 um counts and allele counts, the bins in CalicoST's layout (`ovarian_format_bb`) and cnamaste's sample sheet. Clonalscope's references; the sample's cell types (`celltyping_ovarian`). |
| `run_methods` | Numbat, Clonalscope and CalicoST per case and sample ([Baselines](#baselines)); cnamaste ([cnamaste](#cnamaste)); the runtime tables (`sim_runtime`, `vishd_runtime`). |
| `analysis` | The four steps of `scripts/analysis/` ([Analysis step](#analysis-step)). |

`-n` lists every job, and `--forcerun <rule>` reruns one step. Logs: [Run layout](#run-layout).

### Baselines

Each baseline runs on both cases, the CalicoST simulations (`zenodo_sim`) and the ovarian sample (`ovarian`), a rule
per step and sample, in its own environment, running a program of `scripts/baselines/`:

| Rule | Program | Writes |
|---|---|---|
| `sim_10x_matrix` | `common/h5ad_to_10x.py` | A simulation's 10x matrix from its h5ad, Ensembl IDs from `PATH_TO_ORDERED_GENES`. |
| `sim_numbat_count_mat` | `numbat/build_count_mat.R` | A simulation's Numbat count matrix. |
| `sim_numbat_allele_counts` | | The deposit's `numbat_allele_counts.tsv.gz`, as the deposit ran Numbat. |
| `ovarian_numbat_aggregate_bins` | `numbat/aggregate_bins.R` | The 16 um matrix and UGP's 16 um pileup summed into 32 um bins. |
| `ovarian_numbat_allele_counts` | `numbat/build_allele_counts.R` | Allele counts from that pileup and UGP's phasing (`numbat:::preprocess_allele`). |
| `numbat_run`, `numbat_postprocess` | `numbat/run.R`, `numbat/postprocess.py` | `run_numbat()`, then `<sample>.RNA.cells.tsv` and `cluster_states.tsv`. |
| `clonalscope_refs` | `clonalscope/build_refs_clonalscope.sh` | 200 kb bins, arm segments and `data-raw/` at Clonalscope's commit, shared by both cases. |
| `clonalscope_aggregate_matrix` | `common/aggregate_10x_matrix.R` | A coarser bin's matrix, the 16 um matrix summed. |
| `clonalscope_run`, `clonalscope_decode_states` | `clonalscope/tutorial.R`, `clonalscope/decode_states.R` | The Visium tumor_spots_tutorial, then per-cluster copies. |
| `calicost_sim` | `calicost/repro.sh` | CalicoST on one simulation with the deposit's configfiles. |
| `calicost_stage_input`, `calicost_ovarian` | `calicost/stage_input.sh`, `calicost/ovarian_config.py`, `calicost/run.sh` | UGP's formatted bins staged; one ovarian run's `configuration_cna`, run under a cgroup `MemoryMax`. |

| Method | Tutorial followed | Run |
|---|---|---|
| Numbat | [numbat vignette](https://kharchenkolab.github.io/numbat/articles/numbat.html) (scRNA mode, `ref_hca`), [spatial vignette](https://kharchenkolab.github.io/numbat/articles/spatial-rna.html) | `run_numbat()` with `genome = "hg38"`; the ovarian sample on 32 um bins. |
| Clonalscope | [Visium tumor_spots_tutorial](https://github.com/seasoncloud/Clonalscope/blob/identify_normal_cells_noWGS/samples/V11Y04-378-A1/tumor_spots_tutorial/README.md) (no paired DNA), call for call in `clonalscope/tutorial.R` | `FindNormalReference` (pca), arm segments, `RunCovCluster` and `MalignantAssignment` per sample; `RunCovCluster`'s `result_final` to per-cluster relative copies `Uest`, `TCN = round(ploidy * Uest)`. A run that left no `Cov_obj.rds` is skipped, as its log says. One core per sample: `BayesNonparCluster` is a single-threaded Gibbs sampler. |
| CalicoST | [CalicoST docs](https://calicost.readthedocs.io/), its README command `OMP_NUM_THREADS=1 python <CalicoST>/src/calicost/calicost_main.py -c configuration_cna` | Simulations: the deposit's own `configfile0..4`, only input, output and resource paths remapped; `parse_input.py` once, then the initializations in parallel. Ovarian: UGP's `bb-016um-formatted` staged (`unique_snp_ids.npy` re-saved for CalicoST's numpy 1.24), then a run per `runs` entry. CalicoST's scripts read their paths from `$CONFIGS`, `<OUTDIR>/configs/baselines/` (`common/config.sh`). |

> [!IMPORTANT]
> The ovarian sample's Numbat runs on 32 um bins: `run_numbat()` calls `stats::hclust` on a spot x spot distance, which
> stops above 65536 items, and the 16 um slide exceeds that.

> [!WARNING]
> `calicost/run.sh` needs systemd user scopes with the memory controller; a run over `MemoryMax` is killed (status
> `rc=137 (killed: OOM ...)`).

Diagnostics, run by hand with their paths as arguments: `calicost/compare_repro.py` (a CalicoST rerun against the
deposit's), `numbat/run_stats.py` (Numbat input sizes and step timings), `common/extract_deposit_results.py` (the
deposit's CalicoST, Numbat and STARCH labels and profiles).

### cnamaste

| Rule | Runs |
|---|---|
| `cnamaste_sim` | `run_cnamaste` per simulation, solver and seed, as cnamaste's `batch_zenodo_run.sh` runs each job: `-o hmrf.random_seed`, `hmrf.label_solver`, `hmm.initializer` and `ab_decode.decoder` (`lattice` -> `hmm_ab`, `calicost` -> `calicost`), into `nomixing_cnamaste_related/<sim>/<solver>/run_<id>/`. |
| `cnamaste_sim_plot` | `run_cnamaste_plot` of each run. |
| `cnamaste_sim_qa` | `run_sim_qa` of `qa_solver` against CalicoST and each other solver (`--one_shot --overwrite --calicost_ploidy`). |
| `cnamaste_ovarian`, `cnamaste_ovarian_plot` | `run_cnamaste` on the sample at cnamaste's defaults (`hmrf.label_solver: expansion`, alpha expansion), then `run_cnamaste_plot`. |
| `cnamaste_ovarian_binned_data` | `scripts/baselines/cnamaste/write_binned_data.py`: the run's `binned_data.npz`, which the cnamaste page reads and `cleanup3` no longer writes. |

### Analysis step

Each rule runs a function of `scripts/analysis/` in `CNAMASTE_ANALYSIS_ENV`, with `CNAMASTE_ROOT` set to
`CNAMASTE_ANALYSIS_REPO` and its `python/` and the repository on `PYTHONPATH`; each stops if `cnamaste` is imported from
anywhere else (`options.get_cnamaste_root`).

| Rule | Runs | Writes |
|---|---|---|
| `analysis_options` | | `<OUTDIR>/configs/analysis_options.yaml`, the options from the benchmark yaml and `ANALYSIS`; logs the analysis checkout's commit. |
| `analysis_simulation_metrics` | `simulation.write_simulation_metrics` | `tsvs/<simulation_dataset>.tsv.gz`, `.median.tsv.gz`, `tsvs/gene_level_cnp/`. |
| `analysis_real_data_metrics` | `real_data.write_real_data_metrics` | `tsvs/<real_data_dataset>.tsv.gz`, `.json`, `tsvs/<sample>.*.tsv.gz`. |
| `analysis_tables` | `tables.write_tables` | `text/` ([LaTeX tables](#latex-tables)). |
| `analysis_figures` | `figures.write_figures` | `plots/`, `pdf-pages/`, `tsvs/diagnosis.tsv.gz` ([Figures](#figures)). |

The two metric steps run in parallel, and the tables and figures read their tables. Both read `benchmark_yaml` and
`config_yaml`'s regions on chromosomes 1 to 22 and the genes inside them, compute every table, and write them only
then.

`write_simulation_metrics`, per simulation of `groundtruth_dir` in name order, once on each run's own segments and once
at the genes:

1. Read the truth (`cnamaste.truth.read_truth`) and each `methods` run with its method's reader. A Numbat run that
   stopped finding no CNV is a call of none, every spot clone 0 and `neu` on every true segment; a run with no output
   is left unscored, NaN but for its runtime.
2. Score the run against the truth ([Against the truth](#against-the-truth-simulations)); at the genes, write the
   simulation's `gene_level_cnp`.
3. Take each run's runtime and memory from `runtime` by `simulation_id` and `method`; merge the two levels into one wide
   table (`simulation_metrics.merge_levels`), refused where their runs, spots or runtime differ.

`write_real_data_metrics`, for `sample`:

1. Read the in-tissue spots, their cell types and each `runs` run ([Method inputs](#method-inputs)); give each run's
   clones to the spots at `side_um` (`real_data_metrics.get_run_spots`) and its profile to the genes
   (`gene_level.get_gene_states`).
2. `genome_composition` ([Genome composition](#genome-composition)).
3. Each run against `reference`: `cna_confusion`, clone ARI and CNA ARI; `clone_ari` and `cnp_similarity` of every two
   runs ([Between runs](#between-runs-real-data), [CNP similarity](#cnp-similarity)).
4. `segment_level_cnp` and the copy-number error of the runs giving copies
   ([Copy-number error](#copy-number-error), `real_data_metrics.get_copy_number_tables`); with `sweep_runs`,
   `cn_error.sweep`.
5. Each run's spot coverage, false-positive counts and runtime row; every per-run table merged into the wide table
   ([Tables](#tables)).

> [!NOTE]
> The pipeline passes no sweep run (`sweep_runs`, `calicost_sweep`, `false_positive_sweep`): diagnosis.pdf has no
> sweep pages, its pages 1 and 2 hold Numbat and cnamaste, and no `cn_error.sweep` is written. The sweep's runs are in
> `<OUTDIR>/configs/benchmark.yaml`: `calicost` (`msu20`), `calicost_016um_msu10`, `_msu15`, `_msu25` and `_msu50`.

---

## Metrics

### Method inputs

`scripts/metrics/read_method_input.py` reads a run with its method's reader (`get_method_reader`). Every method's normal
clone is clone 0, and a coarser bin's clone is given to each spot at `side_um` inside it.

| Method | Labels | Profile | Called per segment and clone |
|---|---|---|---|
| truth | `truth_clone_labels.tsv` | `truth_acn_profile.tsv` | AB pair |
| cnamaste | `clone_labels_int.tsv.gz`, the integer clones | `cnv_<ploidy>_seglevel.tsv`, else `cnv_seglevel.tsv` | AB pair |
| CalicoST | `clone_labels.tsv` of the best rectangle by `total_llf` | `cnv_seglevel.tsv`, unconstrained ploidy | AB pair |
| Numbat | `clone_post_<i>.tsv`, the last iteration | `bulk_clones_final.tsv.gz` | `cnv_state` |
| Clonalscope | `<sample>.<assay>.cells.tsv` of a Visium tutorial export's `postprocess/`, its normal clusters first | `<sample>.<assay>.profile.tsv` | total copies `TCN` |
| cell type | `spot_annotation.tsv` | | |

- cnamaste: `run_cnamaste` merges clones whose (A,B) agree into an integer clone named by a member, and keeps every
  clone in its seglevel table; an integer clone's copies are those of the clone naming it
  (`read_method_input.get_labelled_profile`). Its `clone_ari` is therefore of its integer clones, as `run_sim_qa`'s
  `int_clone_ari`.
- Clonalscope: its clusters are kept as called, its normal clusters (`malignant` N) each a clone, on the simulations and
  in the real-data comparisons; the spot counts and the figures merge them into one clone 0
  (`real_data_metrics.get_normal_clones`). A run that halts writes no export.
- CalicoST: a run whose unconstrained table leaves out a labelled clone is read at a fallback ploidy
  ([CalicoST clones missing from the unconstrained table](#calicost-clones-missing-from-the-unconstrained-table)).
- The H&E image and the bins' places come from a Space Ranger directory of the sample, at any bin side.

### CNA states and fields

A CNA state is one of Numbat's six, its `loh` written `cnloh`, by an AB pair's balance and total copies, each pair as
Numbat's [file descriptions](https://kharchenkolab.github.io/numbat/articles/descriptions.html) name it:

|            | total 2     | above 2        | below 2    |
| ---------- | ----------- | -------------- | ---------- |
| balanced   | neu (1:1)   | bamp (2:2)     | bdel (0:0) |
| unbalanced | cnloh (2:0) | amp (2:1, 3:1) | del (1:0)  |

A state is `neu` exactly where its pair is (1, 1). Each method gives, per clone, some of five fields:

| Method | `ascn` | `Uest` | `tcn` | `state` | `tcn_state` |
| --- | --- | --- | --- | --- | --- |
| truth, cnamaste, CalicoST | `A\|B`, as the table phases it | | `A + B` | of the phase-free pair | of `A + B` |
| Numbat | the state's pair: `1\|1`, `1\|0`, `2\|0`, `0\|0`; `unknown` for `amp`, `bamp` | | its total, or `unknown` | `cnv_state` | of the state |
| Clonalscope | | the cluster's mean coverage fold change, `AssignCluster`'s `Uest` (`REL_CN`) | `TCN`, `round(2 * Uest)` | `del` (< 2), `neu` (2), `amp` (> 2) | its `state` |

A field is `uncalled` where the run calls nothing there. `tcn_state` is the state of total copies alone, `del` (< 2),
`neu` (2) or `amp` (> 2), as copytyping-benchmark's `analysis/scripts/cn_transform.py` (1204da9) reads total copies:
`cnloh` is `neu`, `bdel` `del`, `bamp` `amp`.

### Genes and joint segments

- **Genes**: the GTF's `gene` features wholly inside a region of the bed on chromosomes 1 to 22, where the truth and
  every method's input lie; `START` 0-based, `END` exclusive. A gene takes, of a run's segments and of the truth's, the
  one overlapping it most, the first on a tie (`gene_level.get_max_overlap`); a gene no segment overlaps is `uncalled`.
- **Joint segments**: the intervals between consecutive breakpoints of any run, kept where some run has a segment
  (`scripts/metrics/segment_level.py`), so each lies inside at most one segment of each run, as copytyping-benchmark's
  `analysis/scripts/cn_error.py` (1204da9) joins competing runs. Only runs whose profile gives copies are joined
  (`cna_states.get_profile_copies`): CalicoST's and cnamaste's AB pairs and Clonalscope's total copies; Numbat's states
  give none.

### Against the truth (simulations)

A cell is one spot a run labels on one row of its profile that a true segment overlaps: one of its own segments
(`qa.metrics`), its `(segment)` columns, or one gene, its `(gene)` columns. Each row counts once, whatever its length.
A row takes, per true clone, the copies of the true segment overlapping it most (`segment_true_ab_pairs`); a spot's
true AB pair is its true clone's and its call is its fitted clone's, so no clone is matched to another. AB pairs are
compared phase-free (`phase_free`). Rows the run calls no state on in some clone, and spots it does not label, are not
scored.

| Column | Panel | Definition | cnamaste's `qa.metrics` |
|---|---|---|---|
| `sensitivity` | CNA sensitivity | 1 - (truly altered cells called `neu`) / (truly altered cells) | `sensitivity` |
| `specificity` | Normal specificity | 1 - (truly `neu` cells called altered) / (truly `neu` cells) | `specificity` |
| `exact_altered` | N-CNA recall | (truly altered cells called their true AB pair) / (truly altered cells); NaN for a run of states | `exact_altered` |
| `state_altered` | CNA recall | (truly altered cells called their true CNA state) / (truly altered cells) | |
| `clone_ari` | Clone ARI | Adjusted Rand index of the run's clones against the true ones, over the spots it labels | `clone_ari` |

Truly altered is a true pair other than (1, 1). A rate is NaN where its denominator is 0.

> [!NOTE]
> `sensitivity`, `specificity` and `state_altered` are read from one table of cells by `(true CNA state, called CNA
> state)` (`supervised_metrics.get_state_metrics`). As a state is `neu` exactly where its pair is (1, 1), the first two
> count the cells `qa.metrics` counts on AB pairs (`tests/test_simulation_metrics.py`). `state_altered` counts a cell
> only where the two states agree, so `amp` called `cnloh` is wrong. The called state is the called pair's (cnamaste,
> CalicoST), the `cnv_state` (Numbat), or that of the total copies (Clonalscope). Numbat's `amp` covers 2:1 and 3:1,
> so its state can match a true (1, 3) where a pair cannot. Clonalscope calls a true `cnloh` `neu`, missed, and a true
> `bamp` or `bdel` `amp` or `del`, detected but not its state.

> [!IMPORTANT]
> On a run's own segments, a true CNA smaller than the run's segment around it drops out: the segment takes the copies
> of the true segment overlapping it most. A method with coarse segments, such as Numbat's whole chromosomes outside
> its calls or Clonalscope's arms, is scored on fewer and larger CNAs; where no CNA remains, sensitivity and CNA recall
> are NaN. At the genes, every true CNA a gene lies in is scored.

### Between runs (real data)

On the in-tissue spots at `side_um` both runs label, at the finest field both give, their `level`
(`cna_states.get_comparison_field`): `ascn`, AB pairs phase-free as `(major,minor)`, where both give pairs (cnamaste,
CalicoST); `tcn_state` where either gives total copies alone (Clonalscope); else `state` (Numbat).

| Function | Definition |
|---|---|
| `compute_cna_confusion` | Cells by `(state of A, state of B)` at the `level`, `uncalled` included: per shared spot, the genes at each pair of its two clones' states, summed over the spots. |
| `compute_clone_ARI` | Adjusted Rand index of the two runs' clones over the shared spots. |
| `compute_pairwise_clone_ARI` | `compute_clone_ARI` of every two runs, each pair over its own shared spots, a run against itself 1. |
| `compute_cna_ARI` | Adjusted Rand index of the two runs' states over the cells both call a state on: the confusion as a contingency table, `uncalled` left out ([Hubert and Arabie 1985](https://doi.org/10.1007/BF01908075), as scikit-learn's `adjusted_rand_score`). |
| `compute_state_jaccard` | Of A's state `s` and B's state `t`, the cells where A calls `s` and B calls `t` over those where A calls `s` or B calls `t`, over the cells both call ([Jaccard 1912](https://doi.org/10.1111/j.1469-8137.1912.tb05611.x)); NaN where neither holds such a cell. Unlike the ARI, it compares states by name: a run calling each `del` of the other `cnloh` scores 0 at `(del, del)`. |

### CNP similarity

The CNP similarity of two runs compares their tumor clones' copy-number profiles (CNPs) at the genes, at their CNA
`state`, no spot and no truth. Let $\mathcal{G}$ be the [genes](#genes-and-joint-segments), and $s^a_i(g)$ the state of
run $a$'s clone $i$ at gene $g$, or $\varnothing$ where the run calls none there. For tumor clones $i \in K_a$ and
$j \in K_b$, the altered genes are those both call and at least one calls other than `neu`:

```math
\mathcal{A}_{ij} = \{\, g \in \mathcal{G} : s^a_i(g) \neq \varnothing,\ s^b_j(g) \neq \varnothing,\
\big(s^a_i(g), s^b_j(g)\big) \neq (\texttt{neu}, \texttt{neu}) \,\}
```

```math
J_{ij} = \frac{1}{|\mathcal{A}_{ij}|} \sum_{g \in \mathcal{A}_{ij}} \mathbb{1}\left[ s^a_i(g) = s^b_j(g) \right],
\qquad S(a, b) = \max_{i \in K_a,\ j \in K_b} J_{ij}
```

- $K_a$ is run $a$'s clones but its normal ones, clone 0 and Clonalscope's normal clusters
  (`real_data_metrics.get_tumor_cnps`).
- A gene both call `neu` is left out; one called `neu` by one and altered by the other is a mismatch.
- $J_{ij}$ is NaN where $\mathcal{A}_{ij}$ is empty, and the maximum skips it; $S(a, b)$ is NaN where every $J_{ij}$
  is. A tie goes to the first pair in clone order. $S(a, b) = S(b, a)$, and $S(a, a) = 1$ where some tumor clone of
  $a$ alters a gene.
- $J_{ij}$ is the Jaccard index of the two clones' altered genes among those both call, a gene both alter counted in
  the intersection only where their states match.
- `unsupervised_metrics.get_altered_genes` is $\mathcal{A}_{ij}$, `compute_cnp_similarity` $J_{ij}$, and
  `compute_pairwise_cnp_similarity` $S(a, b)$ for every two runs.

> [!IMPORTANT]
> Clonalscope's state is that of its total copies, `del`, `neu` or `amp`, so a gene another run calls `cnloh`, `bamp`
> or `bdel` never matches it.

### Genome composition

`compute_genome_composition`: per clone, the share of the genome in each CNA state and `uncalled`, summing to 1. The
genome is the regions' bases, a segment weighing its bases inside them, or the genes, 1 each. Its `all` row weighs each
clone by its spots: the share of the run's cells, a labelled spot at a gene, in each state.

### Copy-number error

Whether a run's copies reconstruct the counts of its clones, as copytyping-benchmark's `analysis/scripts/cn_error.py`
(1204da9) scores them, with copytyping's pseudobulks (e0341e4, `inference/model_utils.py`), in
`scripts/metrics/cn_error.py`. It scores the runs whose profile gives copies: CalicoST's and cnamaste's AB pairs, RDR
and BAF; Clonalscope's total copies, RDR alone. Numbat's states give no copies.

**Bins.** UGP's bins at `side_um`, `bb/MSR<msr>/<assay>/` (`read_method_input.read_bin_counts`), inside an autosomal
region of the bed. UGP bins the SNPs adaptively, each bin at least `min_snp_reads` SNP reads, and sums into each bin its
SNPs' phased A and B counts and the UMI of each gene it overlaps most; a bin with no SNP has no BAF.

**Copies per bin.** A bin takes the joint segment overlapping it most and, from each run, the copies of the segment
holding it, per clone `k`: `A`, `B` and `C = A + B`; uncalled where the run has no segment there.

**Clone pseudobulks.** A clone's counts are summed over its spots at `side_um`:

    X[k,u] = sum_{n in k} X[u,n]      the clone's UMI in bin u
    T[k]   = sum_u X[k,u]             its UMI over the bins
    b[k,u], d[k,u]                    its B-allele and total allele counts in bin u
    lambda_u = sum_{n in ref} X[u,n] / sum_{n in ref} T_n

`lambda_u`, the baseline, is bin `u`'s share of the UMI of the spots of a normal cell type, neither `Epithelial` nor
`Low confidence` (`read_method_input.get_normal_type`); it sums to 1.

**Observed and expected.**

    RDR_obs[k,u] = X[k,u] / (lambda_u T[k])          RDR_exp[k,u] = C[u,k] / sum_u' lambda_u' C[u',k]
    BAF_obs[k,u] = b[k,u] / d[k,u]                   BAF_exp[k,u] = B[u,k] / (A[u,k] + B[u,k])

`RDR_exp`'s sum runs over the bins the run calls for the clone. The diploid prediction, every clone 1:1, gives
`RDR_exp` 1 and `BAF_exp` 0.5.

**Phasing.** UGP's counts are phased, a bin's B allele one haplotype that may switch between bins; a run's `A|B`
carries its own phase. A bin is read in the orientation, `BAF_exp` or `1 - BAF_exp`, with the smaller squared error
summed over the run's clones weighed by their spots: one orientation for every clone (`cn_error.get_baf_flips`,
`#flipped-bins`).

**Bins scored.** Where every run competing on the channel calls every clone, every observed and expected value is
finite, and some run calls a CNA in some clone with spots, copies not 1:1 or, given total copies alone, not 2
(`cn_error.get_altered`). The bins scored therefore depend on the runs compared.

| Column | Definition, `n_k` the clone's spots |
|---|---|
| `mse_rdr` | `sum_{k,u} n_k (log2 RDR_obs - log2 RDR_exp)^2 / sum_{k,u} n_k`, every clone. |
| `mse_baf` | `sum_{k,u} n_k (BAF_obs - BAF_exp)^2 / sum_{k,u} n_k`, `BAF_exp` in the bin's orientation. |
| `null_rdr`, `null_baf` | The same error of the diploid prediction. |
| `ss_rdr`, `ss_baf` | Skill score, `1 - mse / null`: 1 a perfect reconstruction, 0 the diploid prediction's, below 0 worse. |

MSE weighs a clone by its spots, the normal clone among them. `#tumor-clones` counts a run's clones but clone 0 and
Clonalscope's normal clusters (`real_data_metrics.get_normal_clones`).

> [!NOTE]
> `RDR_exp` is a pure tumor's: a clone's spots holding normal cells read an RDR nearer 1 than its copies give, which
> raises `mse_rdr` and lowers `ss_rdr`. No purity is modelled, as copytyping-benchmark's `cn_error.py:719-721` calls
> `pseudobulk_rdr_matrix` without its `purity`.

### CalicoST clones missing from the unconstrained table

A CalicoST clone can be in `clone_labels.tsv` yet absent from `cnv_seglevel.tsv`, its unconstrained integer copies
(CalicoST c1abcae, [`calicost_main.py`](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/calicost_main.py),
[`find_integer_copynumber.py`](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/find_integer_copynumber.py)):

- The labels are the final fit's `new_assignment` (`calicost_main.py:388`). Every clone merge,
  `similarity_components_rdrbaf_neymanpearson` and `merge_by_minspots`, runs before that fit is saved (`:299`), so no
  merge removes a labelled clone.
- The unconstrained copies (`hill_climbing_integer_copynumber_fixdiploid`, `:334`, `:337`) need the clone's diploid
  balanced state: a state holding at least 10% of its bins, 2% on the retry, its BAF within 0.05 of 0.5
  (`find_integer_copynumber.py:73-81`). Where neither attempt finds one, the bare `except` skips the clone
  (`calicost_main.py:338-340`) and logs nothing; the run's log lacks only its `max med ploidy = None, clone <s>` line
  (`:342`).
- The diploid, triploid and tetraploid tables (`hill_climbing_integer_copynumber_oneclone`, `:331`) need no such state,
  so they hold the clone.

A run's `rdrbaf_final_nstates<n>_smp.npz` shows which: `pred_cnv[:, s]` and `new_p_binom[:, s]` give the share of clone
`s`'s bins in each state and its BAF. `read_method_input.read_calicost` therefore reads the `calicost_ploidy` table
where it holds every labelled tumour clone, else the first of `calicost_fallback_ploidy` that holds them all, else the
`calicost_ploidy` table, logged (`get_calicost_ploidy`). A run read at a fallback names it in its title, e.g.
`CalicoST (16 um, min-snp-umi: 10, diploid)`.

---

## Outputs

### Run layout

| Path | Contents |
|---|---|
| `<OUTDIR>/setup/` | `<KEY>.ok` per checked repository or environment; `inputs.ok`. |
| `<OUTDIR>/configs/` | The configs the steps read: CalicoST's scripts' `baselines/{datasets,metadata,calicost.<case>}.yaml`, `ugp/`, `cnamaste.<case>.yaml`, `benchmark.yaml`, `config.yaml`, `analysis_options.yaml`. |
| `<OUTDIR>/references/clonalscope/` | Clonalscope's references. |
| `<OUTDIR>/zenodo_sim/{inputs,numbat,clonalscope,calicost,cnamaste}/` | The simulations' inputs and runs. |
| `<OUTDIR>/ovarian/{inputs,ugp,celltyping,numbat,clonalscope,calicost,cnamaste}/` | The sample's. |
| `<OUTDIR>/zenodo_sim/simulation.runtime.tsv`, `<OUTDIR>/ovarian/runtime.<sample>.tsv` | Each run's runtime and memory (`scripts/metrics/runtime.py`): a Snakemake `benchmark:`'s `s` / 60 and `max_pss` / 1024, CalicoST's `/usr/bin/time -v`; NA where not recorded. |
| `<OUTDIR>/benchmark_results/` | The analysis' `tsvs/`, `text/`, `plots/` and `pdf-pages/`. |
| `<OUTDIR>/logs/<step>/<rule>[.<wildcards>].log` | Each job's log, its runtime and memory (`benchmark:`) in the `.tsv` beside it. |
| `<OUTDIR>/<case>/<method>/{logs,bench}/` | A method's per-sample run logs and benchmarks; CalicoST's ovarian runs in `<OUTDIR>/logs/baselines/`. |

### Tables

Under `<OUTDIR>/benchmark_results/`, names at their defaults. Tables are tab-separated and gzipped, the gzip header's
time 0, so equal tables give equal bytes (`scripts/metrics/outputs.py`).

| File | Contents |
|---|---|
| `tsvs/calicost-sim.tsv.gz` | The simulations' wide table, a row per simulation and method (below). |
| `tsvs/calicost-sim.median.tsv.gz` | A row per simulated class and method: `sim_id` (the simulation's id less `_random<r>`), `method`, `#replicates`, each metric's median over the class' replicates, NaN skipped (`plots.simulation.get_class_medians`). |
| `tsvs/gene_level_cnp/<sample_id>.gene_level_cnp.tsv.gz` | Per simulation, a row per gene: `gene_id`, `gene`, `CHR`, `START`, `END`; `truth clone<k> ascn` and `state` per true clone, `normal` clone 0; `<method> clone<k> <field>` per run, clone and field. |
| `tsvs/ovarian.tsv.gz` | The sample's wide table, a row per run (below). |
| `tsvs/ovarian.json` | Its summary (`real_data_metrics.get_summary`): `dataset`, `sample`, `reference`, then per run `run`, `method`, `resolution`, its `metrics` (`SUMMARY_METRICS`) and their `provenance` (`PROVENANCE`) counts, NaN written `null`. |
| `tsvs/<sample>.gene_level_cnp.tsv.gz` | A row per gene: each run's fields per clone. |
| `tsvs/<sample>.segment_level_cnp.tsv.gz` | A row per joint segment: `CHR`, `START`, `END`, then each copy-giving run's fields per clone, `uncalled` where it has no segment. |
| `tsvs/<sample>.genome_composition.tsv.gz` | Per run and clone, and its `all` row, the shares by state. |
| `tsvs/<sample>.cna_confusion.tsv.gz` | Per run, the cells by `(reference state, run state)` at its `level`, long. |
| `tsvs/<sample>.clone_ari.tsv.gz` | Every ordered pair of runs: `run_a`, `run_b`, `#shared-spots`, `clone_ari`. |
| `tsvs/<sample>.cnp_similarity.tsv.gz` | Every ordered pair of runs: `run_a`, `run_b`, `clone_a`, `clone_b` (the pair the maximum picks, NA where $S$ is NaN), `#altered-genes`, `cnp_similarity`. |
| `tsvs/<sample>.cn_error.sweep.tsv.gz` | With `sweep_runs` only: the reference and sweep runs' copy-number error, scored as one panel on the bins they share. |
| `tsvs/diagnosis.tsv.gz` | The numbers behind diagnosis.pdf's UMI page, a row per run and UMI measure. |

The simulations' wide table, each metric per level (`simulation_metrics.merge_levels`):

| Columns | Meaning |
|---|---|
| `sample_id` | `numcnas<a.b>_cnasize<s>_ploidy<p>_random<r>` (`plots.simulation.SAMPLE_ID`). |
| `method`, `solver` | A key of the benchmark yaml's `simulation`; a cnamaste solver's `solver`. |
| `#spots`, `#called`, `spot-coverage-rate` | The truth's spots, those the run labels, and `#called` / `#spots`. |
| `<metric> (segment)`, `<metric> (gene)` | [Against the truth](#against-the-truth-simulations), on the run's own segments and at the genes. |
| `Runtime [mins]`, `Memory [GB]` | Its `runtime` row. |

The sample's wide table, a row per `runs` in order, each column once (`real_data_metrics.merge_run_tables`): a column
two tables hold must agree on every run, and a count, `#<name>`, is an integer, NA where a table has no row of the run.

| Columns | From |
|---|---|
| `sample`, `run`, `method`, `resolution` | The run; `resolution` its bins' side, `<side>um`. |
| `#in-tissue`, `#called`, `spot-coverage-rate`, `#normal-spots`, `#tumor-spots` | Its calls on the in-tissue spots at `side_um`, Clonalscope's normal clusters one clone 0 (`get_spot_table`, `qc_metrics.get_coverage`). |
| `#normal-spots-<type>`, `#tumor-spots-<type>`, `#false-positive-tumor-spots`, `#normal-type-spots`, `false-positive-rate` | Its calls per broad cell type, `Epithelial` first; its tumor calls on a normal type, neither `Epithelial` nor `Low confidence`; its calls on such a type; the former over the latter (`qc_metrics.get_false_positive_counts`). |
| `neu` ... `uncalled` | Its genome composition's `all` row. |
| `level`, `#shared-spots`, `clone_ari`, `cna_ari`, `reference` | Against the reference run, NA for the reference. |
| `#clones`, `#tumor-clones`, `#rdr-bins`, `#baf-bins`, `#flipped-bins`, `mse_*`, `null_*`, `ss_*` | Its copy-number error, the runs giving copies. |
| `setting`, `Runtime [mins]`, `Memory [GB]` | Its row of the sample's `runtime` table: its method's, at the side of the bins it labels, its `setting` a word of the run directory's name or its parent's; a run matching no row, or several, stops the step (`get_runtime_table`). |

The false-positive rate is FP / N, the spots of a normal type the negatives, N = TN + FP (James et al., *An
Introduction to Statistical Learning*, 1st ed., sec. 4.4.3, Table 4.7, https://doi.org/10.1007/978-1-4614-7138-7); an
epithelial, low-confidence or untyped spot is in neither count.

### LaTeX tables

`tables.write_tables` writes, under `text/`:

- `table.tex` (`latex_tables.write_latex_table`): from the sample's wide table, a booktabs tabular in a `table` float
  (`[!h]`, centred, a `0.85\linewidth` minipage), captioned and labelled `tab:data_runtime` (`latex_tables.CAPTION`,
  `LABEL`). A row per run in method order, cnamaste, CalicoST, Numbat, Clonalscope (`plots.style.get_ordered_runs`):
  `Method`, `Resolution` [$\mu$m], `Spot coverage rate` and `False-positive rate` to 0.001, `Runtime [mins]` to 0.1,
  `Memory [GB]` to 0.01, `$-$` where not recorded (`latex_tables.get_table_columns`). Needs `\usepackage{booktabs}`.
- `metric.tex` (`metric_text.write_metric_text`): a section per dataset, a paragraph per metric with its columns and
  its definition as computed, each followed by a `% shown in: <page> (<panels>)` line per page drawing it. The step stops
  where a scored metric column has no paragraph (`metric_text.get_missing_columns`). Needs `\usepackage{amsmath,amssymb}`.
- Each typeset on letter pages by `pdflatex` (`latex_tables.get_article_tex`), `table.pdf` and `metric.pdf`, no date or
  ID in the PDF, so equal tables give equal bytes; pdflatex's logs `table.log` and `metric.log`.

### Figures

`figures.write_figures` writes three PDFs under `plots/` and each page alone under `pdf-pages/`, cropped to its content
with 0.02 in of white round it (`plots.pages.SAVE_OPTIONS`): a paper figure's page under its name
(`plots.pages.PAGE_NAMES`), any other as `<stem>-P<k>.pdf`. A page whose condition below fails is left out.

- Every page shows each panel's letter alone, `(a)`, `(b)`, ...; each PDF's last page names the panels, per page
  `Page <k>` and its title, then `(x) <title>` per panel (`layout.plot_panel_titles`). Text sizes:
  `layout.set_text_sizes`.
- Runs are drawn in method order, cnamaste, CalicoST, Numbat, Clonalscope, each method one colour on every page, and
  named by method, bin side and `msu<k>` alone (`plots.style.get_run_label`); CalicoST's `msu<k>`,
  [`min_snpumi_perspot`](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/arg_parse.py#L30), is
  written `min-snp-umi: <k>` (`get_setting_label`).
- Clonalscope's normal clusters are drawn as one clone 0, its states `del`, `neu`, `amp` as Numbat's `DEL`, `NEU`,
  `AMP`; without `clonalscope_dir` it is a placeholder.
- The simulation metrics and solver pages leave out the classes (6,3, 50M) and (3,3, 50M)
  (`plots.pages.DROPPED_CLASSES`); diagnosis.pdf draws them in full.

**combined.pdf**

| Page | Saved as | Draws |
|---|---|---|
| 1 | `main-fig-zenodo-sim-result` | The simulations' metrics at the genes, a box per simulated class and method: (a) CNA sensitivity, (b) clone ARI, (c) CNA recall, (d) N-CNA recall (`plots.pages.plot_simulation_pages`). |
| 2 | `main-fig-cnamaste-result` | With `binned_data.npz`, cnamaste's own page (`plots.cnamaste_page.plot_cnamaste_page`, after `cnamaste.plot.layout.plot_combined_page`): (a) the H&E slide; (b) the integer clones in place, keyed with each one's share of the in-tissue spots; (c) onward, a letter per clone, its pseudobulk RDR then BAF per bin, with cnamaste's expected lines at its normal fraction (`cnamaste_normal_fraction`, `hmm_ab_decoder.mixture`); last, the copy-number profile per integer clone. Unlike cnamaste's page, the RDR track is 0 to 4 (`TRACK_RDR_YLIM`) and a (1,1) bin's dot grey. |
| 3 | `main-fig-cnp-compare` | (a) Each method's clone profile on the cnamaste run's bins (`layout.plot_profile_page`, `plots.pages.get_bin_profile`), only the chromosomes where some clone holds a CNA, each row named with its clone's share of the in-tissue spots: copies as cnamaste's `plot_state_profile` draws them, (1,1) white; Numbat's and Clonalscope's states in Numbat's colours, CNLoH's and DEL's swapped; keyed by CNA state and by `(A,B)` state grouped by BAF. (b) CNP similarity and (c) clone ARI of every two runs (`plots.run_comparison.get_run_comparison_panels`); (d) each run's RDR and BAF skill scores (`plots.cn_error.draw_skill_score_panel`). |
| 4 | `combined-P<k>` | The panel names. |

**diagnosis.pdf**

| Page | Draws |
|---|---|
| 1 | The clones of Numbat, of each `false_positive_sweep` run and of cnamaste, each with its false-positive tumor spots, a tumor call on a spot of a normal type, in a third colour. |
| 2 | With a `bb_dir` at `side_um`: each run's false-positive spots' total UMI and SNP-covering UMI against its other called in-tissue spots', as densities, with each group's spots and median (`plots.umi`). |
| 3 | The copy-number error, each run per channel, RDR and BAF: (a) the MSE on a log axis, the diploid prediction's a black tick across its bar, (b) the skill scores (`plots.cn_error`). |
| | With `calicost_sweep`: each sweep run's clones, then its profile, named `CalicoST (<side> um, <label>)`. |
| | The simulation metrics page, (a) normal specificity and (b) spot coverage rate first, three rows of two, and the solver page, every class (`plots.pages.plot_simulation_diagnosis`). |
| last | The panel names. |

**appendix.pdf**

| Saved as | Draws |
|---|---|
| `suppl-fig-umi` | With `bb_dir`: a row per `bb_dir`, its per-spot total UMI (its `bb` AnnData) and SNP-covering UMI (A + B of its `multi_snp` allele matrices) (`plots.umi.plot_umi_page`). |
| `suppl-fig-ovarian-label` | The H&E slide, each spot's cell type and each method's clones, each panel keyed inside with each label's share of the in-tissue spots (`layout.plot_annotation_page`). |
| `suppl-fig-cna-composition` | Each run's proportion of cells in `del`, `amp`, `cnloh` and `bamp`, its `all` row (`plots.run_comparison.plot_genome_composition`). |
| `suppl-fig-cna-jaccard` | Each run's state Jaccard index against the reference at its `level`, a row per reference state, a column per run state (`plots.run_comparison.plot_cna_jaccard`). |
| `appendix-P<k>` | The simulations' metrics at the genes, as combined.pdf's page 1. |
| `suppl-fig-zenodo-sim-runtime` | The simulations' runtime and memory per method, log axes. |
| `suppl-fig-zenodo-sim-ablation` | Where the table holds cnamaste's solvers, `run_sim_qa`'s solver comparison (`qa.plot.plot_validation`, `ddc060b`): (a) CNA sensitivity, (b) 1 - normal specificity on a log axis, a run at specificity 1 drawn at a floor ticked 0.00 (`plots.simulation.get_complement`), (c) N-CNA recall, (d) clone ARI of the integer clones, a box per class and solver (`plots.simulation.plot_solver_ablation`). |
| `appendix-P<k>` | The panel names. |

---

## Code

| Module | Holds |
|---|---|
| `scripts/analysis/` | The analysis step: `options` (the options and `get_cnamaste_root`), `simulation`, `real_data`, `tables`, `figures`. |
| `scripts/metrics/read_method_input.py` | The readers: each method's run, the spots, their cell types, and UGP's counts per bin and spot. |
| `scripts/metrics/cna_states.py` | The CNA states, and each method's fields per segment and clone. |
| `scripts/metrics/gene_level.py`, `segment_level.py` | The genes inside the region bed and the runs' joint segments, and each run's profile on them. |
| `scripts/metrics/supervised_metrics.py` | A run against the truth; `qa.metrics`' functions copied, on a run's labels and profile. |
| `scripts/metrics/unsupervised_metrics.py` | Two runs against each other. |
| `scripts/metrics/qc_metrics.py` | Genome composition, spot coverage and false-positive counts. |
| `scripts/metrics/cn_error.py`, `pseudobulk.py` | The copy-number error; each clone's pseudobulk RDR and BAF per bin beside its expected ones (`pseudobulk.py`, imported by no module). |
| `scripts/metrics/simulation_metrics.py`, `real_data_metrics.py` | The simulations' and the sample's wide tables. |
| `scripts/metrics/latex_tables.py`, `metric_text.py`, `outputs.py` | The LaTeX table, the metric text, and the outputs' paths and writers. |
| `scripts/metrics/runtime.py` | The runtime tables. |
| `scripts/plots/style.py`, `layout.py`, `pages.py` | Colours, names, order and panel rows; the spatial and profile pages on cnplot's genome axis; each figure written to its PDF and each page alone. |
| `scripts/plots/simulation.py`, `run_comparison.py`, `cn_error.py`, `umi.py`, `cnamaste_page.py` | The simulation pages, the run comparisons, the copy-number error page, the UMI pages, cnamaste's own page. |
| `scripts/platforms/visium.py` | Visium HD square bins: barcodes read as grid places, placed on the image, given to finer bins. |
| `scripts/baselines/` | The baselines' programs ([Baselines](#baselines)). |
| `scripts/celltyping/` | The cell types ([celltyping.md](celltyping.md)). |

---

## Test

`tests/` holds the tests; they need no `tmp/input_data/`, and `pytest.ini` puts the repository on the import path. Run
from the repository, in the cnamaste analysis environment:

```sh
UV_PROJECT_ENVIRONMENT=/path/to/envs/cnamaste_analysis_env uv run --no-sync --project /path/to/repos/cna_maste_analysis pytest
```

The code is formatted with ruff (`ruff.toml`) and snakefmt, both at 120 columns:

```sh
uvx ruff@0.16.8 format . && uvx ruff@0.16.8 check .
uvx snakefmt@2.1.0 --line-length 120 run_cna_maste_reproducibility.smk
```
