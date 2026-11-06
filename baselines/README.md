# baselines: Clonalscope, Numbat and CalicoST on Visium HD 3prime ovarian (H1-N3KXRGQ) and the CalicoST simulations

One command runs a method on a case, `bash scripts/run_baseline.sh <method> <case>`; the cases are `ovarian` (Visium HD 3prime ovarian H1-N3KXRGQ) and `calicost-sim` (90 CalicoST simulations, [doi:10.5281/zenodo.14175627](https://doi.org/10.5281/zenodo.14175627)).

## Layout
```
configs/datasets.yaml           input paths per dataset, shared by every method
configs/metadata.yaml           reference files, conda env per method, CalicoST clone
configs/<method>.<case>.yaml    samples, method arguments and output dir of one method on one case
configs/nondefault_params.yaml  every non-default method argument, its package default and its tutorial source
envs/<method>.yaml              conda recipe
scripts/                        entrypoints: run_baseline.sh, status_report.sh, evaluation
scripts/common/                 helpers shared by the methods
scripts/<method>/               a method's own scripts
logs/                           <method>_<case>.{status,<stage>.log,resources.tsv}, calicost_<tag>.*
```
Paths in `configs/` are on the ragr cluster. `scripts/common/config.sh` holds the conda root, the one path read before the YAML; its base python (with pyyaml) reads the rest. A snakefile reads three configs merged: `--configfile configs/datasets.yaml configs/metadata.yaml configs/<method>.<case>.yaml`. No script reads `configs/nondefault_params.yaml`; it records the arguments that `scripts/numbat/run.R` and `scripts/clonalscope/tutorial.R` set. Scripts run from any directory and write `logs/` here.

Script names below are relative to `scripts/`.

| Script | Role |
| --- | --- |
| `run_baseline.sh` | entrypoint: a method's snakemake stages on one case, with status and resource logs |
| `common/config.sh` | conda root and `get`, the YAML reader of the shell scripts |
| `common/resource_sampler.sh` | per-minute RSS and CPU of a process tree |
| `common/samples.smk` | dataset, samples and 10x matrix dirs, included by every snakefile |
| `common/prepare_sim_inputs.smk`, `common/h5ad_to_10x.py` | simulation h5ad to 10x matrix dir, read by Clonalscope and Numbat |
| `common/aggregate_10x_matrix.R` | a 16um 10x matrix summed into coarser bins |

Install commands run from `baselines/` in bash, after `source scripts/common/config.sh`.

## Clonalscope

- Method: https://github.com/seasoncloud/Clonalscope
- Commit: `07dd0c576b6754d18dda86ec8905d8a9be5e77b0` (main; package version 1.1.0)
- Tutorial: [Visium tumor_spots_tutorial](https://github.com/seasoncloud/Clonalscope/blob/identify_normal_cells_noWGS/samples/V11Y04-378-A1/tumor_spots_tutorial/README.md) (no paired DNA), followed call for call in `scripts/clonalscope/tutorial.R`

Stages: `clonalscope/run.smk` (FindNormalReference(pca), arm segments, RunCovCluster, MalignantAssignment; an `aggregate` sample first sums the 16um matrix), then `clonalscope/decode_states.smk` (RunCovCluster `result_final` to per-cluster relative CN `Uest`, `TCN = round(ploidy * Uest)`, cells, segments and coverage TSVs). Non-default arguments and their sources: `configs/nondefault_params.yaml` (`clonalscope.visium_tutorial`).

```sh
E=$(get 'c["software"]["envs"]["clonalscope"]' configs/metadata.yaml)
$CONDA/bin/mamba env create -p $E -f envs/clonalscope.yaml
$E/bin/R -e 'remotes::install_github("seasoncloud/Clonalscope@07dd0c576b6754d18dda86ec8905d8a9be5e77b0", upgrade="never")'

# NB: data-raw/ at the same commit, into references.clonalscope_data_raw
D=$(get 'c["references"]["clonalscope_data_raw"]' configs/metadata.yaml)
git clone https://github.com/seasoncloud/Clonalscope /n/fs/scratch/rl6004/Clonalscope
git -C /n/fs/scratch/rl6004/Clonalscope checkout 07dd0c576b6754d18dda86ec8905d8a9be5e77b0
mkdir -p $D && cp /n/fs/scratch/rl6004/Clonalscope/data-raw/{sizes.cellranger-GRCh38-1.0.0.txt,cyclegenes.rds,hg38.genes.bed,cytoarm_table_hg38.txt} $D/
```

Run; finished samples are skipped, `-n` is a dry run, `CORES` sets parallel samples:
```sh
bash scripts/run_baseline.sh clonalscope ovarian        # H1-N3KXRGQ 16um and 32um
bash scripts/run_baseline.sh clonalscope calicost-sim   # 10x inputs from h5ad, then every simulation
```
Outputs: `<outdir>/<sample>/visium_tutorial/{raw_output,postprocess}/`; logs `logs/clonalscope_{ovarian,calicost-sim}.*`.

> [!NOTE]
> `BayesNonparCluster` is a single-threaded Gibbs sampler, so each sample is one single-core job.

## Numbat

- Method: https://github.com/kharchenkolab/numbat
- Commit: CRAN 1.5.2 via bioconda `r-numbat`; GitHub has no `v1.5.2` tag, the commit setting `Version: 1.5.2` is `54d8268f3d14bd6afaa6784eeb839db7479a8947`
- Tutorials: [numbat vignette](https://kharchenkolab.github.io/numbat/articles/numbat.html) (scRNA mode, `ref_hca`), [spatial vignette](https://kharchenkolab.github.io/numbat/articles/spatial-rna.html) (`max_entropy = 0.8` for Visium)

Stages: preprocess (ovarian, `numbat/preprocess_ovarian.smk`: 16um counts and UGP pileup summed to 32um bins, allele counts from UGP phasing via `numbat:::preprocess_allele`; simulations, `numbat/preprocess_sim.smk`: the deposit's `numbat_allele_counts.tsv.gz`), then `numbat/run.smk` (`run_numbat()`), then `numbat/postprocess.smk` (`cells.tsv`, `cluster_states.tsv`). Non-default arguments and their sources: `configs/nondefault_params.yaml` (`numbat`).

```sh
$CONDA/bin/mamba env create -p $(get 'c["software"]["envs"]["numbat"]' configs/metadata.yaml) -f envs/numbat.yaml
```

Run; finished samples are skipped, `-n` is a dry run, `CORES` sets snakemake cores:
```sh
bash scripts/run_baseline.sh numbat ovarian        # H1-N3KXRGQ at 32um
bash scripts/run_baseline.sh numbat calicost-sim   # 10x inputs from h5ad, then every simulation
```
Outputs: `<outdir>/<sample>/{prep,max_entropy_08/run_numbat}/`, `<outdir>/postprocess/max_entropy_08/<sample>/`; logs `logs/numbat_{ovarian,calicost-sim}.*`.

> [!IMPORTANT]
> Ovarian runs on 32um bins: `run_numbat()` calls `stats::hclust` on a spot x spot distance, which stops above 65536 items, and the 16um slide exceeds that.

## CalicoST

- Method: https://github.com/raphael-group/CalicoST
- Commit: `c1abcae3e3657e01e547ee4529e3b9d039221453` (main), cloned at `software.calicost_repo` of `configs/metadata.yaml`
- Tutorial: [CalicoST docs](https://calicost.readthedocs.io/) and README command `OMP_NUM_THREADS=1 python <CalicoST>/src/calicost/calicost_main.py -c configuration_cna`

Ovarian (`calicost/ovarian.sh`): UGP `bb-<XXX>um-formatted` dirs are staged (`calicost/stage_input.sh`; `unique_snp_ids.npy` re-saved for numpy 1.24, CalicoST's pin), then each run of `configs/calicost.ovarian.yaml` gets its `configuration_cna` (`calicost/ovarian_config.py`: CalicoST's default configuration with paths filled; `msu<N>` sets `min_snpumi_perspot = N`) and runs under a cgroup `MemoryMax` (`calicost/run.sh`).
Simulations (`calicost/sim.sh`): the deposit's own `configfile0..4` per simulation, only input, output and resource paths remapped; `parse_input.py` once, then the five initializations in parallel (`calicost/repro.sh`).

```sh
C=$(get 'c["software"]["calicost_repo"]' configs/metadata.yaml)
E=$(get 'c["software"]["envs"]["calicost"]' configs/metadata.yaml)
git clone https://github.com/raphael-group/CalicoST $C
git -C $C checkout c1abcae3e3657e01e547ee4529e3b9d039221453
$CONDA/bin/mamba env create -p $E -f envs/calicost.yaml
$E/bin/pip install -e $C
```

Run; finished runs are skipped:
```sh
bash scripts/run_baseline.sh calicost ovarian        # stage inputs, then each run in order
bash scripts/run_baseline.sh calicost calicost-sim   # max_jobs simulations at a time
```
Outputs: ovarian `<outdir>/<tag>/`, logs `logs/calicost_<tag>.*` (with its `configuration_cna`); simulations `<out_root>/<sim>/clone3_rectangle<r>_w1.0/`, logs `<out_root>/logs/`.

> [!WARNING]
> `calicost/run.sh` needs systemd user scopes with the memory controller; a run over `MemoryMax` is killed (status `rc=137 (killed: OOM ...)`).

## Evaluation

Each script takes its paths as arguments; its header gives the usage.

| Script | Output |
| --- | --- |
| `vishd_runtime.py` | runtime and peak memory of every ovarian run |
| `sim_runtime.py` | runtime and peak memory per simulation, with the deposit's Numbat runtime |
| `eval_sim_ari.py` | clone ARI of Numbat and Clonalscope per simulation |
| `compare_clone_labels.py` | ARI of a coarse-bin labelling against a 16um labelling |
| `extract_deposit_results.py` | the deposit's CalicoST, Numbat and STARCH label and profile files |
| `calicost/compare_repro.py` | rerun CalicoST against the deposit's CalicoST, per initialization |
| `numbat/run_stats.py` | Numbat input sizes, step timings and resources per simulation |
| `status_report.sh` | one-screen status of every run in `logs/` |

> [!NOTE]
> Simulation 10x inputs (`inputs/<sim>/filtered_feature_bc_matrix/`) are shared by Clonalscope and Numbat; Ensembl IDs come from Cell Ranger GRCh38-2020-A `genes/ordered_genes.tsv`, which Clonalscope needs to join its gene BED.
