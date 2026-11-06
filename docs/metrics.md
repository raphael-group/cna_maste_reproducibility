# Simulation metrics

Each method's run of a simulation is scored against the simulation's truth by `cnamaste.qa.metrics`, the metrics
`run_sim_qa` tables for cnamaste and CalicoST, plus a CNA recall for methods that call CNA states instead of
allele copies. `scripts/run_benchmark.py`'s simulation step (`write_simulation_metrics`) writes the table, and its
figure step draws it in `combined.pdf` and `appendix.pdf` (see [Run](#run)). A sample's runs are compared with one
another by their clone profiles, the [CNP similarity](#cnp-similarity-real-data).

## Cells

Every count metric counts **cells**: one spot on one segment of the run's own copy-number table
(`qa.metrics`, module docstring).

- A segment is one row of the run's table, and counts once whatever its length.
- A segment's true copies, per true clone, are those of the true segment overlapping it most
  (`segment_true_ab_pairs`). A segment that no true segment overlaps is not scored.
- A spot's true call on a segment is its true clone's there. Its called state is its fitted clone's there, so every
  labelled spot is scored and no clone is matched to another.
- An AB pair is a segment's allele copies `(A, B)`, compared phase-free: `(A, B)` and `(B, A)` are the same
  (`phase_free`).

> [!IMPORTANT]
> Scoring on the run's own segments means a true CNA smaller than the run's segment around it drops out of the truth:
> the segment takes the copies of the true segment overlapping it most. A method with coarse segments, such as
> Numbat's whole-chromosome segments outside its calls or Clonalscope's arm segments, is then scored on fewer and
> larger CNAs. Where no CNA remains in any scored segment, sensitivity and CNA recall are NaN.

> [!NOTE]
> A spot the run does not label is left out of every metric. `#called` is the number of spots the run labels.

## Metrics

| Column | Panel | Definition | Source |
| --- | --- | --- | --- |
| `sensitivity` | CNA sensitivity | 1 - (truly altered cells called (1, 1)) / (truly altered cells) | `qa.metrics.sensitivity` |
| `specificity` | Normal specificity | 1 - (truly (1, 1) cells called altered) / (truly (1, 1) cells) | `qa.metrics.specificity` |
| `exact_altered` | N-CNA recall | (truly altered cells called their true AB pair) / (truly altered cells) | `qa.metrics.exact_altered` |
| `state_altered` | CNA recall | (truly altered cells called their true CNA state) / (truly altered cells) | `metrics/supervised_metrics.get_state_metrics` |
| `clone_ari` | Clone ARI | Adjusted Rand index of the fitted clones against the true ones, over the labelled spots | `qa.metrics.clone_ari` |
| `int_clone_ari` | (N-) Clone ARI | As `clone_ari`, of cnamaste's integer clones (`clone_labels_int.tsv`) | `qa.metrics.int_clone_ari` |

Truly altered means a true pair other than (1, 1). Each rate is NaN where its denominator is 0. `run_sim_qa` draws
`int_clone_ari` for cnamaste and `clone_ari` for CalicoST (`qa.plot.PANELS`); Numbat and Clonalscope have no integer
merge, so the simulation table carries `clone_ari`.

### CNA states

A CNA state is one of Numbat's six `cnv_state` values, its `loh` written `cnloh`. A pair's state follows from whether it is balanced (its two
alleles' copies are equal) and how its total copies compare with 2:

|            | total 2   | above 2        | below 2    |
| ---------- | --------- | -------------- | ---------- |
| balanced   | neu (1:1) | bamp (2:2)     | bdel (0:0) |
| unbalanced | cnloh (2:0) | amp (2:1, 3:1) | del (1:0)  |

Each pair is named as in Numbat's file descriptions (https://kharchenkolab.github.io/numbat/articles/descriptions.html):

- "CNV state 1:1 (neutral)"
- "CNV state 2:0 (CNLoH)"
- "CNV state 1:0 (deletion)"
- "CNV state 2:1 (amplification)"
- "amplification states (2:1, 3:1)"
- "CNV state 2:2 (biallelic amplification)"
- "CNV state 0:0 (biallelic deletion)"

A state is `neu` exactly where its pair is (1, 1). The state-based sensitivity and specificity therefore count the
same cells as `qa.metrics`' (`tests/test_simulation_metrics.py`).

## Methods

| Method | Read by | Called per segment and clone | Metrics computed |
| --- | --- | --- | --- |
| cnamaste | `run_sim_qa` (`qa.metrics.read_estimate`) | AB pair | all |
| CalicoST | `read_method_input.read_calicost`: a run's `clone_labels.tsv`, `cnv_seglevel.tsv`, its unconstrained copies (`calicost.read.read_estimate`) | AB pair | all |
| Numbat | `read_method_input.read_numbat`: `clone_post_<i>.tsv`, `bulk_clones_final.tsv.gz` | state, `cnv_state` | all but `exact_altered` |
| Clonalscope | `read_method_input.read_clonalscope`: the Visium tutorial run's `cells.tsv`, `profile.tsv` | total copies `TCN`, as `del` (< 2), `neu` (2), `amp` (> 2) | all but `exact_altered` |

- **Numbat**
  - Numbat calls states, not copies, so `exact_altered` is NaN.
  - Its `amp` covers both 2:1 and 3:1, so its state can match a true (1, 3) where its pair cannot.
  - A Numbat run that stops with "No CNV remains after filtering ..." (`run_numbat`, R/main.R) writes no clones. It is
    scored as a call of none: every spot in clone 0, `neu` on every true segment.
- **Clonalscope**
  - Clonalscope's export gives total copies alone, so a true `cnloh`, `bamp` or `bdel` cell is never its state.
  - A true `cnloh` cell, total 2, counts as missed.
  - Its total copies are read as copytyping-benchmark's `analysis/scripts/cn_transform.py` (1204da9) reads them:
    "DEL (-1) = total CN < 2, NEU (0) = total CN 2, AMP (1) = total CN > 2".
  - Its clusters are kept as called. Its normal clusters (`malignant` N) come first, and it may call several.
  - A Clonalscope run that halts writes no export; its row is NaN but for its runtime.

## Table

`output/tsvs/calicost-sim.tsv.gz`, one row per simulation and method, each run scored on its own segments and at the
genes, a metric's column per level, `<metric> (segment)` and `<metric> (gene)` (`simulation_metrics.merge_levels`):

| Column | Meaning |
| --- | --- |
| `sample_id` | The simulation, `numcnas<a.b>_cnasize<s>_ploidy<p>_random<r>` (`plots.simulation.SAMPLE_ID`). |
| `method` | A method of the `simulation` entry: `numbat`, `clonalscope`, `calicost`. |
| `#spots` | The truth's spots. |
| `#called` | The spots the run labels. |
| `spot-coverage-rate` | `#called` / `#spots`: the truth's spots the run labels with any clone, normal or tumor; one value at both levels. |
| `sensitivity (<level>)` ... `clone_ari (<level>)` | As above, on the run's own segments, `segment`, or at the genes, `gene`. |
| `Runtime [mins]`, `Memory [GB]` | From `input_data/simulation_data/simulation.runtime.tsv`: the run rule's Snakemake benchmark (`run_numbat`, Clonalscope's `run_tutorial`), its `s` / 60 and `max_pss` / 1024. |

The figures, pages of `output/plots/combined.pdf` and `output/plots/appendix.pdf`
(`run_benchmark.py --simulation_dataset calicost-sim`), draw each column in its own panel, a box per simulated class
and method (`scripts/plots/simulation.py`): a page of the metrics at the genes (`simulation_metrics.get_level_table`),
combined.pdf's page 1, two rows of two panels, `sensitivity` and `clone_ari`, `state_altered` and `exact_altered`, the
classes (6,3, 50M) and (3,3, 50M) left out (`plots.pages.DROPPED_CLASSES`); then the appendix's page of
`Runtime [mins]` and `Memory [GB]` on log axes. `diagnosis.pdf` draws the metrics page in full, every class, three
rows of two panels, `specificity` and `spot-coverage-rate` first. A panel with no method's values reads "not
computed".

## Run

`run_benchmark.sh` runs the benchmark, this table its first step, `cnamaste` imported from `$CNAMASTE_ROOT` (README,
[Run](../README.md#run)):

```sh
CNAMASTE_ROOT=/path/to/cna_maste bash <benchmark>/scripts/run_benchmark.sh
```

Each run is scored on its own segments and at the genes, in one run; `--simulation_dataset` names the table (
[pipeline.md, Metrics](pipeline.md#metrics)).

The inputs are the `simulation` entry of `--benchmark_yaml`, `configs/benchmark.local.yaml` by default, under
`input_data/simulation_data/` (gitignored):

| Key | Contents |
| --- | --- |
| `groundtruth_dir` | `groundtruth/`, per simulation `truth_clone_labels.tsv`, `truth_acn_profile.tsv` |
| `runtime` | `simulation.runtime.tsv`, each run's runtime and memory |
| `numbat.dir` | `numbat/{sample_id}/max_entropy_08/run_numbat`, Numbat's output |
| `clonalscope.dir` | `clonalscope/postprocess/visium_tutorial/{sample_id}/scRNA`, Clonalscope's export |
| `calicost.dir` | `method_results/calicost/{sample_id}`, CalicoST's best rectangle by `total_llf` (`best_init.tsv` beside) |

## CNP similarity (real data)

The CNP similarity of two runs of a sample compares their clones' copy-number profiles (CNPs) at the genes, no spot
and no truth: of their most similar pair of tumor clones, the share of their abnormal region, the genes either clone
calls altered, where the two clones' CNA states match.

### Definition

Let $\mathcal{G}$ be the genes inside the autosomal regions ([Genes](pipeline.md#genes)), and $s^a_i(g)$ the CNA state
([CNA states](#cna-states)) of run $a$'s clone $i$ at gene $g$, the state of its segment overlapping $g$ most, or
$\varnothing$ where the run calls none there (`uncalled`). For a tumor clone $i \in K_a$ of run $a$ and $j \in K_b$ of
run $b$, the altered genes are those both call and at least one calls other than `neu`:

```math
\mathcal{A}_{ij} = \{\, g \in \mathcal{G} : s^a_i(g) \neq \varnothing,\ s^b_j(g) \neq \varnothing,\
\big(s^a_i(g), s^b_j(g)\big) \neq (\texttt{neu}, \texttt{neu}) \,\}
```

the two clones' agreement is the share of the altered genes where their states match:

```math
J_{ij} = \frac{1}{|\mathcal{A}_{ij}|} \sum_{g \in \mathcal{A}_{ij}} \mathbb{1}\left[ s^a_i(g) = s^b_j(g) \right]
```

and the runs' CNP similarity is that of their most similar pair of tumor clones:

```math
S(a, b) = \max_{i \in K_a,\ j \in K_b} J_{ij}
```

- A gene both call `neu` is left out of the numerator and the denominator alike; a gene one calls `neu` and the
  other altered counts as a mismatch.
- $K_a$ is run $a$'s tumor clones: its clones but its normal ones, clone 0 and Clonalscope's normal clusters
  (`real_data_metrics.get_normal_clones`, `get_tumor_cnps`): the CNPs compared are the tumor's.
- $J_{ij}$ is NaN where $\mathcal{A}_{ij}$ is empty, and the maximum skips it; $S(a, b)$ is NaN where every
  $J_{ij}$ is, a run with no tumor clone among them. A tie goes to the first pair in clone order.
- $S(a, b) = S(b, a)$, and $S(a, a) = 1$ where some tumor clone of $a$ alters a gene.
- `unsupervised_metrics.get_altered_genes` is $\mathcal{A}_{ij}$, `compute_cnp_similarity` $J_{ij}$, and
  `compute_pairwise_cnp_similarity` $S(a, b)$ for every two runs.

> [!NOTE]
> $J_{ij}$ is the Jaccard index ([Jaccard 1912](https://doi.org/10.1111/j.1469-8137.1912.tb05611.x)) of the two
> clones' altered genes among those both call, their intersection over their union, a gene both alter counted in the
> intersection only where their states match.

> [!IMPORTANT]
> Clonalscope's CNA state is that of its total copies, `del`, `neu` or `amp` ([Methods](#methods)), so a gene another
> run calls `cnloh`, `bamp` or `bdel` never matches it.

### Table and figure

`run_benchmark.py`'s real-data step (`write_real_data_metrics`) writes `output/tsvs/<sample>.cnp_similarity.tsv.gz`, a
row per ordered pair of runs:

| Column | Meaning |
| --- | --- |
| `run_a`, `run_b` | The two runs, `--run` keys, in `--run` order. |
| `clone_a`, `clone_b` | The pair $(i, j)$ the maximum picks, NA where $S$ is NaN. |
| `#altered-genes` | $\lvert\mathcal{A}_{ij}\rvert$ of that pair. |
| `cnp_similarity` | $S(a, b)$. |

Its figure step draws it as `combined.pdf` page 3 (c), beside the clone ARI's
(`plots.run_comparison.get_run_comparison_panels`): a heatmap, a row and a column per run in method order, $S(a, b)$
written in its cell, `-` where NaN, its colour from white at 0 to blue at 1. The run commands are
[pipeline.md's Steps](pipeline.md#steps).
