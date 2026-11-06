# Pipeline: the steps, the input readers and the metrics

How `run_benchmark.py` scores the simulations, compares a sample's runs and writes the LaTeX tables; the figures' pages are [figures.md](figures.md), the simulation metrics and the CNP similarity [metrics.md](metrics.md).

## Simulations

`run_benchmark.py`'s first step (`write_simulation_metrics`) scores each method's runs of the simulations against their
truth by `cnamaste.qa.metrics`, the metrics `run_sim_qa` tables for cnamaste and CalicoST, plus a CNA recall
([Metrics](#metrics)), each run on its own segments and at the genes, and writes one wide table,
`tsvs/<--simulation_dataset>.tsv.gz`, `calicost-sim` by default: a row per simulation and method, `#spots`, `#called`,
`spot-coverage-rate`, each metric per level, `<metric> (segment)` and `<metric> (gene)`
(`simulation_metrics.merge_levels`), then its runtime and memory; combined.pdf's page 1 and the appendix draw it ([Appendix](figures.md#appendix)).
Beside it, `tsvs/<--simulation_dataset>.median.tsv.gz`: a row per simulated class and method, `sim_id` (the
simulation's id less `_random<r>`), `method`, `#replicates`, then each metric's median over the class' replicates, NaN
skipped (`simulation.get_class_medians`).
The runs are the `simulation` entry of
`--benchmark_yaml` (`configs/benchmark.local.yaml`): `groundtruth_dir`, `runtime`, and per method a `dir` template of
`{sample_id}`, Numbat, Clonalscope and CalicoST in `input_data/simulation_data/`. cnamaste's are
`nomixing_cnamaste_related/`, as its `batch_zenodo_run.sh` writes them, a directory per solver,
`<label_solver>-<initializer>-<decoder>`, holding its `run_<id>`: the `dir` ends `run_*`, the one run it matches
(`simulation_metrics.get_run_dir`). `cnamaste` is the solver `expansion-lattice-lattice`, and `cnamaste-<solver>` each
other solver, each with its `solver`, the table's `solver` column; a key's method is the key before its first `_` or `-`
(`read_method_input.get_run_method`), so every solver is read as cnamaste. cnamaste's runtime in `runtime` is each run's
own, its `cnamaste_<id>.log`'s "Done in"; its memory is not recorded. [metrics.md](metrics.md) states how each
method's output is read.

## Metrics

`scripts/` holds one entry script, `run_benchmark.py`, and `run_benchmark.sh`, which runs it on H1-N3KXRGQ
([README, Run](../README.md#run)): its steps write the metric tables, scoring the simulations' runs against their truth
(`write_simulation_metrics`) and comparing a sample's runs with one another (`write_real_data_metrics`), then the
LaTeX tables, then the three PDFs from those tables (`write_figures`). Its components are `scripts/plots/`, the
figures, `scripts/metrics/`, the metrics, and `scripts/platforms/`, the Visium HD bins:

| Module | Holds |
| --- | --- |
| `plots/style.py` | A method's colour and fill, a run's name and order, the keys, the panel letters, a page of rows of panels laid out tight. |
| `plots/simulation.py` | The simulation pages: per simulated class, a box per method or cnamaste solver; `qa.plot`'s box helpers copied. |
| `plots/run_comparison.py` | A sample's runs compared: cells per altered state, state Jaccard index, every two runs' clone ARI and CNP similarity. |
| `plots/cn_error.py` | The copy-number error page: mean squared error and skill scores. |
| `plots/umi.py` | Per-spot UMI histograms, and the false-positive tumor spots' UMI against the other called spots'. |
| `plots/cnamaste_page.py` | cnamaste's own page, and its inputs: the fit on the integer clones, the per-spot bins, the normal fractions. |
| `plots/pages.py` | A method's title, a profile on the cnamaste run's bins, the metric pages, and each figure written to its PDF and each page alone. |
| `plots/layout.py` | The spatial pages, H&E, cell types and each method's clones, and the profile page on cnplot's genome axis; the panel letters, titles and text sizes. |

| Module | Holds |
| --- | --- |
| `read_method_input.py` | The readers: each method's run, the spots, their cell types, and UGP's counts per bin and spot. |
| `outputs.py` | The outputs' paths under `--outdir`, and the table and summary writers. |
| `cna_states.py` | The CNA states, and each method's fields per segment and clone. |
| `gene_level.py` | The genes inside the region bed, and each run's profile at them. |
| `segment_level.py` | The runs' joint segments, the union of their breakpoints, and each run's profile on them. |
| `supervised_metrics.py` | A run against the truth; `qa.metrics`' functions copied, on a run's labels and profile. |
| `unsupervised_metrics.py` | Two runs against each other. |
| `qc_metrics.py` | A run's genome by CNA state; its spot coverage and false-positive counts. |
| `cn_error.py` | A run's copies against its clones' pseudobulk RDR and BAF. |
| `simulation_metrics.py` | The simulations' runs scored against their truth, with their runtime, both levels merged into one wide table (`simulation`). |
| `real_data_metrics.py` | A sample's runs compared with the reference run, with their runtime, merged into one wide table, and its summary (`real_data`). |
| `latex_tables.py` | The real-data wide table's LaTeX table, and its PDF. |
| `pseudobulk.py` | Each clone's pseudobulk RDR and BAF per bin beside the RDR and BAF its decoded `(A, B)` is expected to read; no module imports it. |

| Module | Holds |
| --- | --- |
| `platforms/visium.py` | Visium HD square bins: barcodes read as grid places, placed on the image, given to finer bins. |

### CNA states and fields

A CNA state is one of Numbat's six, its `loh` written `cnloh`, by an AB pair's balance and total copies, each pair as
Numbat's [file descriptions](https://kharchenkolab.github.io/numbat/articles/descriptions.html) name it:

|            | total 2     | above 2        | below 2    |
| ---------- | ----------- | -------------- | ---------- |
| balanced   | neu (1:1)   | bamp (2:2)     | bdel (0:0) |
| unbalanced | cnloh (2:0) | amp (2:1, 3:1) | del (1:0)  |

Each method gives, per clone, some of five fields:

| Method | `ascn` | `Uest` | `tcn` | `state` | `tcn_state` |
| --- | --- | --- | --- | --- | --- |
| truth, cnamaste, CalicoST | `A\|B`, as the table phases it | | `A + B` | of the phase-free pair | of `A + B` |
| Numbat | the state's pair: `1\|1`, `1\|0`, `2\|0`, `0\|0`; `unknown` for `amp`, `bamp` | | its total, or `unknown` | `cnv_state` | of the state |
| Clonalscope | | the cluster's mean coverage fold change, `AssignCluster`'s `Uest` (`REL_CN`) | `TCN`, `round(2 * Uest)` | `del` (< 2), `neu` (2), `amp` (> 2) | its `state` |

A field is `uncalled` where the run calls nothing there. `tcn_state` is the state of total copies alone, `del`
(< 2), `neu` (2) or `amp` (> 2), as copytyping-benchmark's `analysis/scripts/cn_transform.py` (1204da9) reads
total copies: `cnloh` is `neu`, `bdel` `del`, `bamp` `amp`.

### Genes

With `--config_yaml` (`configs/config.local.yaml` by default) naming `reference.bed_file`, the region bed, and
`reference.refdata.gene_list.gtf_file`, the reference's genes:

- The genes are the GTF's `gene` features wholly inside a region of the bed on chromosomes 1 to 22, where the truth and
  every method's input lie; `START` 0-based, `END` exclusive, as BED's.
- A gene takes, of a run's segments and of the truth's, the one that overlaps it most, the first on a tie; a gene no
  segment overlaps is `uncalled`.
- `<sample>.gene_level_cnp.tsv.gz`, a row per gene, in `tsvs/`, a simulation's in `tsvs/gene_level_cnp/`: `gene_id`, `gene`, `CHR`, `START`, `END`; then `truth clone<k> ascn` and
  `truth clone<k> state` per true clone, `normal` clone 0, on simulations; then `<method> clone<k> <field>` per run,
  clone and field above.

### Joint segments

The runs' joint segments are the intervals between consecutive breakpoints of any run, kept where some run has a segment
(`scripts/metrics/segment_level.py`), so each lies inside at most one segment of each run, as copytyping-benchmark's
`analysis/scripts/cn_error.py` (1204da9) joins competing runs. `write_real_data_metrics` joins the runs whose
profile gives copies (`metrics.cna_states.get_profile_copies`): CalicoST's and cnamaste's AB pairs and Clonalscope's
total copies; Numbat's states give none. `tsvs/<sample>.segment_level_cnp.tsv.gz`, a row per joint segment: `CHR`,
`START`, `END`; then `<method> clone<k> <field>` per such run, clone and field, each joint segment the run's segment
overlapping it most, `uncalled` where none does.

### Against the truth (simulations)

A cell is one spot a run labels on one row of its profile that a true segment overlaps: one of its own segments
(`qa.metrics`), its `(segment)` columns, or one gene, its `(gene)` columns. Each row counts once, whatever its
length. A cell's true AB pair is its true clone's, its called AB pair or state its fitted clone's, so no clone is
matched to another. Rows the run calls no state on in some clone are not scored.

| Column | Definition |
| --- | --- |
| `sensitivity` | 1 - (truly altered cells called `neu`) / (truly altered cells) |
| `specificity` | 1 - (truly `neu` cells called altered) / (truly `neu` cells) |
| `exact_altered` | (truly altered cells called their true AB pair, phase-free) / (truly altered cells); NaN for a run of states |
| `state_altered` | (truly altered cells called their true CNA state) / (truly altered cells) |
| `clone_ari` | Adjusted Rand index of the run's clones against the true ones, over the spots it labels |

Truly altered is a true pair other than (1, 1). A rate is NaN where its denominator is 0.

> [!NOTE]
> `sensitivity`, `specificity` and `state_altered` are read from one table of cells by `(true CNA state, called CNA
> state)`, the six states of [CNA states and fields](#cna-states-and-fields)
> (`supervised_metrics.get_state_metrics`). `state_altered` counts a truly altered cell only where the two states agree, so `amp` called `cnloh` is wrong;
> `sensitivity` and `specificity` use `neu` against the other five alone, and as a state is `neu` exactly where its
> pair is (1, 1), they count the cells `qa.metrics` counts on AB pairs. The true state is the true pair's. The called
> state is the called pair's (cnamaste, CalicoST), the `cnv_state` (Numbat), or that of the total copies (Clonalscope,
> `del`, `neu` or `amp`). Clonalscope's call on a true `cnloh` is `neu`, so the cell is missed; on a true `bamp` it is
> `amp`, detected but not its state; on a true `bdel` it is `del`, likewise. Its `state_altered` is therefore below 1
> wherever the truth holds `cnloh`, `bamp` or `bdel`.

> [!IMPORTANT]
> On a run's own segments, a true CNA smaller than the segment around it drops out: the segment takes the copies of the
> true segment overlapping it most. At the genes, every true CNA a gene lies in is scored.

### Between runs (real data)

On the in-tissue spots at `--side_um` both runs label, a coarser bin's clone given to every spot inside it, at the
finest field both runs give, their `level` (`cna_states.get_comparison_field`): `ascn`, AB pairs phase-free as
`(major,minor)`, where both give pairs (cnamaste, CalicoST); `tcn_state` where either gives total copies alone
(Clonalscope); else `state` (Numbat):

| Function | Definition |
| --- | --- |
| `compute_cna_confusion` | Cells by `(state of A, state of B)` at the `level`, `uncalled` included: per shared spot, the genes at each pair of its two clones' states, summed over the spots. |
| `compute_clone_ARI` | Adjusted Rand index of the two runs' clones over the shared spots. |
| `compute_pairwise_clone_ARI` | `compute_clone_ARI` of every two runs, each pair over its own shared spots, a run against itself 1. |
| `compute_cna_ARI` | Adjusted Rand index of the two runs' states over the cells both call a state on: that of the confusion as a contingency table, `uncalled` left out ([Hubert and Arabie 1985](https://doi.org/10.1007/BF01908075), as scikit-learn's `adjusted_rand_score` computes it). |
| `compute_state_jaccard` | Jaccard index of each two states' cells, one state per run, over the cells both call: of A's state `s` and B's state `t`, the cells where A calls `s` and B calls `t` over those where A calls `s` or B calls `t` ([Jaccard 1912](https://doi.org/10.1111/j.1469-8137.1912.tb05611.x)); NaN where neither state holds such a cell. Unlike the ARI, it compares states by name: a run calling each `del` of the other `cnloh` scores 0 at `(del, del)`. |

The wide table holds, per run, its `reference`, `--reference` (cnamaste by default), its `level`, `#shared-spots`,
`clone_ari` and `cna_ari` ([Wide table and summary](#wide-table-and-summary-real-data));
`tsvs/<sample>.cna_confusion.tsv.gz` each run's confusion at its `level`, long, and `tsvs/<sample>.clone_ari.tsv.gz`
every ordered pair of runs, `run_a`, `run_b`, `#shared-spots` and `clone_ari`.
Clonalscope's clusters are compared as called, each normal cluster its own clone, as on the simulations.

The CNP similarity compares two runs' tumor clones' profiles alone, at the genes, no spot, at their CNA `state`: of
their most similar pair of tumor clones, the share of the genes both call and at least one calls other than `neu`
where the two states match (`compute_pairwise_cnp_similarity`,
[CNP similarity](metrics.md#cnp-similarity-real-data)). `tsvs/<sample>.cnp_similarity.tsv.gz` holds every ordered
pair of runs, `run_a`, `run_b`, `clone_a`, `clone_b`, `#altered-genes` and `cnp_similarity`.

### Genome composition

`compute_genome_composition`: per clone, the share of the genome in each CNA state and `uncalled`, summing to 1. The
genome is the regions' bases, a segment weighing its bases inside them, or the genes, 1 each. Its `all` row weighs each clone by its spots: the share of the run's cells, a labelled spot at a gene, in each
state. `tsvs/<sample>.genome_composition.tsv.gz` holds it per run and clone.

### Copy-number error

Whether a run's copies reconstruct the counts of its clones, as copytyping-benchmark's `analysis/scripts/cn_error.py`
(1204da9) scores them, with copytyping's pseudobulks (e0341e4, `inference/model_utils.py`), in
`scripts/metrics/cn_error.py`. It scores the runs whose profile gives copies: CalicoST's and cnamaste's AB pairs, RDR
and BAF; Clonalscope's total copies, `round(2 * REL_CN)`, RDR alone, its BAF `N/A`. Numbat's states give no copies.

**Bins.** RDR and BAF are scored on UGP's bins at `--side_um`, `bb/MSR<--msr>/<--assay>/`, `--msr` 300 by default
(`read_method_input.read_bin_counts`). UGP bins the SNPs adaptively, each bin at least `min_snp_reads` SNP reads, and
sums into each bin its SNPs' phased A and B counts and the UMI of each gene it overlaps most; a bin with no SNP has no
BAF. The bins inside an autosomal region of the region bed are kept.

**Copies per bin.** A bin takes the joint segment overlapping it most ([Joint segments](#joint-segments)) and, from
each run, the copies of the segment holding that joint segment, per clone `k`: `A`, `B` and `C = A + B`; uncalled where
the run has no segment there.

**Clone pseudobulks.** A clone's counts are summed over its spots, the run's clones on the in-tissue spots at
`--side_um`, a coarser bin's clone given to the spots inside it:

    X[k,u] = sum_{n in k} X[u,n]      the clone's UMI in bin u
    T[k]   = sum_u X[k,u]             its UMI over the bins
    b[k,u], d[k,u]                    its B-allele and total allele counts in bin u
    lambda_u = sum_{n in ref} X[u,n] / sum_{n in ref} T_n

`lambda_u`, the baseline, is bin `u`'s share of the UMI of the spots of a normal cell type, neither `Epithelial` nor
`Low confidence` (`read_method_input.get_normal_type`); it sums to 1.

**Observed and expected.**

    RDR_obs[k,u] = X[k,u] / (lambda_u T[k])          RDR_exp[k,u] = C[u,k] / sum_u' lambda_u' C[u',k]
    BAF_obs[k,u] = b[k,u] / d[k,u]                   BAF_exp[k,u] = B[u,k] / (A[u,k] + B[u,k])
    pi[k,u] = lambda_u RDR_exp[k,u]                  the bin's expected share of the clone's UMI

`RDR_exp`'s sum runs over the bins the run calls for the clone. The diploid prediction, every clone 1:1, gives
`RDR_exp` 1, as `lambda` sums to 1, `pi = lambda_u`, and `BAF_exp` 0.5.

**Phasing.** UGP's counts are phased, a bin's B allele one haplotype, which may switch between bins; a run's `A|B`
carries its own phase. A bin is read in the orientation, `BAF_exp` or `1 - BAF_exp`, with the smaller squared error
summed over the run's clones, each weighed by its spots: one orientation for every clone, as a bin has one phase
(`cn_error.get_baf_flips`, `#flipped-bins`). Both BAF metrics read the bin in that orientation.

**Bins scored.** A bin is scored where every run competing on the channel calls every clone and every observed and
expected value is finite, and where some run calls a CNA in some clone with spots: its copies not 1:1, or, given total
copies alone, not 2 (`cn_error.get_altered`). A bin every run calls neutral throughout, where every run's prediction is
the diploid one, is not scored. The bins scored therefore depend on the runs compared.

**Metrics**, per run, on the bins scored, with `n_k` the clone's spots:

| Column | Definition |
| --- | --- |
| `mse_rdr` | `sum_{k,u} n_k (log2 RDR_obs - log2 RDR_exp)^2 / sum_{k,u} n_k`, every clone. |
| `mse_baf` | `sum_{k,u} n_k (BAF_obs - BAF_exp)^2 / sum_{k,u} n_k`, `BAF_exp` in the bin's orientation, every clone. |
| `null_rdr`, `null_baf` | The same error of the diploid prediction: `log2 RDR_exp` 0, `BAF_exp` 0.5. |
| `ss_rdr`, `ss_baf` | Skill score, `1 - mse / null`: 1 a perfect reconstruction, 0 the diploid prediction's, below 0 worse. |

MSE weighs a clone by its spots, the normal clone among them. `#tumor-clones` counts a run's clones but its normal
ones, the clones of the spots its reader labels clone 0 with the normal clusters merged: clone 0, and Clonalscope's
normal clusters (`real_data_metrics.get_normal_clones`).

> [!NOTE]
> `RDR_exp` is a pure tumor's: a clone's spots that hold normal cells read an RDR nearer 1 than its copies give, which
> raises `mse_rdr` and lowers `ss_rdr`. No purity is modelled, as copytyping-benchmark's `cn_error.py:719-721` calls
> `pseudobulk_rdr_matrix` without its `purity`.

**Outputs.** In the wide table, per run giving copies, `#clones`, `#tumor-clones`, then per channel its bins scored
(`#rdr-bins`; `#baf-bins`, `#flipped-bins`) and the columns above; its spots are its `#called`. With `--sweep_run
<key>`, once per run, `tsvs/<sample>.cn_error.sweep.tsv.gz` holds the same columns, its spots `#spots`, for the
reference run and the sweep runs, scored as one panel on the bins they share; with `--calicost_sweep`, the figure step
draws it in the diagnosis, each sweep run a shade of its method's colour.

### Wide table and summary (real data)

`write_real_data_metrics` computes every table first, then merges the per-run ones into one wide table, a row per
`--run` in order, each column once (`real_data_metrics.merge_run_tables`): a column two tables hold, e.g. `#called`,
must agree on every run, a count, `#<name>`, is an integer, NA where a table has no row of the run. It writes the
detail tables, `tsvs/<sample>.<name>.tsv.gz`, the wide table, `tsvs/<--real_data_dataset>.tsv.gz`, `ovarian` by
default, and its summary, `tsvs/<--real_data_dataset>.json`, only then. The wide table's columns, in order:

| Columns | From |
| --- | --- |
| `sample`, `run`, `method`, `resolution` | The run; `resolution` the side of the bins it ran on, `<side>um`. |
| `#in-tissue`, `#called`, `spot-coverage-rate`, `#normal-spots`, `#tumor-spots` | Its tumor calls on the in-tissue spots at `--side_um`, Clonalscope's normal clusters one clone 0 (`real_data_metrics.get_spot_table`, `qc_metrics.get_coverage`): the spots it labels, `#called` / `#in-tissue`, those of clone 0 and of another clone. |
| `#normal-spots-<type>`, `#tumor-spots-<type>`, `#false-positive-tumor-spots`, `#normal-type-spots`, `false-positive-rate` | Its false-positive counts (`qc_metrics.get_false_positive_counts`): its calls on the spots of each broad cell type, `Epithelial` first, `Low_confidence` last; its tumor calls on a type neither `Epithelial` nor `Low confidence`; its calls, normal or tumor, on such a type; the former over the latter. |
| `neu` ... `uncalled` | Its genome composition's `all` row ([Genome composition](#genome-composition)). |
| `level`, `#shared-spots`, `clone_ari`, `cna_ari`, `reference` | Against the reference run ([Between runs](#between-runs-real-data)), NA for the reference. |
| `#clones` ... `ss_baf` | Its copy-number error, the runs giving copies ([Copy-number error](#copy-number-error)). |
| `setting`, `Runtime [mins]`, `Memory [GB]` | Its row of the sample's `runtime` table (`real_data_metrics.get_runtime_table`). |

The false-positive rate is FP / N, the spots of a normal type the negatives, N = TN + FP (James et al., *An
Introduction to Statistical Learning*, 1st ed., sec. 4.4.3, Table 4.7, https://doi.org/10.1007/978-1-4614-7138-7); an
epithelial, low-confidence or untyped spot is in neither count. A run's runtime row is its method's, at the side of the bins it labels, its `setting` a word of the run
directory's name or its parent's; a run that matches no row, or several, stops the script.

The summary, `tsvs/<--real_data_dataset>.json` (`real_data_metrics.get_summary`), holds `dataset`, `sample` and
`reference`, then per run its `run`, `method` and `resolution`, its `metrics` and their `provenance`, NaN or NA written
`null`:

- `metrics` (`real_data_metrics.SUMMARY_METRICS`): `spot-coverage-rate`, `false-positive-rate`, `clone_ari`,
  `cna_ari`, the shares `neu` ... `uncalled`, `mse` and `ss` of RDR and BAF, `Runtime [mins]`, `Memory [GB]`.
- `provenance` (`real_data_metrics.PROVENANCE`): `setting`, `#in-tissue`, `#called`, `#normal-spots`,
  `#tumor-spots`, `#false-positive-tumor-spots`, `#normal-type-spots`, `level`, `#shared-spots`, `#clones`,
  `#tumor-clones`, `#rdr-bins`, `#baf-bins`, `#flipped-bins`.

The counts per cell type and the diploid prediction's `null_*` stay in the wide table alone; the per-gene,
per-segment, per-state and per-pair tables, and `cn_error.sweep`, are detail tables.

### Steps

`run_benchmark.py`'s two metric steps ([README, Run](../README.md#run)) both read `--benchmark_yaml`, `configs/benchmark.local.yaml` by
default, and `--config_yaml`'s region bed, its regions on chromosomes 1 to 22, and the genes inside them
([Genes](#genes)).

`write_simulation_metrics` (`simulation_metrics.get_simulation_metrics`), once per level, on the runs' own
segments, then at the genes, per simulation of `groundtruth_dir`, in name order:

1. Read the truth (`truth.read_truth`) and each `--method`'s run with its method's reader
   (`read_method_input.get_method_reader`), every key of the `simulation` entry with a `dir` by default. A Numbat run
   that stopped finding no CNV is scored as a call of none, every spot clone 0 and `neu` on every true segment; a run
   with no output is left unscored, NaN but for its runtime.
2. Score the run against the truth (`supervised_metrics.compute_supervised_metrics`) on its own segments, or at the
   genes, each gene its segment overlapping it most: the columns of [Against the truth](#against-the-truth-simulations),
   `#spots`, the truth's spots, `#called`, those the run labels, and `spot-coverage-rate`, `#called` / `#spots`.
3. At the genes, write `tsvs/gene_level_cnp/<sample_id>.gene_level_cnp.tsv.gz`, the truth's and each run's fields per
   gene, as each simulation is scored: a profile, not a metric table.

Each run takes its runtime and memory from `runtime` by `simulation_id` and `method`; the two levels are merged into
one wide table (`simulation_metrics.merge_levels`), refused where their runs, spots or runtime differ, and written
once, `tsvs/<--simulation_dataset>.tsv.gz`, `calicost-sim` by default.

`write_real_data_metrics`, for `--sample <sample>`:

1. Read the sample's in-tissue spots, `ugp.spots`, their cell types, `cell_type.labels`, and each `--run`, `numbat`,
   `clonalscope`, `calicost` and `cnamaste` by default, with its method's reader, the method the key before its first
   `_` or `-`, cnamaste's copies those of the clones its integer clones name: run_cnamaste merges clones whose (A,B)
   agree into an integer clone named by a member, whose copies it calls, and keeps every clone in its seglevel table
   (`read_method_input.get_labelled_profile`). A run's clones are given to the in-tissue spots at `--side_um`, a
   coarser bin's clone to each spot inside it (`real_data_metrics.get_run_spots`), and its profile to the genes
   (`gene_level.get_gene_states`).
2. `genome_composition`: per run and clone, and its `all` row, the genes' shares by CNA state and `uncalled`
   ([Genome composition](#genome-composition)).
3. Each run against `--reference`, cnamaste by default, at their `level`, the finest field both give: the confusion over
   the cells of the spots both label, `cna_confusion`, then clone ARI and CNA ARI
   ([Between runs](#between-runs-real-data)); `clone_ari`, the clone ARI of every two runs; `cnp_similarity`, the
   CNP similarity of every two runs, of their tumor clones.
4. `segment_level_cnp` and the copy-number error: the runs whose profile gives copies, on their joint segments
   ([Joint segments](#joint-segments)), scored at UGP's bins at `--side_um`, `MSR<--msr>/<--assay>/`, inside the
   regions, against the baseline of the spots of a normal cell type ([Copy-number error](#copy-number-error),
   `real_data_metrics.get_copy_number_tables`).
5. With `--sweep_run <key>`, once per run: each sweep run read as in step 1, and scored with the `--reference` run as
   one panel, on the bins they share, `cn_error.sweep`.
6. Each run's spot coverage and false-positive counts, and its runtime row; every per-run table merged into the wide table,
   and the tables, the wide table and its summary written ([Wide table and summary](#wide-table-and-summary-real-data)).

`write_figures` draws the tables in the appendix ([Appendix](figures.md#appendix)) and, with `--calicost_sweep`,
`cn_error.sweep` in the diagnosis.

## CalicoST clones missing from the unconstrained table

A CalicoST clone can be in `clone_labels.tsv` yet absent from `cnv_seglevel.tsv`, its unconstrained integer copies
(CalicoST c1abcae,
[`calicost_main.py`](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/calicost_main.py),
[`find_integer_copynumber.py`](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/find_integer_copynumber.py)):

- The labels are the final fit's `new_assignment` (`calicost_main.py:388`). Every clone merge,
  `similarity_components_rdrbaf_neymanpearson` and `merge_by_minspots`, runs before that fit is saved (`:299`), so no
  merge removes a labelled clone.
- The unconstrained copies (`hill_climbing_integer_copynumber_fixdiploid`, `:334`, `:337`) need the clone's diploid
  balanced state: a state holding at least 10% of its bins, 2% on the retry, its BAF within 0.05 of 0.5
  (`find_integer_copynumber.py:73-81`). Where neither attempt finds one, the bare `except` skips the clone
  (`calicost_main.py:338-340`) and logs nothing; the run's log lacks only its `max med ploidy = None, clone <s>` line
  (`:342`).
- The diploid, triploid and tetraploid tables (`hill_climbing_integer_copynumber_oneclone`, `:331`) need no such
  state, so they hold the clone.

A run's fit shows which: `rdrbaf_final_nstates<n>_smp.npz`'s `pred_cnv[:, s]` and `new_p_binom[:, s]` give the share
of clone `s`'s bins in each state and its BAF. `read_method_input.read_calicost` therefore reads a run's
`--calicost_ploidy` table where it holds every labelled tumour clone; else the first of `--calicost_fallback_ploidy`,
`diploid` then `tetraploid` by default, that holds them all; else the `--calicost_ploidy` table, logged
(`get_calicost_ploidy`). A run read at a fallback names it in its title, e.g.
`CalicoST (16 um, min-snp-umi: 10, diploid)`.

## LaTeX tables

`run_benchmark.py`'s third step (`metrics.latex_tables.write_latex_table`) writes the benchmark's LaTeX table from
the sample's wide table, `tsvs/<--real_data_dataset>.tsv.gz`, `ovarian` by default
([Wide table and summary](#wide-table-and-summary-real-data)), into `text/table.tex` (`--table <name>` for
`text/<name>.tex`), after a comment naming the sample: a booktabs tabular in a `table` float placed `[!h]`, centred, in
a `minipage` of `0.85\linewidth`, its caption and label, `tab:data_runtime`, under it (`latex_tables.CAPTION`,
`LABEL`). A row per run in method order, cnamaste, CalicoST, Numbat, Clonalscope (`plots.style.get_ordered_runs`, as
the figures order them), its columns (`latex_tables.get_table_columns`), each name over its unit or its second line,
every column but `Method` centred:

- `Method`, then `Resolution`, the side of its bins, `[$\mu$m]`;
- `Spot coverage rate`, `spot-coverage-rate`, #called / #in-tissue at `--side_um`, a coarser run's bin clone given to
  each spot inside it, and `False-positive rate`, `false-positive-rate`, each to 0.001;
- `Runtime [mins]` to 0.1 and `Memory [GB]` to 0.01, `$-$` where not recorded.

The including document needs `\usepackage{booktabs}`. Beside it, `text/table.pdf` typesets the float on letter pages
by `pdflatex` (`latex_tables.get_article_tex`), no date or ID in the PDF, so equal tables give equal bytes; pdflatex's
log is `text/table.log`.

The same step writes the metric text, `text/metric.tex` (`--metric_text <name>` for `text/<name>.tex`), and its PDF,
`text/metric.pdf`, letter pages typeset as the tables are (`metrics.metric_text.write_metric_text`): a section per
dataset, simulations then real data, each opening with the cells it scores, then a paragraph per metric, its name, its
table columns and its definition as the code computes it, each followed by a `% shown in: <page> (<panels>)` line per
page drawing it, a paper figure's page by its name (`plots.pages.PAGE_NAMES`), a diagnosis page as `diagnosis-P<k>` of
`run_benchmark.sh`'s run, or `text/table.tex`. The step stops where a metric column the tables score, the simulations'
per level and the real data's summary metrics, has no paragraph (`metric_text.get_missing_columns`). The including
document needs `\usepackage{amsmath,amssymb}`.
