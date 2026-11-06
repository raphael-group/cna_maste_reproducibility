# Figures: the pages of combined.pdf, diagnosis.pdf and appendix.pdf

What each page of the benchmark's three PDFs draws, all written by `run_benchmark.py`'s last step under `<outdir>/plots/`, each page alone under `<outdir>/pdf-pages/` ([README, Output](../README.md#output)).

The benchmark's figures are three PDFs under `plots/`, all written by `run_benchmark.py`'s last step: `combined.pdf`,
`diagnosis.pdf` and `appendix.pdf` ([Appendix](#appendix)); each run also writes every page of each alone, every
page cropped to its content, 0.02 in of white round it, as the lab's figure rules save and LaTeX places a figure
(`plots.pages.SAVE_OPTIONS`; `docs/figure-rule.md` MP-3, LF-3): `pdf-pages/<stem>-P<k>.pdf`, e.g.
`pdf-pages/combined-P4.pdf`, a paper figure's page under its own name (`plots.pages.PAGE_NAMES`): combined.pdf's pages 1
to 3 `main-fig-zenodo-sim-result.pdf`, `main-fig-cnamaste-result.pdf` and `main-fig-cnp-compare.pdf`; the appendix's
UMI, spatial labels, cells per altered state and state Jaccard index pages `suppl-fig-umi.pdf`,
`suppl-fig-ovarian-label.pdf`, `suppl-fig-cna-composition.pdf` and `suppl-fig-cna-jaccard.pdf`; the simulations'
runtime and memory, and cnamaste solver pages `suppl-fig-zenodo-sim-runtime.pdf` and
`suppl-fig-zenodo-sim-ablation.pdf`. Its first two steps write the metric tables, one wide table per dataset and its
detail tables, which the appendix draws
([Metrics](pipeline.md#metrics)). Its figure step reads each method's run once; its pages are drawn with
`scripts.plots.layout.plot_annotation_page`, each method named by the side of the bins it ran on, read from
its barcodes.

> [!IMPORTANT]
> Every page of the three PDFs shows each panel's letter alone, `(a)`, `(b)`, ..., no title beside it. Each PDF's last
> page names the panels: per page with a titled panel, `Page <k>` and the page's own title, then a line per panel,
> `(x) <title>` (`layout.plot_panel_titles`). A panel's title is drawn as its letter's label, not as text on the page
> (`layout.draw_panel_letter`, `plots.style.add_panel_letters`, `plot_annotation_page`), and a new page follows suit.
> Every text is the page's size, `[tool.cnamaste.figures]`'s `font_size`, but axis labels, a colorbar's among them,
> 1.25 times it, the largest text, so no key, tick label or letter is larger (`layout.set_text_sizes`).

## combined.pdf

- Page 1: the simulations' metrics at the genes, (a) CNA sensitivity, (b) clone ARI, (c) CNA recall and (d) N-CNA
  recall, a box per simulated class and method (`plots.pages.plot_simulation_pages`), as the [Appendix](#appendix)
  describes; the key right of (d), centred on its axes.
- Page 2, where the cnamaste run wrote its per-spot bins, `binned_data.npz` (`run.write_binned_data`), beside its
  `rdrbaf_final_nstates<n>_smp.npz` in the `--cnamaste_labels` directory: cnamaste's own page,
  `cnamaste.plot.layout.plot_combined_page`, broken into panels in two columns
  (`plots.cnamaste_page.plot_cnamaste_page`). Left: (a) the H&E slide, (b) cnamaste's clones, coloured as the
  appendix's spatial panels, keyed inside (b)'s top left, a clone a row, with each clone's share of the in-tissue spots, the run logging
  any spot the key covers; both framed on cnamaste's spots and the slide's origin, the x axis ending at the spots' right
  edge. Right, 0.1 in past them and their axis labels: cnamaste's own scatter tracks
  (`plot_genomic.plot_clones_genomic`), a letter per clone, (c) Clone 0 and (d) Clone I, each clone's pseudobulk RDR
  then BAF per bin, the normal clone's BAF alone, its RDR 1 by construction; each clone named once over its first track,
  its letter on that line, a tumour clone's states' shares of its bins keyed right of its name, read row by row, each
  track's `RDR` or `BAF` 1 pt from its tick labels, an RDR track 0 to 4 (`plots.cnamaste_page.TRACK_RDR_YLIM`, where
  cnamaste's default is 6), a bin above it not drawn; (e) under them cnamaste's copy-number profile, a row per integer
  clone, on its bins, as cnamaste's `plot_copy_number_profile.plot_state_profile` draws it, (1,1) white
  (`layout.draw_copy_profile`), the chromosomes named under it. (b), the tracks and the profile draw the integer
  clones, `clone_labels_int.tsv`'s, as cnamaste's page does given `int_groups`: each named by and drawn at its first
  clone, at its normal fraction, its tracks pooling its spots (`plots.cnamaste_page.get_integer_fit`). Each column's
  letters from the left of its axis labels. A (1,1) bin's dot is grey, where cnamaste's tracks draw it steel blue
  (`palette.TRACK_NORMAL`). The tracks' black lines are cnamaste's own expected lines, as its page draws them: each
  clone's decoded (A,B) mixed with normal (1,1) cells at its normal fraction, `--cnamaste_normal_fraction`
  (`hmm_ab_decoder.mixture`), RDR over the clone's Z_c; the normal clone, at 1, reads BAF 1/2 throughout. run_cnamaste
  writes the normal fractions to no file, so they are given. The tracks fill the left column's height, each at least 0.4
  in.
- Page 3 (`layout.plot_profile_page`), in two columns. Left, (a): each method's clone profile, cnamaste, CalicoST,
  Numbat, Clonalscope, top to bottom, its method's name over it, drawn as page 2's (e) on the cnamaste run's bins
  (`cnv_seglevel.tsv`'s rows, `cnamaste.plot.layout.GenomicAxis`): each bin takes the state of the run's segment
  overlapping it most (`plots.pages.get_bin_profile`, `metrics.gene_level.get_max_overlap`); CalicoST's and cnamaste's
  integer copies as cnamaste's own `plot_copy_number_profile.plot_state_profile` draws them
  (`layout.draw_copy_profile`), a bin no segment overlaps drawn (1,1) as it draws one; Numbat's CNV state and
  Clonalscope's total-copy state, in Numbat's colours, CNLoH's and DEL's swapped, and names, as it draws copies
  (`layout.draw_state_profile`), such a bin left blank; each run of one state within a chromosome one segment, a CNA's
  ends outlined, no gap between chromosomes and no centromere line. Only the chromosomes where a drawn clone of any
  method holds a CNA, a state other than `NEU` or `(1,1)`, are drawn, each whole; cnamaste's chromosome names under the
  last profile, `chr` at the left. Each clone's row is named in two lines, its name over its share of the in-tissue
  spots as the appendix's spatial keys count it, e.g. `Clone I` over `(<k>%)`, the names left-aligned in one column for every method;
  the rows fill the column to the right column's height, each at least 0.2 in and a name's height. Two keys under one
  name, `CNA` (`layout.KEY_NAME`): left to right `NEU`, `DEL`, `AMP`, `CNLoH` (`run_benchmark.CNA_KEY_ORDER`), the
  states drawn; under it, unnamed, the `(A,B)` states in cnamaste's profile colours (`cnamaste.palette.COPY_STATES`,
  a `B > A` state lighter, through `plot_copy_number_profile.state_fill`), (1,1) white as Numbat's `NEU`, grouped by
  BAF, min(A, B) / (A + B), ascending, so a state and its phase flip share a group, each group in cnamaste's order,
  `(0,0)` a group first. Both keys are laid out as cnplot's copy-number key (`_plot_cnv_legend`): the name in the first
  key's left margin, the swatches from the profiles' left edge, one swatch width for both keys, each state's name under
  its swatch, untitled, the groups packed into rows of the profiles' width, each row justified to the key's widest; the
  keys under the chromosome names, at the page's bottom. Right (`plots.run_comparison.get_run_comparison_panels`): (b)
  the CNP similarity and (c) the clone ARI of
  every two runs, side by side, each a heatmap of 0.24 in a run, a row and a column per run in the methods' order, its
  value in each cell, the runs named by their methods under each and left of (b); (d) under them, each run's RDR and BAF
  skill score, `1 - MSE / MSE_{(1,1)}` against the diploid, (1,1), prediction (`plots.cn_error.draw_skill_score_panel`),
  as diagnosis.pdf's copy-number error page draws them, its axis up to 0.75 unless a score passes it, the runs keyed
  from its left axis, under its tick labels, in one row where they fit; (d) grows where the left column is taller.
- Page 4: the panels of pages 1 to 3 named, CalicoST with its label; (a) of page 3 names each profile's run.

Clonalscope is `--clonalscope_dir`, its normal clusters, `malignant` N, drawn as one clone 0 so its tumor spots are
its other clones, as every method's are; its states are `del`, `neu` and `amp` by total copies, drawn as Numbat's
`DEL`, `NEU` and `AMP`. Without it, Clonalscope is a placeholder on page 3 and the appendix's spatial page. CalicoST is `--calicost_dir`,
`msu20` in `run_benchmark.sh`'s run, named with its label. A label `msu<k>`, CalicoST's
`min_snpumi_perspot` ([`arg_parse.py`](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/arg_parse.py#L30)),
is written `min-snp-umi: <k>` in every figure (`plots.style.get_setting_label`); tables keep the label. A figure names a
run by its method, the side of its bins and its `msu<k>` alone, whatever else its key holds, e.g. a run's id
(`plots.style.get_run_label`).

## diagnosis.pdf

`diagnosis.pdf`, page 1: the clones of Numbat, of each sweep run named by `--false_positive_sweep <label>` and of
cnamaste, one panel per row, each with its false-positive tumor spots in a third colour: a tumor call on a spot of a
normal type, neither `Epithelial` nor `Low confidence` (`scripts/metrics/qc_metrics.py`'s numerator). A clone is one
colour on every page. Page 2, with a `--bb_dir` at `--side_um`: the text column wide, a row per run, 2.3 in each
(`plots.style.get_panel_rows`), its false-positive spots' total UMI and SNP-covering UMI against its other called
in-tissue spots', as densities, with each group's spots and median; the numbers in `tsvs/diagnosis.tsv.gz`. The other
called spots hold every cell type, epithelial ones among them. Then a page of the
copy-number error ([Copy-number error](pipeline.md#copy-number-error)), each run in method order: (a) the mean squared error on a
log axis and (b) the skill scores against the diploid prediction, `1 - MSE / MSE_{(1,1)}`, each a group per channel, RDR
and BAF, named by the bins scored, the diploid prediction's error on the run's pseudobulk a black tick across its bar,
no bar where a run has no allele copies; side by side, each panel's axes 1.8 in tall
(`plots.cn_error.CN_ERROR_AXES_HEIGHT`), the key under them (`plots.style.get_panel_rows`). Last, the page naming
every panel.

`diagnosis.pdf`, before its simulation pages ([Appendix](#appendix)) and its last page, with
`--calicost_sweep <label>=<dir>` once per run of a CalicoST parameter
sweep: each run's clones, then each run's profile, each clone's row named with its spots as its key names it, every
panel lettered, named `CalicoST (<side> um, <label>)`, e.g.
`CalicoST (16 um, min-snp-umi: 10)`, then the table read where it is the fallback's. A tumour clone that the table
read leaves out is named beside its profile's name
([CalicoST clones missing from the unconstrained table](pipeline.md#calicost-clones-missing-from-the-unconstrained-table)).

> [!NOTE]
> `run_benchmark.sh` passes no sweep run, no `--sweep_run`, `--calicost_sweep` or `--false_positive_sweep`: its
> diagnosis.pdf has no sweep pages, its pages 1 and 2 hold Numbat and cnamaste, and it writes no `cn_error.sweep`. The
> sweep's runs stay in `configs/benchmark.local.yaml`: `calicost_016um` (`H1-N3KXRGQ_016um_msu20`),
> `calicost_016um_msu10`, `_msu15`, `_msu25`, `_msu50` and `calicost_032um`.

## Appendix

`run_benchmark.py`'s figure step writes the appendix, one figure per page, `appendix.pdf` by default:

- Page 1, with `--bb_dir <side_um>=<dir>`: the text column wide, one row per `--bb_dir`, 2.3 in each
  (`plots.style.get_panel_rows`), in the order given and lettered (a), (b), named by sample, bin side and in-tissue
  spots, its per-spot total UMI (its `bb` AnnData) left and SNP-covering UMI (A + B of its `multi_snp` allele matrices)
  right (`plots.umi.plot_umi_page`). UGP's 16 um bins are in `bb/`, its 32 um bins in `bb-032um/`.
- Spatial panels, three a row, in the methods' order (`layout.plot_annotation_page`): the H&E slide, each spot's cell
  type (not epithelial, then epithelial, low-confidence spots left out) and cnamaste's clones; CalicoST's, Numbat's and
  Clonalscope's. Each panel is keyed inside its top left, as combined.pdf page 2's (b) (`layout.draw_inside_key`),
  untitled, a label a row; a key that covers spots takes a label more a row, read row by row, while it fits the panel,
  the layout covering fewest kept and logged; each label with its in-tissue spots' share of the in-tissue spots, so a
  key need not sum to 100%. Each panel's letter is over it, from the left of its y label, or of its ticks. The page is
  at most `page_size("full")` tall, the text block less its caption, so it prints unscaled: the panels narrow, the
  grid centred, where it would run past. A panel's spots are stored 5 pixels a side (`layout.SPOT_IMAGE_REPEAT`), as
  on combined.pdf's page 2.
- From the sample's tables (`--real_data_dataset`), each run in its method's fill, the simulation pages' box colour, and
  named by its method and the side of the bins it ran on, e.g. `Numbat (32 um)`: a page of each run's proportion of
  cells in `del`, `amp`, `cnloh` and `bamp`, its `all` row, a group per state, each bar's value to 0.001 above it, the
  runs in method order (`plots.run_comparison.plot_genome_composition`, `plots.style.get_ordered_runs`); each run's spot
  coverage rate, false-positive rate, runtime and memory are a LaTeX table ([LaTeX tables](pipeline.md#latex-tables)). Then a page of
  each run's state Jaccard index against the reference's at the run's `level` ([Between runs](pipeline.md#between-runs-real-data)), a
  panel per run in method order, named by the `level`: a row per state of the reference, a column per state of the run,
  the states either calls on a cell both call, each entry the Jaccard index of the cells the reference calls its row's
  state and the run its column's, over the cells both call, its diagonal the agreement on each state
  (`plots.run_comparison.plot_cna_jaccard`); the run's states turned a quarter under it; the panels in rows of up to
  ten states, top-aligned 0.1 in apart past their labels, the rows' first frames on one left edge, the colour bar right
  of the first row; a first-row cell 0.24 in a side, as combined.pdf's page 3 (b) and (c); a later row's panels one
  side, so the row with its labels is as wide as the first row with its colour bar (`run_comparison.get_fill_side`);
  the page as wide and tall as its content. The CNP similarity and clone ARI of every two runs are combined.pdf's
  page 3 (b) and (c), the copy-number error a page of diagnosis.pdf.
- From the simulations' wide table (`--simulation_dataset`): a page of its metrics at the genes, combined.pdf's
  page 1 (`plots.pages.plot_simulation_pages`) (`plots.pages.SIMULATION_LEVEL`, `simulation_metrics.get_level_table`; the table keeps both levels), a
  box per simulated class and method, two rows of two panels (`plots.pages.SIMULATION_COLUMNS`): (a) CNA sensitivity
  and (b) clone ARI, (c) CNA recall and (d) N-CNA recall; then one page of the runtime and memory of the methods with
  any; then, where the table holds cnamaste's solvers, a page of them
  (`plots.simulation.plot_solver_ablation`): cnamaste `cleanup3`'s `run_sim_qa` solver comparison
  (`qa.plot.plot_validation`, `ddc060b`) on one page, (a) CNA
  sensitivity, (b) 1 - normal specificity on a log axis, (c) N-CNA recall, (d) clone ARI of the integer clones, a box
  per simulated class and solver, `expansion-lattice-lattice` first in cnamaste's blue, the others in `SOLVER_COLOURS`.
  In (b), a run at specificity 1 is drawn at a floor, the power of 10 a decade under the lowest positive value's, ticked
  0.00 (`plots.simulation.get_complement`), its point drawn as the others; the boxes count it at the floor; the axis is
  ticked 10^-3, 10^-2, 10^-1 and 1.00 above it (`plots.simulation.FLOOR_TICKS`). The metrics
  and solvers pages leave out the classes (6,3, 50M) and (3,3, 50M) (`plots.pages.DROPPED_CLASSES`); diagnosis.pdf
  draws both pages in full, every class, the metrics page with (a) normal specificity and (b) spot coverage rate first,
  three rows of two panels (`plots.pages.plot_simulation_diagnosis`). Each page is the text column wide, laid out
  tight (`plots.style.get_panel_rows`, `docs/figure-rule.md` FP-5). The metrics pages name the classes under the last
  row, no x-axis label, the key 0.1 in right of the last panel, centred on its axes, the rows 0.16 in apart past their labels, a panel's axes 1.43 in
  tall (`plots.simulation.METRIC_AXES_HEIGHT`, `plots.style.set_axes_height`); the runtime and solver pages are 2.3 in
  a row, the key under the panels.
- Every appendix page is untitled; the last page names every panel, a simulation page's panels with the level they are
  scored at.
