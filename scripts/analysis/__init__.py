"""The benchmark's analysis, four steps the pipeline runs in order on one options yaml (`options.read_options`):

- `simulation.write_simulation_metrics` -- the simulations' runs scored against their truth;
- `real_data.write_real_data_metrics` -- a Visium HD sample's runs compared;
- `tables.write_tables` -- the LaTeX table and the metric text;
- `figures.write_figures` -- the three PDFs, from the tables the steps before wrote.

Each metric step computes every metric, merges the per-run ones into one wide table, and writes its tables only then.
`run_cna_maste_reproducibility.smk` runs each step as a rule, in the analysis environment (`pyproject.toml`, cnplot
included), `pdflatex` with `geometry`, `amsmath`, `amssymb` and `booktabs` for the tables. The analysis imports no
cnamaste: what it takes from cnamaste is ported in `scripts/cnamaste_port/`.

Notes/References:
    The metrics are in the benchmark's `docs/reference.md` (Metrics). The regions are the bed's
    on chromosomes 1 to 22, where the truth and every method's input lie. The spots compared are UGP's in-tissue spots
    at `side_um` a run labels. The copy-number error's baseline is the spots of a normal cell type
    (`metrics.read_method_input.get_normal_type`); its `#tumor-clones` each run's clones but its normal ones
    (`metrics.real_data_metrics.get_normal_clones`).
    Numbat's outputs: https://kharchenkolab.github.io/numbat/articles/descriptions.html
    cnplot: https://github.com/raphael-group/cnplot
"""
