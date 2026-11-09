"""Step 2 of the analysis: a Visium HD sample's runs compared with one another, its cell types and its bins.

Runpeng Luo (2026-11-09), from `scripts/run_benchmark.py` (2026-11-04)

Outputs (under outdir):
    tsvs/<real_data_dataset>.tsv.gz: one row per run (`metrics.real_data_metrics.merge_run_tables`): `sample`, `run`,
        `method`, `resolution`; its spot coverage and false-positive counts per broad cell type
        (`metrics.real_data_metrics.get_spot_table`); its genome composition's `all` row, the share of its cells per
        CNA state and `uncalled` (`metrics.qc_metrics.compute_genome_composition`); the reference run, the field
        compared, `level`, their shared spots, clone ARI and CNA ARI (`metrics.real_data_metrics.compare_runs`), NaN
        for the reference; its RDR and BAF error at the bins, the runs giving copies
        (`metrics.cn_error.get_cn_error_table`); its runtime row's `setting`, runtime and memory
        (`metrics.real_data_metrics.get_runtime_table`).
    tsvs/<real_data_dataset>.json: its summary metrics per run, and their provenance
        (`metrics.real_data_metrics.get_summary`).
    tsvs/<sample>.<table>.tsv.gz: the sample's detail tables: `gene_level_cnp`, per gene, each run's fields;
        `genome_composition`, per run and clone, its genes' shares by CNA state and `uncalled`;
        `clone_ari` and `cnp_similarity`, per ordered pair of runs
        (`metrics.unsupervised_metrics.compute_pairwise_clone_ARI`, `compute_pairwise_cnp_similarity`), a run against
        itself 1; `cna_confusion`, per run, the cells by `(reference state, run state)` at its `level`
        (`metrics.cna_states.get_comparison_field`); `segment_level_cnp`, per joint segment of the runs giving copies
        (`metrics.segment_level.get_segment_level_table`); with `sweep_runs`, `cn_error.sweep`, the copy-number error
        of the reference run and the sweep runs, scored as one panel on the bins they share.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scripts.analysis.options import RUNS, AnalysisOptions
from scripts.metrics.gene_level import (
    get_autosomal_regions,
    get_gene_level_table,
    get_gene_states,
    get_region_genes,
    read_genes,
    read_reference_files,
    read_regions,
)
from scripts.metrics.outputs import get_json_path, get_table_path, write_json, write_table
from scripts.metrics.qc_metrics import compute_genome_composition
from scripts.metrics.read_method_input import get_method_reader, get_run_method, read_broad_cell_types, read_spots
from scripts.metrics.real_data_metrics import (
    compare_runs,
    get_copy_number_tables,
    get_run_name,
    get_run_spots,
    get_runtime_table,
    get_spot_table,
    get_summary,
    get_tumor_cnps,
    merge_run_tables,
    read_sample_entry,
)
from scripts.metrics.unsupervised_metrics import compute_pairwise_clone_ARI, compute_pairwise_cnp_similarity

from cnamaste.config import start_time
from cnamaste.logger import get_logger

logger = get_logger("scripts.analysis.real_data", start_time=start_time)


def write_real_data_metrics(options: AnalysisOptions) -> None:
    """Read the sample's runs, compute every table, merge the per-run ones into one wide table, and write them all,
    with the wide table's summary.

    Args:
        options: the analysis' options, its reference among its runs (`options.get_argument_error`).
    """
    logger.info(f"Step 2 of 4: {options.sample}'s runs compared, tsvs/{options.real_data_dataset}.tsv.gz.")
    keys = options.runs or list(RUNS)
    entry = read_sample_entry(options.benchmark_yaml, options.sample)
    spots = read_spots(Path(entry["ugp"]["spots"]))
    broad = read_broad_cell_types(Path(entry["cell_type"]["labels"]))
    bed_file, gtf_file = read_reference_files(options.config_yaml)
    regions = get_autosomal_regions(read_regions(bed_file))
    genes = get_region_genes(read_genes(gtf_file), regions)

    runs, labels, states = {}, {}, {}
    for key in keys:
        run_dir = Path(entry[key]["dir"])
        runs[key] = get_method_reader(get_run_method(key), integer_profile=True)(run_dir)
        labels[key] = get_run_spots(key, runs[key], run_dir, spots, options.side_um)
        states[key] = get_gene_states(runs[key].profile, genes)

    composition = pd.concat(
        [
            compute_genome_composition(
                states[key].to_numpy(dtype=str), np.ones(len(genes)), float(len(genes)), labels[key]
            ).assign(**get_run_name(key))
            for key in keys
        ],
        ignore_index=True,
    )
    composition = composition[["run", "method", *composition.columns[:-2]]]
    agreement, cells = compare_runs(labels, runs, genes, options.reference)
    cnps = {key: get_tumor_cnps(key, runs[key], Path(entry[key]["dir"]), states[key]) for key in keys}
    segment_level, errors = get_copy_number_tables(
        runs, labels, entry, regions, broad, options.side_um, options.msr, options.assay
    )
    # NB a detail table's float format: None writes full precision.
    details: dict[str, tuple[pd.DataFrame, str | None]] = {
        "gene_level_cnp": (get_gene_level_table(genes, runs), None),
        "genome_composition": (composition, None),
        "cna_confusion": (cells, None),
        "clone_ari": (compute_pairwise_clone_ARI(labels), "%.6g"),
        "cnp_similarity": (compute_pairwise_cnp_similarity(cnps), "%.6g"),
        "segment_level_cnp": (segment_level, None),
    }
    if options.sweep_runs:
        for key in options.sweep_runs:
            run_dir = Path(entry[key]["dir"])
            runs[key] = get_method_reader(get_run_method(key), integer_profile=True)(run_dir)
            labels[key] = get_run_spots(key, runs[key], run_dir, spots, options.side_um)
        panel = {key: runs[key] for key in [options.reference, *options.sweep_runs]}
        sweep_errors = get_copy_number_tables(
            panel, labels, entry, regions, broad, options.side_um, options.msr, options.assay
        )[1]
        details["cn_error.sweep"] = (sweep_errors, "%.6g")

    # NB a composition's and a copy error's spots are the run's called spots, checked equal as they merge.
    shares = composition[composition["clone"].astype(str) == "all"].drop(columns="clone")
    wide = merge_run_tables(
        [
            get_spot_table(keys, entry, spots, broad, options.side_um),
            shares.rename(columns={"#spots": "#called"}),
            agreement,
            errors.rename(columns={"#spots": "#called"}),
            get_runtime_table(keys, entry, runs),
        ],
        keys,
    )
    wide.insert(0, "sample", options.sample)
    wide["reference"] = options.reference

    for name, (table, float_format) in details.items():
        write_table(table, get_table_path(options.outdir, f"{options.sample}.{name}"), float_format=float_format)
    write_table(wide, get_table_path(options.outdir, options.real_data_dataset))
    write_json(get_summary(wide, options.real_data_dataset), get_json_path(options.outdir, options.real_data_dataset))
