"""Step 1 of the analysis: each method's runs of the CalicoST simulations scored against their truth.

Runpeng Luo (2026-11-09), from `scripts/run_benchmark.py` (2026-11-04)

Outputs (under outdir):
    tsvs/<simulation_dataset>.tsv.gz: one row per simulation and method, each run scored on its own segments and at
        the genes, a metric's column per level, `<metric> (segment)` and `<metric> (gene)`
        (`metrics.simulation_metrics.merge_levels`).
    tsvs/<simulation_dataset>.median.tsv.gz: one row per simulated class and method, `sim_id`, `method`,
        `#replicates`, then each metric's median over the class' replicates, NaN skipped
        (`plots.simulation.get_class_medians`).
    tsvs/gene_level_cnp/<sample_id>.gene_level_cnp.tsv.gz: per gene, the truth's and each run's fields
        (`metrics.gene_level.get_gene_level_table`), written per simulation as it is scored.
"""

from __future__ import annotations

from scripts.analysis.options import AnalysisOptions
from scripts.cnamaste_port.config import start_time
from scripts.cnamaste_port.logger import get_logger
from scripts.metrics.gene_level import (
    get_autosomal_regions,
    get_region_genes,
    read_genes,
    read_reference_files,
    read_regions,
)
from scripts.metrics.outputs import get_table_dir, get_table_path, write_table
from scripts.metrics.simulation_metrics import (
    get_simulation_methods,
    get_simulation_metrics,
    merge_levels,
    read_simulation_config,
)
from scripts.plots.simulation import get_class_medians

logger = get_logger("scripts.analysis.simulation", start_time=start_time)


def write_simulation_metrics(options: AnalysisOptions) -> None:
    """Score the simulations' runs on their own segments and at the genes, merge the two into one wide table, and
    write it and each method's median per simulated class (`plots.simulation.get_class_medians`).

    Args:
        options: the analysis' options.
    """
    logger.info(f"Step 1 of 4: the simulations' metrics, tsvs/{options.simulation_dataset}.tsv.gz.")
    simulation = read_simulation_config(options.benchmark_yaml)
    methods = options.methods or get_simulation_methods(simulation)
    bed_file, gtf_file = read_reference_files(options.config_yaml)
    regions = get_autosomal_regions(read_regions(bed_file))
    genes = get_region_genes(read_genes(gtf_file), regions)
    gene_level_dir = get_table_dir(options.outdir) / "gene_level_cnp"
    tables = {
        "segment": get_simulation_metrics(simulation, methods),
        "gene": get_simulation_metrics(simulation, methods, genes, gene_level_dir),
    }
    wide = merge_levels(tables)
    write_table(wide, get_table_path(options.outdir, options.simulation_dataset))
    medians = get_table_path(options.outdir, f"{options.simulation_dataset}.median")
    write_table(get_class_medians(wide), medians, float_format="%.6g")
