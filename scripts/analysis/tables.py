"""Step 3 of the analysis: the sample's LaTeX table and every metric's definition as LaTeX text, each typeset.

Runpeng Luo (2026-11-09), from `scripts/run_benchmark.py` (2026-11-04)

Outputs (under outdir):
    text/<table>.tex, text/<table>.pdf: the runs' spot coverage rate, false-positive rate, runtime and memory, a row per
        run in method order, a booktabs tabular in a captioned, labelled `table` float, typeset on letter pages, and
        pdflatex's log, `text/<table>.log` (`metrics.latex_tables`).
    text/<metric_text>.tex, text/<metric_text>.pdf: each metric's definition as it is computed, a paragraph per
        metric, each followed by a `% shown in:` line per page and panel drawing it, typeset on letter pages, and
        pdflatex's log (`metrics.metric_text`).
"""

from __future__ import annotations

from scripts.analysis.options import AnalysisOptions
from scripts.cnamaste_port.config import start_time
from scripts.cnamaste_port.logger import get_logger
from scripts.metrics.latex_tables import write_latex_table
from scripts.metrics.metric_text import write_metric_text

logger = get_logger("scripts.analysis.tables", start_time=start_time)


def write_tables(options: AnalysisOptions) -> None:
    """Write and typeset the sample's LaTeX table, from step 2's wide table, and the metric text.

    Args:
        options: the analysis' options.
    """
    logger.info(
        f"Step 3 of 4: the LaTeX tables and metric text, text/{options.table}.tex and text/{options.metric_text}.tex."
    )
    write_latex_table(options.outdir, options.real_data_dataset, options.table)
    write_metric_text(options.outdir, options.metric_text)
