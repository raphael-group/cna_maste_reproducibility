"""Each method's runs on the simulations scored against their truth (`metrics.supervised_metrics`), with their runtime.

- `get_metric_columns` -- the table's columns;
- `get_no_call` -- a run that calls no CNV;
- `read_simulation_config`, `get_simulation_methods` -- the simulations' entry of a benchmark yaml, and its methods;
- `get_run_dir`, `read_simulation_run` -- a method's run of one simulation, Numbat's run that found no CNV as a call of
  none;
- `get_run_metrics` -- one run's metrics and spot coverage, on its segments or at the genes;
- `read_runtime` -- each run's runtime and memory;
- `get_simulation_metrics` -- every simulation's run of each method, one row each, and its gene-level table;
- `get_levels`, `get_level_column`, `merge_levels`, `get_level_table` -- the levels, each level's table merged into one wide table, a metric's
  column per level, and one level's table read back from it.

On its segments, a run's cells are its spots on its own segments (`qa.metrics`); at the genes, its spots at each gene
inside the regions, the gene taking the segment, and the true segment, that overlaps it most
(`metrics.gene_level.get_gene_profile`). `spot-coverage-rate` is `#called` / `#spots`, the truth's spots the run
labels with any clone, normal or tumor, as the real data's (`metrics.qc_metrics.get_coverage`), one value at both
levels.

A method's key in the simulations' entry is a method, or a method then `-` and a variant, read by the method's reader
(`read_method_input.get_run_method`), e.g. `cnamaste-opengm_icm-lattice-lattice`; a key's `solver`, cnamaste's
`<label_solver>-<initializer>-<decoder>` (`run_sim_qa`), is its rows' `solver`, NaN for a key without one.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from scripts.metrics.gene_level import (
    get_gene_level_table,
    get_gene_profile,
    get_truth_input,
)
from scripts.metrics.outputs import write_table
from scripts.metrics.read_method_input import (
    MethodInput,
    get_method_reader,
    get_run_method,
    get_yaml_paths,
)
from scripts.metrics.supervised_metrics import (
    compute_supervised_metrics,
    get_supervised_metrics,
)

from cnamaste.config import start_time
from cnamaste.logger import get_logger
from cnamaste.truth import Truth, read_truth

logger = get_logger("simulation_metrics", start_time=start_time)

RUNTIME = ("Runtime [mins]", "Memory [GB]")
"""A run's runtime and memory, as `simulation.runtime.tsv` names them."""

LEVELS = ("segment", "gene")
"""The rows a run is scored on: its own segments, or the genes inside the regions."""

NUMBAT_NO_CNV = "No CNV remains"
"""The start of the message Numbat's `run_numbat` logs and returns where it stops on finding no CNV (R/main.R)."""


def get_runtime_columns() -> tuple[str, str]:
    """Get a run's runtime and memory columns, `Runtime [mins]` and `Memory [GB]`."""
    return RUNTIME


def get_metric_columns() -> tuple[str, ...]:
    """Get the table's columns: the run, the level it is scored at, its spots and spot coverage rate, the supervised
    metrics, then `RUNTIME`."""
    return (
        "sample_id",
        "method",
        "solver",
        "level",
        "#spots",
        "#called",
        "spot-coverage-rate",
        *get_supervised_metrics(),
        *RUNTIME,
    )


def get_no_call(truth: Truth) -> MethodInput:
    """Get a run that calls no CNV: every truth spot in clone 0, `neu` on every true segment.

    Args:
        truth: the simulation's truth.

    Returns:
        The `MethodInput`.
    """
    profile = truth.profile[["chr", "start", "end"]].set_axis(["CHR", "START", "END"], axis=1).assign(clone0="neu")
    return MethodInput(pd.Series(0, index=truth.barcodes), profile)


def read_simulation_config(benchmark_yaml: Path) -> dict[str, Any]:
    """Read the simulations' entry of a benchmark yaml.

    Args:
        benchmark_yaml: a benchmark yaml, e.g. `configs/benchmark.local.yaml`, its `simulation`: `groundtruth_dir`, a
            directory per simulation; `runtime`, a `simulation.runtime.tsv`; per method, `dir`, a template of a run's
            directory in its `{sample_id}`.

    Returns:
        The entry, each path absolute (`read_method_input.get_yaml_paths`).
    """
    with open(benchmark_yaml) as handle:
        return get_yaml_paths(dict(yaml.safe_load(handle)["simulation"]))


def get_simulation_methods(simulation: dict[str, Any]) -> list[str]:
    """Get the methods with runs in a simulations' entry: its keys holding a `dir`, in order.

    Args:
        simulation: `read_simulation_config`'s entry.

    Returns:
        The methods.
    """
    return [key for key, value in simulation.items() if isinstance(value, dict) and "dir" in value]


def get_run_dir(simulation: dict[str, Any], method: str, sample_id: str) -> Path:
    """Get a method's run directory of one simulation: its `dir` template filled with `sample_id`; a template whose last
    part is a pattern, e.g. `run_*`, cnamaste's `run_<id>` named by each run's configuration, the one directory it
    matches.

    Args:
        simulation: `read_simulation_config`'s entry.
        method: a method of `get_simulation_methods(simulation)`.
        sample_id: the simulation.

    Returns:
        The directory; where the pattern matches none, the template filled, so the run reads as having no output.

    Raises:
        ValueError: where the pattern matches more than one directory.
    """
    path = Path(simulation[method]["dir"].format(sample_id=sample_id))
    if not any(c in path.name for c in "*?["):
        return path
    found = sorted(p for p in path.parent.glob(path.name) if p.is_dir())
    if len(found) > 1:
        raise ValueError(f"{sample_id} {method}: {len(found)} runs match {path}, one wanted")
    return found[0] if found else path


def read_simulation_run(simulation: dict[str, Any], method: str, sample_id: str, truth: Truth) -> MethodInput | None:
    """Read a method's run of one simulation (`read_method_input.get_method_reader`); a Numbat run that stopped finding
    no CNV (`NUMBAT_NO_CNV`) as a call of none (`get_no_call`).

    Args:
        simulation: `read_simulation_config`'s entry.
        method: the method.
        sample_id: the simulation.
        truth: its truth.

    Returns:
        The run, or None where it wrote no output, logged.
    """
    run_dir = get_run_dir(simulation, method, sample_id)
    log = run_dir / "log.txt"
    stopped = not any(run_dir.glob("clone_post_*.tsv")) and log.is_file() and NUMBAT_NO_CNV in log.read_text()
    if method == "numbat" and stopped:
        logger.info(f"{sample_id} numbat: {NUMBAT_NO_CNV}, scored as a call of none.")
        return get_no_call(truth)
    try:
        return get_method_reader(get_run_method(method))(run_dir)
    except FileNotFoundError as error:
        logger.warning(f"{sample_id} {method}: no output, left unscored ({error}).")
        return None


def get_run_metrics(run: MethodInput, truth: Truth, genes: pd.DataFrame | None = None) -> dict[str, Any]:
    """Get one run's spots and metrics against its truth, on its segments or at the genes.

    Args:
        run: the run.
        truth: its simulation's truth.
        genes: the genes inside the autosomal regions (`metrics.gene_level.get_region_genes`), or None to score the
            run's segments.

    Returns:
        `#spots`, the truth's, `#called`, those the run labels, `spot-coverage-rate`, `#called` / `#spots`, then the
        supervised metrics.
    """
    if genes is None:
        found = compute_supervised_metrics(run, truth)
    else:
        found = compute_supervised_metrics(MethodInput(run.labels, get_gene_profile(run.profile, genes)), truth)
    called = int(run.labels.reindex(truth.barcodes).notna().sum())
    return {
        "#spots": len(truth.barcodes),
        "#called": called,
        "spot-coverage-rate": called / len(truth.barcodes),
        **found,
    }


def read_runtime(path: Path) -> pd.DataFrame:
    """Read each run's runtime and memory.

    Args:
        path: a `simulation.runtime.tsv`, `simulation_id`, `method`, then `RUNTIME`.

    Returns:
        `sample_id`, `method`, then `RUNTIME`, one row per run.
    """
    return pd.read_csv(path, sep="\t").rename(columns={"simulation_id": "sample_id"})[["sample_id", "method", *RUNTIME]]


def get_simulation_metrics(
    simulation: dict[str, Any],
    methods: list[str],
    genes: pd.DataFrame | None = None,
    gene_level_dir: Path | None = None,
) -> pd.DataFrame:
    """Score every simulation's run of each method, add its runtime, and write each simulation's gene-level table.

    Args:
        simulation: `read_simulation_config`'s entry: each simulation's truth (`truth.read_truth`), the methods' runs
            (`get_run_dir`) and their runtime.
        methods: the methods, of `get_simulation_methods(simulation)`.
        genes: the genes inside the autosomal regions, to score each run at the genes, or None for its segments.
        gene_level_dir: where `<sample_id>.gene_level_cnp.tsv.gz` is written per simulation
            (`metrics.gene_level.get_gene_level_table`), or None; only with `genes`.

    Returns:
        One row per simulation and method, `get_metric_columns()`, its `level` `gene` with `genes`, else `segment`; a
        run with no output NaN but for its runtime.
    """
    rows = []
    for truth_dir in sorted(p for p in Path(simulation["groundtruth_dir"]).iterdir() if p.is_dir()):
        truth = read_truth(truth_dir)
        runs = {}
        for method in methods:
            run = read_simulation_run(simulation, method, truth_dir.name, truth)
            found = get_run_metrics(run, truth, genes) if run is not None else {}
            rows.append(
                {"sample_id": truth_dir.name, "method": method, "solver": simulation[method].get("solver"), **found}
            )
            if run is not None:
                runs[method] = run
        if genes is not None and gene_level_dir is not None:
            path = Path(gene_level_dir) / f"{truth_dir.name}.gene_level_cnp.tsv.gz"
            write_table(get_gene_level_table(genes, runs, get_truth_input(truth)), path)
    scored = pd.DataFrame(rows).merge(read_runtime(Path(simulation["runtime"])), how="left")
    scored["level"] = "segment" if genes is None else "gene"
    logger.info(f"Scored {len(scored)} runs of {len(scored) // max(len(methods), 1)} simulations.")
    return scored.reindex(columns=list(get_metric_columns())).astype({"#spots": "Int64", "#called": "Int64"})


def get_levels() -> tuple[str, ...]:
    """Get the rows a run is scored on, `LEVELS`: `segment`, its own segments, and `gene`, the genes."""
    return LEVELS


def get_level_column(column: str, level: str) -> str:
    """Get a metric's column at a level in the wide table, `<column> (<level>)`, e.g. `sensitivity (gene)`.

    Args:
        column: the metric, of `get_supervised_metrics()`.
        level: one of `LEVELS`.

    Returns:
        The column.
    """
    return f"{column} ({level})"


def get_level_metrics() -> tuple[str, ...]:
    """Get the metrics scored per level, `get_supervised_metrics()`."""
    return get_supervised_metrics()


def merge_levels(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Merge each level's table into one wide table, a row per simulation and method: `sample_id`, `method`, `solver`,
    `#spots`, `#called`, `spot-coverage-rate`, each level's metrics (`get_level_column`), then `RUNTIME`.

    Args:
        tables: each level's table (`get_simulation_metrics`), by level, its rows in one order.

    Returns:
        The wide table.

    Raises:
        ValueError: where the levels' runs, spots or runtime differ.
    """
    shared = ["sample_id", "method", "solver", "#spots", "#called", "spot-coverage-rate", *RUNTIME]
    first = next(iter(tables.values()))[shared].reset_index(drop=True)
    scored = []
    for level, table in tables.items():
        if not table[shared].reset_index(drop=True).equals(first):
            raise ValueError(f"the {level} table's runs, spots or runtime differ from the first level's")
        metrics = table[list(get_level_metrics())].reset_index(drop=True)
        scored.append(metrics.rename(columns=lambda column, level=level: get_level_column(column, level)))
    return pd.concat([first[shared[:6]], *scored, first[list(RUNTIME)]], axis=1)


def get_level_table(wide: pd.DataFrame, level: str) -> pd.DataFrame:
    """Get one level's table from the wide table (`merge_levels`), as `get_simulation_metrics` returns it.

    Args:
        wide: the wide table.
        level: one of `LEVELS`.

    Returns:
        `get_metric_columns()`, its `level` `level`.
    """
    columns = {get_level_column(column, level): column for column in get_level_metrics()}
    table = wide.rename(columns=columns).assign(level=level)
    return table.reindex(columns=list(get_metric_columns()))
