"""A sample's runs compared without a truth (`run_benchmark.write_real_data_metrics`): each run's clones and CNA states against
a reference run's, its calls per broad cell type, and its copies against its clones' pseudobulk counts.

- `get_key_label`, `read_sample_entry` -- a run key's label, and a sample's entry of the benchmark yaml;
- `get_run_spots`, `get_run_name` -- a run's clone per in-tissue spot, and its `run` and `method` columns;
- `compare_runs` -- each run against the reference run: its clone and CNA ARI, and its cells by state, at the finest
  field both give;
- `get_tumor_spot_labels`, `get_spot_table` -- each run's tumor calls on the spots, and its spot coverage and
  false-positive counts;
- `get_normal_clones`, `get_tumor_cnps` -- a run's normal clones, and its tumor clones' CNA states per gene;
- `get_copy_number_tables` -- the runs' joint segments and copy errors;
- `read_runtime_table`, `get_run_runtime`, `get_runtime_table` -- each run's runtime and memory;
- `merge_run_tables`, `get_summary` -- the per-run tables merged into one wide table, and its summary metrics and
  their provenance per run, `SUMMARY_METRICS` and `PROVENANCE`.

The spots compared are UGP's in-tissue spots at the spots' side that a run labels; a run on coarser bins gives each
bin's clone to every spot inside it (`read_method_input.get_spot_labels`).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from scripts.metrics.cn_error import get_cn_error_table
from scripts.metrics.cna_states import (
    get_cna_states,
    get_comparison_field,
    get_field_order,
    get_phase_free_pairs,
    get_profile_copies,
    get_uncalled,
)
from scripts.metrics.gene_level import get_gene_states, get_inside_regions
from scripts.metrics.outputs import get_json_value
from scripts.metrics.qc_metrics import (
    get_coverage,
    get_false_positive_column,
    get_false_positive_counts,
)
from scripts.metrics.read_method_input import (
    BinCounts,
    MethodInput,
    Run,
    SpotLabels,
    get_bin_side,
    get_method_reader,
    get_normal_type,
    get_run_method,
    get_spot_labels,
    get_yaml_paths,
    read_bin_counts,
)
from scripts.metrics.segment_level import get_joint_segments, get_segment_level_table
from scripts.metrics.simulation_metrics import get_runtime_columns
from scripts.metrics.unsupervised_metrics import (
    compute_clone_ARI,
    compute_cna_ARI,
    compute_cna_confusion,
    get_shared_spots,
)

from cnamaste.config import start_time
from cnamaste.logger import get_logger

logger = get_logger("real_data_metrics", start_time=start_time)

RESOLUTION = "resolution [um]"
"""The runtime table's column of the side, in um, of the bins a run ran on."""

SUMMARY_METRICS = (
    "spot-coverage-rate",
    get_false_positive_column(),
    "clone_ari",
    "cna_ari",
    *get_cna_states(),
    get_uncalled(),
    *(f"{name}_{channel}" for channel in ("rdr", "baf") for name in ("mse", "ss")),
    *get_runtime_columns(),
)
"""A run's summary metrics, the wide table's columns a summary keeps (`get_summary`)."""


def get_summary_metrics() -> tuple[str, ...]:
    """Get a run's summary metrics, `SUMMARY_METRICS`: its spot coverage and false-positive rates, clone and CNA ARI,
    genome composition, copy-number error and runtime columns."""
    return SUMMARY_METRICS


PROVENANCE = (
    "setting",
    "#in-tissue",
    "#called",
    "#normal-spots",
    "#tumor-spots",
    "#false-positive-tumor-spots",
    "#normal-type-spots",
    "level",
    "#shared-spots",
    "#clones",
    "#tumor-clones",
    "#rdr-bins",
    "#baf-bins",
    "#flipped-bins",
)
"""What a run's summary metrics were computed on: its runtime row's setting, its spots, the field and spots it is
compared on, its clones and the bins scored (`get_summary`)."""


def get_key_label(key: str) -> str:
    """Get a run key's label: the key after its method and the `_` or `-` after it, e.g. `016um_msu10` of
    `calicost_016um_msu10`, `ac88605f` of `cnamaste-ac88605f`, empty for a method's own key.

    Args:
        key: the run key.

    Returns:
        The label.
    """
    return key.removeprefix(get_run_method(key))[1:]


def read_sample_entry(benchmark_yaml: Path, sample: str) -> dict[str, Any]:
    """Read a sample's entry of the benchmark yaml, under whichever platform holds it.

    Args:
        benchmark_yaml: the benchmark yaml, `<platform>: <sample>: <run key>: {dir, labels, profile}`.
        sample: the sample.

    Returns:
        The sample's runs by key, `ugp` among them, each path absolute (`read_method_input.get_yaml_paths`).

    Raises:
        KeyError: where no platform holds the sample.
    """
    with open(benchmark_yaml) as handle:
        platforms = yaml.safe_load(handle)
    for samples in platforms.values():
        if sample in samples:
            return get_yaml_paths(dict(samples[sample]))
    raise KeyError(f"{sample} is under no platform of {benchmark_yaml}")


def get_run_spots(key: str, run: MethodInput, run_dir: Path, spots: pd.Index, side_um: int) -> pd.Series:
    """Get a run's clone per in-tissue spot of side `side_um`.

    Args:
        key: the run key.
        run: the run.
        run_dir: its directory.
        spots: the in-tissue spots at `side_um`.
        side_um: the spots' side in um.

    Returns:
        The clone per in-tissue spot the run labels.
    """
    labels = get_spot_labels(Run(get_run_method(key), get_key_label(key), run_dir), run.labels, side_um).labels
    labels = labels[labels.index.isin(spots)].astype(int)
    logger.info(f"{key}: {len(labels):_} of {len(spots):_} in-tissue spots labelled.")
    return labels


def get_run_name(key: str) -> dict[str, str]:
    """Get a run's `run` and `method` columns.

    Args:
        key: the run key.

    Returns:
        `run`, the key, and `method`, its method.
    """
    return {"run": key, "method": get_run_method(key)}


def compare_runs(
    labels: dict[str, pd.Series], runs: dict[str, MethodInput], genes: pd.DataFrame, reference: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare each run with the reference run on the spots both label, at the genes, at the finest field both give
    (`cna_states.get_comparison_field`): phase-free AB pairs, CNA states, or states of total copies.

    Args:
        labels: each run's clone per in-tissue spot, by key.
        runs: each run by key.
        genes: the genes compared (`gene_level.get_region_genes`).
        reference: the reference run's key.

    Returns:
        The agreement, a row per run: `run`, `method`, `reference`, `level`, the field compared, `#shared-spots`,
        `clone_ari`, `cna_ari`; and the confusions long, `run`, `method`, `reference`, `level`, `reference_state`,
        `run_state`, `cells`, each confusion's states in table order (`cna_states.get_field_order`).
    """
    agreement, long = [], []
    for key, run_labels in labels.items():
        if key == reference:
            continue
        level = get_comparison_field(runs[reference].profile, runs[key].profile)
        states = {name: get_gene_states(runs[name].profile, genes, level) for name in (reference, key)}
        if level == "ascn":
            states = {name: get_phase_free_pairs(table) for name, table in states.items()}
        values = get_field_order(level, list(states.values()))
        confusion = compute_cna_confusion(labels[reference], states[reference], run_labels, states[key], values)
        agreement.append(
            {
                **get_run_name(key),
                "reference": reference,
                "level": level,
                "#shared-spots": len(get_shared_spots(labels[reference], run_labels)),
                "clone_ari": compute_clone_ARI(labels[reference], run_labels),
                "cna_ari": compute_cna_ARI(confusion),
            }
        )
        cells = (
            confusion.stack().rename("cells").reset_index().set_axis(["reference_state", "run_state", "cells"], axis=1)
        )
        cells = cells.assign(**get_run_name(key), reference=reference, level=level)
        long.append(cells[["run", "method", "reference", "level", "reference_state", "run_state", "cells"]])
    return pd.DataFrame(agreement), pd.concat(long, ignore_index=True)


def get_tumor_spot_labels(keys: list[str], entry: dict[str, Any], side_um: int) -> list[SpotLabels]:
    """Get each run's clone per spot of side `side_um`, Clonalscope's normal clusters one clone 0
    (`read_method_input.get_method_reader`'s `merge_normal`), so a clone other than 0 is a tumor call.

    Args:
        keys: the run keys.
        entry: the sample's runs (`read_sample_entry`).
        side_um: the spots' side in um.

    Returns:
        Each run's `SpotLabels`, in `keys`' order.
    """
    labels = []
    for key in keys:
        method, run_dir = get_run_method(key), Path(entry[key]["dir"])
        tumor = get_method_reader(method, merge_normal=True)(run_dir)
        run = Run(method, get_key_label(key), run_dir)
        labels.append(get_spot_labels(run, tumor.labels, side_um))
    return labels


def get_spot_table(
    keys: list[str], entry: dict[str, Any], spots: pd.Index, broad: pd.Series, side_um: int
) -> pd.DataFrame:
    """Get each run's coverage of the in-tissue spots (`qc_metrics.get_coverage`) and its false-positive counts
    (`qc_metrics.get_false_positive_counts`), on its tumor calls (`get_tumor_spot_labels`).

    Args:
        keys: the run keys.
        entry: the sample's runs (`read_sample_entry`).
        spots: the in-tissue spots at `side_um`.
        broad: each spot's broad cell type, per barcode.
        side_um: the spots' side in um.

    Returns:
        A row per run, in `keys`' order: `run`, `method`, `resolution`, then `get_coverage`'s columns and
        `get_false_positive_counts`', each but its `method` and `resolution`.
    """
    labels = get_tumor_spot_labels(keys, entry, side_um)
    coverage = get_coverage(labels, spots).drop(columns="method")
    counts = get_false_positive_counts(labels, spots, broad).drop(columns=["method", "resolution"])
    return pd.concat([pd.DataFrame([get_run_name(key) for key in keys]), coverage, counts], axis=1)


def get_normal_clones(key: str, run: MethodInput, run_dir: Path) -> list[int]:
    """Get a run's normal clones: those of the spots its method's reader labels clone 0 with the normal clusters
    merged (`read_method_input.get_method_reader`, `merge_normal`), clone 0, or Clonalscope's normal clusters.

    Args:
        key: the run key.
        run: the run, its clusters as called.
        run_dir: its directory.

    Returns:
        The normal clones, sorted.
    """
    merged = get_method_reader(get_run_method(key), merge_normal=True)(run_dir).labels
    normal = merged.index[merged.to_numpy() == 0]
    return sorted(int(k) for k in set(run.labels[run.labels.index.isin(normal)]))


def get_tumor_cnps(key: str, run: MethodInput, run_dir: Path, states: pd.DataFrame) -> pd.DataFrame:
    """Get a run's tumor clones' CNPs: its CNA states per gene, its normal clones (`get_normal_clones`) left out.

    Args:
        key: the run key.
        run: the run, its clusters as called.
        run_dir: its directory.
        states: its CNA state per gene, `clone<k>` per clone (`gene_level.get_gene_states`).

    Returns:
        `states`' columns of the run's tumor clones, in order.
    """
    normal = {f"clone{k}" for k in get_normal_clones(key, run, run_dir)}
    tumor = [column for column in states.columns if column not in normal]
    logger.info(f"{key}: CNPs of tumor clones {tumor}, normal {sorted(normal)} left out.")
    return states[tumor]


def get_copy_number_tables(
    runs: dict[str, MethodInput],
    labels: dict[str, pd.Series],
    entry: dict[str, Any],
    regions: pd.DataFrame,
    broad: pd.Series,
    side_um: int,
    msr: int,
    assay: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Get the runs giving copies (`metrics.cna_states.get_profile_copies`) on their joint segments, and their copy
    errors at UGP's bins inside the autosomal regions.

    Args:
        runs: each run by key.
        labels: each run's clone per in-tissue spot, by key.
        entry: the sample's entry of the benchmark yaml, `ugp.bb_dir` per side.
        regions: the autosomal regions.
        broad: the broad cell type per spot (`read_method_input.read_broad_cell_types`), its normal types the baseline.
        side_um: the spots' side in um, the `bb` directory read.
        msr: UGP's `min_snp_reads` of the bins.
        assay: UGP's assay type.

    Returns:
        The segment-level table (`metrics.segment_level.get_segment_level_table`) and the error table
        (`metrics.cn_error.get_cn_error_table`), its `method` and `resolution`, `<side>um` of the bins the run ran on
        (`read_method_input.get_bin_side`), after `run`.
    """
    copy_runs = {key: run for key, run in runs.items() if get_profile_copies(run.profile) is not None}
    normal = {key: get_normal_clones(key, run, Path(entry[key]["dir"])) for key, run in copy_runs.items()}
    logger.info(f"Normal clones, left out of `#tumor-clones`: {normal}.")
    segments = get_joint_segments([run.profile for run in copy_runs.values()])
    bins = read_bin_counts(Path(entry["ugp"]["bb_dir"][side_um]) / f"MSR{msr}" / assay)
    inside = get_inside_regions(bins.bins, regions)
    logger.info(f"Bins inside the autosomal regions: {int(inside.sum()):_} of {len(inside):_}.")
    bins = BinCounts(
        bins.bins[inside].reset_index(drop=True),
        bins.x_counts[inside],
        bins.b_counts[inside],
        bins.totals[inside],
        bins.spots,
    )
    reference = broad.index[get_normal_type(broad).to_numpy()]
    logger.info(f"Baseline: {len(reference):_} spots of a normal cell type.")
    errors = get_cn_error_table(copy_runs, labels, normal, segments, bins, reference)
    errors.insert(1, "method", errors["run"].map(get_run_method))
    errors.insert(2, "resolution", [f"{get_bin_side(runs[key].labels)}um" for key in errors["run"]])
    return get_segment_level_table(segments, copy_runs), errors


def read_runtime_table(path: Path) -> pd.DataFrame:
    """Read a sample's runtime table.

    Args:
        path: the table, tab-separated, `method`, `setting`, `RESOLUTION`, then `get_runtime_columns`, per run.

    Returns:
        The table, those columns, NA where not recorded.

    Raises:
        ValueError: where a column is missing.
    """
    columns = ["method", "setting", RESOLUTION, *get_runtime_columns()]
    table = pd.read_csv(path, sep="\t")
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise ValueError(f"{path} has no {', '.join(missing)} column")
    return table[columns]


def get_run_runtime(key: str, run_dir: Path, side_um: int, runtime: pd.DataFrame) -> dict[str, Any]:
    """Get a run's row of the runtime table: its method's, at `side_um`, its `setting` a word of the run directory's
    name or its parent's, set apart by neither a letter, a digit nor a `.`.

    Args:
        key: the run key.
        run_dir: its directory.
        side_um: the side, in um, of the bins it labels.
        runtime: the sample's runtime table (`read_runtime_table`).

    Returns:
        `run`, the key, then the row's columns.

    Raises:
        ValueError: where not exactly one row matches.
    """
    method = get_run_method(key)
    names = f"{run_dir.parent.name}/{run_dir.name}"
    in_names = (
        runtime["setting"]
        .astype(str)
        .map(lambda setting: re.search(rf"(?<![A-Za-z0-9.]){re.escape(setting)}(?![A-Za-z0-9.])", names) is not None)
    )
    found = runtime[(runtime["method"] == method) & (runtime[RESOLUTION] == side_um) & in_names]
    if len(found) != 1:
        raise ValueError(f"{key}: {len(found)} rows of the runtime table are {method} at {side_um} um in {names}")
    row = found.iloc[0].to_dict()
    logger.info(f"{key}: {row['setting']} at {side_um} um, from {names}.")
    return {"run": key, **row}


def get_runtime_table(keys: list[str], entry: dict[str, Any], runs: dict[str, MethodInput]) -> pd.DataFrame:
    """Get each run's row of the sample's runtime table (`get_run_runtime`), at the side of the bins it labels
    (`read_method_input.get_bin_side`).

    Args:
        keys: the run keys.
        entry: the sample's entry of the benchmark yaml, `runtime` and `dir` per run key.
        runs: each run by key.

    Returns:
        A row per run, in `keys`' order: `run`, `setting`, then `get_runtime_columns()`.
    """
    runtime = read_runtime_table(Path(entry["runtime"]))
    rows = [get_run_runtime(key, Path(entry[key]["dir"]), get_bin_side(runs[key].labels), runtime) for key in keys]
    return pd.DataFrame(rows)[["run", "setting", *get_runtime_columns()]]


def merge_run_tables(tables: list[pd.DataFrame], keys: list[str]) -> pd.DataFrame:
    """Merge per-run tables into one wide table, a row per run in `keys`' order, NaN where a table has no row of it; a
    column two tables hold, e.g. `method`, kept once, from the first; a count, a column `#<name>`, an integer.

    Args:
        tables: each a row per run of some runs, by `run`.
        keys: the run keys.

    Returns:
        The wide table: `run`, then each table's columns in order.

    Raises:
        ValueError: where two tables hold a column with different values on a run both hold.
    """
    wide = pd.DataFrame({"run": keys})
    for table in tables:
        shared = sorted((set(table.columns) & set(wide.columns)) - {"run"})
        both = wide.set_index("run")[shared].join(table.set_index("run")[shared], how="inner", rsuffix=" (next)")
        for column in shared:
            held, next_held = both[column], both[f"{column} (next)"]
            if not ((held == next_held) | (held.isna() & next_held.isna())).all():
                raise ValueError(f"{column} differs between two per-run tables")
        wide = wide.merge(table.drop(columns=shared), on="run", how="left", validate="one_to_one")
    counts = [column for column in wide.columns if column.startswith("#")]
    return wide.astype(dict.fromkeys(counts, "Int64"))


def get_summary(wide: pd.DataFrame, dataset: str) -> dict[str, Any]:
    """Get a dataset's summary metrics: its `dataset`, `sample` and `reference`, then per run its `run`, `method` and
    `resolution`, its `metrics` (`SUMMARY_METRICS`) and their `provenance` (`PROVENANCE`), each a column of the wide
    table it has, NaN written None.

    Args:
        wide: the wide table, a row per run (`merge_run_tables`), with `sample` and `reference`.
        dataset: the dataset's name, e.g. `ovarian`.

    Returns:
        The summary.
    """
    runs = [
        {
            **{column: get_json_value(row[column]) for column in ("run", "method", "resolution")},
            "metrics": {column: get_json_value(row[column]) for column in SUMMARY_METRICS if column in row},
            "provenance": {column: get_json_value(row[column]) for column in PROVENANCE if column in row},
        }
        for row in wide.to_dict("records")
    ]
    sample, reference = (str(wide[column].dropna().iloc[0]) for column in ("sample", "reference"))
    return {"dataset": dataset, "sample": sample, "reference": reference, "runs": runs}
