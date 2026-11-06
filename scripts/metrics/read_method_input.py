"""Readers of the benchmark's inputs: each method's run, its spot labels and clone profile, and the sample's spots.

- `read_cnamaste`, `read_calicost`, `read_numbat`, `read_clonalscope` -- one run directory each, its `MethodInput`;
  `get_labelled_profile` -- a cnamaste profile's clones its integer clones name;
- `get_calicost_run`, `get_calicost_ploidy` -- a CalicoST directory's run, and the table it is read at, its
  fallbacks `get_calicost_fallbacks`;
- `get_method_reader`, `get_method_name`, `get_run_method` -- a method's reader, how a figure or table names it,
  and a run key's method;
- `get_bin_side`, `get_spot_labels` -- the side of the bins a run ran on, and its clone per spot of one side;
- `get_raw_barcodes`, `read_spots` -- UGP keys' Space Ranger barcodes, every in-tissue spot;
- `read_bin_counts` -- UGP's UMI and phased allele counts per bin and spot;
- `get_spot_sums`, `get_spot_umis` -- per-spot sums, and a UGP `bb` directory's total and SNP-covering UMI per
  spot;
- `read_broad_cell_types` -- each spot's broad cell type;
- `get_epithelial`, `get_normal_type` -- whether a spot's cell type is epithelial, or of a normal type;
- `get_benchmark_dir`, `get_benchmark_path`, `get_yaml_paths` -- the benchmark's directory, and its yamls' paths,
  relative to it or absolute.

Every method's normal clone is clone 0: cnamaste's by `hmrf.reindex_clones`; CalicoST's by `reorder_results`, "re-order
clones such that normal clones are always clone 0" (CalicoST c1abcae, `calicost_main.py:296`); Numbat's clone 1, "always
the normal cells" (https://kharchenkolab.github.io/numbat/articles/results.html), less one. Clonalscope's normal
clusters come first, clone 0 among them, and it may call several, or they are one clone 0 (`read_clonalscope`).
"""

from __future__ import annotations

import re
from collections.abc import Callable, Sequence
from functools import partial
from pathlib import Path
from typing import Any, NamedTuple

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scripts.platforms.visium import get_visiumhd_bins, get_visiumhd_finer_labels

from cnamaste.annotation import load_celltype
from cnamaste.calicost.read import RUN, best_run, read_estimate as read_calicost_estimate, seglevel as calicost_seglevel
from cnamaste.config import start_time
from cnamaste.logger import get_logger
from cnamaste.qa.metrics import read_labels
from cnamaste.qa.plot import NAMES as VALIDATION_NAMES
from cnamaste.results import get_output_path

logger = get_logger(__name__, start_time=start_time)


class MethodInput(NamedTuple):
    """One method's run, its spot labels and clone profile.

    Attributes:
        labels: the clone per barcode, the normal clone 0, at the side of the bins it ran on.
        profile: `CHR`, `START`, `END` per segment, then per clone `clone<k> A` and `clone<k> B` integer copies, or
            `clone<k>` a state.
    """

    labels: pd.Series
    profile: pd.DataFrame


class BinCounts(NamedTuple):
    """UMI and phased allele counts per bin and spot (`read_bin_counts`).

    Attributes:
        bins: `CHR`, `START`, `END` per row of the matrices.
        x_counts: `(bins, spots)` UMI.
        b_counts: `(bins, spots)` B-allele counts.
        totals: `(bins, spots)` A- and B-allele counts.
        spots: the Space Ranger barcode of each column.
    """

    bins: pd.DataFrame
    x_counts: sp.csr_matrix
    b_counts: sp.csr_matrix
    totals: sp.csr_matrix
    spots: pd.Index


CLONE_FIELD = re.compile(r"clone\d+ ")
"""A profile's column of one clone's field, `clone<k> <field>`, e.g. `clone1 A`."""

CALICOST_FALLBACKS = ("diploid", "tetraploid")
"""The tables tried, in order, where a CalicoST run's table leaves out a labelled clone (`get_calicost_ploidy`)."""

NAMES = {**VALIDATION_NAMES, "numbat": r"$\mathtt{Numbat}$", "clonalscope": r"$\mathtt{Clonalscope}$"}
"""How a figure names a method: as the validation figure's key does (`qa.plot.NAMES`), and Numbat and Clonalscope."""

EPITHELIAL = "Epithelial"
"""The broad cell type read as epithelial, a `cell_type`'s class before its `(<fine>)`."""

LOW_CONFIDENCE = "Low confidence"
"""The broad class of a spot the cell typing could not call, left out."""

NUMBAT_STATE_COLOURS = {
    "neu": "#FFFFFF",
    "bamp": "#FA8072",
    "bdel": "#87CEEB",
    "amp": "#CD4F39",
    "loh": "#4169E1",
    "del": "#698B22",
}
"""Each CNV state's colour on the profiles: Numbat's clone profile's, `cnv_heatmap`'s fill (R/vis.R), `neu` white, the
rest its `cnv_colors`, salmon, skyblue, tomato3, olivedrab4 and royalblue, as R gives them in hex, in `cnv_colors`'
order, with `loh` and `del` swapped: `loh` royalblue, `del` olivedrab4."""

NUMBAT_STATE_LABELS = {"neu": "NEU", "bamp": "BAMP", "bdel": "BDEL", "amp": "AMP", "loh": "CNLoH", "del": "DEL"}
"""Each CNV state's name in Numbat's key, `cnv_labels` (R/vis.R): upper case, `loh` written CNLoH."""

BENCHMARK = Path(__file__).resolve().parents[2]
"""The benchmark's directory, `scripts/metrics/`' grandparent: a relative path of its yamls is under it."""

YAML_PATHS = ("dir", "labels", "profile", "spots", "spaceranger_dir", "bb_dir", "runtime", "groundtruth_dir")
"""The benchmark yaml's keys whose values are paths: the value, or each value of the mapping it holds, e.g. `bb_dir`'s
per bin side."""


def get_benchmark_dir() -> Path:
    """Get the benchmark's directory, `BENCHMARK`: it holds `scripts/`, `configs/`, `data/` and `input_data/`."""
    return BENCHMARK


def get_benchmark_path(path: str) -> str:
    """Get a path of the benchmark's yamls: as given where absolute, else under the benchmark's directory.

    Args:
        path: the path, absolute or relative to `BENCHMARK`; a template's `{sample_id}` and `run_*` are kept.

    Returns:
        The absolute path.
    """
    return str(BENCHMARK / path)


def get_yaml_paths(entry: dict[str, Any]) -> dict[str, Any]:
    """Get a benchmark yaml entry, its paths under the benchmark's directory (`get_benchmark_path`) at any depth.

    A path is the string value of a `YAML_PATHS` key, or each string value of the mapping it holds; every other value
    is kept, e.g. a run's `labels: null` or a cnamaste run's `solver`.

    Args:
        entry: the entry, e.g. a sample's or the `simulation` entry.

    Returns:
        The entry, its paths absolute.
    """
    resolved: dict[str, Any] = {}
    for key, value in entry.items():
        if isinstance(value, dict) and key in YAML_PATHS:
            resolved[key] = {k: get_benchmark_path(v) if isinstance(v, str) else v for k, v in value.items()}
        elif isinstance(value, dict):
            resolved[key] = get_yaml_paths(value)
        elif isinstance(value, str) and key in YAML_PATHS:
            resolved[key] = get_benchmark_path(value)
        else:
            resolved[key] = value
    return resolved


def get_run_method(key: str) -> str:
    """Get a run key's method: the key before its first `_` or `-`, e.g. `calicost` of `calicost_016um_msu10`,
    `cnamaste` of `cnamaste-ac88605f`.

    Args:
        key: the run key.

    Returns:
        The method.
    """
    return re.split(r"[_-]", key, maxsplit=1)[0]


def get_method_name(method: str) -> str:
    """Get how a figure or table names a method (`NAMES`), mathtext that LaTeX also reads; the method as given where
    `NAMES` has none."""
    return NAMES.get(method, method)


def read_cnamaste(cnamaste_dir: Path, integer_profile: bool = False) -> MethodInput:
    """Read a cnamaste run: its integer clones, `clone_labels_int.tsv[.gz]`, and its segment-level copies.

    The copies are its `cnv_<ploidy>_seglevel.tsv` (`results.get_output_path`), else its `cnv_seglevel.tsv`, every
    clone's, or with `integer_profile` the clones its integer clones name alone (`get_labelled_profile`).

    Args:
        cnamaste_dir: a `run_cnamaste` output directory.
        integer_profile: whether the copies are those of the integer clones alone.

    Returns:
        The `MethodInput`.

    Raises:
        FileNotFoundError: where it holds no labels, or no or several segment tables.
    """
    cnamaste_dir = Path(cnamaste_dir)
    labels = read_labels(cnamaste_dir / "clone_labels_int.tsv").astype(int)
    try:
        seglevel = get_output_path(cnamaste_dir, "seglevel")
    except FileNotFoundError:
        seglevel = cnamaste_dir / "cnv_seglevel.tsv"
        if not seglevel.is_file():
            raise
    logger.info(f"Read cnamaste {cnamaste_dir}: {len(labels):_} spots, clones {sorted(set(labels))}, {seglevel.name}.")
    profile = pd.read_csv(seglevel, sep="\t")
    return MethodInput(labels, get_labelled_profile(profile, labels) if integer_profile else profile)


def get_labelled_profile(profile: pd.DataFrame, labels: pd.Series) -> pd.DataFrame:
    """Get a cnamaste profile's columns of the clones its integer clones name, and its other columns.

    run_cnamaste merges clones whose `(A, B)` agree into integer clones (`hmm_ab_decoder.integer_clone_groups`), and
    "Each names its group by a member clone, whose copies it calls" (`qa.metrics.merged_at_barcode`); its
    `cnv_<ploidy>_seglevel.tsv` keeps every clone, so a clone merged into another's group holds copies no spot carries.

    Args:
        profile: `clone<k> <field>` per clone, and other columns, e.g. a `cnv_<ploidy>_seglevel.tsv`.
        labels: the integer clone per spot, `clone_labels_int.tsv`.

    Returns:
        `profile` without the columns of a clone `labels` does not carry.
    """
    carried = {f"clone{k}" for k in pd.unique(labels.astype(int))}
    dropped = sorted({c.split(" ", 1)[0] for c in profile.columns if CLONE_FIELD.match(c)} - carried)
    if dropped:
        logger.info(f"cnamaste clones {dropped} merged into others' integer clones, their copies left out.")
    return profile[[c for c in profile.columns if not CLONE_FIELD.match(c) or c.split(" ", 1)[0] in carried]]


def get_calicost_run(calicost_dir: Path) -> Path:
    """Get the CalicoST run a directory names: itself where it is one run, a `clone*_rectangle*_w*` directory or one
    holding `clone_labels.tsv`, else its best rectangle by `total_llf` (`calicost.read.best_run`).

    Args:
        calicost_dir: a CalicoST sample directory or one of its runs.

    Returns:
        The run's directory.

    Raises:
        FileNotFoundError: where it is no run and none of its runs holds a final fit.
    """
    calicost_dir = Path(calicost_dir)
    one_run = RUN.fullmatch(calicost_dir.name) or (calicost_dir / "clone_labels.tsv").is_file()
    return calicost_dir if one_run else best_run(calicost_dir)[1]


def get_calicost_table_clones(run: Path, ploidy: str) -> set[int] | None:
    """Get the clones a CalicoST run's table at `ploidy` holds, its `clone<k> A` columns.

    Args:
        run: the run's directory.
        ploidy: one of `CALICOST_PLOIDIES`.

    Returns:
        The clones, or None where the run has no table at `ploidy`.
    """
    try:
        table = calicost_seglevel(run, ploidy)
    except FileNotFoundError:
        return None
    return {int(m.group(1)) for c in table.columns if (m := re.fullmatch(r"clone(\d+) A", str(c)))}


def get_calicost_fallbacks() -> tuple[str, ...]:
    """Get the tables tried, in order, where a CalicoST run's table leaves out a labelled clone: `diploid`, then
    `tetraploid`."""
    return CALICOST_FALLBACKS


def get_calicost_ploidy(run: Path, ploidy: str = "unconstrained", fallbacks: Sequence[str] = CALICOST_FALLBACKS) -> str:
    """Get the table a CalicoST run is read at: `ploidy`'s where it holds every tumour clone of `clone_labels.tsv`,
    else the first of `fallbacks` that does, else `ploidy`'s.

    CalicoST's unconstrained table leaves out a clone with no diploid balanced state, a state holding 10% of its bins,
    2% on the retry, its BAF within 0.05 of 0.5: `find_diploid_balanced_state` raises and a bare `except` skips the
    clone (CalicoST c1abcae, `find_integer_copynumber.py:73-81`, `calicost_main.py:334-340`). Its diploid, triploid and
    tetraploid tables (`hill_climbing_integer_copynumber_oneclone`, `calicost_main.py:331`) need no such state.

    Args:
        run: the run's directory (`get_calicost_run`).
        ploidy: one of `CALICOST_PLOIDIES`, the table read where it holds every labelled clone, or where no fallback
            does.
        fallbacks: `CALICOST_PLOIDIES` tried in order; empty for none.

    Returns:
        `ploidy` or one of `fallbacks`.
    """
    labelled = set(pd.read_csv(run / "clone_labels.tsv", sep="\t")["clone_label"].astype(int)) - {0}
    if labelled <= (get_calicost_table_clones(run, ploidy) or set()):
        return ploidy
    for fallback in fallbacks:
        if labelled <= (get_calicost_table_clones(run, fallback) or set()):
            return fallback
    return ploidy


def read_calicost(
    calicost_dir: Path, ploidy: str = "unconstrained", fallbacks: Sequence[str] = CALICOST_FALLBACKS
) -> MethodInput:
    """Read a CalicoST run: its `clone_labels.tsv` and a `cnv<_ploidy>_seglevel.tsv` (`calicost.read.read_estimate`),
    `ploidy`'s, or the first of `fallbacks` holding every labelled tumour clone where `ploidy`'s leaves one out
    (`get_calicost_ploidy`).

    Args:
        calicost_dir: a CalicoST sample directory or one of its runs (`get_calicost_run`).
        ploidy: one of `CALICOST_PLOIDIES`; `unconstrained` is its `cnv_seglevel.tsv`, which `plots/acn_genome.pdf`
            draws.
        fallbacks: `CALICOST_PLOIDIES` tried in order; empty to read `ploidy`'s table alone.

    Returns:
        The `MethodInput`.

    Raises:
        FileNotFoundError: where no run holds a final fit, or the run no table at the ploidy chosen.
    """
    run = get_calicost_run(calicost_dir)
    chosen = get_calicost_ploidy(run, ploidy, fallbacks)
    estimate = read_calicost_estimate(run, chosen)
    labelled = set(estimate.int_clones.astype(int)) - {0}
    if chosen != ploidy:
        logger.info(f"CalicoST {run}: its {ploidy} table leaves out a labelled clone, read at {chosen}.")
    elif not labelled <= (get_calicost_table_clones(run, ploidy) or set()):
        logger.warning(f"CalicoST {run}: none of {', '.join([ploidy, *fallbacks])} holds every labelled clone.")
    logger.info(f"Read CalicoST {run}: {len(estimate.int_clones):_} spots, its {chosen} copies.")
    return MethodInput(estimate.int_clones, estimate.seglevel)


def get_numbat_iteration(numbat_dir: Path) -> int:
    """Get the last iteration of phylogeny optimization a Numbat output directory holds.

    "The file names are post-fixed with the `i`th iteration of phylogeny optimization" (descriptions.html).

    Args:
        numbat_dir: a Numbat output directory.

    Returns:
        The largest `i` of its `clone_post_{i}.tsv`.

    Raises:
        FileNotFoundError: where it holds none.
    """
    found = [
        int(m.group(1)) for p in Path(numbat_dir).iterdir() if (m := re.fullmatch(r"clone_post_(\d+)\.tsv", p.name))
    ]
    if not found:
        raise FileNotFoundError(f"no clone_post_<i>.tsv in {numbat_dir}")
    return max(found)


def read_numbat_labels(numbat_dir: Path) -> pd.Series:
    """Read a Numbat run's clone per barcode, its normal clone 0.

    At its `get_numbat_iteration` `i`: `clone_post_{i}.tsv`'s `clone_opt`, "Maximum a posteriori clone assignment",
    less one, by `cell`. Numbat's "clone 1 is always the normal cells" (results.html), so the normal clone is clone 0,
    as cnamaste's is. Numbat's own tumor call, `compartment_opt`, is `p_cnv > 0.5`, `p_cnv` the posterior summed over
    the tumor clones (`get_clone_post`, R/main.R): a barcode whose most probable clone is the normal one, at below 0.5,
    is a tumor call there and clone 0 here.

    Args:
        numbat_dir: a Numbat output directory.

    Returns:
        The clone per barcode, Numbat's barcodes.
    """
    i = get_numbat_iteration(numbat_dir)
    post = pd.read_csv(Path(numbat_dir) / f"clone_post_{i}.tsv", sep="\t")
    labels = pd.Series(post["clone_opt"].astype(int).to_numpy() - 1, index=post["cell"].astype(str).to_numpy())
    logger.info(f"Read Numbat iteration {i} labels: {len(labels):_} barcodes, clones {sorted(set(labels))}.")
    return labels


def read_numbat_profile(numbat_dir: Path) -> pd.DataFrame:
    """Read a Numbat run's clone profile: each clone's CNV state per consensus segment.

    From `bulk_clones_final.tsv.gz`, "Clone-level pseudobulk profiles based on final cell lineage tree"
    (descriptions.html): each clone's `cnv_state` per consensus segment (`seg_cons`), a clone its `sample`, numbered as
    `read_numbat_labels` numbers them.

    Args:
        numbat_dir: a Numbat output directory.

    Returns:
        `CHR`, `START`, `END` per segment, in genome order, then `clone<k>` per clone, its state.

    Raises:
        ValueError: where a clone carries several states on a segment, or a state has no `NUMBAT_STATE_COLOURS`.
    """
    columns = ["sample", "CHROM", "seg_cons", "seg_start", "seg_end", "cnv_state"]
    bulk = pd.read_csv(Path(numbat_dir) / "bulk_clones_final.tsv.gz", sep="\t", usecols=columns, low_memory=False)
    segment = ["CHROM", "seg_cons", "seg_start", "seg_end"]
    if (bulk.groupby(["sample", *segment])["cnv_state"].nunique() > 1).any():
        raise ValueError("bulk_clones_final.tsv.gz: a clone carries several states on one segment")
    unknown = set(bulk["cnv_state"]) - set(NUMBAT_STATE_COLOURS)
    if unknown:
        raise ValueError(f"bulk_clones_final.tsv.gz: states {sorted(unknown)} not in NUMBAT_STATE_COLOURS")

    states = bulk.drop_duplicates(["sample", *segment]).pivot(index=segment, columns="sample", values="cnv_state")
    states = states.rename(columns=lambda k: f"clone{int(k) - 1}").rename_axis(columns=None).reset_index()
    states = states.sort_values(["CHROM", "seg_start"]).drop(columns="seg_cons")
    logger.info(f"Read Numbat clone profile: {len(states)} segments, {states.shape[1] - 3} clones.")
    return states.rename(columns={"CHROM": "CHR", "seg_start": "START", "seg_end": "END"}).reset_index(drop=True)


def read_numbat(numbat_dir: Path) -> MethodInput:
    """Read a Numbat run: its clones (`read_numbat_labels`) and profile (`read_numbat_profile`).

    Args:
        numbat_dir: a Numbat output directory.

    Returns:
        The `MethodInput`.
    """
    return MethodInput(read_numbat_labels(numbat_dir), read_numbat_profile(numbat_dir))


def get_total_copy_state(total: pd.Series) -> pd.Series:
    """Get total copies' CNA states: `del` below 2, `neu` at 2, `amp` above, as copytyping-benchmark's
    `analysis/scripts/cn_transform.py` (1204da9) reads Clonalscope's total copies.

    Args:
        total: total copies, NA where none is called.

    Returns:
        The state per value, NA where `total` is.
    """
    states = pd.Series(pd.NA, index=total.index, dtype="object")
    states[total < 2] = "del"
    states[total == 2] = "neu"
    states[total > 2] = "amp"
    return states


def read_clonalscope(clonalscope_dir: Path, merge_normal: bool = False) -> MethodInput:
    """Read a Clonalscope run's export: its clusters, and each cluster's total copies per segment as a CNA state.

    The export, of Clonalscope's Visium tutorial run (copytyping-benchmark's `export_tutorial`), holds one
    `<sample>.<assay>.cells.tsv`, `BARCODE`, `malignant` and `cluster` per spot, and one `<sample>.<assay>.profile.tsv`,
    per cluster `CLUSTER` and `MALIGNANT`, then per segment `#CHR`, `START`, `END`, `REL_CN`, its `Uest`, and `TCN`,
    its total copies, read as states by `get_total_copy_state`. `REL_CN` is the cluster's column of the export's
    `cluster_states.tsv`, `t(AssignCluster Uest)` (copytyping-benchmark 1204da9, `analysis/scripts/cn_transform.py`),
    and `TCN` is `round(2 * REL_CN)`, as `three_state_clonalscope` decodes scRNA. Clusters are numbered normal first,
    `malignant` N, then tumor, each in cluster order: clone 0 is a normal cluster, and Clonalscope may call several,
    clones 0 to m - 1. `malignant` is Clonalscope's `annot`, "T: malignant cells and N: non-malingnat cells"
    (Clonalscope 07dd0c5, `R/AssignCluster.R:14`).

    With `merge_normal`, every normal cluster is clone 0 and the tumor clusters are clones 1 to t, in cluster order;
    the profile holds the tumor clones alone, Clonalscope profiling each normal cluster on its own.

    Args:
        clonalscope_dir: a Clonalscope export directory, e.g. `postprocess/visium_tutorial/<sample>/scRNA`.
        merge_normal: whether the normal clusters are one clone, clone 0.

    Returns:
        The `MethodInput`, its profile per cluster `clone<k>`, its state, `clone<k> Uest` and `clone<k> tcn`.

    Raises:
        FileNotFoundError: where the directory holds no `*.cells.tsv` or `*.profile.tsv`, or several.
    """
    found = {kind: sorted(Path(clonalscope_dir).glob(f"*.{kind}.tsv")) for kind in ("cells", "profile")}
    for kind, paths in found.items():
        if len(paths) != 1:
            raise FileNotFoundError(f"{len(paths)} *.{kind}.tsv in {clonalscope_dir}, not one")
    cells = pd.read_csv(found["cells"][0], sep="\t")
    profile = pd.read_csv(found["profile"][0], sep="\t")

    clusters = profile.drop_duplicates("CLUSTER").sort_values(["MALIGNANT", "CLUSTER"])
    normal = clusters.loc[clusters["MALIGNANT"] == "N", "CLUSTER"].tolist()
    if merge_normal:
        tumor = clusters.loc[clusters["MALIGNANT"] != "N", "CLUSTER"]
        clone = {**dict.fromkeys(normal, 0), **{cluster: k for k, cluster in enumerate(tumor, start=1)}}
        profile = profile[profile["MALIGNANT"] != "N"]
    else:
        clone = {cluster: k for k, cluster in enumerate(clusters["CLUSTER"])}
    labels = pd.Series(cells["cluster"].map(clone).to_numpy(), index=cells["BARCODE"].astype(str).to_numpy())
    if labels.isna().any():
        raise ValueError(f"{found['cells'][0]}: clusters {sorted(set(cells['cluster']) - set(clone))} have no profile")

    states = profile.assign(state=get_total_copy_state(profile["TCN"]), clone=profile["CLUSTER"].map(clone))
    states = states.rename(columns={"REL_CN": "Uest", "TCN": "tcn"})
    states = states.pivot(index=["#CHR", "START", "END"], columns="clone", values=["state", "Uest", "tcn"])
    states.columns = [f"clone{k}" if field == "state" else f"clone{k} {field}" for field, k in states.columns]
    order = [f"clone{k}{field}" for k in sorted(set(profile["CLUSTER"].map(clone))) for field in ("", " Uest", " tcn")]
    numbers = {c: "float64" for c in order if c.endswith((" Uest", " tcn"))}
    states = states[order].astype(numbers).reset_index()
    # NB the export writes some coordinates in scientific notation, e.g. `1e+08`.
    states = states.rename(columns={"#CHR": "CHR"}).astype({"START": "int64", "END": "int64"})
    states = states.sort_values(["CHR", "START"]).reset_index(drop=True)
    logger.info(
        f"Read Clonalscope {clonalscope_dir}: {len(labels):_} spots, {len(clusters)} clusters ({len(normal)} normal"
        f"{', as clone 0' if merge_normal else ''}), {len(states)} segments."
    )
    return MethodInput(labels.astype(int), states)


def get_method_reader(
    method: str, merge_normal: bool = False, integer_profile: bool = False
) -> Callable[[Path], MethodInput]:
    """Get a method's reader of a run directory: `read_numbat`, `read_calicost` (unconstrained), `read_cnamaste` or
    `read_clonalscope`, its normal clusters one clone 0 with `merge_normal`, cnamaste's copies its integer clones'
    alone with `integer_profile`.

    Args:
        method: `numbat`, `calicost`, `cnamaste` or `clonalscope`.
        merge_normal: whether Clonalscope's normal clusters are one clone 0; every other method has one normal clone.
        integer_profile: whether cnamaste's copies are its integer clones' alone (`read_cnamaste`).

    Returns:
        The reader.

    Raises:
        ValueError: where the method has none.
    """
    readers: dict[str, Callable[[Path], MethodInput]] = {
        "numbat": read_numbat,
        "calicost": read_calicost,
        "cnamaste": partial(read_cnamaste, integer_profile=integer_profile),
        "clonalscope": partial(read_clonalscope, merge_normal=merge_normal),
    }
    if method not in readers:
        raise ValueError(f"no reader of {method!r} runs: one of {', '.join(readers)}")
    return readers[method]


class Run(NamedTuple):
    """One method's run.

    Attributes:
        method: a key of `NAMES`: `numbat`, `calicost`, `clonalscope` or `cnamaste`.
        label: the run's label, empty where the method has one run.
        directory: the run's directory.
    """

    method: str
    label: str
    directory: Path


def get_run_argument(text: str) -> tuple[str, Path]:
    """Get a run's label and directory from its command-line value, `[<label>=]<dir>`.

    Args:
        text: the value.

    Returns:
        The label, empty where none is given, and the directory.
    """
    label, sep, directory = text.partition("=")
    return (label, Path(directory)) if sep and "/" not in label else ("", Path(text))


def get_run_name(run: Run) -> str:
    """Get a run's name in a log or a table: its method, then its label where it has one.

    Args:
        run: the run.

    Returns:
        The name.
    """
    return f"{run.method} {run.label}".strip()


class SpotLabels(NamedTuple):
    """One run's clones on the spots: a run is its name (`get_run_name`) at its side.

    Attributes:
        run: the run.
        labels: the clone per spot of the side asked for.
        run_side_um: the side in um of the bins the run ran on.
    """

    run: Run
    labels: pd.Series
    run_side_um: int


def get_bin_side(labels: pd.Series) -> int:
    """Get the side in um of the bins a run ran on, from its first labelled barcode.

    Args:
        labels: its clone per barcode, at the bins it ran on, e.g. its `MethodInput.labels`.

    Returns:
        The side in um.
    """
    return int(get_visiumhd_bins(labels.dropna().index)[0, 0])


def get_spot_labels(run: Run, labels: pd.Series, side_um: int) -> SpotLabels:
    """Get a run's clone per spot of side `side_um`, from its clone per bin of the bins it ran on.

    A coarser bin's clone is given to every spot inside it (`platforms.visium.get_visiumhd_finer_labels`);
    a spot of no clone is left out.

    Args:
        run: the run.
        labels: its clone per barcode, at the bins it ran on, e.g. its `MethodInput.labels`.
        side_um: the spots' side in um.

    Returns:
        The run's `SpotLabels`.

    Raises:
        ValueError: where the run's bins are of several sides, or `side_um` does not divide theirs.
    """
    labels = labels.dropna()
    return SpotLabels(run, get_visiumhd_finer_labels(labels, side_um), get_bin_side(labels))


def get_raw_barcodes(keys: pd.Series) -> pd.Index:
    """Get each UGP key's Space Ranger barcode: `raw` of `{raw}_{dataset_id}_{assay_type}`.

    Read as UGP reads it, its last two `_` bounding the other two fields (`io_utils.read_barcodes_by_dataset`,
    Universal-Genotyping-Pipeline c5c0b9b).

    Args:
        keys: the UGP keys, e.g. a `barcodes.tsv.gz`'s lines or a `bb` AnnData's `obs_names`.

    Returns:
        The Space Ranger barcodes, in `keys`' order.

    Raises:
        ValueError: where a key has fewer than three fields.
    """
    parts = pd.Series(keys, dtype=str).str.rsplit("_", n=2, expand=True)
    if parts.shape[1] != 3 or parts.isna().any().any():
        raise ValueError("every UGP key is {raw}_{dataset_id}_{assay_type}")
    return pd.Index(parts[0])


def read_spots(barcodes_file: Path) -> pd.Index:
    """Read every in-tissue spot from a UGP `barcodes.tsv.gz`, e.g. `bb/<MSR>/<assay>/barcodes.tsv.gz`.

    One `{raw}_{dataset_id}_{assay_type}` per line, read back to its Space Ranger barcode (`get_raw_barcodes`). UGP
    keeps Space Ranger's filtered barcodes, the bins under tissue, and orders a `bb` directory's allele matrices'
    columns as this file's lines (`combine_counts_nonbulk.py`, Universal-Genotyping-Pipeline c5c0b9b).

    Args:
        barcodes_file: the barcode file, plain or `.gz`, no header.

    Returns:
        The spots' Space Ranger barcodes, in file order.

    Raises:
        ValueError: where a line has fewer than three fields.
    """
    spots = get_raw_barcodes(pd.read_csv(barcodes_file, sep="\t", header=None, dtype=str).iloc[:, 0])
    logger.info(f"Read {len(spots):_} in-tissue spots from {barcodes_file}.")
    return spots


def read_bin_counts(bin_dir: Path) -> BinCounts:
    """Read UGP's counts per bin and spot, a `bb/MSR<min_snp_reads>/<assay>/` directory: the bins, `bb.tsv.gz`; per bin
    and spot, its UMI, `bb.Xcount.npz`, and its B-allele and total allele counts, `bb.Ballele.npz` and `bb.Tallele.npz`;
    the spots, `barcodes.tsv.gz` (`read_spots`), in the matrices' column order.

    UGP bins the SNPs adaptively, each bin at least `min_snp_reads` SNP reads (`build_adaptive_bins`,
    `combine_counts_nonbulk.py:367-375`), and sums into each bin its SNPs' phased A and B counts (`sum_features_to_bbs`,
    `combine_counts_nonbulk.py:423-425`; the phase, `phase_and_concat_nonbulk.py:197-202`) and the UMI of each gene it
    overlaps most (`feature_utils.sum_umis_to_bins`), Universal-Genotyping-Pipeline c5c0b9b. A bin with no SNP has
    allele counts 0.

    Args:
        bin_dir: the directory, e.g. `bb/MSR300/VISIUM3prime`.

    Returns:
        The counts.

    Raises:
        ValueError: where a matrix's shape is not the bins by the spots.
    """
    bin_dir = Path(bin_dir)
    bins = pd.read_csv(bin_dir / "bb.tsv.gz", sep="\t", usecols=["#CHR", "START", "END"]).rename(
        columns={"#CHR": "CHR"}
    )
    x_counts, b_counts, totals = (
        sp.csr_matrix(sp.load_npz(bin_dir / f"bb.{n}.npz")) for n in ("Xcount", "Ballele", "Tallele")
    )
    spots = read_spots(bin_dir / "barcodes.tsv.gz")
    for matrix in (x_counts, b_counts, totals):
        if matrix.shape != (len(bins), len(spots)):
            raise ValueError(f"{bin_dir}: a matrix of {matrix.shape}, not {len(bins)} bins by {len(spots)} spots")
    logger.info(
        f"Read {bin_dir}: {len(bins):_} bins, {len(spots):_} spots, {int(x_counts.sum()):_} UMI, "
        f"{int(totals.sum()):_} allele counts."
    )
    return BinCounts(bins, x_counts, b_counts, totals, spots)


def get_spot_sums(counts: Any, spot_axis: int = 0) -> np.ndarray:
    """Sum `counts` over every axis but the spots'.

    Args:
        counts: a 2-D count matrix, sparse or dense.
        spot_axis: the axis indexing spots.

    Returns:
        The per-spot sums, int64 for integer counts.
    """
    sums = np.asarray(counts.sum(axis=1 - spot_axis)).ravel()
    return sums.astype(np.int64) if np.all(np.mod(sums, 1) == 0) else sums


def get_spot_umis(bb_dir: Path, assay: str, barcodes_file: Path) -> pd.DataFrame:
    """Get each in-tissue spot's total UMI and SNP-covering UMI from a UGP `bb` directory.

    The total UMI is the spot's sum over `<assay>.h5ad`'s genes, matched by its key; the SNP-covering UMI is A + B of
    `multi_snp/<assay>/bb.{A,B}allele.npz`, whose columns are `barcodes_file`'s lines in order
    (`read_method_input.read_spots`).

    Args:
        bb_dir: a UGP `bb` directory.
        assay: UGP's assay type, the stem of the AnnData and the allele matrices' directory.
        barcodes_file: the directory's `barcodes.tsv.gz`, `<MSR>/<assay>/barcodes.tsv.gz`.

    Returns:
        `total_umi` and `snp_umi` per Space Ranger barcode, `barcodes_file`'s order.

    Raises:
        ValueError: where the AnnData, the allele matrices and `barcodes_file` hold different spots.
    """
    spots = read_spots(barcodes_file)
    logger.info(f"Reading {bb_dir / f'{assay}.h5ad'}.")
    adata = anndata.read_h5ad(bb_dir / f"{assay}.h5ad")
    total = pd.Series(get_spot_sums(adata.X, spot_axis=0), index=get_raw_barcodes(pd.Series(adata.obs_names)))
    snp_dir = bb_dir / "multi_snp" / assay
    a_allele, b_allele = sp.load_npz(snp_dir / "bb.Aallele.npz"), sp.load_npz(snp_dir / "bb.Ballele.npz")
    if a_allele.shape[1] != len(spots) or not total.index.sort_values().equals(spots.sort_values()):
        raise ValueError(f"{bb_dir}: its AnnData, allele matrices and {barcodes_file} hold different spots")
    snp = get_spot_sums(a_allele, spot_axis=1) + get_spot_sums(b_allele, spot_axis=1)
    return pd.DataFrame({"total_umi": total.reindex(spots).to_numpy(), "snp_umi": snp}, index=spots)


def read_broad_cell_types(celltype_file: Path) -> pd.Series:
    """Read each spot's broad cell type: its `cell_type` before the ` (<fine>)` (`annotation.load_celltype`).

    Args:
        celltype_file: per-spot cell types, `barcode` and `cell_type`, `<broad> (<fine>)`.

    Returns:
        The broad cell type per barcode, e.g. `EPITHELIAL`, `LOW_CONFIDENCE`.
    """
    celltype = load_celltype(celltype_file)
    return pd.Series(celltype["cell_type"].str.split(" (", regex=False).str[0].to_numpy(), index=celltype["barcode"])


def get_normal_type(broad: pd.Series) -> pd.Series:
    """Get whether each spot's broad cell type is a normal type: neither `EPITHELIAL` nor `LOW_CONFIDENCE`.

    The cell types give normal truth only: an epithelial spot may be tumor or normal, a spot of another type is normal.

    Args:
        broad: the broad cell type per barcode (`read_broad_cell_types`).

    Returns:
        True or False per barcode.
    """
    return ~broad.isin([EPITHELIAL, LOW_CONFIDENCE])


def get_epithelial(broad: pd.Series) -> pd.Series:
    """Get whether each spot's broad cell type is `EPITHELIAL`, a `LOW_CONFIDENCE` spot left out.

    Args:
        broad: the broad cell type per barcode (`read_broad_cell_types`).

    Returns:
        True or False per barcode.
    """
    epithelial = broad[broad != LOW_CONFIDENCE] == EPITHELIAL
    n_low = int((broad == LOW_CONFIDENCE).sum())
    logger.info(
        f"Cell types: {int(epithelial.sum()):_} epithelial, {int((~epithelial).sum()):_} not, "
        f"{n_low:_} {LOW_CONFIDENCE.lower()} left out."
    )
    return epithelial
