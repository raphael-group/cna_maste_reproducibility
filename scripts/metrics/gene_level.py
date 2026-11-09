"""The reference's genes inside the region bed, and each run's clone profile at them: a gene takes the segment that
overlaps it most.

- `read_reference_files` -- the region bed and the genes' GTF a benchmark config names;
- `read_genes` -- the reference's genes, from its GTF;
- `read_regions` -- the region bed;
- `get_inside_regions`, `get_region_genes` -- whether each interval lies inside an autosomal region, and the genes
  that do;
- `get_max_overlap` -- per interval, the target interval that overlaps it most;
- `get_gene_profile` -- a run's profile at the genes it calls, a row per gene, for the segment-level metrics;
- `get_truth_input` -- a simulation's truth as a run;
- `get_interval_fields`, `get_gene_states` -- a run's fields per interval (a gene, a joint segment), or its states per
  gene, `uncalled` where no segment overlaps it;
- `get_gene_level_table` -- the truth's and each run's fields per gene, a `<sample>.gene_level_cnp.tsv`.

Coordinates are 0-based and half-open, `[START, END)`, as BED's: a GTF's 1-based, inclusive `[start, end]` is
`[start - 1, end)`. Chromosomes are compared without their `chr` prefix.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scripts.metrics.cna_states import get_field_table, get_uncalled
from scripts.metrics.read_method_input import MethodInput, get_benchmark_path

from cnamaste.config import start_time
from cnamaste.logger import get_logger
from cnamaste.truth import Truth

logger = get_logger("gene_level", start_time=start_time)

AUTOSOMES = tuple(str(k) for k in range(1, 23))
"""The chromosomes kept, without `chr`: the simulations' truth and every method's input are autosomal."""

GENE_COLUMNS = ("gene_id", "gene", "CHR", "START", "END")
"""A gene's identity and interval, as cnamaste names its gene (`cnv_genelevel.tsv`) and segment columns."""


def get_chromosomes(values: pd.Series) -> np.ndarray:
    """Get chromosome names without their `chr` prefix, as strings.

    Args:
        values: chromosome names or numbers, e.g. `chr1` or `1`.

    Returns:
        The names, e.g. `1`.
    """
    return pd.Series(values).astype(str).str.removeprefix("chr").to_numpy()


def read_reference_files(config_yaml: Path) -> tuple[Path, Path]:
    """Read the region bed and the genes' GTF a benchmark config names.

    Args:
        config_yaml: a benchmark config, e.g. the pipeline's `<OUTDIR>/configs/config.yaml`: `reference.bed_file`,
            the region bed, and `reference.refdata.gene_list.gtf_file`, the reference's genes, each relative to the
            benchmark's directory or absolute (`read_method_input.get_benchmark_path`).

    Returns:
        The bed and the GTF.
    """
    with open(config_yaml) as handle:
        reference = yaml.safe_load(handle)["reference"]
    bed_file, gtf_file = reference["bed_file"], reference["refdata"]["gene_list"]["gtf_file"]
    return Path(get_benchmark_path(bed_file)), Path(get_benchmark_path(gtf_file))


def read_genes(gtf_file: Path) -> pd.DataFrame:
    """Read the reference's genes: each `gene` feature of its GTF.

    Args:
        gtf_file: a GTF, plain or `.gz`, e.g. Cell Ranger's `genes/genes.gtf.gz`.

    Returns:
        `GENE_COLUMNS`, a row per gene, in the GTF's order, `START` 0-based.
    """
    rows = []
    with gzip.open(gtf_file, "rt") if str(gtf_file).endswith(".gz") else open(gtf_file) as handle:
        for line in handle:
            if "\tgene\t" not in line or line.startswith("#"):
                continue
            chrom, _, feature, start, end, _, _, _, attributes = line.rstrip("\n").split("\t", 8)
            if feature != "gene":
                continue
            gene_id = re.search(r'gene_id "([^"]+)"', attributes)
            name = re.search(r'gene_name "([^"]+)"', attributes)
            rows.append(
                (gene_id.group(1) if gene_id else "", name.group(1) if name else "", chrom, int(start) - 1, int(end))
            )
    genes = pd.DataFrame(rows, columns=list(GENE_COLUMNS))
    genes["CHR"] = get_chromosomes(genes["CHR"])
    logger.info(f"Read {len(genes):_} genes from {gtf_file}.")
    return genes


def read_regions(bed_file: Path) -> pd.DataFrame:
    """Read the region bed: a region per line, `chrom`, `start`, `end`, then its name.

    Args:
        bed_file: the bed, e.g. `data/hg38.regions.bed`, a chromosome arm per line.

    Returns:
        `CHR` (no `chr`), `START`, `END`, a row per region.
    """
    regions = pd.read_csv(bed_file, sep="\t", header=None, comment="#", usecols=[0, 1, 2])
    regions.columns = ["CHR", "START", "END"]
    regions["CHR"] = get_chromosomes(regions["CHR"])
    return regions


def get_autosomal_regions(regions: pd.DataFrame) -> pd.DataFrame:
    """Get the regions on `AUTOSOMES`.

    Args:
        regions: `read_regions`' table.

    Returns:
        Its rows on chromosomes 1 to 22.
    """
    return regions[regions["CHR"].isin(AUTOSOMES)].reset_index(drop=True)


def get_inside_regions(intervals: pd.DataFrame, regions: pd.DataFrame) -> np.ndarray:
    """Get whether each interval lies wholly inside one region on `AUTOSOMES`.

    Args:
        intervals: `CHR`, `START`, `END` per interval, `CHR` with or without `chr`.
        regions: `read_regions`' table.

    Returns:
        True or False per interval.
    """
    regions = get_autosomal_regions(regions)
    chroms = get_chromosomes(intervals["CHR"])
    starts, ends = intervals["START"].to_numpy(), intervals["END"].to_numpy()
    inside = np.zeros(len(intervals), dtype=bool)
    for chrom, start, end in regions[["CHR", "START", "END"]].itertuples(index=False):
        inside |= (chroms == chrom) & (starts >= start) & (ends <= end)
    return inside


def get_region_genes(genes: pd.DataFrame, regions: pd.DataFrame) -> pd.DataFrame:
    """Get the genes inside a region, on `AUTOSOMES` (`get_inside_regions`).

    Args:
        genes: `read_genes`' table.
        regions: `read_regions`' table.

    Returns:
        The genes each wholly inside one region, indexed by `gene_id`, in genome order.
    """
    regions = get_autosomal_regions(regions)
    kept = (
        genes[get_inside_regions(genes, regions)]
        .assign(order=lambda d: d["CHR"].astype(int))
        .sort_values(["order", "START", "END"])
    )
    logger.info(f"Kept {len(kept):_} of {len(genes):_} genes inside {len(regions)} autosomal regions.")
    return kept.drop(columns="order").set_index("gene_id", drop=False).rename_axis(None)


def get_max_overlap(query: pd.DataFrame, target: pd.DataFrame) -> np.ndarray:
    """Get, per query interval, the target interval that overlaps it most.

    Ties go to the first in `target`'s order, as `qa.metrics.segment_true_ab_pairs` breaks them.

    Args:
        query: `CHR`, `START`, `END` per interval.
        target: `CHR`, `START`, `END` per interval.

    Returns:
        Per query interval, the target's row position, -1 where no target interval overlaps it by a base or more.
    """
    q_chrom, t_chrom = get_chromosomes(query["CHR"]), get_chromosomes(target["CHR"])
    q_lo, q_hi = query["START"].to_numpy(dtype=float), query["END"].to_numpy(dtype=float)
    t_lo, t_hi = target["START"].to_numpy(dtype=float), target["END"].to_numpy(dtype=float)
    best = np.full(len(query), -1, dtype=np.int64)
    for chrom in np.unique(q_chrom):
        qi, ti = np.flatnonzero(q_chrom == chrom), np.flatnonzero(t_chrom == chrom)
        if not ti.size:
            continue
        shared = np.minimum(q_hi[qi, None], t_hi[None, ti]) - np.maximum(q_lo[qi, None], t_lo[None, ti])
        found = shared.max(axis=1) > 0
        best[qi[found]] = ti[np.argmax(shared[found], axis=1)]
    return best


def get_gene_profile(profile: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    """Get a run's profile at the genes it calls: each gene the row of the segment overlapping it most.

    The metrics then score genes as they score segments, each gene once.

    Args:
        profile: a `MethodInput.profile`.
        genes: `get_region_genes`' table.

    Returns:
        The profile's columns, a row per gene some segment overlaps, `CHR`, `START` and `END` the gene's, indexed by
        `gene_id`.
    """
    best = get_max_overlap(genes, profile)
    called = best >= 0
    rows = profile.iloc[best[called]].reset_index(drop=True)
    rows[["CHR", "START", "END"]] = genes.loc[called, ["CHR", "START", "END"]].to_numpy()
    return rows.set_axis(genes.index[called])


def get_truth_input(truth: Truth) -> MethodInput:
    """Get a simulation's truth as a run: its clone per spot and AB pairs per true segment.

    Args:
        truth: the truth, its clones numbered as `truth.read_truth` numbers them, `normal` 0.

    Returns:
        The `MethodInput`, `clone<k> A`, `clone<k> B` per true clone `k`.
    """
    profile = truth.profile
    columns = {"CHR": profile["chr"], "START": profile["start"], "END": profile["end"]}
    for k, clone in enumerate(truth.clones):
        columns[f"clone{k} A"] = profile[f"{clone}_A_copy"]
        columns[f"clone{k} B"] = profile[f"{clone}_B_copy"]
    return MethodInput(pd.Series(truth.labels, index=truth.barcodes), pd.DataFrame(columns))


def get_interval_fields(profile: pd.DataFrame, intervals: pd.DataFrame) -> pd.DataFrame:
    """Get a run's fields per interval (`cna_states.get_field_table`), each interval its segment overlapping it most.

    Args:
        profile: a `MethodInput.profile`.
        intervals: `CHR`, `START`, `END` per interval, e.g. `get_region_genes`' table or the joint segments.

    Returns:
        `clone<k> <field>` per clone and field, a row per interval, `uncalled` where no segment overlaps it.
    """
    best = get_max_overlap(intervals, profile)
    fields = get_field_table(profile).iloc[np.maximum(best, 0)].set_axis(intervals.index)
    fields.loc[best < 0, :] = get_uncalled()
    return fields


def get_gene_states(profile: pd.DataFrame, genes: pd.DataFrame, field: str = "state") -> pd.DataFrame:
    """Get a run's call per gene and clone at one field (`cna_states.get_profile_fields`), its CNA state by default.

    Args:
        profile: a `MethodInput.profile`.
        genes: `get_region_genes`' table.
        field: the field, e.g. `state`, `tcn_state` or `ascn`.

    Returns:
        `clone<k>` per clone, a row per gene, the field's value or `uncalled`.
    """
    fields = get_interval_fields(profile, genes)
    states = fields[[c for c in fields.columns if c.endswith(f" {field}")]]
    return states.rename(columns=lambda c: c.removesuffix(f" {field}"))


def get_gene_level_table(
    genes: pd.DataFrame, runs: dict[str, MethodInput], truth: MethodInput | None = None
) -> pd.DataFrame:
    """Get the truth's and each run's fields per gene: a `<sample>.gene_level_cnp.tsv`.

    Args:
        genes: `get_region_genes`' table.
        runs: each run by its name, e.g. `numbat`.
        truth: the truth as a run (`get_truth_input`), or None.

    Returns:
        `GENE_COLUMNS`, then `truth clone<k> ascn` and `truth clone<k> state` per true clone, then `<name> clone<k>
        <field>` per run, clone and field (`cna_states.get_profile_fields`).
    """
    parts = [genes[list(GENE_COLUMNS)].reset_index(drop=True)]
    if truth is not None:
        fields = get_interval_fields(truth.profile, genes)
        fields = fields[[c for c in fields.columns if c.endswith((" ascn", " state"))]]
        parts.append(fields.add_prefix("truth ").reset_index(drop=True))
    for name, run in runs.items():
        parts.append(get_interval_fields(run.profile, genes).add_prefix(f"{name} ").reset_index(drop=True))
    return pd.concat(parts, axis=1)
