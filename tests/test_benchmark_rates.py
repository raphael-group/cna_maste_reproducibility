"""The benchmark's rates (`scripts`): a run's coverage of the in-tissue spots
(`metrics.qc_metrics.get_coverage`) and its calls per broad cell type and its false-positive rate
(`metrics.qc_metrics.get_false_positive_counts`), on its labels given to 16 um spots.

Referees: the counts made by hand on six spots; UGP's barcode key `{raw}_{dataset_id}_{assay_type}`
(Universal-Genotyping-Pipeline `io_utils.read_barcodes_by_dataset`).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scripts.metrics.qc_metrics import (
    get_coverage,
    get_false_positive_column,
    get_false_positive_counts,
)
from scripts.metrics.read_method_input import Run, SpotLabels, read_spots
from scripts.platforms.visium import get_visiumhd_finer_labels


def test_rates_count_in_tissue_spots_by_cell_type(tmp_path: Path) -> None:
    """Six in-tissue 16 um spots read back from UGP's keys, four in the 32 um bin (0, 0). A 32 um run labels that bin
    clone 1 and an off-tissue bin: four tumor calls, the stromal and immune ones false positives, the epithelial and
    low-confidence ones in neither count; its two normal-type calls both false positives. A 16 um run calls the
    epithelial spot and one stromal spot tumor, the other stromal spot normal: one false positive of two normal-type
    calls. A normal-type spot a run does not label is in neither count."""
    raw = [f"s_016um_{r:05d}_{c:05d}-1" for r, c in ((0, 0), (0, 1), (1, 0), (1, 1), (0, 2), (0, 3))]
    (tmp_path / "barcodes.tsv").write_text("".join(f"{b}_U1_VISIUM3prime\n" for b in raw))
    coarse = pd.Series([1, 2], index=["s_032um_00000_00000-1", "s_032um_00005_00005-1"])
    fine = pd.Series([2, 0, 1], index=[raw[0], raw[4], raw[5]])
    labels = [
        SpotLabels(Run("numbat", "", Path("coarse")), get_visiumhd_finer_labels(coarse, 16), 32),
        SpotLabels(Run("cnamaste", "", Path("fine")), get_visiumhd_finer_labels(fine, 16), 16),
    ]
    broad = pd.Series(["Epithelial", "Stromal", "Immune", "Low confidence", "Stromal", "Stromal"], index=raw)

    spots = read_spots(tmp_path / "barcodes.tsv")
    assert spots.tolist() == raw

    coverage = get_coverage(labels, spots)
    assert coverage[["resolution", "#called", "#normal-spots", "#tumor-spots"]].to_numpy().tolist() == [
        ["32um", 4, 0, 4],
        ["16um", 3, 1, 2],
    ]
    assert coverage["spot-coverage-rate"].tolist() == [4 / 6, 3 / 6]

    counts = get_false_positive_counts(labels, spots, broad)
    types = ["Epithelial", "Immune", "Stromal", "Low_confidence"]
    assert list(counts.columns[2:-3]) == [f"#{call}-spots-{t}" for t in types for call in ("normal", "tumor")]
    assert counts.iloc[:, 2:-3].to_numpy().tolist() == [[0, 1, 0, 1, 0, 1, 0, 1], [0, 1, 0, 0, 1, 1, 0, 0]]
    assert counts["#false-positive-tumor-spots"].tolist() == [2, 1]
    assert counts["#normal-type-spots"].tolist() == [2, 2]
    assert counts[get_false_positive_column()].tolist() == [2 / 2, 1 / 2]
