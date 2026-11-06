"""The simulations' state metrics (`scripts/metrics`): CNA states of AB pairs, and the cells
scored by state, on segments and at genes.

Referees: Numbat's file descriptions, each state's pair (https://kharchenkolab.github.io/numbat/articles/descriptions.html);
`cnamaste.qa.metrics.sensitivity`, `specificity` and `exact_altered` on the same run; the cells counted by hand on four
spots.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scripts.metrics.cna_states import get_ab_states
from scripts.metrics.gene_level import get_gene_profile
from scripts.metrics.read_method_input import MethodInput
from scripts.metrics.simulation_metrics import (
    get_level_table,
    get_metric_columns,
    get_run_dir,
    merge_levels,
)
from scripts.metrics.supervised_metrics import (
    compute_supervised_metrics,
    fitted_at_barcode,
    get_state_confusion,
    get_state_metrics,
)

from cnamaste.qa.metrics import RunEstimate, exact_altered, phase_free, sensitivity, specificity
from cnamaste.truth import Truth


def test_ab_states_are_numbats_pairs_either_phase() -> None:
    """Numbat's pairs, 1:1 neutral, 1:0 deletion, 2:0 CNLoH, 2:1 and 3:1 amplification, 2:2 and 0:0 biallelic, read
    as its states in either phase."""
    pairs = {
        (1, 1): "neu",
        (1, 0): "del",
        (2, 0): "cnloh",
        (2, 1): "amp",
        (3, 1): "amp",
        (2, 2): "bamp",
        (0, 0): "bdel",
    }
    for (a, b), state in pairs.items():
        for code in (a * 1_000 + b, b * 1_000 + a):
            assert get_ab_states(phase_free(np.array([code])))[0] == state


def test_state_metrics_count_qa_cells_and_states() -> None:
    """Two true segments, the first a (2, 1) gain in the tumor clone. Spots: two normal, two tumor; the run's clone 1
    calls (3, 1) on the first segment, so one normal spot in it is a false positive, one tumor spot in it its true
    state (amp) but not its pair, and the tumor spot in clone 0 a miss. Of 2 truly altered cells, 1 missed and 1 its
    state; of 6 truly neutral cells, 1 called altered. Sensitivity, specificity and N-CNA recall are qa.metrics', from
    AB pairs or from Numbat-like states; at one gene inside each segment, the genes score as the segments."""
    truth = Truth(
        barcodes=np.array(["n1", "n2", "t1", "t2"]),
        coords=np.zeros((4, 2)),
        labels=np.array([0, 0, 1, 1]),
        clones=("normal", "clone_1"),
        profile=pd.DataFrame(
            {
                "chr": [1, 1],
                "start": [0, 100],
                "end": [100, 200],
                "normal_A_copy": [1, 1],
                "normal_B_copy": [1, 1],
                "clone_1_A_copy": [2, 1],
                "clone_1_B_copy": [1, 1],
            }
        ),
    )
    segments = {"CHR": [1, 1], "START": [0, 100], "END": [100, 200]}
    pairs = pd.DataFrame({**segments, "clone0 A": [1, 1], "clone0 B": [1, 1], "clone1 A": [3, 1], "clone1 B": [1, 1]})
    states = pd.DataFrame({**segments, "clone0": ["neu", "neu"], "clone1": ["amp", "neu"]})
    labels = pd.Series([0, 1, 1, 0], index=truth.barcodes)
    genes = pd.DataFrame({"CHR": ["1", "1"], "START": [10, 150], "END": [20, 160]}, index=["g1", "g2"])

    expected = {"sensitivity": 0.5, "specificity": 5 / 6, "state_altered": 0.5}
    for profile in (pairs, states):
        assert get_state_metrics(get_state_confusion(profile, truth, fitted_at_barcode(labels, truth))) == expected
        at_genes = compute_supervised_metrics(MethodInput(labels, get_gene_profile(profile, genes)), truth)
        on_segments = compute_supervised_metrics(MethodInput(labels, profile), truth)
        assert pd.Series(at_genes).equals(pd.Series(on_segments))
    estimate = RunEstimate(clones=labels, int_clones=labels, seglevel=pairs)
    assert (sensitivity(estimate, truth), specificity(estimate, truth)) == (0.5, 5 / 6)
    assert compute_supervised_metrics(MethodInput(labels, pairs), truth)["exact_altered"] == exact_altered(
        estimate, truth
    )


def test_levels_merge_into_one_wide_table_and_back() -> None:
    """Each level's table, read back from the wide table, is itself; a level whose runs differ is refused."""
    columns = list(get_metric_columns())
    tables = {
        level: pd.DataFrame(
            [
                [
                    "sim1",
                    "cnamaste",
                    "expansion-lattice-lattice",
                    level,
                    10,
                    8,
                    0.8,
                    *([value] * (len(columns) - 9)),
                    2.0,
                    np.nan,
                ]
            ],
            columns=columns,
        )
        for level, value in (("segment", 0.5), ("gene", 0.25))
    }
    wide = merge_levels(tables)
    assert list(wide.columns[:7]) == [
        "sample_id",
        "method",
        "solver",
        "#spots",
        "#called",
        "spot-coverage-rate",
        "sensitivity (segment)",
    ]
    for level, table in tables.items():
        pd.testing.assert_frame_equal(get_level_table(wide, level), table, check_dtype=False)
    tables["gene"].loc[0, "#called"] = 9
    with pytest.raises(ValueError, match="gene"):
        merge_levels(tables)


def test_run_dir_pattern_finds_the_one_run(tmp_path: Path) -> None:
    """A `run_*` template finds the simulation's one cnamaste run, named by its configuration; two runs are refused."""
    simulation = {"cnamaste": {"dir": str(tmp_path / "{sample_id}" / "expansion-lattice-lattice" / "run_*")}}
    run = tmp_path / "sim1" / "expansion-lattice-lattice" / "run_64a6bcbc"
    run.mkdir(parents=True)
    assert get_run_dir(simulation, "cnamaste", "sim1") == run
    (run.parent / "run_20439582").mkdir()
    with pytest.raises(ValueError, match="2 runs"):
        get_run_dir(simulation, "cnamaste", "sim1")
