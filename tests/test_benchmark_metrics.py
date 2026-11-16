"""The benchmark's gene-level profiles and truth-free metrics (`scripts/metrics`).

Referees: genes, segments, spots and shares counted by hand; scikit-learn's `adjusted_rand_score` on the cells
realized one by one; counts built from copies by copytyping's RDR and BAF expectations, so the copies they were built
from reconstruct them exactly.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from scripts.metrics.cn_error import get_cn_error_table
from scripts.metrics.cna_states import (
    get_comparison_field,
    get_field_order,
    get_phase_free_pairs,
    get_profile_copies,
    get_uncalled,
)
from scripts.metrics.gene_level import get_gene_states, get_interval_fields
from scripts.metrics.qc_metrics import compute_genome_composition, get_region_weights
from scripts.metrics.read_method_input import (
    BinCounts,
    MethodInput,
    get_labelled_profile,
)
from scripts.metrics.real_data_metrics import merge_run_tables
from scripts.metrics.segment_level import get_joint_segments, get_segment_level_table
from scripts.metrics.unsupervised_metrics import (
    compute_clone_ARI,
    compute_cna_ARI,
    compute_cna_confusion,
    compute_cnp_similarity,
    compute_pairwise_clone_ARI,
    compute_pairwise_cnp_similarity,
    compute_state_jaccard,
)
from sklearn.metrics import adjusted_rand_score

GENES = pd.DataFrame({"CHR": ["1", "1", "2"], "START": [80, 300, 0], "END": [120, 310, 10]}, index=["g1", "g2", "g3"])
"""Three genes: g1 straddling base 100 of chromosome 1, g2 at base 300, g3 on chromosome 2."""


def test_gene_fields_take_the_segment_overlapping_most() -> None:
    """A gene takes the segment that overlaps it most, the first on a tie (g1, 20 bases each), and is `uncalled` where
    none does (g3); Numbat's `amp` has no one pair, its `loh` is `cnloh` 2|0; AB pairs give `A|B` and `A + B`; every
    method's `tcn_state` is its state of total copies, `cnloh` at 2 `neu`."""
    segments = {"CHR": ["chr1", "chr1"], "START": [0, 100], "END": [100, 400]}
    numbat = pd.DataFrame({**segments, "clone0": ["amp", "loh"]})
    pairs = pd.DataFrame({**segments, "clone0 A": [2, 0], "clone0 B": [1, 2]})

    fields = get_interval_fields(numbat, GENES)
    assert fields["clone0 state"].tolist() == ["amp", "cnloh", "uncalled"]
    assert fields["clone0 ascn"].tolist() == ["unknown", "2|0", "uncalled"]
    assert fields["clone0 tcn"].tolist() == ["unknown", "2", "uncalled"]
    assert fields["clone0 tcn_state"].tolist() == ["amp", "neu", "uncalled"]
    fields = get_interval_fields(pairs, GENES)
    assert fields[["clone0 ascn", "clone0 tcn", "clone0 state", "clone0 tcn_state"]].to_numpy().tolist() == [
        ["2|1", "3", "amp", "amp"],
        ["0|2", "2", "cnloh", "neu"],
        ["uncalled", "uncalled", "uncalled", "uncalled"],
    ]


def test_cna_confusion_and_ari_count_shared_spot_genes() -> None:
    """Spots s1-s3 both runs label, s4 run A alone. Per shared spot, its clones' state at each gene: s1 and s2 (A
    clone 0, B clone 0) neu/neu at g1-g3; s3 (A clone 1, B clone 0) del/neu at g1, amp/neu at g2, uncalled/neu at g3.
    The CNA ARI is scikit-learn's over the cells both call, realized; the clone ARI over the shared spots, every two
    runs' alike, a run against itself 1. Of the 8 cells both call, B's 8 neu: the state Jaccard index of A's neu and
    B's neu is 6 / 8, of A's del and B's neu 1 / 8, of A's del and B's del 0 / 1; a state neither calls NaN. With equal
    calls, it is 1 on the diagonal and 0 off it."""
    labels_a = pd.Series([0, 0, 1, 1], index=["s1", "s2", "s3", "s4"])
    labels_b = pd.Series([0, 0, 0], index=["s1", "s2", "s3"])
    states_a = pd.DataFrame({"clone0": ["neu", "neu", "neu"], "clone1": ["del", "amp", "uncalled"]}, index=GENES.index)
    states_b = pd.DataFrame({"clone0": ["neu", "neu", "neu"]}, index=GENES.index)

    confusion = compute_cna_confusion(labels_a, states_a, labels_b, states_b)
    assert int(confusion.to_numpy().sum()) == 3 * 3
    assert (confusion.loc["neu", "neu"], confusion.loc["del", "neu"], confusion.loc["amp", "neu"]) == (6, 1, 1)
    assert confusion.loc["uncalled", "neu"] == 1

    cells_a = ["neu"] * 6 + ["del", "amp"]
    cells_b = ["neu"] * 8
    assert compute_cna_ARI(confusion) == adjusted_rand_score(cells_a, cells_b)
    jaccard = compute_state_jaccard(confusion)
    assert get_uncalled() not in jaccard.index
    assert (jaccard.loc["neu", "neu"], jaccard.loc["del", "neu"], jaccard.loc["del", "del"]) == (6 / 8, 1 / 8, 0.0)
    assert np.isnan(jaccard.loc["cnloh", "cnloh"])
    assert compute_clone_ARI(labels_a, labels_b) == adjusted_rand_score([0, 0, 1], [0, 0, 0])
    pairwise = compute_pairwise_clone_ARI({"a": labels_a, "b": labels_b}).set_index(["run_a", "run_b"])
    assert pairwise.loc[("a", "b"), "clone_ari"] == pairwise.loc[("b", "a"), "clone_ari"]
    assert pairwise.loc[("a", "b"), "clone_ari"] == compute_clone_ARI(labels_a, labels_b)
    assert pairwise.loc[[("a", "a"), ("b", "b")], "clone_ari"].tolist() == [1.0, 1.0]
    assert pairwise["#shared-spots"].tolist() == [4, 3, 3, 3]

    states_b = pd.DataFrame({"clone0": ["neu", "neu", "neu"], "clone1": ["del", "amp", "neu"]}, index=GENES.index)
    labels_b = pd.Series([0, 0, 1], index=["s1", "s2", "s3"])
    confusion = compute_cna_confusion(labels_a, states_a, labels_b, states_b)
    assert compute_cna_ARI(confusion) == adjusted_rand_score(cells_a, cells_a) == 1.0
    jaccard = compute_state_jaccard(confusion).loc[["neu", "del", "amp"], ["neu", "del", "amp"]].to_numpy()
    assert (jaccard == np.eye(3)).all()


def test_cnamaste_profile_keeps_its_integer_clones() -> None:
    """Clones 0-3 fitted, 2 merged into 0's integer clone and 3 into 1's: the labels carry 0 and 1, so the profile keeps
    clone 0's and clone 1's fields, as `qa.metrics.merged_at_barcode` names a group by a member whose copies it calls,
    and every column that is no clone's."""
    profile = pd.DataFrame(
        {"CHR": [1], "START": [0], "END": [10], **{f"clone{k} {f}": [1] for k in range(4) for f in ("A", "B")}}
    )
    kept = get_labelled_profile(profile, pd.Series([0, 1, 1, 0], index=["s1", "s2", "s3", "s4"]))
    assert list(kept.columns) == ["CHR", "START", "END", "clone0 A", "clone0 B", "clone1 A", "clone1 B"]


def test_cnp_similarity_takes_the_most_similar_clone_pair() -> None:
    """At the genes both clones call, at least one not `neu`: A1 against B2 scores g1-g2, A1's g4 uncalled and g3 `neu`
    in both, and matches g1, 1/2; A2 against B2 scores g1, g2, g4 and matches g1, 1/3; B1 calls nothing and two `neu`
    clones alter nothing, NaN. The runs' similarity is A1 against B2's, B against A its clones swapped; a run against
    itself 1; a run with no clone NaN against every run."""
    genes = ["g1", "g2", "g3", "g4"]
    states_a = pd.DataFrame(
        {"clone1": ["amp", "del", "neu", "uncalled"], "clone2": ["amp", "del", "neu", "cnloh"]}, index=genes
    )
    states_b = pd.DataFrame({"clone1": ["uncalled"] * 4, "clone2": ["amp", "neu", "neu", "neu"]}, index=genes)
    assert compute_cnp_similarity(states_a["clone1"], states_b["clone2"]) == 1 / 2
    assert compute_cnp_similarity(states_a["clone2"], states_b["clone2"]) == 1 / 3
    assert np.isnan(compute_cnp_similarity(states_a["clone1"], states_b["clone1"]))
    neutral = pd.Series(["neu"] * 4, index=genes)
    assert np.isnan(compute_cnp_similarity(neutral, neutral))

    cnps = {"a": states_a, "b": states_b, "c": states_a[[]]}
    pairwise = compute_pairwise_cnp_similarity(cnps).set_index(["run_a", "run_b"])
    best = ["clone_a", "clone_b", "#altered-genes", "cnp_similarity"]
    assert pairwise.loc[("a", "b"), best].tolist() == [1, 2, 2, 1 / 2]
    assert pairwise.loc[("b", "a"), best].tolist() == [2, 1, 2, 1 / 2]
    assert pairwise.loc[[("a", "a"), ("b", "b")], "cnp_similarity"].tolist() == [1.0, 1.0]
    assert pairwise.xs("c", level="run_a")["cnp_similarity"].isna().all()


def test_runs_compare_at_the_finest_field_both_give() -> None:
    """Two AB-pair runs compare at their pairs, phase-free: A's 2|1 and B's 1|2 at g1 agree, A's 0|2 and B's 1|1 at g2
    do not; a run of total copies alone (Clonalscope) compares at `tcn_state`, Numbat's states at `state`."""
    segments = {"CHR": ["1", "1"], "START": [0, 200], "END": [200, 400]}
    pairs_a = pd.DataFrame({**segments, "clone0 A": [2, 0], "clone0 B": [1, 2]})
    pairs_b = pd.DataFrame({**segments, "clone0 A": [1, 1], "clone0 B": [2, 1]})
    totals = pd.DataFrame({**segments, "clone0": ["amp", "neu"], "clone0 Uest": [1.5, 1.0], "clone0 tcn": [3, 2]})
    numbat = pd.DataFrame({**segments, "clone0": ["amp", "neu"]})
    assert get_comparison_field(pairs_a, pairs_b) == "ascn"
    assert get_comparison_field(pairs_a, totals) == "tcn_state"
    assert get_comparison_field(pairs_a, numbat) == "state"

    states = [get_phase_free_pairs(get_gene_states(run, GENES, "ascn")) for run in (pairs_a, pairs_b)]
    values = get_field_order("ascn", states)
    assert values == ["(1,1)", "(2,0)", "(2,1)"]
    labels = pd.Series([0], index=["s1"])
    confusion = compute_cna_confusion(labels, states[0], labels, states[1], values)
    assert (confusion.loc["(2,1)", "(2,1)"], confusion.loc["(2,0)", "(1,1)"]) == (1, 1)
    assert confusion.loc["uncalled", "uncalled"] == 1


def test_genome_composition_shares_sum_to_one() -> None:
    """Regions of 300 bases; segments of 100 neu, 50 del outside them and 100 amp in clone 1 inside: per clone, shares
    of 300 summing to 1 with `uncalled`; `all` weighs clone 0 by 3 spots and clone 1 by 1; at genes, 1 per gene, g1
    `neu`, g2 and g3 `uncalled`."""
    regions = pd.DataFrame({"CHR": ["1", "2"], "START": [0, 0], "END": [200, 100]})
    profile = pd.DataFrame({"CHR": [1, 1, 3], "START": [0, 100, 0], "END": [100, 200, 50]})
    states = np.array([["neu", "neu"], ["neu", "amp"], ["del", "del"]])
    weights = get_region_weights(profile, regions)
    assert weights.tolist() == [100, 100, 0]

    table = compute_genome_composition(states, weights, 300.0, pd.Series([0, 0, 0, 1]))
    assert table["clone"].tolist() == ["clone0", "clone1", "all"]
    assert table[["neu", "amp", "uncalled"]].round(6).to_numpy().tolist() == [
        [round(2 / 3, 6), 0.0, round(1 / 3, 6)],
        [round(1 / 3, 6), round(1 / 3, 6), round(1 / 3, 6)],
        [round(7 / 12, 6), round(1 / 12, 6), round(1 / 3, 6)],
    ]
    shares = ["neu", "del", "amp", "cnloh", "bamp", "bdel", "uncalled"]
    assert np.allclose(table[shares].sum(axis=1), 1.0)

    genes = get_gene_states(pd.DataFrame({**profile, "clone0": ["neu", "amp", "del"]}), GENES)
    table = compute_genome_composition(genes.to_numpy(dtype=str), np.ones(len(GENES)), float(len(GENES)))
    assert np.allclose(table.loc[0, ["neu", "amp", "uncalled"]].to_numpy(dtype=float), [1 / 3, 0, 2 / 3])


def test_joint_segments_split_at_every_runs_breakpoint() -> None:
    """`pairs` breaks chr1 at 35 and 55; `partial` at 55, and calls 100-120 beyond the others' 80: the joint segments
    are the intervals between consecutive breakpoints some run's segment holds, the gap at 80-100 none. Each run takes,
    per joint segment, its segment holding it, `uncalled` where it has none. Numbat's states give no copies."""
    labels = pd.Series([0], index=["s1"])
    pairs = pd.DataFrame(
        {"CHR": 1, "START": [0, 35, 55], "END": [35, 55, 80], "clone0 A": [1, 3, 1], "clone0 B": [1, 1, 0]}
    )
    partial = pd.DataFrame(
        {"CHR": "chr1", "START": [0, 100], "END": [55, 120], "clone0": "neu", "clone0 Uest": 1.0, "clone0 tcn": 2}
    )
    runs = {"pairs": MethodInput(labels, pairs), "partial": MethodInput(labels, partial)}

    segments = get_joint_segments([run.profile for run in runs.values()])
    assert segments[["START", "END"]].to_numpy().tolist() == [[0, 35], [35, 55], [55, 80], [100, 120]]
    table = get_segment_level_table(segments, runs)
    assert table["pairs clone0 ascn"].tolist() == ["1|1", "3|1", "1|0", "uncalled"]
    assert table["partial clone0 state"].tolist() == ["neu", "neu", "uncalled", "neu"]

    copies = get_profile_copies(pairs)
    assert copies is not None and copies.total[:, 0].tolist() == [2.0, 4.0, 1.0]
    assert get_profile_copies(pd.DataFrame({"CHR": [1], "START": [0], "END": [80], "clone0": ["amp"]})) is None


def test_cn_error_scores_copies_against_their_clones_pseudobulk() -> None:
    """Bins u1-u4 on chr1, spots n1, n2 of normal clone 0 (10 UMI a bin, the baseline 1/4 each) and t1, t2 of tumor
    clone 1 (20, 20, 40, 10 UMI: copies 2, 2, 4, 1 times 90 UMI over their 9). `pairs` calls those copies, 3:1 at u3:
    its log2 RDR is exact, MSE 0, SS 1. `flat` calls 2 throughout, the diploid prediction, SS 0. `partial` leaves u4
    uncalled, so no run is scored there. u1 and u2, neutral in every run, are not scored: u3 alone is. u3 reads BAF
    0.75 in the tumor against 3:1's 0.25, scored mirrored, flipped, 0 error; u4 has no SNP. The two total-copy runs
    have no BAF."""
    pairs = pd.DataFrame(
        {
            "CHR": 1,
            "START": [0, 35, 55],
            "END": [35, 55, 80],
            "clone0 A": [1, 1, 1],
            "clone0 B": [1, 1, 1],
            "clone1 A": [1, 3, 1],
            "clone1 B": [1, 1, 0],
        }
    )
    flat = pd.DataFrame({"CHR": "chr1", "START": [0], "END": [80], "clone0 tcn": [2.0], "clone1 tcn": [2.0]})
    partial = pd.DataFrame({"CHR": "1", "START": [0, 100], "END": [55, 120], "clone0 tcn": 2.0, "clone1 tcn": 2.0})
    spots = pd.Index(["n1", "n2", "t1", "t2"])
    labels = pd.Series([0, 0, 1, 1], index=spots)
    runs = {
        key: MethodInput(labels, profile) for key, profile in zip(("pairs", "flat", "partial"), (pairs, flat, partial))
    }
    segments = get_joint_segments([run.profile for run in runs.values()])

    bins = pd.DataFrame({"CHR": "chr1", "START": [0, 20, 40, 60], "END": [10, 30, 50, 70]})
    umi = sp.csr_matrix([[10, 10, 20, 20], [10, 10, 20, 20], [10, 10, 40, 40], [10, 10, 10, 10]])
    b_counts = sp.csr_matrix([[2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 3, 3], [0, 0, 0, 0]])
    totals = sp.csr_matrix([[4, 4, 4, 4], [4, 4, 4, 4], [4, 4, 4, 4], [0, 0, 0, 0]])

    table = get_cn_error_table(
        runs,
        {k: labels for k in runs},
        {k: [0] for k in runs},
        segments,
        BinCounts(bins, umi, b_counts, totals, spots),
        spots[:2],
    )
    table = table.set_index("run")
    assert table["#rdr-bins"].tolist() == [1, 1, 1]
    assert table.loc["pairs", "mse_rdr"] == 0.0 and table.loc["pairs", "ss_rdr"] == 1.0
    assert table.loc["flat", "mse_rdr"] == table.loc["flat", "null_rdr"] > 0 and table.loc["flat", "ss_rdr"] == 0.0
    assert table.loc["pairs", ["#baf-bins", "#flipped-bins", "mse_baf", "ss_baf"]].tolist() == [1, 1, 0.0, 1.0]
    assert table.loc[["flat", "partial"], "mse_baf"].isna().all()


def test_per_run_tables_merge_in_run_order_and_agree() -> None:
    """Per-run tables merge into a row per run in the keys' order, a run a table lacks NaN, a column two tables hold
    kept once; a run whose shared column differs between tables is refused."""
    spot = pd.DataFrame({"run": ["b", "a"], "method": ["m", "n"], "#called": [5, 7]})
    errors = pd.DataFrame({"run": ["a"], "method": ["n"], "#called": [7], "mse_rdr": [0.1]})
    wide = merge_run_tables([spot, errors], ["a", "b"])
    assert wide.columns.tolist() == ["run", "method", "#called", "mse_rdr"]
    assert wide["#called"].tolist() == [7, 5]
    assert wide["mse_rdr"].isna().tolist() == [False, True]
    with pytest.raises(ValueError, match="#called"):
        merge_run_tables([spot, errors.assign(**{"#called": [6]})], ["a", "b"])


def test_metric_text_defines_every_metric_scored() -> None:
    """Every metric column the simulation and real-data tables score has a paragraph of its dataset in the metric text,
    and every paragraph names a page drawing it, so the text cannot drift from the metrics."""
    pytest.importorskip("cnplot", reason="cnplot is in the analysis environment: uv sync (pyproject.toml)")
    from scripts.metrics.metric_text import get_metric_texts, get_missing_columns

    texts = get_metric_texts()

    assert get_missing_columns(texts) == []
    assert all(text.shown for text in texts)
    assert get_missing_columns(tuple(text for text in texts if text.name != "Skill score")) == [
        "real_data: ss_rdr",
        "real_data: ss_baf",
    ]
