"""Two runs compared without a truth, on the spots both label, or on their clones' profiles (CNPs) alone.

- `get_shared_spots` -- the spots two runs both label;
- `compute_clone_ARI` -- the adjusted Rand index of their clones;
- `compute_pairwise_clone_ARI` -- that of every two runs' clones;
- `compute_cna_confusion` -- cells by `(state of A, state of B)`, a cell one shared spot at one gene;
- `get_contingency_ari` -- the adjusted Rand index of a contingency table;
- `compute_cna_ARI` -- the adjusted Rand index of the two runs' states over the cells both call;
- `compute_state_jaccard` -- the Jaccard index of each two states' cells, one state per run, over the cells both call;
- `get_altered_genes` -- the genes two clones' CNPs both call a state on, at least one other than `neu`;
- `compute_cnp_similarity` -- the share of those genes where the two CNPs' states match;
- `compute_pairwise_cnp_similarity` -- that of every two runs' most similar clones.

The adjusted Rand index is Hubert and Arabie's (1985, https://doi.org/10.1007/BF01908075), as scikit-learn's
`adjusted_rand_score` computes it from the pair confusion matrix (`sklearn.metrics.cluster._supervised`,
scikit-learn 1.9.1). It is a function of the contingency table alone, so the CNA ARI is that of the CNA confusion: the
two runs' partitions of the cells by state. Its cells are the shared spots at every gene, as `qa.metrics` counts a
run's cells against the truth, without realizing them. The ARI is unchanged where one run renames its states, so it
credits a run that calls every `del` of the other `cnloh`; the state Jaccard index matches states by name: of run A's
state `s` and run B's state `t`, `|A_s & B_t| / |A_s | B_t|` over the cells both call (Jaccard 1912), its diagonal the
agreement on each state.

The CNP similarity compares clone profiles alone, no spot: of two runs' clones `i` and `j`, at the genes both call and
at least one calls other than `neu`, `A_ij`, the share where their CNA states match, `J_ij`; of the two runs, the
highest `J_ij` over their clone pairs (`docs/reference.md`, CNP similarity). `J_ij` is the Jaccard index of the two
clones' altered genes (Jaccard 1912, https://doi.org/10.1111/j.1469-8137.1912.tb05611.x), a gene both alter shared only
where their states match.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scripts.metrics.cna_states import get_cna_states, get_uncalled
from sklearn.metrics import adjusted_rand_score


def get_shared_spots(labels_a: pd.Series, labels_b: pd.Series) -> pd.Index:
    """Get the spots two runs both label.

    Args:
        labels_a: run A's clone per spot, NaN where none.
        labels_b: run B's clone per spot, NaN where none.

    Returns:
        The spots, in `labels_a`'s order.
    """
    return labels_a.dropna().index.intersection(labels_b.dropna().index, sort=False)


def compute_clone_ARI(labels_a: pd.Series, labels_b: pd.Series) -> float:
    """Get the adjusted Rand index of two runs' clones over the spots both label.

    Args:
        labels_a: run A's clone per spot.
        labels_b: run B's clone per spot.

    Returns:
        The ARI, NaN where no spot is shared.
    """
    spots = get_shared_spots(labels_a, labels_b)
    if spots.empty:
        return float("nan")
    return float(adjusted_rand_score(labels_a[spots].to_numpy(), labels_b[spots].to_numpy()))


def compute_pairwise_clone_ARI(labels: dict[str, pd.Series]) -> pd.DataFrame:
    """Get the adjusted Rand index of every two runs' clones, each pair over the spots both label
    (`compute_clone_ARI`); a run against itself is 1, and A against B is B against A.

    Args:
        labels: each run's clone per spot, by key, in order.

    Returns:
        A row per ordered pair of runs, `run_a` then `run_b` in `labels`' order: `run_a`, `run_b`, `#shared-spots`,
        `clone_ari`.
    """
    keys = list(labels)
    found = {}
    for i, key_a in enumerate(keys):
        for key_b in keys[i:]:
            shared = len(get_shared_spots(labels[key_a], labels[key_b]))
            found[key_a, key_b] = found[key_b, key_a] = (shared, compute_clone_ARI(labels[key_a], labels[key_b]))
    rows = [(key_a, key_b, *found[key_a, key_b]) for key_a in keys for key_b in keys]
    return pd.DataFrame(rows, columns=["run_a", "run_b", "#shared-spots", "clone_ari"])


def compute_cna_confusion(
    labels_a: pd.Series,
    states_a: pd.DataFrame,
    labels_b: pd.Series,
    states_b: pd.DataFrame,
    values: list[str] | None = None,
) -> pd.DataFrame:
    """Get cells by `(state of A, state of B)`: per spot both label, the genes at each pair of its clones' states.

    Spots of clone `i` in A and `j` in B add their count times the genes at `(A's clone i state, B's clone j state)`,
    summed over the clone pairs. A state is a value of the field the runs are compared at
    (`cna_states.get_comparison_field`): a CNA state, a state of total copies, or a phase-free AB pair.

    Args:
        labels_a: run A's clone per spot.
        states_a: run A's state per gene, `clone<k>` per clone (`gene_level.get_gene_states`).
        labels_b: run B's clone per spot.
        states_b: run B's state per gene, the genes of `states_a`, in its order.
        values: the field's values in table order (`cna_states.get_field_order`), or None for `get_cna_states()`.

    Returns:
        Cells, A's states by row and B's by column, each `values` then `uncalled`.

    Raises:
        ValueError: where the two runs' genes differ, or a state is not among `values` and `uncalled`.
    """
    if not states_a.index.equals(states_b.index):
        raise ValueError("the two runs' states are on different genes")
    order = [*(get_cna_states() if values is None else values), get_uncalled()]
    found = set(np.ravel(states_a.to_numpy())) | set(np.ravel(states_b.to_numpy()))
    if found - set(order):
        raise ValueError(f"states {sorted(found - set(order))} not among {order}")
    spots = get_shared_spots(labels_a, labels_b)
    pairs = pd.DataFrame({"a": labels_a[spots].astype(int).to_numpy(), "b": labels_b[spots].astype(int).to_numpy()})
    confusion = pd.DataFrame(0, index=pd.Index(order, name="A"), columns=pd.Index(order, name="B"), dtype=np.int64)
    for (i, j), n in pairs.value_counts().items():
        genes = pd.crosstab(states_a[f"clone{i}"], states_b[f"clone{j}"])
        confusion = confusion.add(genes.reindex(index=order, columns=order, fill_value=0) * int(n), fill_value=0)
    return confusion.astype(np.int64)


def get_contingency_ari(table: np.ndarray) -> float:
    """Get the adjusted Rand index of the two partitions a contingency table counts.

    scikit-learn's `adjusted_rand_score` from its `pair_confusion_matrix`, the pair counts from the table, in Python
    integers: `n` cells give `n^2`-sized sums.

    Args:
        table: `(rows, columns)` counts, the items of partition A's part `i` and partition B's part `j`.

    Returns:
        The ARI, 1 where the two partitions agree on every pair, an empty table among them.
    """
    counts = [[int(x) for x in row] for row in np.asarray(table)]
    n = sum(map(sum, counts))
    rows = [sum(row) for row in counts]
    columns = [sum(column) for column in zip(*counts, strict=True)] if counts else []
    squares = sum(x * x for row in counts for x in row)
    tp = squares - n
    fp = sum(x * columns[j] for row in counts for j, x in enumerate(row)) - squares
    fn = sum(x * rows[i] for i, row in enumerate(counts) for x in row) - squares
    tn = n * n - fp - fn - squares
    if fn == 0 and fp == 0:
        return 1.0
    return 2.0 * (tp * tn - fn * fp) / ((tp + fn) * (fn + tn) + (tp + fp) * (fp + tn))


def compute_cna_ARI(confusion: pd.DataFrame) -> float:
    """Get the adjusted Rand index of two runs' states over the cells both call a state on.

    Args:
        confusion: `compute_cna_confusion`'s cells.

    Returns:
        The ARI of its rows and columns, `uncalled` left out.
    """
    states = [state for state in confusion.index if state != get_uncalled()]
    return get_contingency_ari(confusion.loc[states, states].to_numpy())


def compute_state_jaccard(confusion: pd.DataFrame) -> pd.DataFrame:
    """Get the Jaccard index of each two states' cells over the cells both runs call: of run A's state `s` and run B's
    state `t`, the cells where A calls `s` and B calls `t`, over those where A calls `s` or B calls `t`.

    Args:
        confusion: `compute_cna_confusion`'s cells.

    Returns:
        The index per A's state (rows) and B's state (columns), `uncalled` left out, NaN where neither state holds a
        cell both runs call.
    """
    states = [state for state in confusion.index if state != get_uncalled()]
    cells = confusion.loc[states, states].to_numpy(dtype=float)
    union = cells.sum(axis=1)[:, None] + cells.sum(axis=0)[None, :] - cells
    index = np.divide(cells, union, out=np.full_like(cells, np.nan), where=union > 0)
    return pd.DataFrame(index, index=confusion.loc[states, states].index, columns=confusion.loc[states, states].columns)


def get_altered_genes(cnp_a: pd.Series, cnp_b: pd.Series) -> np.ndarray:
    """Get whether each gene is one two clones' CNPs both call a state on, at least one of them other than `neu`.

    Args:
        cnp_a: a clone's CNA state per gene, `uncalled` where none (`gene_level.get_gene_states`).
        cnp_b: another clone's, on the genes of `cnp_a`, in its order.

    Returns:
        True or False per gene, in `cnp_a`'s order.

    Raises:
        ValueError: where the two CNPs' genes differ.
    """
    if not cnp_a.index.equals(cnp_b.index):
        raise ValueError("the two CNPs are on different genes")
    a, b = cnp_a.to_numpy(), cnp_b.to_numpy()
    called = (a != get_uncalled()) & (b != get_uncalled())
    return np.asarray(called & ((a != "neu") | (b != "neu")), dtype=bool)


def compute_cnp_similarity(cnp_a: pd.Series, cnp_b: pd.Series) -> float:
    """Get the CNP similarity of two clones: of the genes both call a state on, at least one other than `neu`
    (`get_altered_genes`), the share where their CNA states match.

    Args:
        cnp_a: a clone's CNA state per gene, `uncalled` where none (`gene_level.get_gene_states`).
        cnp_b: another clone's, on the genes of `cnp_a`, in its order.

    Returns:
        The share, NaN where no such gene.
    """
    altered = get_altered_genes(cnp_a, cnp_b)
    if not altered.any():
        return float("nan")
    return float(np.mean(cnp_a.to_numpy()[altered] == cnp_b.to_numpy()[altered]))


def compute_pairwise_cnp_similarity(cnps: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Get the CNP similarity of every two runs: that of their most similar pair of clones (`compute_cnp_similarity`),
    the first in clone order on a tie; a run against itself 1, and A against B is B against A, its clones swapped.

    Args:
        cnps: each run's CNPs compared, its tumor clones' (`real_data_metrics.get_tumor_cnps`), `clone<k>` per clone,
            by key, in order, every run on the same genes.

    Returns:
        A row per ordered pair of runs, `run_a` then `run_b` in `cnps`' order: `run_a`, `run_b`, `clone_a`, `clone_b`,
        the most similar pair's clones, `#altered-genes`, the genes it scores (`get_altered_genes`), and
        `cnp_similarity`; NA clones, 0 genes and NaN where no pair of clones has a gene to score.
    """
    keys = list(cnps)
    found = {}
    for i, key_a in enumerate(keys):
        for key_b in keys[i:]:
            pairs = [
                (compute_cnp_similarity(cnps[key_a][clone_a], cnps[key_b][clone_b]), clone_a, clone_b)
                for clone_a in cnps[key_a].columns
                for clone_b in cnps[key_b].columns
            ]
            scored = [pair for pair in pairs if not np.isnan(pair[0])]
            if not scored:
                found[key_a, key_b] = found[key_b, key_a] = (pd.NA, pd.NA, 0, float("nan"))
                continue
            similarity, clone_a, clone_b = max(scored, key=lambda pair: pair[0])
            altered = int(get_altered_genes(cnps[key_a][clone_a], cnps[key_b][clone_b]).sum())
            k_a, k_b = (int(str(clone).removeprefix("clone")) for clone in (clone_a, clone_b))
            found[key_a, key_b] = (k_a, k_b, altered, similarity)
            found[key_b, key_a] = (k_b, k_a, altered, similarity)
    rows = [(key_a, key_b, *found[key_a, key_b]) for key_a in keys for key_b in keys]
    columns = ["run_a", "run_b", "clone_a", "clone_b", "#altered-genes", "cnp_similarity"]
    return pd.DataFrame(rows, columns=columns).astype({"clone_a": "Int64", "clone_b": "Int64"})
