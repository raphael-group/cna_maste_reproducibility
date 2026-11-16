"""Copy-number reconstruction error: each run's clone copies against the pseudobulk counts of its clones' spots, as
copytyping-benchmark's `analysis/scripts/cn_error.py` (1204da9) scores them, at UGP's bins for RDR and BAF.

- `get_baseline` -- each bin's read-depth baseline from the reference spots;
- `get_clone_spots`, `get_pseudobulk` -- each count column's clone, and per clone a count matrix's sum over its spots;
- `get_observed_rdr`, `get_expected_rdr` -- per clone and bin, its pseudobulk's log2 RDR, and its copies';
- `get_observed_baf`, `get_expected_baf` -- per clone and bin, its pseudobulk's BAF, and its copies' `B / (A + B)`;
- `get_bin_copies`, `get_altered` -- a run's copies per bin, those on its joint segment, and whether it calls a CNA
  there;
- `score_rdr`, `score_baf`, `get_baf_flips`, `get_skill_score` -- a run's mean squared error, the diploid
  prediction's, each bin's orientation, and its skill score;
- `get_cn_error_table` -- each run's errors on the bins every run calls and some run calls a CNA on.

The pseudobulks, with `X` UMI per bin and spot, `T_n` spot `n`'s UMI over the bins, `B` and `D` B-allele and total
allele counts per bin and spot (`read_method_input.read_bin_counts`), `C`, `A`, `B` a run's total and allele copies
(copytyping e0341e4 `inference/model_utils.py`: `compute_baseline_proportions`, `pseudobulk_rdr_matrix`,
`pseudobulk_baf_matrix`, `clone_rdr_gk`):

    lambda_u     = sum_{n in ref} X[u,n] / sum_{n in ref} T_n
    RDR_obs[k,u] = sum_{n in k} X[u,n] / (lambda_u * sum_{n in k} T_n)
    RDR_exp[k,u] = C[u,k] / sum_u' lambda_u' C[u',k]
    BAF_obs[k,u] = sum_{n in k} B[u,n] / sum_{n in k} D[u,n]
    BAF_exp[k,u] = B[u,k] / (A[u,k] + B[u,k])

RDR is compared in log2. Each mean squared error weighs a clone by its spots `n_k`, over the bins `U` scored:

    mse_rdr = sum_{u in U} sum_k n_k (RDR_obs - RDR_exp)^2 / sum_{u in U} sum_k n_k
    mse_baf = sum_{u in U} min(sum_k n_k (BAF_obs - BAF_exp)^2, sum_k n_k (BAF_obs - (1 - BAF_exp))^2)
              / sum_{u in U} sum_k n_k

BAF is compared mirrored, "a method's A|B ordering carries no phase, so each bin is scored against whichever of exp /
(1 - exp) fits better. The flip is chosen PER BIN but JOINTLY over the clones -- a bin has a single phase"
(cn_error.py). The null is the diploid prediction's error, every clone 1:1 throughout: its RDR `2 / sum_u lambda_u 2`
is 1, log2 0, as `lambda` sums to 1, and its BAF 0.5. The skill score `ss = 1 - mse / null`: 1 a perfect
reconstruction, 0 the diploid prediction's, below 0 worse. A bin is scored where every run competing on the channel
calls every clone and every clone's observed and expected values are finite, the bins shared by the panel, and where
some run of the panel calls a CNA in some clone (`get_altered`): a bin every run calls neutral throughout, where the
diploid prediction is every run's, is left out.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scripts.cnamaste_port.config import start_time
from scripts.cnamaste_port.logger import get_logger
from scripts.metrics.cna_states import ProfileCopies, get_profile_copies
from scripts.metrics.gene_level import get_max_overlap
from scripts.metrics.read_method_input import BinCounts, MethodInput

logger = get_logger("cn_error", start_time=start_time)

CN_ERROR_COLUMNS = (
    "run",
    "#clones",
    "#tumor-clones",
    "#spots",
    "#rdr-bins",
    "mse_rdr",
    "null_rdr",
    "ss_rdr",
    "#baf-bins",
    "#flipped-bins",
    "mse_baf",
    "null_baf",
    "ss_baf",
)
"""A run's row of the error table: its clones, tumor clones and spots, then per channel its bins scored, its error, the
diploid prediction's error `null_<channel>` and its skill score."""


def get_cn_error_columns() -> tuple[str, ...]:
    """Get the error table's columns: `run`, `#clones`, `#tumor-clones`, `#spots`, then per channel, RDR then BAF, its
    bins scored (`#rdr-bins`; `#baf-bins`, `#flipped-bins`), `mse_`, `null_` and `ss_`."""
    return CN_ERROR_COLUMNS


def get_baseline(counts: sp.csr_matrix, reference: np.ndarray) -> np.ndarray:
    """Get each bin's read-depth baseline: its UMI over the reference spots' UMI (`compute_baseline_proportions`).

    Args:
        counts: `(bins, spots)` UMI.
        reference: True per reference spot, `(spots,)`.

    Returns:
        `(bins,)` baselines, summing to 1.

    Raises:
        ValueError: where no spot is a reference spot.
    """
    if not reference.any():
        raise ValueError("no reference spots")
    reference_umi = np.asarray(counts[:, reference].sum(axis=1)).ravel().astype(float)
    return reference_umi / reference_umi.sum()


def get_clone_spots(labels: pd.Series, spots: pd.Index, clones: list[int]) -> np.ndarray:
    """Get each count column's clone among `clones`, -1 where the run labels none.

    Args:
        labels: the run's clone per spot.
        spots: the count matrix's spots, its columns.
        clones: the clones the run's profile gives copies of.

    Returns:
        `(spots,)` clone numbers.

    Raises:
        ValueError: where the run labels a spot with a clone its profile has no copies of.
    """
    found = labels.reindex(spots)
    missing = set(found.dropna().astype(int)) - set(clones)
    if missing:
        raise ValueError(f"clones {sorted(missing)} label spots but have no copies")
    return found.fillna(-1).astype(int).to_numpy()


def get_pseudobulk(counts: sp.csr_matrix, clone_spots: np.ndarray, clones: list[int]) -> np.ndarray:
    """Get per clone a count matrix's sum over its spots.

    Args:
        counts: `(bins, spots)` counts.
        clone_spots: each spot's clone, -1 where none (`get_clone_spots`).
        clones: the clones, the output's column order.

    Returns:
        `(bins, clones)` sums.
    """
    members = sp.csr_matrix((clone_spots[:, None] == np.asarray(clones)[None, :]).astype(float))
    return np.asarray((counts @ members).todense(), dtype=float)


def get_observed_rdr(pooled: np.ndarray, baseline: np.ndarray) -> np.ndarray:
    """Get per clone and bin its pseudobulk's log2 RDR (`pseudobulk_rdr_matrix`), its UMI over the bin's baseline
    times its UMI over the bins; NaN where its count or baseline is 0.

    Args:
        pooled: `(bins, clones)` pseudobulk UMI (`get_pseudobulk`).
        baseline: `(bins,)` baselines (`get_baseline`).

    Returns:
        `(clones, bins)` log2 RDR.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        rdr = pooled / (baseline[:, None] * pooled.sum(axis=0)[None, :])
        return np.where(rdr > 0, np.log2(rdr), np.nan).T


def get_expected_rdr(total: np.ndarray, baseline: np.ndarray) -> np.ndarray:
    """Get per clone and bin its copies' log2 RDR (`clone_rdr_gk`), the denominator over the bins it calls; NaN where
    it calls none or no copies.

    Args:
        total: `(bins, clones)` total copies, NaN where uncalled.
        baseline: `(bins,)` baselines.

    Returns:
        `(clones, bins)` log2 RDR.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        denominator = np.nansum(baseline[:, None] * total, axis=0)
        rdr = np.where(denominator > 0, total / denominator, np.nan)
        return np.where(rdr > 0, np.log2(rdr), np.nan).T


def get_observed_baf(b_counts: np.ndarray, totals: np.ndarray) -> np.ndarray:
    """Get per clone and bin its pseudobulk's BAF (`pseudobulk_baf_matrix`), NaN where it has no allele counts.

    Args:
        b_counts: `(bins, clones)` pseudobulk B-allele counts (`get_pseudobulk`).
        totals: `(bins, clones)` pseudobulk allele counts.

    Returns:
        `(clones, bins)` BAF.
    """
    return np.divide(b_counts, totals, out=np.full(totals.shape, np.nan), where=totals > 0).T


def get_expected_baf(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Get per clone and bin its copies' BAF, `B / (A + B)`, NaN where `A + B` is 0 or uncalled.

    Args:
        a: `(bins, clones)` A copies.
        b: `(bins, clones)` B copies.

    Returns:
        `(clones, bins)` BAF.
    """
    total = a + b
    return np.divide(b, total, out=np.full(total.shape, np.nan), where=total > 0).T


def get_bin_copies(
    bins: pd.DataFrame, segments: pd.DataFrame, profile: pd.DataFrame, copies: ProfileCopies
) -> ProfileCopies:
    """Get a run's copies per bin: those of its segment holding the bin's joint segment, the one overlapping the bin
    most (`gene_level.get_max_overlap`), NaN where no joint segment or no segment of the run overlaps it.

    Args:
        bins: `CHR`, `START`, `END` per bin.
        segments: the joint segments (`segment_level.get_joint_segments`).
        profile: the run's `MethodInput.profile`.
        copies: its copies per segment (`cna_states.get_profile_copies`).

    Returns:
        Its copies, `(bins, clones)` each.
    """
    joint = get_max_overlap(bins, segments)
    rows = get_max_overlap(segments, profile)
    row = np.where(joint >= 0, rows[np.maximum(joint, 0)], -1)
    called = (row >= 0)[:, None]

    def at_bins(values: np.ndarray) -> np.ndarray:
        return np.where(called, values[np.maximum(row, 0)], np.nan)

    return ProfileCopies(copies.clones, at_bins(copies.a), at_bins(copies.b), at_bins(copies.total))


def get_altered(a: np.ndarray, b: np.ndarray, total: np.ndarray) -> np.ndarray:
    """Get whether a run calls a CNA in some clone in each bin: a clone's copies not 1:1, or, given total copies
    alone, not 2; a clone it calls nothing for none.

    Args:
        a: `(bins, clones)` A copies, NaN where none.
        b: `(bins, clones)` B copies, NaN where none.
        total: `(bins, clones)` total copies, NaN where uncalled.

    Returns:
        `(bins,)` True or False.
    """
    allelic = np.isfinite(a) & np.isfinite(b)
    neutral = np.where(allelic, (a == 1) & (b == 1), total == 2)
    return (np.isfinite(total) & ~neutral).any(axis=1)


def get_skill_score(mse: float, null: float) -> float:
    """Get the skill score `1 - mse / null`, NaN where either is not finite or the null is not positive.

    Args:
        mse: the run's mean squared error.
        null: the diploid prediction's.

    Returns:
        The score.
    """
    if not (np.isfinite(mse) and np.isfinite(null)) or null <= 0:
        return np.nan
    return float(1.0 - mse / null)


def score_rdr(observed: np.ndarray, expected: np.ndarray, sizes: np.ndarray, keep: np.ndarray) -> dict[str, float]:
    """Score a run's log2 RDR on the kept bins: its mean squared error, the diploid prediction's (log2 RDR 0), and its
    skill score.

    Args:
        observed: `(clones, bins)` observed log2 RDR.
        expected: `(clones, bins)` expected log2 RDR.
        sizes: `(clones,)` spots per clone, each clone's weight.
        keep: `(bins,)` the bins scored.

    Returns:
        `#rdr-bins`, `mse_rdr`, `null_rdr`, `ss_rdr`.
    """
    ok = keep[None, :] & np.isfinite(observed) & np.isfinite(expected)
    weights = np.where(ok, sizes[:, None], 0.0)
    total = weights.sum()
    if total <= 0:
        return {"#rdr-bins": 0, "mse_rdr": np.nan, "null_rdr": np.nan, "ss_rdr": np.nan}
    mse = float((weights * np.where(ok, observed - expected, 0.0) ** 2).sum() / total)
    null = float((weights * np.where(ok, observed, 0.0) ** 2).sum() / total)
    return {
        "#rdr-bins": int(ok.any(axis=0).sum()),
        "mse_rdr": mse,
        "null_rdr": null,
        "ss_rdr": get_skill_score(mse, null),
    }


def get_baf_flips(observed: np.ndarray, expected: np.ndarray, sizes: np.ndarray, keep: np.ndarray) -> np.ndarray:
    """Get each bin's orientation: flipped, its expected BAF read as `1 - BAF`, where that fits the observed BAF
    better, its squared error weighed by spots summed over the clones, one orientation for every clone.

    Args:
        observed: `(clones, bins)` observed BAF.
        expected: `(clones, bins)` expected BAF.
        sizes: `(clones,)` spots per clone.
        keep: `(bins,)` the bins scored.

    Returns:
        `(bins,)` True where flipped.
    """
    ok = keep[None, :] & np.isfinite(observed) & np.isfinite(expected)
    weights = np.where(ok, sizes[:, None], 0.0)
    as_is = (weights * np.where(ok, observed - expected, 0.0) ** 2).sum(axis=0)
    flipped = (weights * np.where(ok, observed - (1.0 - expected), 0.0) ** 2).sum(axis=0)
    return ok.any(axis=0) & (flipped < as_is)


def score_baf(observed: np.ndarray, expected: np.ndarray, sizes: np.ndarray, keep: np.ndarray) -> dict[str, float]:
    """Score a run's BAF on the kept bins, mirrored per bin jointly over its clones: its mean squared error, the
    diploid prediction's (BAF 0.5), its skill score, and the bins scored flipped.

    Args:
        observed: `(clones, bins)` observed BAF.
        expected: `(clones, bins)` expected BAF.
        sizes: `(clones,)` spots per clone, each clone's weight.
        keep: `(bins,)` the bins scored.

    Returns:
        `#baf-bins`, `#flipped-bins`, `mse_baf`, `null_baf`, `ss_baf`.
    """
    ok = keep[None, :] & np.isfinite(observed) & np.isfinite(expected)
    weights = np.where(ok, sizes[:, None], 0.0)
    total = weights.sum()
    if total <= 0:
        return {"#baf-bins": 0, "#flipped-bins": 0, "mse_baf": np.nan, "null_baf": np.nan, "ss_baf": np.nan}
    flip = get_baf_flips(observed, expected, sizes, keep)
    scored = ok.any(axis=0)
    mse = float((weights * np.where(ok, observed - np.where(flip, 1.0 - expected, expected), 0.0) ** 2).sum() / total)
    null = float((weights * np.where(ok, observed - 0.5, 0.0) ** 2).sum() / total)
    return {
        "#baf-bins": int(scored.sum()),
        "#flipped-bins": int(flip.sum()),
        "mse_baf": mse,
        "null_baf": null,
        "ss_baf": get_skill_score(mse, null),
    }


def get_cn_error_table(
    runs: dict[str, MethodInput],
    labels: dict[str, pd.Series],
    normal: dict[str, list[int]],
    segments: pd.DataFrame,
    bins: BinCounts,
    reference: pd.Index,
) -> pd.DataFrame:
    """Get each run's RDR and BAF error at the bins, on the bins every run competing on the channel calls, finite for
    each of its clones, and some run calls a CNA on (`get_altered`); a run with no allele copies (Clonalscope) has no
    BAF error.

    Args:
        runs: each run by its key, its profile giving copies (`cna_states.get_profile_copies`).
        labels: each run's clone per spot of the counts' side, by key.
        normal: each run's normal clones, by key, left out of its `#tumor-clones`.
        segments: the runs' joint segments (`segment_level.get_joint_segments`).
        bins: the UMI and allele counts per bin and spot (`read_method_input.read_bin_counts`), its bins those scored.
        reference: the reference spots, whose UMI give the baseline.

    Returns:
        `get_cn_error_columns()`, a row per run in `runs`' order, BAF's NaN and 0 for a run with no allele copies.
    """
    baseline = get_baseline(bins.x_counts, bins.spots.isin(reference))
    rdr, baf, altered = {}, {}, []
    for key, run in runs.items():
        copies = get_profile_copies(run.profile)
        if copies is None:
            raise ValueError(f"{key}: its profile gives no copies")
        clone_spots = get_clone_spots(labels[key], bins.spots, copies.clones)
        clones = [k for k in copies.clones if (clone_spots == k).any()]
        columns = [copies.clones.index(k) for k in clones]
        at_bins = get_bin_copies(bins.bins, segments, run.profile, copies)
        a, b, total = at_bins.a[:, columns], at_bins.b[:, columns], at_bins.total[:, columns]
        altered.append(get_altered(a, b, total))
        sizes = np.array([(clone_spots == k).sum() for k in clones], dtype=float)
        tumor = np.array([k not in normal[key] for k in clones])
        pooled = get_pseudobulk(bins.x_counts, clone_spots, clones)
        rdr[key] = {
            "observed": get_observed_rdr(pooled, baseline),
            "expected": get_expected_rdr(total, baseline),
            "sizes": sizes,
            "tumor": tumor,
        }
        if np.isfinite(copies.a).any():
            b_counts = get_pseudobulk(bins.b_counts, clone_spots, clones)
            totals = get_pseudobulk(bins.totals, clone_spots, clones)
            baf[key] = {
                "observed": get_observed_baf(b_counts, totals),
                "expected": get_expected_baf(a, b),
                "sizes": sizes,
            }

    def get_shared(channel: dict[str, dict[str, np.ndarray]]) -> np.ndarray:
        if not channel:
            return np.zeros(len(bins.bins), dtype=bool)
        found = [(np.isfinite(v["observed"]) & np.isfinite(v["expected"])).all(axis=0) for v in channel.values()]
        return np.logical_and.reduce(found)

    shared_rdr, shared_baf = get_shared(rdr), get_shared(baf)
    altered_any = np.logical_or.reduce(altered)
    keep_rdr, keep_baf = shared_rdr & altered_any, shared_baf & altered_any
    logger.info(
        f"Bins scored, shared and with a CNA, of {len(bins.bins):_}: RDR {int(keep_rdr.sum()):_} of "
        f"{int(shared_rdr.sum()):_} shared for {len(rdr)} runs, BAF {int(keep_baf.sum()):_} of "
        f"{int(shared_baf.sum()):_} shared for {len(baf)} runs."
    )

    rows = []
    for key, found in rdr.items():
        row: dict[str, object] = {
            "run": key,
            "#clones": len(found["sizes"]),
            "#tumor-clones": int(found["tumor"].sum()),
            "#spots": int(found["sizes"].sum()),
        }
        row |= score_rdr(found["observed"], found["expected"], found["sizes"], keep_rdr)
        if key in baf:
            alleles = baf[key]
            row |= score_baf(alleles["observed"], alleles["expected"], alleles["sizes"], keep_baf)
        else:
            row |= {"#baf-bins": 0, "#flipped-bins": 0}
        rows.append(row)
    return pd.DataFrame(rows, columns=list(CN_ERROR_COLUMNS))
