"""Each clone's pseudobulk RDR and BAF per bin beside the RDR and BAF its decoded `(A, B)` is expected to read."""

from typing import Any

import numpy as np
import pandas as pd

from cnamaste.hmm_ab_decoder import mixture
from cnamaste.pseudobulk import merge_pseudobulk_by_index_mix
from cnamaste.spatio_genomic_counts import from_npz


def get_clone_expectation(
    seglevel: pd.DataFrame, clone: int, log_mu_shift: float | None = None, normal_fraction: float | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """A clone's expected RDR and BAF per bin, on its pseudobulk dots' scale, as `plot_clones_genomic` draws them.

    With `normal_fraction` f_c, the decoded `(A, B)` mixed with normal cells at f_c (`hmm_ab_decoder.mixture`):
    RDR `(f_c + (1 - f_c)(A + B) / 2) / Z_c`, BAF `(f_c + (1 - f_c) B) / (2 f_c + (1 - f_c)(A + B))`. Without it, the
    fitted state's, `exp(logmu) / Z_c` and `p`. `Z_c = exp(log_mu_shift)`, 1 where the fit has none.

    Args: `seglevel` the run's `cnv_<ploidy>_seglevel.tsv`, `clone<k> A`, `B`, `logmu`, `p` per bin; `clone` k;
        `log_mu_shift` the fit's `new_log_mu_shift` of the clone, or None; `normal_fraction` its decoded f_c, or None.
    Returns: the expected RDR and BAF, (n_bins,) each.
    """
    z = 1.0 if log_mu_shift is None else float(np.exp(log_mu_shift))
    if normal_fraction is None:
        return np.exp(seglevel[f"clone{clone} logmu"].to_numpy()) / z, seglevel[f"clone{clone} p"].to_numpy()
    copies = seglevel[[f"clone{clone} A", f"clone{clone} B"]].to_numpy()
    depth, b_fraction = mixture(copies, float(normal_fraction))
    return depth / z, b_fraction


def get_clone_pseudobulk(
    binned: Any,
    seglevel: pd.DataFrame,
    clones: pd.Series,
    log_mu_shift: np.ndarray | None = None,
    normal_fractions: dict[int, float] | None = None,
) -> pd.DataFrame:
    """Each clone's pseudobulk RDR and BAF per bin, as `plot_clones_genomic` draws them, beside its decoded `(A, B)`
    and the RDR and BAF it is expected to read (`get_clone_expectation`).

    A clone's RDR is its read depth over its normal baseline and its BAF its B count over its SNP-covering depth, each
    summed over its spots (`merge_pseudobulk_by_index_mix`); a bin with no depth reads NaN.

    Args: `binned` an open `binned_data.npz` (`utils.write_binned_data`); `seglevel` the run's
        `cnv_<ploidy>_seglevel.tsv`, a row per bin, `clone<k> A`, `B`, `logmu`, `p` per clone; `clones` clone per
        barcode, the run's `clone_labels.tsv`, whose clones `seglevel` decodes; `log_mu_shift` the fit's
        `new_log_mu_shift`, one per clone in ascending order, or None; `normal_fractions` clone -> its decoded f_c,
        a clone absent taking its fitted states.
    Returns: `CHR`, `START`, `END`, `clone`, `RDR`, `BAF`, `A`, `B`, `exp_RDR`, `exp_BAF`, a row per clone and bin,
        clones ascending.
    Raises: `ValueError` where `seglevel` has not a row per bin, a binned spot has no clone, or `log_mu_shift` has not
        one value per clone.
    """
    barcodes = np.asarray(binned["barcodes"]).astype(str)
    single_X = from_npz(binned, "single_X")
    if single_X.shape[0] != len(seglevel):
        raise ValueError(f"{single_X.shape[0]} bins against {len(seglevel)} seglevel rows")
    labels = clones.reindex(barcodes)
    if labels.isna().any():
        raise ValueError(f"{int(labels.isna().sum())} binned spots have no clone")
    labels = labels.astype(int).to_numpy()
    ids = np.unique(labels)
    shifts = None if log_mu_shift is None else np.asarray(log_mu_shift, dtype=np.float64).reshape(-1)
    if shifts is not None and len(shifts) != len(ids):
        raise ValueError(f"{len(shifts)} log_mu_shift values for {len(ids)} clones")
    fractions = normal_fractions or {}
    expected = [
        get_clone_expectation(seglevel, k, None if shifts is None else shifts[i], fractions.get(int(k)))
        for i, k in enumerate(ids)
    ]
    X, base_nb_mean, total_bb_RD, _ = merge_pseudobulk_by_index_mix(
        single_X,
        from_npz(binned, "single_base_nb_mean"),
        from_npz(binned, "single_total_bb_RD"),
        [np.flatnonzero(labels == k) for k in ids],
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        rdr, baf = X[:, 0, :] / base_nb_mean, X[:, 1, :] / total_bb_RD
    bins = seglevel[["CHR", "START", "END"]].reset_index(drop=True)
    return pd.concat(
        [
            bins.assign(
                clone=k,
                RDR=rdr[:, i],
                BAF=baf[:, i],
                A=seglevel[f"clone{k} A"].to_numpy(),
                B=seglevel[f"clone{k} B"].to_numpy(),
                exp_RDR=expected[i][0],
                exp_BAF=expected[i][1],
            )
            for i, k in enumerate(ids)
        ],
        ignore_index=True,
    )
