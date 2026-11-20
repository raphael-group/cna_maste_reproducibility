"""Ported from cna_maste `python/cnamaste/hmm_ab_decoder.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim, with their imports pointed at `scripts.cnamaste_port`.

(A,B) decoder for (allele-specific) integer copy numbers given a K-state (mu, p) hmm fit.

The reference is `cna-maste-paper` (81813c5), section `sec:ab_decoder`. Equations are cited by their label, not by
their number, which the sections before them move: `eqn:nb_emission` and `eqn:bb_emission` the emission,
`eqn:ab_admixture` the normal admixture, `eqn:ab_hmm` the scored path, `eqn:ab_coupling` the clone coupling.

The names here are this module's and are not the paper's symbols. The correspondence, once, so a reader of either
can follow the other:

    paper       here                              paper    here
    theta       normal_fraction, f_c              eta_s    specificity_prior
    mu, p       `mixture`'s two returns           eta_c    clone_coupling
    Z_m         Z, the exp of log_mu_shift        eta_t    the price `_mean_viterbi` bisects
    G, M        n_obs, n_clones                   P        ploidy
    T_bar_m     the decoded path's mean A + B     K        the lattice's pairs

Known departures from the reference, each stated again where it is made: `_mean_viterbi` bisects the ploidy price
where the paper steps it; `NORMAL_FRACTIONS` is a grid the paper does not name; `fold` scores the minor count,
which the paper does not do, and is reached only by `hmm_initialize`.
"""

import numpy as np
from scripts.cnamaste_port.config.defaults import CONFIG_DEFAULTS

# TODO cnamaste.config.defaults.int_copy_num if necessary, otherwise via cnamaste.config, which should be defaulted if not specified.
INT_COPY_NUM = CONFIG_DEFAULTS["int_copy_num"]
"""The `int_copy_num` keys' defaults, the configuration's own: what a decoder takes where its caller passes none.

Every keyword of `_hmm_ab_ploidy_decoder` that names an `int_copy_num` key takes its default from here, so the
decoder states no default of its own and the two cannot drift. A run passes the loaded configuration's values
(`run_cnamaste`), so these bind only where a caller passes none."""


def lattice_pairs(
    *,
    params="mp",
    fold=False,
    max_allele_copy=INT_COPY_NUM["max_allele_copy"],
    max_total_copy=INT_COPY_NUM["max_total_copy"],
):
    """The `(A, B)` states the hmm_ab decoders decode over, under the same `params` and `fold`.

    `sec:ab_decoder`'s lattice: every `(A, B)` with `A >= B`, `0 < A + B <= max_total_copy` and each at most
    `max_allele_copy`, which at the defaults 6 and 5 is the paper's set exactly. Without read depth
    (`"m" not in params`) states with one `p` are one, represented by the lowest total; without folding each class also
    enters at its reflection. One function, so the decoder, `hmm_initialize.lattice_init` and the GMM seed agree.
    `fold` is a departure: the paper scores both `(A, B)` and `(B, A)`, the phased counts, which is `fold=False`.

    Args: `params` the fit's channels; `fold` keep `A >= B` only; `max_allele_copy`, `max_total_copy` bounds.
    Returns: the lattice's `(A, B)` pairs.
    """
    pairs = np.array([(a, b) for a in range(max_allele_copy + 1) for b in range(a + 1) if 0 < a + b <= max_total_copy])

    if "m" not in params:
        fraction = pairs[:, 1] / pairs.sum(axis=1)
        _, first = np.unique(np.round(fraction, 9), return_index=True)
        pairs = pairs[np.sort(first)]

    # TODO `fold` is reached only by `hmm_initialize` (`only_minor=True`) and has no counterpart in the reference,
    #      which scores both `(A, B)` and `(B, A)`. Keep it as an initialization-only option with its own
    #      measurement, or retire it.
    if not fold:
        flipped = pairs[pairs[:, 0] != pairs[:, 1]][:, ::-1]
        pairs = np.vstack([pairs, flipped])

    return pairs


def mixture(copies, normal_fraction):
    """The read depth and B-allele fraction of `(A, B)` cells mixed with normal (1, 1) cells at `normal_fraction`.

    `eqn:ab_admixture` with `mu = (A + B) / 2` and `p = B / (A + B)` substituted and the halves cancelled: the depth
    is relative to a diploid normal, `f + (1 - f) (A + B) / 2`; the B fraction `(f + (1 - f) B) / (2 f + (1 - f)
    (A + B))`. A clone's expected count is `base_nb_mean` times the depth over `Z_c = exp(log_mu_shift_c)`.

    Args: `copies` (..., 2) integer `(A, B)`; `normal_fraction` f in [0, 1].
    Returns: the depth and the B fraction, each of shape `copies.shape[:-1]`.
    """
    copies = np.asarray(copies)
    tot = copies.sum(axis=-1)
    fc = normal_fraction
    return fc + (1 - fc) * tot / 2.0, (fc + (1 - fc) * copies[..., 1]) / (2 * fc + (1 - fc) * tot)
