"""Ported from cna_maste `python/cnamaste/config/defaults.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim, with their imports pointed at `scripts.cnamaste_port`.

The configuration's defaults and the values each choice takes: what a configuration that does not carry a key
reads, and what one that does may set.
"""

from __future__ import annotations

from typing import Any

MERGE_STAGES = ("baf", "rdr")
"""The stages that fit clones and so can merge them: the BAF-only search, and the BAF/RDR one."""


INTEGER_DECODERS = ("hmm_ab", "calicost")
"""`int_copy_num.decoder`: the `(A, B)` HMM (the default) or CalicoST's diploid hill climb (copied in `hmm_ab_decoder`).
No other value is taken: `lattice_cluster`, one pair per cluster of bins, is set aside in
`cnamaste.sandbox.decoder.integer_copy`."""


HMM_AB_SOLVERS = ("icm", "expansion", "panclone_viterbi")
"""`int_copy_num.solver`: how the hmm_ab decoders decode the tumour clones coupled by `clone_coupling`
(`hmm_ab_decoder.HMM_AB_SOLVER`). `icm` the default; `expansion` OpenGM's alpha expansion over the whole field, without
`genome_wide_ploidy_constraint`; `panclone_viterbi` the exact joint decode, `O(n K^C)` in time and memory."""


INITIALIZATION_METHODS = ("best_equal_partition",)
"""`hmrf.initialization_method`: how the BAF clones start, `spatial.best_equal_partition`, the spots cut into 3 x 3
rectangles of equal count and grouped (`run_cnamaste`). No other start is taken.

`phasing`'s rectangles (`spatial.initialize_clones`) and the RDR refinement's split
(`spatial.initialize_rectangular_clones`) are their own stages' and not named by this key."""


CONFIG_DEFAULTS: dict[str, dict[str, Any]] = {
    # NB the slide's classes, one per percentile band of its grayscale (`he.get_he_image`), and whether each spot's
    #    class refines the initial partition (`run_cnamaste`), where the slide gives one; and the image the H&E
    #    network reads (`cnamaste.cnn.he`, `run_cnamaste_he`), Space Ranger's hires or lowres.
    "he": {"apply_to_initial_partition": True, "n_components": 4, "resolution": "hires"},
    # NB `write_binned_data`: whether the run writes `binned_data.npz`, the per-spot bins a page is redrawn from after
    #    the run (the sandbox's truth pages). Off by default: it grows with the spots, 1.33 GB at 40,000 Visium HD bins,
    #    and no shipped reader needs it.
    "run": {"plot_field_argmax": False, "write_binned_data": False},
    # NB the HMM's read-depth ceiling (`hmm_nophasing.get_max_rdr_zeroed_baseline`) and BAF deadband about 1/2
    #    (`hmm_nophasing.get_normal_clone`, `get_neutral_state`), the values the code held as `MAX_RDR` and `EPS_BAF`.
    "hmm": {"max_rdr": 5.0, "eps_baf": 0.05},
    "hmrf": {
        "tumorprop_threshold": 0.0,
        "n_clones_rdr": 2,
        # NB whether a stage merges the clone pair that lowers the labelling cost, while one does.
        **{f"merge_clones_{stage}": True for stage in MERGE_STAGES},
        "label_solver": "opengm_expansion",
        "opengm_mincut_solver": "boykov_kolmogorov",
        "initialization_method": INITIALIZATION_METHODS[0],
    },
    "int_copy_num": {
        "rdr": True,
        "decoder": INTEGER_DECODERS[0],
        # NB CalicoST's diploid hill climb's thresholds (`calicost.adapters.calicost_fixdiploid`).
        "nonbalance_bafdist": 1.0,
        "nondiploid_rdrdist": 10.0,
        # NB the lattice decoder's: `hmm_ab_decoder` reads them here.
        "specificity_prior": 0.5,
        # NB the BB concentration cap: flat likelihood above it. Set here alone, not in the shipped configurations.
        "max_tau": 1e6,
        "max_total_copy": 6,
        "max_allele_copy": 5,
        "min_alpha": 0.0,
        "fit_normal": True,
        "apply_logmu_shift": True,
        "recompute_logmu_shift": False,
        # NB the decoder's stay probability and clone coupling, `hmm_ab_decoder._hmm_ab_ploidy_decoder`'s own defaults.
        "t": 1.0 - 1e-6,
        "clone_coupling": 0.0,
        "assumed_ploidy": "diploid",
        "genome_wide_ploidy_constraint": True,
        # NB two clones are one integer clone where their decoded (A, B) agree on at least this share of bins: the
        #    ledger's `clone_ari_int_099`.
        "int_clone_agreement": 0.99,
        "solver": HMM_AB_SOLVERS[0],
    },
}
"""Section -> key -> the value a configuration that does not carry the key takes.

A key is here rather than in a shipped file where every file set it alike and nothing chose it per run: the
default is then one value rather than two copies that can drift. A section the file does not carry is taken whole,
so `config.<section>.<key>` reads a default as it reads a set value. A default also types its key (`typed`): the
file's value and `-o`'s text are read as a bool, int or float as the default is, so a reader takes the value as it is.
"""
