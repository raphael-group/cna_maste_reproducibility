"""Ported from cna_maste `python/cnamaste/qa/plot.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim, with their imports pointed at `scripts.cnamaste_port`.

The cross-sample validation figure: `all_sample_metrics.tsv`'s scores per simulation group, one method at a time.

`plot_validation` draws cnaster's `plot_validation_stats.plot_metrics` from `run_sim_qa`'s table: a 2 x 2 grid, one
metric per panel, a box per simulated class (samples sharing a CNA number and size) with each sample overlaid. A panel whose
metric the table does not carry for the method reads "not computed". `sim_groups` parses the sample ids into groups.

The figure is `combined.png`'s format (`figure_style`): half a page, the text column wide, every text in the stated
serif at its one size, panels lettered (a)-(d).

Departures from cnaster's: matplotlib only (no seaborn), so the boxes and points are drawn by hand; every panel keeps
one x position per group present for the method, so the bottom row's tick labels name every panel's boxes; the
panels are CNA sensitivity and normal specificity over the phase-free N-CNA recall and (N-) Clone ARI, with no log-likelihood
or third row; specificity's axis starts at 0.8 (`FLOORS`); there is no title; the size is written `50M` rather than
`5e7`; groups tied on size and number of CNAs are ordered by name.
"""

from __future__ import annotations

NAMES = {"cnamaste": r"$\mathtt{cn}\alpha$-$\mathtt{ma}\mathit{st}\mathtt{e}$", "calicost": r"$\mathtt{CalicoST}$"}
"""How the key writes a method; a method not here is written as named."""
