"""Ported from cna_maste `python/cnamaste/plot/plot_genomic.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim, with their imports pointed at `scripts.cnamaste_port`.
"""

from typing import Any, Dict, Optional

import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.transforms import Transform, blended_transform_factory, offset_copy
from scripts.cnamaste_port.config.clock import start_time
from scripts.cnamaste_port.figure_style import DPI, column_size, stagger, style
from scripts.cnamaste_port.hmm_ab_decoder import lattice_pairs, mixture
from scripts.cnamaste_port.logger import get_logger
from scripts.cnamaste_port.palette import FLIP_OPACITY, get_full_palette
from scripts.cnamaste_port.pseudobulk import merge_pseudobulk_by_index_mix
from scripts.cnamaste_port.utils import cast_clone_label, get_intervals

NAME_GAP = 4.0
"""Points between a track's y axis furniture (tick labels and axis label) and the clone's name left of it."""


class _LeftOfYAxis(Transform):
    """Move a point's x to `NAME_GAP` left of `ax`'s y axis as drawn, tick labels and axis label included; read at
    every draw, as the figure's fonts are set on the way out (`restyle`).

    Args: `ax` the track.
    """

    input_dims = output_dims = 2
    is_affine = False

    def __init__(self, ax: Any) -> None:
        super().__init__()
        self.ax = ax

    def transform_non_affine(self, values: Any) -> Any:
        figure = self.ax.figure
        edge = self.ax.yaxis.get_tightbbox().x0 - NAME_GAP * figure.dpi / 72.0
        values = np.array(values, dtype=float)
        values[..., 0] = edge
        return values


logger = get_logger(__name__, start_time=start_time)


NORMAL_OPACITY = 0.75


RDR_FLOOR, RDR_YLIM = -0.5, 6.0
"""An RDR track's y range: from just below 0 to `rdr_ylim`, by default `RDR_YLIM`."""


def _create_clone_gridspec(
    n_pairs: int, axes_per_clone: int, base_height: float, sample_list: list = None, top: bool = False
):
    """Create a GridSpec for clone tracks, with vertical gaps between clones, and a panel above them where `top`.

    Args: `n_pairs` clones; `axes_per_clone` tracks each; `base_height` per track; `sample_list` the samples; `top`.
    Returns: the figure, its clone axes, and the top panel's axes or None.
    """
    n_axes_total = axes_per_clone * n_pairs
    fig = plt.figure(figsize=column_size(base_height * (n_pairs + 0.5 * top)), dpi=DPI, facecolor="white")

    height_ratios = [0.5, 0.25] if top else []
    for i in range(n_pairs):
        height_ratios.extend([1] * axes_per_clone)
        if i < n_pairs - 1:
            height_ratios.append(0.25)  # Gap between clone tracks

    gs = gridspec.GridSpec(len(height_ratios), 1, height_ratios=height_ratios, hspace=0)

    top_ax = fig.add_subplot(gs[0, 0]) if top else None
    axes, row = [], 2 if top else 0
    for i in range(n_axes_total):
        axes.append(fig.add_subplot(gs[row, 0]))
        row += 1
        if (i % axes_per_clone == axes_per_clone - 1) and (i < n_axes_total - 1):
            row += 1  # Skip the gap row

    if sample_list is not None:
        fig.suptitle(", ".join(sample_list), x=0.5, y=0.99, ha="center")

    return fig, axes, top_ax


def _format_track_axis(ax, ylabel, ylim, yticks, remove_xticks, n_obs, grid=None):
    """Style an RDR or BAF genomic track: light gray lines at `grid`, `yticks` where it is None.

    Args: `ax` the track; `ylabel`, `ylim`, `yticks`; `remove_xticks`; `n_obs` the x extent; `grid` the lines' levels.
    """
    ax.set_ylabel(ylabel)
    ax.set_ylim(ylim)
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{y:.1f}" for y in yticks])
    ax.set_xlim([0, n_obs])

    if remove_xticks:
        ax.set_xticks([])

    for y in yticks if grid is None else grid:
        ax.axhline(y=y, c="lightgray", linewidth=0.5, zorder=0, gid="grid")


def _draw_chromosome_boundaries(axes, lengths, unique_chrs, chrtext_shift):
    """Draw vertical contig boundaries and label chromosomes on the bottom axis.

    Args: `axes` the tracks; `lengths` bins per chromosome; `unique_chrs` names; `chrtext_shift` the label offset.
    """
    starts = [float(np.sum(lengths[:i])) for i in range(len(lengths))]
    names = [f"chr{c}" for c in unique_chrs]
    middles = [start + float(length) / 2.0 for start, length in zip(starts, lengths)]

    # NB each name level under its chromosome's middle, and on a second line where it would meet its neighbour's:
    #    combined.png's chromosome axis (`layout.GenomicAxis.frame`). Rotated names anchored at the boundary ran
    #    into each other from chr8 on at the text column's width, and a rotated name is read at an angle (T- #46).
    size = style()["font_size"]
    points = (
        axes[-1].get_position().width * axes[-1].figure.get_size_inches()[0] * 72.0 / max(float(np.sum(lengths)), 1.0)
    )
    # NB three lines, not combined.png's two: easy's last chromosomes are short enough that two leave chr20
    #    on chr22 at this width.
    lines = stagger([m * points for m in middles], names, size, lines=3)

    for start_len, middle, name, line in zip(starts, middles, names, lines):
        # Label only on the bottom-most axis
        at = offset_copy(axes[-1].get_xaxis_transform(), axes[-1].figure, y=-2.0 - line * 1.1 * size, units="points")

        axes[-1].text(middle, 0.0, name, transform=at, ha="center", va="top")

        # Draw boundaries across all axes
        for ax in axes:
            ax.axvline(x=start_len, c="black", linewidth=0.5)


def _annotate_clone_stats(
    ax,
    clone_label,
    n_spots,
    n_umis,
    n_snp_umis,
    tumor_prop=None,
    paired_ax=None,
):
    """Annotate the axis with the clone's identity and its spot and UMI statistics.

    The name sits `NAME_GAP` points left of the y axis's tick labels and axis label as drawn (`_LeftOfYAxis`), not
    a fixed distance from the axes, which landed it on the axis label ("BAF").

    Args: `ax` the track; `clone_label` the name drawn; `n_spots`, `n_umis`, `n_snp_umis`; `tumor_prop`; `paired_ax` its partner.
    """
    at = ax.transAxes + _LeftOfYAxis(ax)

    ax.text(
        0.0,
        0.5 if paired_ax is None else 0.0,
        clone_label,
        ha="right",
        va="center",
        rotation="vertical",
        transform=at,
        clip_on=False,
    )

    theta_text = f"$\\hat{{\\theta}}={tumor_prop:.2f}$" if tumor_prop is not None else ""
    stats_text = f"{n_spots:_} spots; {int(n_umis):_} umis; {int(n_snp_umis):_} snp-umis; {theta_text}"

    ax.text(
        0.0,
        1.02,
        stats_text.strip("; "),
        ha="left",
        va="bottom",
        transform=ax.transAxes,
    )


def _state_legend(
    phased: np.ndarray, map_cn: dict[Any, int], default: int, palette: list, normal: bool, drop_unused: bool = True
) -> list:
    """One legend dot per `(larger, smaller)` pair a clone's phased `(A, B)` fold to, in the palette's order, with the
    share of bins it and its mirror hold together: `35% (1, 0)` for 27% (1, 0) and 8% (0, 1). (1, 1) only on the
    normal clone, so a tumour clone's key lists its aberrations alone; only the profile tells a pair from its mirror.
    With `drop_unused`, a pair whose share reads `0%` is left out of the key; its bins are drawn either way.

    Args: `phased` `(n_obs, 2)` integer copies; `map_cn` palette state -> index; `default` the index of a state
        outside the palette; `palette` RGBA per index; `normal` whether the clone is the normal one; `drop_unused`.
    Returns: the legend handles.
    """
    folded = np.sort(phased, axis=1)[:, ::-1]
    pairs, counts = np.unique(folded, axis=0, return_counts=True)
    rows = sorted((map_cn.get((a, b), default), a, b, 100.0 * n / len(phased))
                  for (a, b), n in zip(pairs.tolist(), counts.tolist()) if normal or (a, b) != (1, 1))  # fmt: skip
    if drop_unused:
        rows = [row for row in rows if not reads_zero(row[3])]
    return [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=palette[index], markersize=10, linestyle="None",
               label=f"{share:.0f}% ({a}, {b})")
        for index, a, b, share in rows
    ]  # fmt: skip


def reads_zero(share: float) -> bool:
    """Whether a percentage prints as `0%` at the keys' precision, `:.0f`."""
    return f"{share:.0f}" in ("0", "-0")


BAF_GRID_TOTAL = 4
"""The largest total copy `A + B` whose B fractions `BAF_GRID` draws."""


def lattice_bafs(max_total_copy: int = BAF_GRID_TOTAL) -> np.ndarray:
    """The B fractions of the decoder's lattice pairs up to `max_total_copy`, with no normal admixture, each once.

    At no admixture a pair reads `B / (A + B)`, so equivalent pairs read one value exactly: every balanced pair 1/2, a
    pair and its mirror `p` and `1 - p`, and the pairs of one ratio `B : A` one value; each is drawn once.

    Args: `max_total_copy` the largest `A + B`.
    Returns: the sorted values.
    """
    _, b = mixture(lattice_pairs(max_total_copy=max_total_copy), 0.0)
    return np.unique(np.round(b, 9))


BAF_GRID = lattice_bafs()
"""Every BAF track's light gray lines: 0, 1/4, 1/3, 1/2, 2/3, 3/4 and 1, the pairs up to total copy 4 at no admixture."""


def clone_states(param: np.ndarray, clone: int) -> np.ndarray:
    """One clone's state parameters: `(n_states,)` shared by every clone, or column `clone` of `(n_states, n_clones)`.

    cnamaste's fit shares its states across clones; CalicoST's final fit holds one column per clone, in
    `new_assignment`'s sorted order, as its `pred_cnv` does.

    Args: `param` `(n_states,)` or `(n_states, n_clones)`; `clone` the clone's index in sorted order.
    Returns: `(n_states,)`.
    """
    param = np.asarray(param)
    assert param.ndim in (1, 2), f"the fit's parameters are (n_states,) or (n_states, n_clones), not {param.shape}"
    return param if param.ndim == 1 else param[:, clone]


def plot_clones_genomic(
    lengths: np.ndarray,
    single_X: np.ndarray,
    single_base_nb_mean: np.ndarray,
    single_total_bb_RD: np.ndarray,
    df_cnv: Optional[pd.DataFrame] = None,
    res_combine: Optional[Dict[str, Any]] = None,
    single_tumor_prop: Optional[np.ndarray] = None,
    clone_index: Optional[list] = None,
    sample_list: Optional[list] = None,
    remove_xticks: bool = True,
    rdr_ylim: float = RDR_YLIM,
    chrtext_shift: float = -0.25,
    base_height: float = 3.2,
    pointsize: float = 3.0,
    linewidth: float = 1.0,
    palette_name: str = "tracks",
    plot_baf_errors: str = "beta",
    plot_rdr_errors: str = "poisson",
    known_nb_baseline=None,
    axes: Optional[list] = None,
    genome=None,
    normal_fractions: Optional[Dict[Any, float]] = None,
    phase_switches: Optional[np.ndarray] = None,
    colors: Optional[list] = None,
    labels: Optional[list] = None,
    drop_unused_states: bool = True,
):
    """Plot each clone's aggregated RDR and BAF, with error models, and the expectation each segment is fitted to.

    The black lines are what each segment's counts are expected to read, on the dots' own scale (count over
    `base_nb_mean`, B count over depth), a clone's RDR divided by its `Z_c = exp(new_log_mu_shift_c)` where the fit
    has one. Given `df_cnv` and the clone's normal fraction in `normal_fractions`, they are the integer `(A, B)` mixed
    with normal cells at that fraction (`hmm_ab_decoder.mixture`), as the decoder scores it; otherwise the fitted
    states, `exp(new_log_mu[s])` and `new_p_binom[s]` along the clone's `pred_cnv`. The dashed line is `1 - p`.

    With `df_cnv` None, plots the raw data without integer states. Given `axes`, one per track (RDR then BAF per clone),
    draws into them; given `genome`, a `cnamaste.plot.layout.GenomicAxis`, every track takes its x axis. Given
    `phase_switches`, a panel (a) above the tracks (b) draws the run's inferred phase switches per bin, counted along
    each contig (`phasing.cumulative_phase_switches`); only on a figure this function makes. Given `colors` and
    `labels`, one per group of `clone_index`, the raw tracks draw each group in its colour and name it by its label, as
    is, in place of the clone's numeral; only without a fit. With `drop_unused_states`, a clone's key leaves out the
    states whose share of its bins reads `0%`; the bins are drawn either way.

    Args: `single_X` (n_obs, 2, n_spots); `df_cnv`, `res_combine` the fit; `clone_index` spots per clone; `axes`;
        `normal_fractions` clone label (`clone {id}`) -> the decoder's normal fraction; `phase_switches` (n_obs,);
        `colors`, `labels` per group; `drop_unused_states`.
    Returns: the figure.
    """
    logger.info("Plotting aggregated rdr and baf for clones.")

    # 1. State resolution and hierarchy setup
    if df_cnv is not None:
        assert res_combine is not None, "res_combine required if df_cnv is provided."
        unique_chrs = np.unique(df_cnv.CHR.values)
        final_clone_ids = np.sort(np.unique(res_combine["new_assignment"]))
        clone_index = [np.where(res_combine["new_assignment"] == c)[0] for c in final_clone_ids]
        color_palette, ordered_acn = get_full_palette(palette_name)
        state_colors = [color_palette[c] for c in ordered_acn]
        map_cn = {x: i for i, x in enumerate(ordered_acn)}
        default_idx = map_cn.get((1, 1), 0)

    elif res_combine is not None:
        unique_chrs = 1 + np.arange(len(lengths))
        final_clone_ids = np.sort(np.unique(res_combine["new_assignment"]))
        clone_index = [np.where(res_combine["new_assignment"] == c)[0] for c in final_clone_ids]
        n_states = res_combine["new_p_binom"].shape[0]
        # Pre-compute fallback palette outside the loop
        base_pal = sns.color_palette("deep", n_states)
        palette = [mcolors.to_rgba(color, alpha=1.0) for color in base_pal]

    else:
        assert clone_index is not None, "clone_index must be provided."
        assert lengths is not None, "lengths must be provided."
        unique_chrs = 1 + np.arange(len(lengths))
        final_clone_ids = [str(i) for i in range(len(clone_index))]

    assert (colors is None and labels is None) or res_combine is None, "colors and labels name raw groups, not clones"
    for given in (colors, labels):
        assert given is None or len(given) == len(clone_index), f"{len(given)} for {len(clone_index)} groups"

    assert single_X.shape[0] == np.sum(lengths), "Mismatch in genomic segment defined X and lengths."

    X, base_nb_mean, total_bb_RD, tumor_prop = merge_pseudobulk_by_index_mix(
        single_X,
        single_base_nb_mean,
        single_total_bb_RD,
        clone_index,
        single_tumor_prop,
    )

    n_obs = X.shape[0]
    spots_per_clone = [len(xx) for xx in clone_index]
    nonempty_clones = np.where(np.sum(total_bb_RD, axis=0) > 0)[0]

    assert len(nonempty_clones) == total_bb_RD.shape[1]
    assert np.all(nonempty_clones == np.arange(len(final_clone_ids)))

    if known_nb_baseline is not None:
        base_nb_mean = known_nb_baseline.copy()

    has_rdr = base_nb_mean is not None and np.max(base_nb_mean) > 0
    axes_per_clone = 2 if has_rdr else 1

    top_ax = None
    if axes is None:
        fig, axes, top_ax = _create_clone_gridspec(
            len(nonempty_clones), axes_per_clone, base_height, sample_list, top=phase_switches is not None
        )
        own_figure = True
    else:
        assert phase_switches is None, "the phase-switch panel is drawn on a figure plot_clones_genomic makes"
        assert len(axes) == axes_per_clone * len(nonempty_clones), f"{len(axes)} axes for {len(nonempty_clones)} clones"
        fig, own_figure = axes[0].figure, False
    if top_ax is not None:
        switches = np.asarray(phase_switches)
        assert switches.shape == (n_obs,), f"{switches.shape} phase switches for {n_obs} bins"
        top_ax.plot(
            np.arange(n_obs), switches, drawstyle="steps-mid", color="black", linewidth=linewidth, gid="phase_switches"
        )
        high = max(int(switches.max()), 1)
        _format_track_axis(top_ax, "phase\nswitches", (0, 1.05 * high), [0, high], remove_xticks, n_obs)
        top_ax.set_yticklabels(["0", str(high)])
        # NB lettered as truth_combined.png's panels, in the left margin above each: (a) the switches, (b) the tracks.
        for letter, ax in (("(a)", top_ax), ("(b)", axes[0])):
            fig.text(
                0.01,
                1.0,
                letter,
                ha="left",
                va="bottom",
                transform=blended_transform_factory(fig.transFigure, ax.transAxes),
            )
    x_vals = np.arange(n_obs)

    for s, c in enumerate(nonempty_clones):
        cid = final_clone_ids[c]
        ax_idx = s * axes_per_clone
        ax_rdr = axes[ax_idx] if has_rdr else None
        ax_baf = axes[ax_idx + 1] if has_rdr else axes[ax_idx]

        if df_cnv is not None:
            # NB phased, as cnv_<ploidy>_seglevel holds them: B may exceed A. A pair takes its (larger, smaller)'s colour, the
            #    palette's, and the legend keys it with its mirror.
            phased = df_cnv[[f"clone{cid} A", f"clone{cid} B"]].to_numpy().astype(int)
            state_tuples = pd.Series(zip(phased.max(axis=1), phased.min(axis=1)))
            hue_indices = state_tuples.map(map_cn).fillna(default_idx).astype(int).values

            if palette_name in ("chisel", "tracks"):
                palette = [
                    mcolors.to_rgba(
                        color,
                        alpha=(NORMAL_OPACITY if ordered_acn[i] == (1, 1) else 1.0),
                    )
                    for i, color in enumerate(state_colors)
                ]
            else:
                palette = [mcolors.to_rgba(color, alpha=1.0) for color in state_colors]

            point_colors = np.array(palette)[hue_indices]
            # NB the dots are phased where the key is not: a bin decoded B > A takes its pair's colour at FLIP_OPACITY.
            if palette_name == "tracks":
                point_colors[:, 3] = np.where(phased[:, 1] > phased[:, 0], FLIP_OPACITY, point_colors[:, 3])

        elif res_combine is not None:
            if res_combine["pred_cnv"].ndim == 1 or res_combine["pred_cnv"].shape[1] == 1:
                this_pred = res_combine["pred_cnv"][(c * n_obs) : (c * n_obs + n_obs)].flatten() % n_states
            else:
                this_pred = res_combine["pred_cnv"][:, c] % n_states

            assert len(this_pred) == n_obs, f"Clone {cid} copy states defined for {len(this_pred)}, expected {n_obs}."

            point_colors = np.array(palette)[this_pred]
            unique_hues = np.unique(this_pred)

        else:
            point_colors = "#4C72B0" if colors is None else colors[c]

        # --- RDR Plotting ---
        if has_rdr:
            y_vals_rdr = X[:, 0, c] / base_nb_mean[:, c]

            if plot_rdr_errors == "poisson":
                with np.errstate(divide="ignore", invalid="ignore"):
                    std_err_rdr = np.sqrt(X[:, 0, c]) / base_nb_mean[:, c]
                    std_err_rdr[~np.isfinite(std_err_rdr)] = 0.0

                ax_rdr.errorbar(
                    x_vals,
                    y_vals_rdr,
                    yerr=std_err_rdr,
                    fmt="none",
                    ecolor=point_colors,
                    elinewidth=0.5,
                    zorder=0,
                )

            # Replaced sns.scatterplot with raw matplotlib scatter for speed
            ax_rdr.scatter(
                x_vals,
                y_vals_rdr,
                s=pointsize,
                c=point_colors,
                edgecolors="none",
                linewidth=linewidth,
                zorder=1,
            )

            _format_track_axis(
                ax_rdr,
                "\nRDR",
                [RDR_FLOOR, rdr_ylim],
                np.arange(0, rdr_ylim + 1.0, 1.0),
                remove_xticks,
                n_obs,
            )

        # --- BAF Plotting ---
        baf_vals = X[:, 1, c] / total_bb_RD[:, c]

        if plot_baf_errors == "beta":
            k, n_trials = X[:, 1, c], total_bb_RD[:, c]
            alpha_param, beta_param = k + 1, n_trials - k + 1
            alpha_beta_sum = alpha_param + beta_param
            std_err_baf = np.sqrt((alpha_param * beta_param) / (np.square(alpha_beta_sum) * (alpha_beta_sum + 1)))

            ax_baf.errorbar(
                x_vals,
                baf_vals,
                yerr=std_err_baf,
                fmt="none",
                ecolor=point_colors,
                elinewidth=0.5,
                zorder=0,
            )

        ax_baf.scatter(
            x_vals,
            baf_vals,
            s=pointsize,
            c=point_colors,
            edgecolors="none",
            zorder=1,
        )

        _format_track_axis(
            ax_baf,
            "\nBAF",
            [-0.05, 1.05],
            np.arange(0.0, 1.1, 0.2),
            remove_xticks,
            n_obs,
            grid=BAF_GRID,
        )

        # --- Viterbi Segments ---
        if res_combine is not None:
            # NB the clone's library normalizer: its counts are expected at base_nb_mean * mu / Z_c.
            shift = res_combine["new_log_mu_shift"] if "new_log_mu_shift" in res_combine.keys() else None
            z = 1.0 if shift is None else float(np.exp(np.asarray(shift, dtype=np.float64).reshape(-1)[c]))
            fraction = None if df_cnv is None or normal_fractions is None else normal_fractions.get(f"clone {cid}")

            if fraction is not None:
                # NB the decoded (A, B) per bin, phased as cnv_<ploidy>_seglevel holds them, mixed at the clone's normal fraction.
                copies = np.stack([df_cnv[f"clone{cid} A"].values, df_cnv[f"clone{cid} B"].values], axis=-1)
                depth, b_fraction = mixture(copies, float(fraction))
                segments, _ = get_intervals(copies[:, 0] * (copies[:, 1].max() + 1) + copies[:, 1])
                y_rdrs = [depth[seg[0]] / z for seg in segments]
                y_bafs = [b_fraction[seg[0]] for seg in segments]
            else:
                if res_combine["pred_cnv"].ndim == 1 or res_combine["pred_cnv"].shape[1] == 1:
                    this_pred = (
                        res_combine["pred_cnv"][(c * n_obs) : (c * n_obs + n_obs)].flatten()
                        % res_combine["new_log_mu"].shape[0]
                    )
                else:
                    n_states = res_combine["new_log_mu"].shape[0]
                    this_pred = res_combine["pred_cnv"][:, c] % n_states
                segments, run_states = get_intervals(this_pred)
                exp_log_mu = np.exp(clone_states(res_combine["new_log_mu"], c)) / z
                p_binom_arr = clone_states(res_combine["new_p_binom"], c)
                y_rdrs = [exp_log_mu[lbl] for lbl in run_states]
                y_bafs = [p_binom_arr[lbl] for lbl in run_states]

            rdr_lines, baf_major_lines, baf_minor_lines = [], [], []

            for seg, y_rdr, y_baf in zip(segments, y_rdrs, y_bafs):
                x_start, x_end = seg[0], seg[-1]

                if has_rdr:
                    rdr_lines.append([(x_start, y_rdr), (x_end, y_rdr)])

                baf_major_lines.append([(x_start, y_baf), (x_end, y_baf)])
                baf_minor_lines.append([(x_start, 1.0 - y_baf), (x_end, 1.0 - y_baf)])

            # NB each segment's expectation and its phase flip `1 - p`, black over the dots.
            line = {"colors": "k", "linewidths": 0.5, "zorder": 2}
            if has_rdr and rdr_lines:
                ax_rdr.add_collection(LineCollection(rdr_lines, gid="segments", **line))

            if baf_major_lines:
                ax_baf.add_collection(LineCollection(baf_major_lines, gid="segments", **line))
                ax_baf.add_collection(LineCollection(baf_minor_lines, gid="segments", **line))

        # NB the diploid reference, dashed: RDR 1 and BAF 1/2.
        if has_rdr:
            ax_rdr.axhline(1.0, color="k", linewidth=0.5, linestyle="--", zorder=2, gid="reference")
        ax_baf.axhline(0.5, color="k", linewidth=0.5, linestyle="--", zorder=2, gid="reference")

        # --- Legend & Annotations ---
        if df_cnv is not None or res_combine is not None:
            if df_cnv is not None:
                legend_elements = _state_legend(
                    phased, map_cn, default_idx, palette, normal=cid == 0, drop_unused=drop_unused_states
                )
            else:
                legend_elements = [
                    Line2D([0], [0], marker="o", color="w", markerfacecolor=palette[i], markersize=10, linestyle="None",
                           label=f"{100. * np.mean(this_pred == i):.0f}% ")
                    for i in unique_hues
                    if not (drop_unused_states and reads_zero(100.0 * np.mean(this_pred == i)))
                ]  # fmt: skip

            ax_legend = ax_rdr if has_rdr else ax_baf
            # NB a tumour clone all (1, 1) keys nothing, the normal clone's key holding (1, 1).
            if legend_elements:
                ax_legend.legend(
                    handles=legend_elements,
                    loc="upper right",
                    bbox_to_anchor=(1, 1.25),
                    ncol=len(legend_elements),
                    frameon=False,
                    bbox_transform=ax_legend.transAxes,
                )

        t_prop = tumor_prop[c] if (single_tumor_prop is not None and tumor_prop is not None) else None
        _annotate_clone_stats(
            ax_rdr if has_rdr else ax_baf,
            cast_clone_label(str(cid)) if labels is None else labels[c],
            spots_per_clone[c],
            np.sum(X[:, 0, c]),
            np.sum(total_bb_RD[:, c]),
            t_prop,
            paired_ax=ax_baf if has_rdr else None,
        )

    if genome is None:
        _draw_chromosome_boundaries(
            ([top_ax] if top_ax is not None else []) + list(axes), lengths, unique_chrs, chrtext_shift
        )
    else:
        assert genome.n_bins == n_obs, f"the genomic axis has {genome.n_bins} bins, the tracks {n_obs}"
        for ax in axes:
            genome.frame(ax, labels=ax is axes[-1])

    if own_figure:
        fig.tight_layout()

    return fig
