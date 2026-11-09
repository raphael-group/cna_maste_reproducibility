"""The benchmark's three PDFs' pages: a method's title, a profile on the cnamaste run's bins, the metric pages, and
each figure written to its PDF and each page alone.

- `get_title`, `get_calicost_label` -- a method's title on a page, CalicoST's with its label;
- `get_labelled_states`, `get_bin_profile` -- a state table's states as its key names them, and a profile on bins;
- `plot_simulation_pages` -- each simulation table's metrics, combined.pdf's first pages;
- `plot_metric_pages` -- the appendix's metric figures, a page each, from the run's tables;
- `plot_simulation_diagnosis` -- the simulation pages in full, every metric and simulated class, for diagnosis.pdf;
- `write_pages`, `write_appendix` -- a figure's pages to one PDF and each alone, a paper figure's page under its
  name (`PAGE_NAMES`).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from scripts.metrics.gene_level import get_max_overlap
from scripts.metrics.outputs import get_page_path, get_plot_path, get_table_path
from scripts.metrics.read_method_input import (
    NAMES,
    get_calicost_ploidy,
    get_calicost_run,
)
from scripts.metrics.simulation_metrics import get_level_table
from scripts.platforms.visium import get_visiumhd_bins
from scripts.plots.layout import set_text_sizes
from scripts.plots.run_comparison import plot_cna_jaccard, plot_genome_composition
from scripts.plots.simulation import (
    drop_simulated_classes,
    plot_simulation_metrics,
    plot_simulation_resources,
    plot_solver_ablation,
)
from scripts.plots.style import get_run_labels, get_run_order, get_setting_label

from cnamaste.config import start_time
from cnamaste.figure_style import DPI, page_style, restyle
from cnamaste.logger import get_logger

logger = get_logger("pages", start_time=start_time)


PAGE_NAMES = {
    "cnamaste": "main-fig-cnamaste-result",
    "labels": "suppl-fig-ovarian-label",
    "profiles": "main-fig-cnp-compare",
    "umi": "suppl-fig-umi",
    "composition": "suppl-fig-cna-composition",
    "state_jaccard": "suppl-fig-cna-jaccard",
    "simulation": "main-fig-zenodo-sim-result",
    "runtime": "suppl-fig-zenodo-sim-runtime",
    "ablation": "suppl-fig-zenodo-sim-ablation",
}
"""Each paper figure's page alone, `pdf-pages/<name>.pdf`, set as its figure's label (`Figure.set_label`) and written
in place of `<stem>-P<k>.pdf` (`outputs.get_page_path`): combined.pdf's first simulation table's metrics, cnamaste
page and profiles; the appendix's spatial labels, UMI, cells per altered state, state Jaccard index, and the first
simulation table's runtime and memory, and cnamaste solvers."""


SIMULATION_LEVEL = "gene"
"""The level the simulation pages draw (`metrics.simulation_metrics.get_levels`): the genes; the wide table holds
both."""


SIMULATION_COLUMNS = ("sensitivity", "clone_ari", "state_altered", "exact_altered")
"""The simulation page's metrics; diagnosis.pdf draws every one (`plot_simulation_diagnosis`)."""


DROPPED_CLASSES = (("6.3", 5e7), ("3.3", 5e7))
"""The simulated classes, `(numcnas, cnasize)`, the simulation and solver pages leave out, (6,3, 50M) and
(3,3, 50M); diagnosis.pdf draws every class (`plot_simulation_diagnosis`)."""


PAGE_PAD = 0.02
"""Inches of white kept round a page cropped to its content: the lab saves a figure with `bbox_inches="tight"`, and
LaTeX places a PDF by its bounding box (lab figure rules MP-3, LF-3)."""


SAVE_OPTIONS = {"dpi": DPI, "facecolor": "white", "bbox_inches": "tight", "pad_inches": PAGE_PAD}
"""How every page is saved, to a PDF of pages and alone: images at `DPI`, on white, cropped to its content."""


def get_page_name(page: str) -> str:
    """Get a paper figure's page name, `PAGE_NAMES`': `cnamaste`, `labels`, `profiles`, `umi`, `composition`,
    `state_jaccard`, `simulation`, `runtime` or `ablation`."""
    return PAGE_NAMES[page]


def get_calicost_label(directory: Path, label: str, ploidy: str, fallbacks: list[str]) -> str:
    """Get a CalicoST run's label in its title: `label`, as `style.get_setting_label` writes it, then the table
    it is read at where that is a fallback's.

    Args:
        directory: the run's or its sample's directory.
        label: its sweep label, empty for none.
        ploidy: the table asked for, `calicost_ploidy`.
        fallbacks: `calicost_fallback_ploidy`, in order.

    Returns:
        The label, e.g. `min-snp-umi: 10, diploid`.
    """
    chosen = get_calicost_ploidy(get_calicost_run(directory), ploidy, fallbacks)
    return ", ".join(part for part in (get_setting_label(label), chosen if chosen != ploidy else "") if part)


def get_title(method: str, labels: pd.Series, label: str = "") -> str:
    """Get a method's title on the page: its name, then the side of the bins it ran on and its label, in brackets.

    The name is `read_method_input.NAMES`'; the side is read from its barcodes
    (`platforms.visium.get_visiumhd_bins`), e.g. `CalicoST (16 um, min-snp-umi: 10)`.

    Args:
        method: a key of `NAMES`.
        labels: its clone per barcode, at the bins it ran on.
        label: the run's label, or empty.

    Returns:
        The title.
    """
    sides = [rf"{side} $\mu$m" for side in np.unique(get_visiumhd_bins(labels.index)[:, 0])]
    return f"{NAMES[method]} ({', '.join([*sides, label] if label else sides)})"


def get_labelled_states(profile: pd.DataFrame, labels: dict[str, str]) -> pd.DataFrame:
    """Get a state table's segments and each clone's state, `clone<k>`, as its key names it.

    Args:
        profile: `CHR`, `START`, `END`, then `clone<k>` per clone and any other columns, e.g. `clone<k> Uest`.
        labels: state -> its name in the key, e.g. `NUMBAT_STATE_LABELS`.

    Returns:
        `CHR`, `START`, `END` and `clone<k>` per clone.
    """
    clones = [c for c in profile.columns if c.startswith("clone") and c.removeprefix("clone").isdigit()]
    return profile[["CHR", "START", "END"]].join(profile[clones].apply(lambda states: states.map(labels)))


def get_bin_profile(profile: pd.DataFrame, bins: pd.DataFrame) -> pd.DataFrame:
    """Get a run's profile on bins: each bin its clones' columns of the segment overlapping it most
    (`metrics.gene_level.get_max_overlap`), NaN where none does.

    Args:
        profile: `CHR`, `START`, `END` per segment, then its clones' columns, e.g. `clone<k> A` and `clone<k> B`.
        bins: `CHR`, `START`, `END` per bin, e.g. a cnamaste run's `cnv_seglevel.tsv`.

    Returns:
        `CHR`, `START`, `END` per bin, then the profile's other columns.
    """
    best = get_max_overlap(bins, profile)
    values = profile.drop(columns=["CHR", "START", "END"]).iloc[np.maximum(best, 0)].reset_index(drop=True)
    values = values.where(np.broadcast_to((best >= 0)[:, None], values.shape))
    return pd.concat([bins[["CHR", "START", "END"]].reset_index(drop=True), values], axis=1)


def write_pages(figures: list[Any], outdir: Path, name: str) -> list[Path]:
    """Write `figures`, a page each, to `outputs.get_plot_path`'s `name`, and each page alone to a PDF of its own
    (`outputs.get_page_path`), each cropped to its content (`SAVE_OPTIONS`), under the page's face
    (`figure_style.page_style`), which a text's generic family resolves to when written.

    A `.pdf` `name` holds every page in order; in another format, page k after the first is `<stem>_page<k><suffix>`.
    A page's own PDF is named by its figure's label where it has one (`PAGE_NAMES`).

    Args:
        figures: the pages, in order, each closed once written.
        outdir: the outputs' directory.
        name: the figure's file name, e.g. `combined.pdf`.

    Returns:
        The written paths, the figure's, then each page's PDF.
    """
    output = get_plot_path(outdir, name)
    output.parent.mkdir(parents=True, exist_ok=True)
    pages = [get_page_path(outdir, name, k, figure.get_label()) for k, figure in enumerate(figures, 1)]
    if output.suffix.lower() == ".pdf":
        with page_style(), PdfPages(output) as pdf:
            for page, figure in zip(pages, figures, strict=True):
                pdf.savefig(figure, **SAVE_OPTIONS)
                figure.savefig(page, **SAVE_OPTIONS)
                plt.close(figure)
        return [output, *pages]
    paths = [output, *(output.with_name(f"{output.stem}_page{k}{output.suffix}") for k in range(2, len(figures) + 1))]
    with page_style():
        for path, page, figure in zip(paths, pages, figures, strict=True):
            figure.savefig(page, **SAVE_OPTIONS)
            figure.savefig(path, **SAVE_OPTIONS)
            plt.close(figure)
    return [*paths, *pages]


def plot_simulation_pages(outdir: Path, simulation_metrics: list[str]) -> list[Figure]:
    """Draw each `simulation_metrics` table's `SIMULATION_COLUMNS` at `SIMULATION_LEVEL`, the genes, less
    `DROPPED_CLASSES` (`plot_simulation_metrics`, `read_simulation_levels`), the first page labelled its paper figure's
    name (`PAGE_NAMES`): combined.pdf's first pages.

    Args:
        outdir: the outputs' directory, the tables in its `tsvs/` (`outputs.get_table_path`).
        simulation_metrics: the simulation wide tables' names, e.g. `calicost-sim`.

    Returns:
        A page per table, in order.
    """
    levels = read_simulation_levels(outdir, simulation_metrics)
    pages = [
        plot_simulation_metrics(drop_simulated_classes(level, DROPPED_CLASSES), columns=SIMULATION_COLUMNS)
        for level in levels
    ]
    if pages:
        pages[0].set_label(PAGE_NAMES["simulation"])
    return pages


def plot_metric_pages(
    outdir: Path, sample: str, real_data_metrics: str | None, simulation_metrics: list[str]
) -> list[Figure]:
    """Draw the appendix's metric figures from their tables under `outdir`, a page each (`simulation`,
    `run_comparison`).

    With `real_data_metrics`, the sample's runs against the reference run, `analysis.real_data`'s wide table
    and detail tables: their proportion of cells per altered state (`plot_genome_composition`), then each run's state
    Jaccard index against the reference run's, from their confusion (`plot_cna_jaccard`), the runs in method order
    (`style.get_run_order`), each run named by its method and the side of the bins it ran on, the wide table's
    `resolution` (`style.get_run_labels`), each page labelled its paper figure's name (`PAGE_NAMES`). Then the first
    `simulation_metrics` table's runtime and memory at `SIMULATION_LEVEL` (`plot_simulation_resources`,
    `read_simulation_levels`) and, where it holds two or more cnamaste solvers, their metrics less `DROPPED_CLASSES`
    (`plot_solver_ablation`), `analysis.simulation`'s, each labelled its paper figure's name; the
    tables' metrics are combined.pdf's (`plot_simulation_pages`).

    Args:
        outdir: the outputs' directory, the tables in its `tsvs/` (`outputs.get_table_path`).
        sample: the sample, the real-data tables' prefix.
        real_data_metrics: the real-data wide table's name, e.g. `ovarian`, or None for no real-data pages.
        simulation_metrics: the simulation wide tables' names, e.g. `calicost-sim`.

    Returns:
        The pages, in order.
    """
    pages = []
    if real_data_metrics:
        wide = pd.read_csv(get_table_path(outdir, real_data_metrics), sep="\t")
        composition, cells = (
            pd.read_csv(get_table_path(outdir, f"{sample}.{name}"), sep="\t")
            for name in ("genome_composition", "cna_confusion")
        )
        reference, run_names = str(wide["reference"].iloc[0]), get_run_labels(wide)
        confusions, levels = {}, {}
        compared = set(cells["run"])
        for run in (run for run in get_run_order(wide) if run in compared):
            table = cells[cells["run"] == run]
            order = list(dict.fromkeys(table["reference_state"]))
            pivot = table.pivot(index="reference_state", columns="run_state", values="cells")
            confusions[run], levels[run] = pivot.reindex(index=order, columns=order), str(table["level"].iloc[0])
        named = {
            "composition": plot_genome_composition(composition, sample, run_names),
            "state_jaccard": plot_cna_jaccard(confusions, reference, run_names, levels),
        }
        for name, figure in named.items():
            figure.set_label(PAGE_NAMES[name])
        pages += list(named.values())
    levels = read_simulation_levels(outdir, simulation_metrics)
    if levels:
        named = {"runtime": plot_simulation_resources(levels[0])}
        if levels[0]["solver"].nunique() > 1:
            named["ablation"] = plot_solver_ablation(drop_simulated_classes(levels[0], DROPPED_CLASSES))
        for name, figure in named.items():
            figure.set_label(PAGE_NAMES[name])
        pages += list(named.values())
    return pages


def read_simulation_levels(outdir: Path, simulation_metrics: list[str]) -> list[pd.DataFrame]:
    """Read each simulation wide table at `SIMULATION_LEVEL` (`metrics.simulation_metrics.get_level_table`).

    Args:
        outdir: the outputs' directory, the tables in its `tsvs/` (`outputs.get_table_path`).
        simulation_metrics: the simulation wide tables' names, e.g. `calicost-sim`.

    Returns:
        Each table's rows at the level, in order.
    """
    tables = [pd.read_csv(get_table_path(outdir, name), sep="\t") for name in simulation_metrics]
    return [get_level_table(table, SIMULATION_LEVEL) for table in tables]


def plot_simulation_diagnosis(outdir: Path, simulation_metrics: list[str]) -> list[Figure]:
    """Draw the simulation pages in full, every simulated class: each table's metrics, all of
    `simulation.METRIC_PANELS` (`plot_simulation_metrics`), and, where the first holds two or more cnamaste solvers,
    their metrics (`plot_solver_ablation`).

    Args:
        outdir: the outputs' directory, the tables in its `tsvs/` (`outputs.get_table_path`).
        simulation_metrics: the simulation wide tables' names, e.g. `calicost-sim`.

    Returns:
        The pages, in order, unlabelled.
    """
    levels = read_simulation_levels(outdir, simulation_metrics)
    pages = [plot_simulation_metrics(level) for level in levels]
    if levels and levels[0]["solver"].nunique() > 1:
        pages.append(plot_solver_ablation(levels[0]))
    return pages


def write_appendix(pages: list[Figure], metric_pages: list[Figure], outdir: Path, name: str) -> Path:
    """Write the appendix to one PDF, `outputs.get_plot_path`'s `name`, a page per figure in order, and each page
    alone to a PDF of its own (`outputs.get_page_path`), named by its figure's label where it has one (`PAGE_NAMES`),
    each cropped to its content (`SAVE_OPTIONS`): `pages` in the page's face (`figure_style.restyle`), their axis labels
    the largest text (`layout.set_text_sizes`), then `metric_pages` as `plots` sized their texts.

    Args:
        pages: the pages drawn here, e.g. the UMI page (`plot_umi_page`); each is closed once written.
        metric_pages: the pages drawn at their sizes, e.g. the spatial labels and the metric figures; each is
            closed once written.
        outdir: the outputs' directory.
        name: the appendix's file name, e.g. `appendix.pdf`; its directory is created.

    Returns:
        The appendix's path.
    """
    output = get_plot_path(outdir, name)
    output.parent.mkdir(parents=True, exist_ok=True)
    with page_style() as size, PdfPages(output) as pdf:
        for figure in pages:
            restyle(figure)
            set_text_sizes(figure, size)
        written = []
        for k, figure in enumerate([*pages, *metric_pages], 1):
            written.append(get_page_path(outdir, name, k, figure.get_label()))
            pdf.savefig(figure, **SAVE_OPTIONS)
            figure.savefig(written[-1], **SAVE_OPTIONS)
            plt.close(figure)
    logger.info(f"Wrote {output}: {len(written)} pages, each alone: {', '.join(str(path) for path in written)}.")
    return output
