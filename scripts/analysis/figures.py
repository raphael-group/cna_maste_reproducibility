"""Step 4 of the analysis: the benchmark's three PDFs, from the tables steps 1 and 2 wrote and every run's clones and
copies, read once.

Runpeng Luo (2026-11-09), from `scripts/run_benchmark.py` (2026-11-04)

Outputs (under outdir):
    plots/<combined>: `combined.pdf` by default, `.pdf` or `.png` (`plots.pages.write_pages`). Page 1, the
        simulations' metrics at the genes, less the classes `plots.pages.DROPPED_CLASSES`
        (`plots.pages.plot_simulation_pages`); then, where the cnamaste run wrote its bins, cnamaste's own page
        (`plots.cnamaste_page.plot_cnamaste_page`); then the clone profiles of cnamaste, CalicoST, Numbat and
        Clonalscope on the cnamaste run's bins (`plots.pages.get_bin_profile`), beside the CNP similarity and the
        clone ARI of every two runs and each run's skill scores (`plots.run_comparison.get_run_comparison_panels`).
        Last, the page naming every panel.
    plots/<diagnosis>: `diagnosis.pdf` by default: the clones of Numbat, the `false_positive_sweep` runs and
        cnamaste, each with its false-positive tumor spots; with a `bb_dir` at `side_um`, per run, its false-positive
        spots' total and SNP-covering UMI against its other called spots' (`plots.umi.plot_false_positive_umi`); each
        run's copy-number error (`plots.cn_error.plot_cn_error`), and with `calicost_sweep` the sweep runs' and their
        clones; the simulation pages in full, every metric and simulated class, and cnamaste's solvers
        (`plots.pages.plot_simulation_diagnosis`); last, the page naming every panel.
    tsvs/<diagnosis stem>.tsv.gz: the diagnosis' UMI page's numbers, a row per run and UMI measure
        (`metrics.qc_metrics.get_umi_comparison`).
    plots/<appendix>: `appendix.pdf` by default (`plots.pages.write_appendix`): with `bb_dir`, one row of
        per-spot total and SNP-covering UMI histograms per side (`plots.umi.plot_umi_page`); then the H&E slide and
        cell types, then the clones of cnamaste, CalicoST, Numbat and Clonalscope, three a row
        (`layout.plot_annotation_page`); then the metric figures, a page each (`plots.pages.plot_metric_pages`): the
        runs' cells per altered state, each run's state Jaccard index against the reference run's, the simulations'
        runtime and memory, and cnamaste's solvers, less the classes `plots.pages.DROPPED_CLASSES`. Last, the page
        naming every panel. Its pages are untitled.
    pdf-pages/<stem>-P<k>.pdf: page k of each of the three figures alone, written with the figure; a paper figure's
        page as `<name>.pdf` instead (`plots.pages.PAGE_NAMES`).
Notes/References:
    Every page of the three PDFs shows each panel's letter alone. Each method is named by the side of the bins it ran
    on, read from its barcodes (`plots.pages.get_title`).
    matplotlib's multi-page PDF: https://matplotlib.org/stable/api/backend_pdf_api.html
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scripts.analysis.options import AnalysisOptions
from scripts.metrics.outputs import get_table_path, write_table
from scripts.metrics.qc_metrics import get_false_positive_tumor
from scripts.metrics.read_method_input import (
    NAMES,
    NUMBAT_STATE_COLOURS,
    NUMBAT_STATE_LABELS,
    Run,
    get_epithelial,
    get_labelled_profile,
    get_normal_type,
    get_spot_labels,
    get_spot_umis,
    read_broad_cell_types,
    read_calicost,
    read_clonalscope,
    read_numbat,
    read_spots,
)
from scripts.platforms.visium import get_visiumhd_bin_pixels
from scripts.plots.cn_error import get_copy_number_rows, plot_cn_error
from scripts.plots.cnamaste_page import get_normal_fractions, plot_cnamaste_page, read_clone_tracks
from scripts.plots.layout import (
    FALSE_POSITIVE_COLOUR,
    get_genome_axis,
    plot_annotation_page,
    plot_panel_titles,
    plot_profile_page,
)
from scripts.plots.pages import (
    get_bin_profile,
    get_calicost_label,
    get_labelled_states,
    get_page_name,
    get_title,
    plot_metric_pages,
    plot_simulation_diagnosis,
    plot_simulation_pages,
    write_appendix,
    write_pages,
)
from scripts.plots.run_comparison import get_run_comparison_panels
from scripts.plots.style import get_run_labels
from scripts.plots.umi import plot_false_positive_umi, plot_umi_page

from cnamaste.config import start_time
from cnamaste.he import slide_image_and_meta
from cnamaste.io import get_spatial_positions
from cnamaste.logger import get_logger
from cnamaste.platforms.visium import get_visiumhd_microns_per_pixel, pixel_physical
from cnamaste.qa.metrics import read_labels

logger = get_logger("scripts.analysis.figures", start_time=start_time)

CLONALSCOPE = "clonalscope"
"""The method drawn as a placeholder, named by `NAMES`, where no `clonalscope_dir` is given."""


CNA_KEY_ORDER = ("neu", "del", "bdel", "amp", "bamp", "loh")
"""The `CNA` key's states, left to right, those drawn shown: NEU, DEL, AMP, CNLoH, each balanced state after its
unbalanced one."""


FONT_TYPES = {"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"}
"""Text written as editable text: TrueType fonts embedded in a PDF or PS, in place of matplotlib's Type 3, and text
kept text in an SVG (https://matplotlib.org/stable/users/explain/customizing.html)."""


def write_figures(options: AnalysisOptions) -> None:
    """Read every method's labels and copies once, place the spots, and draw the three PDFs, their text editable
    (`FONT_TYPES`), the metric pages from the run's tables.

    Args:
        options: the analysis' options.

    Raises:
        ValueError: where the runs contradict the options: two `calicost_sweep` runs of one label, or a tumour clone
            of the cnamaste fit with no normal fraction.
    """
    logger.info("Step 4 of 4: the figures, plots/ and pdf-pages/.")
    plt.rcParams.update(FONT_TYPES)

    broad = read_broad_cell_types(options.celltype)
    epithelial = get_epithelial(broad)
    fallbacks = [] if "none" in options.calicost_fallback_ploidy else list(options.calicost_fallback_ploidy)
    calicost_label, calicost_dir = options.calicost_dir
    calicost = read_calicost(calicost_dir, options.calicost_ploidy, fallbacks)
    numbat = read_numbat(options.numbat_dir)
    clonalscope = (
        None if options.clonalscope_dir is None else read_clonalscope(options.clonalscope_dir, merge_normal=True)
    )
    cnamaste_labels = read_labels(options.cnamaste_labels)
    numbat_title = get_title("numbat", numbat.labels)
    clonalscope_title = NAMES[CLONALSCOPE] if clonalscope is None else get_title(CLONALSCOPE, clonalscope.labels)
    calicost_title = get_title(
        "calicost",
        calicost.labels,
        get_calicost_label(calicost_dir, calicost_label, options.calicost_ploidy, fallbacks),
    )
    cnamaste_title = get_title("cnamaste", cnamaste_labels)
    sweep, sweep_runs = {}, {}
    for label, directory in options.calicost_sweep:
        run = read_calicost(directory, options.calicost_ploidy, fallbacks)
        title = get_title(
            "calicost", run.labels, get_calicost_label(directory, label, options.calicost_ploidy, fallbacks)
        )
        sweep[title], sweep_runs[title] = run, Run("calicost", label, directory)
    if len(sweep) != len(options.calicost_sweep):
        raise ValueError(f"calicost_sweep runs need distinct labels: {[label for label, _ in options.calicost_sweep]}")

    # NB each run read once: its clones per spot of side_um, in the appendix's bar order.
    spot_labels = {
        numbat_title: get_spot_labels(Run("numbat", "", options.numbat_dir), numbat.labels, options.side_um),
        **(
            {}
            if clonalscope is None
            else {
                clonalscope_title: get_spot_labels(
                    Run(CLONALSCOPE, "", options.clonalscope_dir), clonalscope.labels, options.side_um
                )
            }
        ),
        calicost_title: get_spot_labels(
            Run("calicost", calicost_label, calicost_dir), calicost.labels, options.side_um
        ),
        **{title: get_spot_labels(sweep_runs[title], run.labels, options.side_um) for title, run in sweep.items()},
        cnamaste_title: get_spot_labels(
            Run("cnamaste", "", options.cnamaste_labels.parent), cnamaste_labels, options.side_um
        ),
    }
    # NB the methods in the figures' order: the spatial grid three a row, H&E, cell type and cnamaste, then CalicoST,
    #    Numbat and Clonalscope; the profiles top to bottom.
    clones = {
        cnamaste_title: spot_labels[cnamaste_title].labels,
        calicost_title: spot_labels[calicost_title].labels,
        numbat_title: spot_labels[numbat_title].labels,
        clonalscope_title: None if clonalscope is None else spot_labels[clonalscope_title].labels,
    }
    sweep_clones = {title: spot_labels[title].labels for title in sweep}

    cnamaste_seglevel = pd.read_csv(options.cnamaste_seglevel, sep="\t")
    profiles = {
        cnamaste_title: get_labelled_profile(cnamaste_seglevel, cnamaste_labels),
        calicost_title: calicost.profile,
        numbat_title: get_labelled_states(numbat.profile, NUMBAT_STATE_LABELS),
        clonalscope_title: None
        if clonalscope is None
        else get_labelled_states(clonalscope.profile, NUMBAT_STATE_LABELS),
    }
    methods = {
        cnamaste_title: NAMES["cnamaste"],
        calicost_title: NAMES["calicost"],
        numbat_title: NAMES["numbat"],
        clonalscope_title: NAMES[CLONALSCOPE],
    }
    state_colours = {NUMBAT_STATE_LABELS[s]: NUMBAT_STATE_COLOURS[s] for s in CNA_KEY_ORDER}

    drawn = [labels for labels in [*clones.values(), *sweep_clones.values()] if labels is not None]
    barcodes = pd.Index(epithelial.index).append([pd.Index(c.index) for c in drawn]).unique()
    positions = get_spatial_positions(str(options.spaceranger_dir))
    microns = get_visiumhd_microns_per_pixel(options.spaceranger_dir)
    xy = pixel_physical(get_visiumhd_bin_pixels(barcodes, positions), microns)
    physical = pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1]}, index=barcodes)
    logger.info(f"Placed {len(physical):_} bins of {options.side_um} um from {options.spaceranger_dir}.")

    he_frame = slide_image_and_meta([options.spaceranger_dir])
    if he_frame is not None:
        he_xy = pixel_physical(he_frame[["x", "y"]].to_numpy(), microns)
        he_frame = he_frame.assign(x=he_xy[:, 0], y=he_xy[:, 1])

    # NB the diagnosis page: Numbat, the chosen sweep runs, cnamaste, each with its false-positive tumor spots.
    fp_clones = {
        numbat_title: clones[numbat_title],
        **{title: sweep_clones[title] for title in sweep if sweep_runs[title].label in options.false_positive_sweep},
        cnamaste_title: clones[cnamaste_title],
    }
    normal_type = broad.index[get_normal_type(broad).to_numpy()]
    false_positive = {title: get_false_positive_tumor(labels, normal_type) for title, labels in fp_clones.items()}
    for title, found in false_positive.items():
        logger.info(f"{title}: {len(found):_} false-positive tumor spots.")

    genome_axis = get_genome_axis(str(options.region_bed), str(options.chrom_sizes))
    # NB every profile on the cnamaste run's bins, as its own page draws its profile.
    bins = cnamaste_seglevel[["CHR", "START", "END"]]
    bin_profiles = {title: None if p is None else get_bin_profile(p, bins) for title, p in profiles.items()}
    for title, table in bin_profiles.items():
        if table is not None:
            uncalled = int(table.drop(columns=["CHR", "START", "END"]).isna().all(axis=1).sum())
            logger.info(f"{title}: {uncalled:_} of {len(bins):_} cnamaste bins overlap no segment.")
    wide = pd.read_csv(get_table_path(options.outdir, options.real_data_dataset), sep="\t")
    errors = get_copy_number_rows(wide)
    pairs = {
        metric: pd.read_csv(get_table_path(options.outdir, f"{options.sample}.{metric}"), sep="\t")
        for metric in ("cnp_similarity", "clone_ari")
    }
    comparison = get_run_comparison_panels(wide, pairs, errors)
    spots = read_spots(options.spots)
    # NB every page's clones, so a clone is one colour on each page.
    clone_ids = sorted(set().union(*(set(c.astype(int)) for c in drawn)))
    pages = plot_simulation_pages(options.outdir, [options.simulation_dataset])
    tracks = read_clone_tracks(options.cnamaste_labels.parent, cnamaste_labels)
    if tracks is None:
        logger.info(f"No binned_data.npz in {options.cnamaste_labels.parent}: combined.pdf has no cnamaste page.")
    else:
        fit_clones = sorted(np.unique(tracks.fit["new_assignment"]).astype(int).tolist())
        fractions = get_normal_fractions(
            options.cnamaste_labels.parent, fit_clones, dict(options.cnamaste_normal_fraction)
        )
        missing = sorted(set(fit_clones) - set(fractions))
        if missing:
            raise ValueError(f"cnamaste_normal_fraction: give the decoder's normal fraction of clone(s) {missing}")
        pages.append(
            plot_cnamaste_page(
                physical,
                he_frame,
                clones[cnamaste_title],
                cnamaste_seglevel,
                tracks,
                fractions,
                cnamaste_title,
                clone_ids=clone_ids,
                in_tissue=spots,
            )
        )
        pages[-1].set_label(get_page_name("cnamaste"))
    # NB the per-spot bins are held for that page alone.
    del tracks
    pages.append(
        plot_profile_page(
            clones,
            bin_profiles,
            methods,
            NUMBAT_STATE_LABELS["neu"],
            state_colours=state_colours,
            clone_ids=clone_ids,
            in_tissue=spots,
            square_panels=comparison[0],
            square_side=comparison[1],
            foot_panel=comparison[2],
        )
    )
    pages[-1].set_label(get_page_name("profiles"))
    pages.append(plot_panel_titles(pages))
    written = write_pages(pages, options.outdir, str(options.combined))
    logger.info(f"Wrote {len(pages)} pages: {', '.join(str(path) for path in written)}.")

    diagnosis = [
        plot_annotation_page(
            physical,
            None,
            fp_clones,
            {},
            genome_axis,
            clone_ids=clone_ids,
            false_positive=false_positive,
            panels_per_row=1,
        )
    ]
    umi_table = None
    side_bb_dir = dict(options.bb_dir).get(options.side_um)
    if side_bb_dir is None:
        logger.info(f"No bb_dir at {options.side_um} um: the diagnosis has no UMI page.")
    else:
        umis = get_spot_umis(side_bb_dir, options.assay, options.spots)
        groups = [(spot_labels[title], false_positive[title]) for title in fp_clones]
        umi_page, umi_table = plot_false_positive_umi(groups, umis, FALSE_POSITIVE_COLOUR, options.sample)
        diagnosis.append(umi_page)
        logger.info(f"False-positive spots' UMI against the other called spots':\n{umi_table.to_string(index=False)}")
    diagnosis.append(plot_cn_error(errors, f"{options.sample}, copy-number error", get_run_labels(errors)))
    if sweep:
        sweep_errors = pd.read_csv(get_table_path(options.outdir, f"{options.sample}.cn_error.sweep"), sep="\t")
        title = f"{options.sample}, copy-number error, CalicoST sweep"
        diagnosis.append(plot_cn_error(sweep_errors, title, get_run_labels(sweep_errors)))
        diagnosis.append(
            plot_annotation_page(
                physical,
                None,
                sweep_clones,
                {title: run.profile for title, run in sweep.items()},
                genome_axis,
                clone_ids=clone_ids,
            )
        )
    diagnosis += plot_simulation_diagnosis(options.outdir, [options.simulation_dataset])
    diagnosis.append(plot_panel_titles(diagnosis))
    written = write_pages(diagnosis, options.outdir, str(options.diagnosis))
    logger.info(f"Wrote {', '.join(str(path) for path in written)}.")
    if umi_table is not None:
        write_table(umi_table, get_table_path(options.outdir, options.diagnosis.stem), float_format="%.6g")

    appendix = [plot_umi_page(options.sample, options.bb_dir, options.assay, options.bins)] if options.bb_dir else []
    for page in appendix:
        page.set_label(get_page_name("umi"))
    labels_page = plot_annotation_page(
        physical,
        epithelial,
        clones,
        {},
        genome_axis,
        he_frame=he_frame,
        clone_ids=clone_ids,
        panels_per_row=3,
        in_tissue=spots,
        keys_inside=True,
        fit_page=True,
    )
    labels_page.set_label(get_page_name("labels"))
    metric_pages = [
        labels_page,
        *plot_metric_pages(options.outdir, options.sample, options.real_data_dataset, [options.simulation_dataset]),
    ]
    metric_pages.append(plot_panel_titles([*appendix, *metric_pages]))
    write_appendix(appendix, metric_pages, options.outdir, str(options.appendix))
