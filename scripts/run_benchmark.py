"""The benchmark, one run: the simulations' runs scored against their truth, a Visium HD sample's runs compared, their
LaTeX tables, and the benchmark's three PDFs.

Author: Runpeng Luo (2026-11-04), from `compute_metrics.py` (2026-11-01), `generate_tables.py` (2026-11-02) and
`plot_combined_pdf.py` (2026-10-30)

Dependencies: a cnamaste repository, `$CNAMASTE_ROOT`, its `python/` on `PYTHONPATH` (`get_cnamaste_root`) and its
`luo` dependency group, cnplot, for the profiles and tracks; `pdflatex` with the `geometry`, `amsmath`, `amssymb` and
`booktabs` packages (TeX Live); `metrics/`, `plots/` and `platforms/`, beside this script, imported as `scripts.*`, the
benchmark's directory on `PYTHONPATH`.

Usage:
    CNAMASTE_ROOT=<cnamaste> PYTHONPATH=<cnamaste>/python:<benchmark> \
        uv run --project <cnamaste> --group luo python <benchmark>/scripts/run_benchmark.py --sample NAME \
        [--benchmark_yaml YAML] [--config_yaml YAML] [--method METHOD ...] [--simulation_dataset calicost-sim] \
        [--run KEY ...] [--reference cnamaste] [--sweep_run KEY ...] [--msr 300] [--real_data_dataset ovarian] \
        [--table table] [--metric_text metric] --spaceranger_dir DIR --celltype FILE --spots FILE [--bb_dir SIDE_UM=DIR ...] \
        --cnamaste_labels FILE --cnamaste_seglevel FILE [--cnamaste_normal_fraction CLONE=F ...] \
        --calicost_dir [LABEL=]DIR [--calicost_sweep LABEL=DIR ...] [--false_positive_sweep LABEL ...] \
        --numbat_dir DIR [--clonalscope_dir DIR] --outdir DIR [--combined combined.pdf] [--diagnosis diagnosis.pdf] \
        [--appendix appendix.pdf]

    `scripts/run_benchmark.sh` sets both and runs it on H1-N3KXRGQ in a new `outdir`. Four steps, in order:
    the simulations' metrics (`write_simulation_metrics`), the sample's runs compared (`write_real_data_metrics`), the
    LaTeX table and metric text (`metrics.latex_tables.write_latex_table`, `metrics.metric_text.write_metric_text`),
    and the figures from those tables (`write_figures`).
    Each metric step computes every metric, merges the per-run ones into one wide table, and writes its tables only
    then.

Inputs:
    benchmark_yaml: the benchmark's `configs/benchmark.local.yaml` by default.
        simulation: its `simulation` entry (`metrics.simulation_metrics.read_simulation_config`):
            groundtruth_dir: a directory per simulation, its `truth_clone_labels.tsv` and `truth_acn_profile.tsv`
                (`truth.read_truth`).
            runtime: a `simulation.runtime.tsv`, `simulation_id`, `method`, `Runtime [mins]` and `Memory [GB]` per run.
            <method>.dir: a run's directory, a template of its `{sample_id}`, read by the method's reader
                (`metrics.read_method_input.get_method_reader`): Numbat's output directory, Clonalscope's Visium
                tutorial export, CalicoST's run (`clone_labels.tsv`, `cnv_seglevel.tsv`, its unconstrained copies),
                cnamaste's run (`clone_labels_int.tsv.gz`, `cnv_<ploidy>_seglevel.tsv`).
        real_data: per sample, its runs, `dir` per run key, its runs' runtime table, `runtime`
            (`metrics.real_data_metrics.read_runtime_table`), UGP's in-tissue spots, `ugp.spots`, its `bb` directory
            per bin side, `ugp.bb_dir`, and the spots' cell types, `cell_type.labels`
            (`metrics.read_method_input.read_broad_cell_types`). The `bb` directory at `side_um` holds UGP's bins,
            `MSR<msr>/<assay>/`, their UMI and phased allele counts per spot
            (`metrics.read_method_input.read_bin_counts`).
    config_yaml: `reference.bed_file`, the region bed, and `reference.refdata.gene_list.gtf_file`, the genes
        (`metrics.gene_level.read_reference_files`), the benchmark's `configs/config.local.yaml` by default.
    spaceranger_dir: a Space Ranger directory of the sample, any bin side: `spatial/tissue_positions.parquet`, which
        places every bin (`platforms.visium.get_visiumhd_bin_pixels`), `scalefactors_json.json` and the
        H&E image.
    celltype: per-spot cell types at `side_um`, `barcode` and `cell_type` (`annotation.load_celltype`),
        `<broad> (<fine>)` (`read_method_input.read_broad_cell_types`): drawn as epithelial or not, low-confidence spots
        left out; a spot of a normal type, neither epithelial nor low-confidence, that a method calls tumor is a
        false-positive tumor spot (`metrics.qc_metrics.get_false_positive_tumor`).
    spots: every in-tissue spot at `side_um`, a UGP `barcodes.tsv.gz` (`read_method_input.read_spots`): the
        `bb_dir` at `side_um`'s own, `<MSR>/<assay>/barcodes.tsv.gz`, its lines the allele matrices' columns.
    bb_dir: a UGP `bb` directory and its bins' side in um, `<side_um>=<dir>`, once per side: its `<assay>.h5ad` and
        `multi_snp/<assay>/bb.{A,B}allele.npz`, (bins, spots); the appendix's UMI page, where given.
    cnamaste_labels: a run's `clone_labels_int.tsv[.gz]`, its directory the run's; cnamaste_seglevel: its
        `cnv_<ploidy>_seglevel.tsv`, the profiles page its integer clones' alone
        (`read_method_input.get_labelled_profile`). With the run's `binned_data.npz` (`run.write_binned_data`) and
        `rdrbaf_final_nstates<n>_smp.npz` beside them (`plots.cnamaste_page.read_clone_tracks`), combined.pdf opens
        with the run's own page, its expected lines at each tumour clone's normal fraction.
    calicost_dir: a CalicoST sample directory, `[<label>=]<dir>`, read at its best rectangle
        (`read_method_input.read_calicost`): its `clone_labels.tsv` and `cnv<_ploidy>_seglevel.tsv` at
        `calicost_ploidy`; its label, e.g. its sweep's, in its title. A `calicost_sweep` run of the same label and
        directory is the same run.
    calicost_sweep: a CalicoST run of a parameter sweep, `<label>=<dir>`, once per run, read as `calicost_dir` is.
    numbat_dir: a Numbat output directory (`read_method_input.read_numbat`): its clones and each clone's CNV state per
        segment, its `bulk_clones_final.tsv.gz`, in Numbat's colours, `loh`'s and `del`'s swapped, and names
        (`NUMBAT_STATE_COLOURS`, `NUMBAT_STATE_LABELS`).
    clonalscope_dir: a Clonalscope Visium tutorial export, its `*.cells.tsv` and `*.profile.tsv`
        (`read_method_input.read_clonalscope`), its normal clusters one clone 0, so its tumor spots are its other
        clones, as every method's are; each tumor clone's total-copy state per segment, `del`, `neu` or `amp`, in
        Numbat's colours and names, keyed with Numbat's. Without it, Clonalscope is a placeholder.
Parameters:
    sample: the sample, a key under the benchmark yaml's platform, e.g. `vishd`; in the figures' titles and tables.
    method: a simulation method scored, once per method, every method of the `simulation` entry by default.
    simulation_dataset: the simulations' wide table's name, `calicost-sim` by default.
    run: a run key of the sample, once per run, `numbat`, `clonalscope`, `calicost` and `cnamaste` by default
        (`RUNS`); its method the key before its first `_` or `-`, e.g. `cnamaste` of `cnamaste-ac88605f`, read by its
        method's reader, Clonalscope's clusters as called, cnamaste's copies those of the clones its integer clones
        name (`metrics.read_method_input.get_labelled_profile`).
    reference: the run every other run is compared with, `cnamaste` by default.
    sweep_run: a run key of the sample, once per run, e.g. a CalicoST parameter sweep's run: the sweep runs and the
        reference run are scored on the copy-number error as one panel, on the bins they share; none by default.
    msr: UGP's `min_snp_reads` of the bins scored, their directory `MSR<msr>`, 300 by default.
    real_data_dataset: the sample's wide table's and its summary's name, `ovarian` by default.
    table: the LaTeX tables' name under `text/`, `table` by default.
    metric_text: the metric text's name under `text/`, `metric` by default.
    cnamaste_normal_fraction: a tumour clone of the run's fit and its decoder's normal fraction, `<clone>=<fraction>`,
        once per tumour clone; by default each integer clone's from the run's log, run_cnamaste's `normal fractions
        per clone` (`plots.cnamaste_page.get_normal_fractions`); the normal clone is held at 1, as the decoder holds
        it (`plots.cnamaste_page.NORMAL_CLONE_FRACTION`).
    false_positive_sweep: the label of a `calicost_sweep` run, once per run, on the diagnosis page; every run of that
        label, at each side.
    side_um: the side of the spots compared and drawn; a run on coarser bins gives each bin's clone to the spots
        inside it (`metrics.read_method_input.get_spot_labels`). 16 by default.
    assay: UGP's assay type, the bins' subdirectory and the stem of a `bb_dir`'s AnnData and allele matrices'
        directory, `VISIUM3prime` by default.
    calicost_ploidy: CalicoST's table, `unconstrained` its `cnv_seglevel.tsv`, which its `plots/acn_genome.pdf` draws.
    calicost_fallback_ploidy: the tables tried, in order, where a CalicoST run's `calicost_ploidy` table leaves out a
        labelled tumour clone: the first holding every labelled clone is read, else `calicost_ploidy`'s; `diploid
        tetraploid` by default, `none` for none (`read_method_input.get_calicost_ploidy`); the run's title names a
        fallback read.
    region_bed, chrom_sizes: the regions and chromosome lengths the diagnosis' sweep profiles' genome axis is drawn
        on, by default the benchmark's `data/` copies of UGP's resources/data/hg38.regions.bed
        and hg38.chrom.sizes: a chromosome's arms, its centromere between them.
    bins: each UMI histogram's bin count.
    outdir: a benchmark run's directory, `run_benchmark.sh`'s the benchmark's `output/benchmark_results/`, emptied
        before each run; tables in its `tsvs/`, LaTeX tables in its `text/`, figures in its `plots/`, each figure's
        pages in its `pdf-pages/`, each existing beforehand (`metrics.outputs`).
Outputs (under outdir):
    tsvs/<simulation_dataset>.tsv.gz: one row per simulation and method, each run scored on its own segments and at
        the genes, a metric's column per level, `<metric> (segment)` and `<metric> (gene)`
        (`metrics.simulation_metrics.merge_levels`).
    tsvs/<simulation_dataset>.median.tsv.gz: one row per simulated class and method, `sim_id`, `method`,
        `#replicates`, then each metric's median over the class' replicates, NaN skipped
        (`plots.simulation.get_class_medians`).
    tsvs/gene_level_cnp/<sample_id>.gene_level_cnp.tsv.gz: per gene, the truth's and each run's fields
        (`metrics.gene_level.get_gene_level_table`), written per simulation as it is scored.
    tsvs/<real_data_dataset>.tsv.gz: one row per run (`metrics.real_data_metrics.merge_run_tables`): `sample`, `run`,
        `method`, `resolution`; its spot coverage and false-positive counts per broad cell type
        (`metrics.real_data_metrics.get_spot_table`); its genome composition's `all` row, the share of its cells per
        CNA state and `uncalled` (`metrics.qc_metrics.compute_genome_composition`); the reference run, the field
        compared, `level`, their shared spots, clone ARI and CNA ARI (`metrics.real_data_metrics.compare_runs`), NaN
        for the reference; its RDR and BAF error at the bins, the runs giving copies
        (`metrics.cn_error.get_cn_error_table`); its runtime row's `setting`, runtime and memory (`metrics.real_data_metrics.get_runtime_table`).
    tsvs/<real_data_dataset>.json: its summary metrics per run, and their provenance
        (`metrics.real_data_metrics.get_summary`).
    tsvs/<sample>.<table>.tsv.gz: the sample's detail tables: `gene_level_cnp`, per gene, each run's fields;
        `genome_composition`, per run and clone, its genes' shares by CNA state and `uncalled`;
        `clone_ari` and `cnp_similarity`, per ordered pair of runs
        (`metrics.unsupervised_metrics.compute_pairwise_clone_ARI`, `compute_pairwise_cnp_similarity`), a run against
        itself 1; `cna_confusion`, per run, the cells by `(reference state, run state)` at its `level`
        (`metrics.cna_states.get_comparison_field`); `segment_level_cnp`, per joint segment of the runs giving copies
        (`metrics.segment_level.get_segment_level_table`); with `sweep_run`, `cn_error.sweep`, the copy-number error
        of the reference run and the sweep runs, scored as one panel on the bins they share.
    text/<table>.tex, text/<table>.pdf: the runs' spot coverage rate, false-positive rate, runtime and memory, a row per
        run in method order, a booktabs tabular in a captioned, labelled `table` float, typeset on letter pages, and
        pdflatex's log, `text/<table>.log` (`metrics.latex_tables`).
    text/<metric_text>.tex, text/<metric_text>.pdf: each metric's definition as it is computed, a paragraph per
        metric, each followed by a `% shown in:` line per page and panel drawing it, typeset on letter pages, and
        pdflatex's log (`metrics.metric_text`).
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
    The metrics are in the benchmark's `docs/pipeline.md` (Metrics) and `docs/metrics.md`. The regions
    are the bed's on chromosomes 1 to 22, where the truth and every method's input lie. The spots compared are UGP's
    in-tissue spots at `side_um` a run labels. The copy-number error's baseline is the spots of a normal cell type
    (`metrics.read_method_input.get_normal_type`); its `#tumor-clones` each run's clones but its normal ones
    (`metrics.real_data_metrics.get_normal_clones`).
    Every page of the three PDFs shows each panel's letter alone. Each method is named by the side of the bins it ran
    on, read from its barcodes (`plots.pages.get_title`).
    Numbat's outputs: https://kharchenkolab.github.io/numbat/articles/descriptions.html
    cnplot: https://github.com/raphael-group/cnplot
    matplotlib's multi-page PDF: https://matplotlib.org/stable/api/backend_pdf_api.html
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scripts.metrics.gene_level import (
    get_autosomal_regions,
    get_gene_level_table,
    get_gene_states,
    get_region_genes,
    read_genes,
    read_reference_files,
    read_regions,
)
from scripts.metrics.latex_tables import write_latex_table
from scripts.metrics.metric_text import write_metric_text
from scripts.metrics.outputs import (
    get_json_path,
    get_table_dir,
    get_table_path,
    write_json,
    write_table,
)
from scripts.metrics.qc_metrics import (
    compute_genome_composition,
    get_false_positive_tumor,
)
from scripts.metrics.read_method_input import (
    NAMES,
    NUMBAT_STATE_COLOURS,
    NUMBAT_STATE_LABELS,
    Run,
    get_benchmark_dir,
    get_calicost_fallbacks,
    get_epithelial,
    get_labelled_profile,
    get_method_reader,
    get_normal_type,
    get_run_argument,
    get_run_method,
    get_spot_labels,
    get_spot_umis,
    read_broad_cell_types,
    read_calicost,
    read_clonalscope,
    read_numbat,
    read_spots,
)
from scripts.metrics.real_data_metrics import (
    compare_runs,
    get_copy_number_tables,
    get_run_name,
    get_run_spots,
    get_runtime_table,
    get_spot_table,
    get_summary,
    get_tumor_cnps,
    merge_run_tables,
    read_sample_entry,
)
from scripts.metrics.simulation_metrics import (
    get_simulation_methods,
    get_simulation_metrics,
    merge_levels,
    read_simulation_config,
)
from scripts.metrics.unsupervised_metrics import (
    compute_pairwise_clone_ARI,
    compute_pairwise_cnp_similarity,
)
from scripts.platforms.visium import get_visiumhd_bin_pixels
from scripts.plots.cn_error import get_copy_number_rows, plot_cn_error
from scripts.plots.cnamaste_page import (
    get_cnamaste_input_mismatch,
    get_normal_fractions,
    plot_cnamaste_page,
    read_clone_tracks,
)
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
from scripts.plots.simulation import get_class_medians
from scripts.plots.style import get_run_labels
from scripts.plots.umi import plot_false_positive_umi, plot_umi_page

import cnamaste
from cnamaste.calicost.read import CALICOST_PLOIDIES
from cnamaste.config import start_time
from cnamaste.he import slide_image_and_meta
from cnamaste.io import get_spatial_positions
from cnamaste.logger import get_logger
from cnamaste.platforms.visium import get_visiumhd_microns_per_pixel, pixel_physical
from cnamaste.qa.metrics import read_labels

logger = get_logger("run_benchmark", start_time=start_time)

BENCHMARK = get_benchmark_dir()
"""The benchmark's directory, this script's parent's: `scripts/`, `configs/` and `data/`."""

DATA = BENCHMARK / "data"
"""The benchmark's `data/`: the genome axis' default files, found from this script whatever the working directory."""

CNAMASTE_ROOT = "CNAMASTE_ROOT"
"""The environment variable naming the cnamaste repository `cnamaste` is imported from (`run_benchmark.sh`)."""

RUNS = ("numbat", "clonalscope", "calicost", "cnamaste")
"""The sample's runs scored where no `--run` is given."""


CLONALSCOPE = "clonalscope"
"""The method drawn as a placeholder, named by `NAMES`, where no `--clonalscope_dir` is given."""


CNA_KEY_ORDER = ("neu", "del", "bdel", "amp", "bamp", "loh")
"""The `CNA` key's states, left to right, those drawn shown: NEU, DEL, AMP, CNLoH, each balanced state after its
unbalanced one."""


FONT_TYPES = {"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"}
"""Text written as editable text: TrueType fonts embedded in a PDF or PS, in place of matplotlib's Type 3, and text
kept text in an SVG (https://matplotlib.org/stable/users/explain/customizing.html)."""


def get_normal_fraction(text: str) -> tuple[int, float]:
    """Read a `--cnamaste_normal_fraction` value, `<clone>=<fraction>`.

    Args:
        text: the value, e.g. `1=0.38`.

    Returns:
        The clone and its normal fraction.

    Raises:
        argparse.ArgumentTypeError: where `text` is not `<int>=<fraction in [0, 1]>`.
    """
    clone, sep, fraction = text.partition("=")
    try:
        value = float(fraction)
    except ValueError:
        value = -1.0
    if not sep or not clone.isdigit() or not 0.0 <= value <= 1.0:
        raise argparse.ArgumentTypeError(f"expected <clone>=<fraction in [0, 1]>, got {text!r}")
    return int(clone), value


def get_bb_dir(text: str) -> tuple[int, Path]:
    """Read a `--bb_dir` value, `<side_um>=<dir>`.

    Args:
        text: the value.

    Returns:
        The bins' side in um and the directory.

    Raises:
        argparse.ArgumentTypeError: where `text` is not `<int>=<dir>`.
    """
    side, sep, directory = text.partition("=")
    if not sep or not side.isdigit() or not directory:
        raise argparse.ArgumentTypeError(f"expected <side_um>=<dir>, got {text!r}")
    return int(side), Path(directory)


def write_simulation_metrics(args: argparse.Namespace) -> None:
    """Score the simulations' runs on their own segments and at the genes, merge the two into one wide table, and
    write it and each method's median per simulated class (`plots.simulation.get_class_medians`).

    Args:
        args: the command line (`parser`).
    """
    simulation = read_simulation_config(args.benchmark_yaml)
    methods = args.methods or get_simulation_methods(simulation)
    bed_file, gtf_file = read_reference_files(args.config_yaml)
    regions = get_autosomal_regions(read_regions(bed_file))
    genes = get_region_genes(read_genes(gtf_file), regions)
    gene_level_dir = get_table_dir(args.outdir) / "gene_level_cnp"
    tables = {
        "segment": get_simulation_metrics(simulation, methods),
        "gene": get_simulation_metrics(simulation, methods, genes, gene_level_dir),
    }
    wide = merge_levels(tables)
    write_table(wide, get_table_path(args.outdir, args.simulation_dataset))
    medians = get_table_path(args.outdir, f"{args.simulation_dataset}.median")
    write_table(get_class_medians(wide), medians, float_format="%.6g")


def write_real_data_metrics(args: argparse.Namespace) -> None:
    """Read the sample's runs, compute every table, merge the per-run ones into one wide table, and write them all,
    with the wide table's summary.

    Args:
        args: the command line (`parser`), its reference among its runs (`get_argument_error`).
    """
    keys = args.runs or list(RUNS)
    entry = read_sample_entry(args.benchmark_yaml, args.sample)
    spots = read_spots(Path(entry["ugp"]["spots"]))
    broad = read_broad_cell_types(Path(entry["cell_type"]["labels"]))
    bed_file, gtf_file = read_reference_files(args.config_yaml)
    regions = get_autosomal_regions(read_regions(bed_file))
    genes = get_region_genes(read_genes(gtf_file), regions)

    runs, labels, states = {}, {}, {}
    for key in keys:
        run_dir = Path(entry[key]["dir"])
        runs[key] = get_method_reader(get_run_method(key), integer_profile=True)(run_dir)
        labels[key] = get_run_spots(key, runs[key], run_dir, spots, args.side_um)
        states[key] = get_gene_states(runs[key].profile, genes)

    composition = pd.concat(
        [
            compute_genome_composition(
                states[key].to_numpy(dtype=str), np.ones(len(genes)), float(len(genes)), labels[key]
            ).assign(**get_run_name(key))
            for key in keys
        ],
        ignore_index=True,
    )
    composition = composition[["run", "method", *composition.columns[:-2]]]
    agreement, cells = compare_runs(labels, runs, genes, args.reference)
    cnps = {key: get_tumor_cnps(key, runs[key], Path(entry[key]["dir"]), states[key]) for key in keys}
    segment_level, errors = get_copy_number_tables(
        runs, labels, entry, regions, broad, args.side_um, args.msr, args.assay
    )
    # NB a detail table's float format: None writes full precision.
    details: dict[str, tuple[pd.DataFrame, str | None]] = {
        "gene_level_cnp": (get_gene_level_table(genes, runs), None),
        "genome_composition": (composition, None),
        "cna_confusion": (cells, None),
        "clone_ari": (compute_pairwise_clone_ARI(labels), "%.6g"),
        "cnp_similarity": (compute_pairwise_cnp_similarity(cnps), "%.6g"),
        "segment_level_cnp": (segment_level, None),
    }
    if args.sweep_runs:
        for key in args.sweep_runs:
            run_dir = Path(entry[key]["dir"])
            runs[key] = get_method_reader(get_run_method(key), integer_profile=True)(run_dir)
            labels[key] = get_run_spots(key, runs[key], run_dir, spots, args.side_um)
        panel = {key: runs[key] for key in [args.reference, *args.sweep_runs]}
        sweep_errors = get_copy_number_tables(panel, labels, entry, regions, broad, args.side_um, args.msr, args.assay)[
            1
        ]
        details["cn_error.sweep"] = (sweep_errors, "%.6g")

    # NB a composition's and a copy error's spots are the run's called spots, checked equal as they merge.
    shares = composition[composition["clone"].astype(str) == "all"].drop(columns="clone")
    wide = merge_run_tables(
        [
            get_spot_table(keys, entry, spots, broad, args.side_um),
            shares.rename(columns={"#spots": "#called"}),
            agreement,
            errors.rename(columns={"#spots": "#called"}),
            get_runtime_table(keys, entry, runs),
        ],
        keys,
    )
    wide.insert(0, "sample", args.sample)
    wide["reference"] = args.reference

    for name, (table, float_format) in details.items():
        write_table(table, get_table_path(args.outdir, f"{args.sample}.{name}"), float_format=float_format)
    write_table(wide, get_table_path(args.outdir, args.real_data_dataset))
    write_json(get_summary(wide, args.real_data_dataset), get_json_path(args.outdir, args.real_data_dataset))


def write_figures(args: argparse.Namespace, command: argparse.ArgumentParser) -> None:
    """Read every method's labels and copies once, place the spots, and draw the three PDFs, their text editable
    (`FONT_TYPES`), the metric pages from the run's tables.

    Args:
        args: the command line (`parser`).
        command: its parser, which reports an input the runs contradict.
    """
    plt.rcParams.update(FONT_TYPES)

    broad = read_broad_cell_types(args.celltype)
    epithelial = get_epithelial(broad)
    fallbacks = [] if "none" in args.calicost_fallback_ploidy else list(args.calicost_fallback_ploidy)
    calicost_label, calicost_dir = args.calicost_dir
    calicost = read_calicost(calicost_dir, args.calicost_ploidy, fallbacks)
    numbat = read_numbat(args.numbat_dir)
    clonalscope = None if args.clonalscope_dir is None else read_clonalscope(args.clonalscope_dir, merge_normal=True)
    cnamaste_labels = read_labels(args.cnamaste_labels)
    numbat_title = get_title("numbat", numbat.labels)
    clonalscope_title = NAMES[CLONALSCOPE] if clonalscope is None else get_title(CLONALSCOPE, clonalscope.labels)
    calicost_title = get_title(
        "calicost", calicost.labels, get_calicost_label(calicost_dir, calicost_label, args.calicost_ploidy, fallbacks)
    )
    cnamaste_title = get_title("cnamaste", cnamaste_labels)
    sweep, sweep_runs = {}, {}
    for label, directory in args.calicost_sweep:
        run = read_calicost(directory, args.calicost_ploidy, fallbacks)
        title = get_title("calicost", run.labels, get_calicost_label(directory, label, args.calicost_ploidy, fallbacks))
        sweep[title], sweep_runs[title] = run, Run("calicost", label, directory)
    if len(sweep) != len(args.calicost_sweep):
        command.error(f"--calicost_sweep runs need distinct labels: {[label for label, _ in args.calicost_sweep]}")

    # NB each run read once: its clones per spot of side_um, in the appendix's bar order.
    spot_labels = {
        numbat_title: get_spot_labels(Run("numbat", "", args.numbat_dir), numbat.labels, args.side_um),
        **(
            {}
            if clonalscope is None
            else {
                clonalscope_title: get_spot_labels(
                    Run(CLONALSCOPE, "", args.clonalscope_dir), clonalscope.labels, args.side_um
                )
            }
        ),
        calicost_title: get_spot_labels(Run("calicost", calicost_label, calicost_dir), calicost.labels, args.side_um),
        **{title: get_spot_labels(sweep_runs[title], run.labels, args.side_um) for title, run in sweep.items()},
        cnamaste_title: get_spot_labels(
            Run("cnamaste", "", args.cnamaste_labels.parent), cnamaste_labels, args.side_um
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

    cnamaste_seglevel = pd.read_csv(args.cnamaste_seglevel, sep="\t")
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
    positions = get_spatial_positions(str(args.spaceranger_dir))
    microns = get_visiumhd_microns_per_pixel(args.spaceranger_dir)
    xy = pixel_physical(get_visiumhd_bin_pixels(barcodes, positions), microns)
    physical = pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1]}, index=barcodes)
    logger.info(f"Placed {len(physical):_} bins of {args.side_um} um from {args.spaceranger_dir}.")

    he_frame = slide_image_and_meta([args.spaceranger_dir])
    if he_frame is not None:
        he_xy = pixel_physical(he_frame[["x", "y"]].to_numpy(), microns)
        he_frame = he_frame.assign(x=he_xy[:, 0], y=he_xy[:, 1])

    # NB the diagnosis page: Numbat, the chosen sweep runs, cnamaste, each with its false-positive tumor spots.
    fp_clones = {
        numbat_title: clones[numbat_title],
        **{title: sweep_clones[title] for title in sweep if sweep_runs[title].label in args.false_positive_sweep},
        cnamaste_title: clones[cnamaste_title],
    }
    normal_type = broad.index[get_normal_type(broad).to_numpy()]
    false_positive = {title: get_false_positive_tumor(labels, normal_type) for title, labels in fp_clones.items()}
    for title, found in false_positive.items():
        logger.info(f"{title}: {len(found):_} false-positive tumor spots.")

    genome_axis = get_genome_axis(str(args.region_bed), str(args.chrom_sizes))
    # NB every profile on the cnamaste run's bins, as its own page draws its profile.
    bins = cnamaste_seglevel[["CHR", "START", "END"]]
    bin_profiles = {title: None if p is None else get_bin_profile(p, bins) for title, p in profiles.items()}
    for title, table in bin_profiles.items():
        if table is not None:
            uncalled = int(table.drop(columns=["CHR", "START", "END"]).isna().all(axis=1).sum())
            logger.info(f"{title}: {uncalled:_} of {len(bins):_} cnamaste bins overlap no segment.")
    wide = pd.read_csv(get_table_path(args.outdir, args.real_data_dataset), sep="\t")
    errors = get_copy_number_rows(wide)
    pairs = {
        metric: pd.read_csv(get_table_path(args.outdir, f"{args.sample}.{metric}"), sep="\t")
        for metric in ("cnp_similarity", "clone_ari")
    }
    comparison = get_run_comparison_panels(wide, pairs, errors)
    spots = read_spots(args.spots)
    # NB every page's clones, so a clone is one colour on each page.
    clone_ids = sorted(set().union(*(set(c.astype(int)) for c in drawn)))
    pages = plot_simulation_pages(args.outdir, [args.simulation_dataset])
    tracks = read_clone_tracks(args.cnamaste_labels.parent, cnamaste_labels)
    if tracks is None:
        logger.info(f"No binned_data.npz in {args.cnamaste_labels.parent}: combined.pdf has no cnamaste page.")
    else:
        fit_clones = sorted(np.unique(tracks.fit["new_assignment"]).astype(int).tolist())
        fractions = get_normal_fractions(args.cnamaste_labels.parent, fit_clones, dict(args.cnamaste_normal_fraction))
        missing = sorted(set(fit_clones) - set(fractions))
        if missing:
            command.error(f"--cnamaste_normal_fraction: give the decoder's normal fraction of clone(s) {missing}")
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
    written = write_pages(pages, args.outdir, str(args.combined))
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
    side_bb_dir = dict(args.bb_dir).get(args.side_um)
    if side_bb_dir is None:
        logger.info(f"No --bb_dir at {args.side_um} um: the diagnosis has no UMI page.")
    else:
        umis = get_spot_umis(side_bb_dir, args.assay, args.spots)
        groups = [(spot_labels[title], false_positive[title]) for title in fp_clones]
        umi_page, umi_table = plot_false_positive_umi(groups, umis, FALSE_POSITIVE_COLOUR, args.sample)
        diagnosis.append(umi_page)
        logger.info(f"False-positive spots' UMI against the other called spots':\n{umi_table.to_string(index=False)}")
    diagnosis.append(plot_cn_error(errors, f"{args.sample}, copy-number error", get_run_labels(errors)))
    if sweep:
        sweep_errors = pd.read_csv(get_table_path(args.outdir, f"{args.sample}.cn_error.sweep"), sep="\t")
        title = f"{args.sample}, copy-number error, CalicoST sweep"
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
    diagnosis += plot_simulation_diagnosis(args.outdir, [args.simulation_dataset])
    diagnosis.append(plot_panel_titles(diagnosis))
    written = write_pages(diagnosis, args.outdir, str(args.diagnosis))
    logger.info(f"Wrote {', '.join(str(path) for path in written)}.")
    if umi_table is not None:
        write_table(umi_table, get_table_path(args.outdir, args.diagnosis.stem), float_format="%.6g")

    appendix = [plot_umi_page(args.sample, args.bb_dir, args.assay, args.bins)] if args.bb_dir else []
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
        *plot_metric_pages(args.outdir, args.sample, args.real_data_dataset, [args.simulation_dataset]),
    ]
    metric_pages.append(plot_panel_titles([*appendix, *metric_pages]))
    write_appendix(appendix, metric_pages, args.outdir, str(args.appendix))


def get_argument_error(args: argparse.Namespace) -> str:
    """Get why the command line's inputs disagree, or "", before any step runs: the cnamaste run's inputs
    (`cnamaste_page.get_cnamaste_input_mismatch`), the reference among the runs, a label per `--calicost_sweep` run,
    and each `--false_positive_sweep` a sweep run's label.

    Args:
        args: the command line (`parser`).

    Returns:
        The error, or "" where there is none.
    """
    mismatch = get_cnamaste_input_mismatch(args.cnamaste_labels, args.cnamaste_seglevel)
    if mismatch:
        return mismatch
    keys = args.runs or list(RUNS)
    if args.reference not in keys:
        return f"--reference {args.reference} is not among the runs {keys}"
    unlabelled = [str(directory) for label, directory in args.calicost_sweep if not label]
    if unlabelled:
        return f"--calicost_sweep {unlabelled}: give each a label, LABEL=DIR"
    unknown = sorted(set(args.false_positive_sweep) - {label for label, _ in args.calicost_sweep})
    if unknown:
        return f"--false_positive_sweep {unknown}: not a --calicost_sweep label"
    return ""


def parser() -> argparse.ArgumentParser:
    """Get the command line's parser: the yamls and the simulations' methods, the sample's runs, the figures' inputs, the
    drawn spots' side, CalicoST's table, the genome axis' files and the outputs."""
    command = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    command.add_argument("--sample", required=True)
    command.add_argument("--benchmark_yaml", type=Path, default=BENCHMARK / "configs" / "benchmark.local.yaml")
    command.add_argument("--config_yaml", type=Path, default=BENCHMARK / "configs" / "config.local.yaml")
    command.add_argument("--method", action="append", dest="methods")
    command.add_argument("--simulation_dataset", default="calicost-sim")
    command.add_argument("--run", action="append", dest="runs")
    command.add_argument("--reference", default="cnamaste")
    command.add_argument("--sweep_run", action="append", dest="sweep_runs", default=[], metavar="KEY")
    command.add_argument("--msr", type=int, default=300)
    command.add_argument("--real_data_dataset", default="ovarian")
    command.add_argument("--table", default="table")
    command.add_argument("--metric_text", default="metric")
    command.add_argument("--spaceranger_dir", type=Path, required=True)
    command.add_argument("--celltype", type=Path, required=True)
    command.add_argument("--spots", type=Path, required=True)
    command.add_argument("--bb_dir", type=get_bb_dir, action="append", default=[])
    command.add_argument("--cnamaste_labels", type=Path, required=True)
    command.add_argument("--cnamaste_seglevel", type=Path, required=True)
    command.add_argument(
        "--cnamaste_normal_fraction", type=get_normal_fraction, action="append", default=[], metavar="CLONE=F"
    )
    command.add_argument("--calicost_dir", type=get_run_argument, required=True, metavar="[LABEL=]DIR")
    command.add_argument("--calicost_sweep", type=get_run_argument, action="append", default=[], metavar="LABEL=DIR")
    command.add_argument("--false_positive_sweep", action="append", default=[], metavar="LABEL")
    command.add_argument("--numbat_dir", type=Path, required=True)
    command.add_argument("--clonalscope_dir", type=Path)
    command.add_argument("--side_um", type=int, default=16)
    command.add_argument("--assay", default="VISIUM3prime")
    command.add_argument("--calicost_ploidy", choices=CALICOST_PLOIDIES, default="unconstrained")
    command.add_argument(
        "--calicost_fallback_ploidy",
        nargs="+",
        choices=(*CALICOST_PLOIDIES, "none"),
        default=list(get_calicost_fallbacks()),
    )
    command.add_argument("--region_bed", type=Path, default=DATA / "hg38.regions.bed")
    command.add_argument("--chrom_sizes", type=Path, default=DATA / "hg38.chrom.sizes")
    command.add_argument("--bins", type=int, default=100)
    command.add_argument("--outdir", type=Path, required=True)
    command.add_argument("--combined", type=Path, default=Path("combined.pdf"))
    command.add_argument("--diagnosis", type=Path, default=Path("diagnosis.pdf"))
    command.add_argument("--appendix", type=Path, default=Path("appendix.pdf"))
    return command


def get_cnamaste_root() -> Path:
    """Get the cnamaste repository `cnamaste` is imported from: `$CNAMASTE_ROOT`, its package `python/cnamaste/`.

    Returns:
        The repository's root.

    Raises:
        RuntimeError: where `CNAMASTE_ROOT` is unset, or `cnamaste` is imported from another directory.
    """
    root = os.environ.get(CNAMASTE_ROOT)
    if not root:
        raise RuntimeError(f"{CNAMASTE_ROOT} is unset: set it to a cnamaste repository's root (README.md, Run).")
    package = Path(root).resolve() / "python" / "cnamaste"
    imported = Path(cnamaste.__file__).resolve().parent
    if imported != package:
        raise RuntimeError(f"cnamaste is imported from {imported}, not {package}: put {package.parent} on PYTHONPATH.")
    return package.parents[1]


def main(argv: list[str] | None = None) -> int:
    """Run the benchmark: score the simulations' runs and the sample's runs, write the real-data LaTeX tables, then
    draw the three PDFs from the tables.

    Args:
        argv: the command line, or None for `sys.argv`.

    Returns:
        0.
    """
    command = parser()
    args = command.parse_args(argv)
    error = get_argument_error(args)
    if error:
        command.error(error)
    logger.info(f"cnamaste imported from {get_cnamaste_root()}.")
    logger.info(f"Step 1 of 4: the simulations' metrics, tsvs/{args.simulation_dataset}.tsv.gz.")
    write_simulation_metrics(args)
    logger.info(f"Step 2 of 4: {args.sample}'s runs compared, tsvs/{args.real_data_dataset}.tsv.gz.")
    write_real_data_metrics(args)
    logger.info(
        f"Step 3 of 4: the LaTeX tables and metric text, text/{args.table}.tex and text/{args.metric_text}.tex."
    )
    write_latex_table(args.outdir, args.real_data_dataset, args.table)
    write_metric_text(args.outdir, args.metric_text)
    logger.info("Step 4 of 4: the figures, plots/ and pdf-pages/.")
    write_figures(args, command)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
