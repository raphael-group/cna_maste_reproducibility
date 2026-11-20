"""The analysis' options: every input, parameter and output name of its four steps, read from one yaml.

Runpeng Luo (2026-11-09), from `scripts/run_benchmark.py`'s command line (2026-11-04)

Dependencies: pyyaml; the analysis environment (`pyproject.toml`). No cnamaste: what the analysis takes from it is
ported in `scripts/cnamaste_port/`.

`run_cna_maste_reproducibility.smk` writes `<OUTDIR>/configs/analysis_options.yaml`, one key per field of
`AnalysisOptions`, from the benchmark yaml's `vishd.<sample>` and the config's `ANALYSIS`; each step reads it
(`read_options`); an option it leaves out keeps its default.

Inputs:
    benchmark_yaml: the runs, the pipeline's `<OUTDIR>/configs/benchmark.yaml`.
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
        (`metrics.gene_level.read_reference_files`), the pipeline's `<OUTDIR>/configs/config.yaml`.
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
    methods: the simulation methods scored, every method of the `simulation` entry by default.
    simulation_dataset: the simulations' wide table's name, `calicost-sim` by default.
    runs: the sample's run keys, `numbat`, `clonalscope`, `calicost` and `cnamaste` by default (`RUNS`); a run's
        method the key before its first `_` or `-`, e.g. `cnamaste` of `cnamaste-ac88605f`, read by its method's
        reader, Clonalscope's clusters as called, cnamaste's copies those of the clones its integer clones name
        (`metrics.read_method_input.get_labelled_profile`).
    reference: the run every other run is compared with, `cnamaste` by default.
    sweep_runs: run keys of the sample, e.g. a CalicoST parameter sweep's runs: the sweep runs and the reference run
        are scored on the copy-number error as one panel, on the bins they share; none by default.
    msr: UGP's `min_snp_reads` of the bins scored, their directory `MSR<msr>`, 300 by default.
    real_data_dataset: the sample's wide table's and its summary's name, `ovarian` by default.
    table: the LaTeX tables' name under `text/`, `table` by default.
    metric_text: the metric text's name under `text/`, `metric` by default.
    cnamaste_normal_fraction: a tumour clone of the run's fit and its decoder's normal fraction, `<clone>=<fraction>`,
        once per tumour clone; by default each integer clone's from the run's log, run_cnamaste's `normal fractions
        per clone` (`plots.cnamaste_page.get_normal_fractions`); the normal clone is held at 1, as the decoder holds
        it (`plots.cnamaste_page.NORMAL_CLONE_FRACTION`).
    false_positive_sweep: the labels of `calicost_sweep` runs on the diagnosis page; every run of that label, at
        each side.
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
    outdir: the analysis' directory, the pipeline's `<OUTDIR>/benchmark_results/`: tables in its `tsvs/`, LaTeX tables
        in its `text/`, figures in its `plots/`, each figure's pages in its `pdf-pages/`, each existing beforehand
        (`metrics.outputs`).
    combined, diagnosis, appendix: the three figures' names under `plots/`.
Notes/References:
    A list option is a yaml list; `bb_dir`, `cnamaste_normal_fraction`, `calicost_dir` and `calicost_sweep` take the
    command line's text forms, e.g. `bb_dir: ["16=<dir>"]`.
"""

from __future__ import annotations

from dataclasses import MISSING, dataclass, field, fields
from pathlib import Path

import yaml
from scripts.cnamaste_port.calicost.read import CALICOST_PLOIDIES
from scripts.metrics.read_method_input import get_benchmark_dir, get_calicost_fallbacks, get_run_argument
from scripts.plots.cnamaste_page import get_cnamaste_input_mismatch

BENCHMARK = get_benchmark_dir()
"""The benchmark's directory: `scripts/`, `configs/` and `data/`."""

DATA = BENCHMARK / "data"
"""The benchmark's `data/`: the genome axis' default files."""

RUNS = ("numbat", "clonalscope", "calicost", "cnamaste")
"""The sample's runs scored where no `runs` is given."""


@dataclass
class AnalysisOptions:
    """The analysis' options, as the module docstring defines them."""

    sample: str
    spaceranger_dir: Path
    celltype: Path
    spots: Path
    cnamaste_labels: Path
    cnamaste_seglevel: Path
    calicost_dir: tuple[str, Path]
    numbat_dir: Path
    outdir: Path
    benchmark_yaml: Path
    config_yaml: Path
    methods: list[str] | None = None
    simulation_dataset: str = "calicost-sim"
    runs: list[str] | None = None
    reference: str = "cnamaste"
    sweep_runs: list[str] = field(default_factory=list)
    msr: int = 300
    real_data_dataset: str = "ovarian"
    table: str = "table"
    metric_text: str = "metric"
    bb_dir: list[tuple[int, Path]] = field(default_factory=list)
    cnamaste_normal_fraction: list[tuple[int, float]] = field(default_factory=list)
    calicost_sweep: list[tuple[str, Path]] = field(default_factory=list)
    false_positive_sweep: list[str] = field(default_factory=list)
    clonalscope_dir: Path | None = None
    side_um: int = 16
    assay: str = "VISIUM3prime"
    calicost_ploidy: str = "unconstrained"
    calicost_fallback_ploidy: list[str] = field(default_factory=lambda: list(get_calicost_fallbacks()))
    region_bed: Path = DATA / "hg38.regions.bed"
    chrom_sizes: Path = DATA / "hg38.chrom.sizes"
    bins: int = 100
    combined: Path = Path("combined.pdf")
    diagnosis: Path = Path("diagnosis.pdf")
    appendix: Path = Path("appendix.pdf")


PATHS = (
    "spaceranger_dir",
    "celltype",
    "spots",
    "cnamaste_labels",
    "cnamaste_seglevel",
    "numbat_dir",
    "outdir",
    "benchmark_yaml",
    "config_yaml",
    "clonalscope_dir",
    "region_bed",
    "chrom_sizes",
    "combined",
    "diagnosis",
    "appendix",
)
"""The options read as paths."""


def get_normal_fraction(text: str) -> tuple[int, float]:
    """Read a `cnamaste_normal_fraction` value, `<clone>=<fraction>`.

    Args:
        text: the value, e.g. `1=0.38`.

    Returns:
        The clone and its normal fraction.

    Raises:
        ValueError: where `text` is not `<int>=<fraction in [0, 1]>`.
    """
    clone, sep, fraction = text.partition("=")
    try:
        value = float(fraction)
    except ValueError:
        value = -1.0
    if not sep or not clone.isdigit() or not 0.0 <= value <= 1.0:
        raise ValueError(f"expected <clone>=<fraction in [0, 1]>, got {text!r}")
    return int(clone), value


def get_bb_dir(text: str) -> tuple[int, Path]:
    """Read a `bb_dir` value, `<side_um>=<dir>`.

    Args:
        text: the value.

    Returns:
        The bins' side in um and the directory.

    Raises:
        ValueError: where `text` is not `<int>=<dir>`.
    """
    side, sep, directory = text.partition("=")
    if not sep or not side.isdigit() or not directory:
        raise ValueError(f"expected <side_um>=<dir>, got {text!r}")
    return int(side), Path(directory)


def get_option_value(name: str, value):
    """Read one option from its yaml value.

    Args:
        name: a field of `AnalysisOptions`.
        value: its yaml value.

    Returns:
        The value as `AnalysisOptions` holds it.

    Raises:
        ValueError: where a CalicoST table is not one of `CALICOST_PLOIDIES`.
    """
    if value is None:
        return None
    if name in PATHS:
        return Path(value)
    if name in ("msr", "side_um", "bins"):
        return int(value)
    if name == "bb_dir":
        return [get_bb_dir(str(v)) for v in value]
    if name == "cnamaste_normal_fraction":
        return [get_normal_fraction(str(v)) for v in value]
    if name == "calicost_dir":
        return get_run_argument(str(value))
    if name == "calicost_sweep":
        return [get_run_argument(str(v)) for v in value]
    if name == "calicost_ploidy" and value not in CALICOST_PLOIDIES:
        raise ValueError(f"calicost_ploidy {value!r} is not one of {CALICOST_PLOIDIES}")
    if name == "calicost_fallback_ploidy" and set(value) - {*CALICOST_PLOIDIES, "none"}:
        raise ValueError(f"calicost_fallback_ploidy {value!r}: each one of {(*CALICOST_PLOIDIES, 'none')}")
    return value


def get_argument_error(options: AnalysisOptions) -> str:
    """Get why the options' inputs disagree, or "", before any step runs: the cnamaste run's inputs
    (`cnamaste_page.get_cnamaste_input_mismatch`), the reference among the runs, a label per `calicost_sweep` run,
    and each `false_positive_sweep` a sweep run's label.

    Args:
        options: the options.

    Returns:
        The error, or "" where there is none.
    """
    mismatch = get_cnamaste_input_mismatch(options.cnamaste_labels, options.cnamaste_seglevel)
    if mismatch:
        return mismatch
    keys = options.runs or list(RUNS)
    if options.reference not in keys:
        return f"reference {options.reference} is not among the runs {keys}"
    unlabelled = [str(directory) for label, directory in options.calicost_sweep if not label]
    if unlabelled:
        return f"calicost_sweep {unlabelled}: give each a label, LABEL=DIR"
    unknown = sorted(set(options.false_positive_sweep) - {label for label, _ in options.calicost_sweep})
    if unknown:
        return f"false_positive_sweep {unknown}: not a calicost_sweep label"
    return ""


def read_options(path: str | Path) -> AnalysisOptions:
    """Read the analysis' options from a yaml, one key per field of `AnalysisOptions`.

    Args:
        path: the yaml.

    Returns:
        The options.

    Raises:
        ValueError: where a key is not an option, a required option is missing, or the options disagree
            (`get_argument_error`).
    """
    with open(path) as handle:
        values = yaml.safe_load(handle) or {}
    names = {f.name for f in fields(AnalysisOptions)}
    unknown = sorted(set(values) - names)
    if unknown:
        raise ValueError(f"{path}: not analysis options: {unknown}")
    required = [
        f.name
        for f in fields(AnalysisOptions)
        if f.name not in values and f.default is MISSING and f.default_factory is MISSING
    ]
    if required:
        raise ValueError(f"{path}: options missing: {required}")
    options = AnalysisOptions(**{name: get_option_value(name, value) for name, value in values.items()})
    error = get_argument_error(options)
    if error:
        raise ValueError(f"{path}: {error}")
    return options
