# The paper's benchmark end to end, its one entrypoint: setup, prepare, run every method, then the analysis.
#
# Runpeng Luo (2026-11-08)
#
# Dependencies:
#   CNA_MASTE_REPRODUCIBILITY_ENV (configs/envs/cna_maste_reproducibility.yaml): snakemake, pyyaml, pandas, pysam, uv,
#   git; conda at CONDA_ROOT; pdflatex for the analysis' LaTeX tables; Linux for the method runs (GNU time, systemd
#   scopes for CalicoST's memory cap).
#
# Usage:
#   snakemake -s run_cna_maste_reproducibility.smk --cores <CORES> --resources mem_mb=<MEM_GB * 1024> [target ...] \
#       [--configfile <yaml>] [--config KEY=VALUE ...]
#           target  # all (default) | setup | prepare | run_methods | analysis, or any rule's output
#
# Inputs
#   configs/cna_maste_reproducibility.yaml: every path, resource, repository, environment and method config; a
#       --configfile given after it overrides its keys.
# Outputs:
#   <OUTDIR>/setup/<KEY>.ok, inputs.ok: a checked repository or environment, the checked inputs.
#   <OUTDIR>/configs/: the configs written for the steps (baselines/, ugp/, cnamaste.<case>.yaml, benchmark.yaml,
#       analysis_options.yaml).
#   <OUTDIR>/zenodo_sim/{inputs,numbat,clonalscope,calicost,cnamaste}/, simulation.runtime.tsv: the simulations' runs.
#   <OUTDIR>/ovarian/{inputs,ugp,celltyping,numbat,clonalscope,calicost,cnamaste}/, runtime.<sample>.tsv: the sample's.
#   <OUTDIR>/benchmark_results/: the analysis' tsvs/, text/, plots/ and pdf-pages/ (docs/reference.md, Outputs).
#   <OUTDIR>/logs/<step>/<rule>[.<wildcards>].log and .tsv: each job's log and its time and memory (benchmark:);
#       a method's per-sample run logs beside its outputs, <OUTDIR>/<case>/<method>/{logs,bench}/; CalicoST's run logs
#       in <OUTDIR>/logs/baselines/.
# Notes/References:
#   Steps: setup (check the inputs, clone a missing repository, build a missing environment), prepare (each method's
#   inputs per simulation; Universal-Genotyping-Pipeline on the sample, then each method's inputs from it; cell types),
#   run_methods (Numbat, Clonalscope and CalicoST per case and sample, scripts/baselines/<method>/; cnamaste per
#   simulation, solver and seed, as cnamaste's batch_zenodo_run.sh runs it, and on the sample), analysis (the four
#   steps of scripts/analysis/, each a rule, on this pipeline's runs). Every step is a rule of this file: no nested
#   workflow.

import contextlib
import json
import os
import sys

import yaml


configfile: os.path.join(workflow.basedir, "configs", "cna_maste_reproducibility.yaml")


REPO = os.path.abspath(workflow.basedir)
# NB: scripts/ as a package, for the rules that run its functions in this environment
sys.path.insert(0, REPO)


def under_repo(path):
    """A config path: absolute as given, relative under this repository."""
    return path if os.path.isabs(path) else os.path.join(REPO, path)


def read_yaml(path):
    with open(under_repo(path)) as f:
        return yaml.safe_load(f)


def write_yaml(path, data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False)


C = config
OUT = under_repo(C["OUTDIR"])
TMP = under_repo(C["TMPDIR"]) if C.get("TMPDIR") else f"{OUT}/tmp"
LOGS = f"{OUT}/logs"
SETUP = f"{OUT}/setup"
CFG = f"{OUT}/configs"
CORES = int(C["CORES"])
MEM_MB = int(C["MEM_GB"]) * 1024
CONDA = C["CONDA_ROOT"]
MAMBA = f"{CONDA}/bin/mamba" if os.path.exists(f"{CONDA}/bin/mamba") else f"{CONDA}/bin/conda"
PY = sys.executable
BIN = os.path.dirname(PY)
BASELINES = f"{REPO}/scripts/baselines"

REPOS = {key: under_repo(C[key]) for key in C["REPOS"]}
ENVS = {key: under_repo(C[key]) for key in C["ENV_RECIPES"]}
CALICOST = read_yaml(C["CALICOST_CONFIG"])
NUMBAT = read_yaml(C["NUMBAT_CONFIG"])
CLONALSCOPE = read_yaml(C["CLONALSCOPE_CONFIG"])
CNAMASTE = C["CNAMASTE"]
ANALYSIS = C["ANALYSIS"]

# NB: simulations: every numcnas* folder of the deposit's simulated_data_related/
ZENODO = C["CALICOST_ZENODO_SIM_ROOT"]
SIM_DIR = f"{ZENODO}/simulated_data_related"
SIMS = sorted(
    d
    for d in (os.listdir(SIM_DIR) if os.path.isdir(SIM_DIR) else [])
    if d.startswith("numcnas") and os.path.isdir(os.path.join(SIM_DIR, d))
)
ZS = f"{OUT}/zenodo_sim"
ZS_INPUTS = f"{ZS}/inputs"
CN = f"{ZS}/cnamaste"
CN_RUNS = f"{CN}/nomixing_cnamaste_related"
SOLVERS = CNAMASTE["zenodo_sim"]["solvers"]
SEEDS = [str(s) for s in CNAMASTE["zenodo_sim"]["random_seeds"]]
QA_SOLVER = CNAMASTE["zenodo_sim"]["qa_solver"]

# NB: the ovarian sample: Space Ranger's 16um bins, and its 16um BAM, Space Ranger's or retagged from the 2um BAM
OUTS = C["OVARIAN_SPACERANGER_ROOT"]
SAMPLE = C["OVARIAN_SAMPLE_ID"]
B16 = f"{OUTS}/binned_outputs/square_016um"
OV = f"{OUT}/ovarian"
OV_INPUTS = f"{OV}/inputs"
BAM16_SPACERANGER = f"{B16}/possorted_genome_bam.bam"
BAM16 = (
    BAM16_SPACERANGER
    if os.path.exists(BAM16_SPACERANGER)
    else f"{OV_INPUTS}/{SAMPLE}/square_016um/possorted_genome_bam.bam"
)
UGP_CONFIG_TEMPLATE = read_yaml(C["UGP_CONFIG_TEMPLATE"])
with open(under_repo(C["UGP_SAMPLE_TEMPLATE"])) as f:
    UGP_SAMPLE_TEMPLATE = json.load(f)
DATASET_ID = UGP_SAMPLE_TEMPLATE["samples"][0]["dataset_id"]
MSR = UGP_CONFIG_TEMPLATE["params_combine_counts"]["min_snp_reads"][0]
UGP_DIR = f"{OV}/ugp"
UGP_PILEUP = f"{UGP_DIR}/pileup/VISIUM3prime_{DATASET_ID}"
UGP_OUTPUTS = {
    "spots": f"{UGP_DIR}/bb/MSR{MSR}/VISIUM3prime/barcodes.tsv.gz",
    "h5ad": f"{UGP_DIR}/bb/VISIUM3prime.h5ad",
    "pileup": f"{UGP_PILEUP}/cellSNP.tag.AD.mtx",
    "phase": f"{UGP_DIR}/phase/chr1.vcf.gz",
}
CELLTYPES = f"{OV}/celltyping/{SAMPLE}/square_016um"
# NB: UGP's 16um bins in CalicoST's input layout, CalicoST's and cnamaste's input
BB16 = f"{OV_INPUTS}/{SAMPLE}/bb-016um-formatted"

# NB: Clonalscope's references, shared by both cases (clonalscope/build_refs_clonalscope.sh at its pinned commit)
CLONALSCOPE_REF = f"{OUT}/references/clonalscope"
CLONALSCOPE_DATA_RAW = f"{CLONALSCOPE_REF}/clonalscope_data-raw_07dd0c5"
CLONALSCOPE_REFS = [
    f"{CLONALSCOPE_REF}/hg38.200kb.bins.bed",
    f"{CLONALSCOPE_REF}/seg_table.arm.tsv",
    *[
        f"{CLONALSCOPE_DATA_RAW}/{name}"
        for name in ("sizes.cellranger-GRCh38-1.0.0.txt", "cyclegenes.rds", "hg38.genes.bed", "cytoarm_table_hg38.txt")
    ],
]
# NB: CalicoST's GRCh38_resources at its commit, its configuration_cna keys; cnamaste's references too
CALICOST_REFERENCES = {
    key: f"{REPOS['CALICOST_REPO']}/GRCh38_resources/{name}"
    for key, name in {
        "geneticmap_file": "genetic_map_GRCh38_merged.tab.gz",
        "hgtable_file": "hgTables_hg38_gencode.txt",
        "filtergenelist_file": "ig_gene_list.txt",
        "filterregion_file": "HLA_regions.bed",
    }.items()
}

# NB: <OUTDIR>/<case> -> the method configs' case
CASES = {"zenodo_sim": "calicost-sim", "ovarian": "ovarian"}
METHODS = ["numbat", "clonalscope", "calicost"]

# NB: each method's samples per case: the simulations, or the ovarian sample at the bin sides it runs on
NUMBAT_SAMPLES = {"zenodo_sim": SIMS, "ovarian": [SAMPLE]}
NUMBAT_PILEUP = f"{OV}/numbat/{SAMPLE}/prep/pileup_{NUMBAT['ovarian']['bin_um']:03d}um"
CLONALSCOPE_AGGREGATE = {f"{SAMPLE}_{b:03d}um": b for b in CLONALSCOPE["ovarian"]["bins_um"] if b != 16}
CLONALSCOPE_SAMPLES = {
    "zenodo_sim": SIMS,
    "ovarian": [SAMPLE if b == 16 else f"{SAMPLE}_{b:03d}um" for b in CLONALSCOPE["ovarian"]["bins_um"]],
}
CALICOST_SIMS = SIMS if CALICOST["calicost-sim"]["sims"] == "all" else list(CALICOST["calicost-sim"]["sims"])
CALICOST_INPUTS = {name.format(sample=SAMPLE): spec for name, spec in CALICOST["ovarian"]["inputs"].items()}
CALICOST_RUNS = {
    tag.format(sample=SAMPLE): {
        key: str(value).format(sample=SAMPLE) if key == "input" else value for key, value in run.items()
    }
    for tag, run in CALICOST["ovarian"]["runs"].items()
}
if any(spec["bin_um"] != 16 for spec in CALICOST_INPUTS.values()):
    raise ValueError(f"{C['CALICOST_CONFIG']}: ovarian inputs are UGP's 16um bins alone, bin_um 16")
CALICOST_INPUT_SOURCES = {name: BB16 for name in CALICOST_INPUTS}
# NB: the ovarian CalicoST run the analysis scores, named in the figures by its tag's last word, e.g. msu20
CALICOST_ANALYSIS_RUN = ANALYSIS["calicost_run"].format(sample=SAMPLE)
if CALICOST_ANALYSIS_RUN not in CALICOST_RUNS:
    raise ValueError(f"ANALYSIS calicost_run {CALICOST_ANALYSIS_RUN}: not a run of {C['CALICOST_CONFIG']}")
CALICOST_LABEL = CALICOST_ANALYSIS_RUN.rsplit("_", 1)[-1]
if ANALYSIS["calicost_sim"] not in ("deposit", "rerun"):
    raise ValueError(f"ANALYSIS calicost_sim {ANALYSIS['calicost_sim']}: deposit or rerun")


def in_env(env):
    """A shell prefix activating the conda environment `env` (a key of ENV_RECIPES)."""
    return f"source {CONDA}/etc/profile.d/conda.sh && set +u && conda activate {ENVS[env]} && set -u && "


IN_CLONALSCOPE = in_env("CLONALSCOPE_ENV")
# NB: what CalicoST's scripts read: their configs, conda, a python with pyyaml and numpy >= 2, and the logs' directory
CALICOST_SHELL = f"CONFIGS={CFG}/baselines CONDA={CONDA} PY={PY} LOGS={LOGS}/baselines"


def clonalscope_matrix(case, sample):
    """The 10x matrix directory Clonalscope reads for `sample` of `case`."""
    if case == "zenodo_sim":
        return f"{ZS_INPUTS}/{sample}/filtered_feature_bc_matrix"
    if sample in CLONALSCOPE_AGGREGATE:
        return f"{OV}/clonalscope/{sample}/input/filtered_feature_bc_matrix"
    return f"{B16}/filtered_feature_bc_matrix"


def baseline_outputs(case, method):
    """A method's last outputs on `case`, one per sample or run."""
    if method == "numbat":
        label = NUMBAT[CASES[case]]["run_label"]
        return [f"{OUT}/{case}/numbat/postprocess/{label}/{s}/{s}.RNA.cells.tsv" for s in NUMBAT_SAMPLES[case]]
    if method == "clonalscope":
        label = CLONALSCOPE[CASES[case]]["run_label"]
        return [f"{OUT}/{case}/clonalscope/{s}/{label}/decode_states.done" for s in CLONALSCOPE_SAMPLES[case]]
    if case == "zenodo_sim":
        return [f"{ZS}/calicost/{sim}/repro.done" for sim in CALICOST_SIMS]
    return [f"{OV}/calicost/{tag}/run.done" for tag in CALICOST_RUNS]


def prepared_outputs():
    """Every method's inputs, both cases."""
    return [
        *CLONALSCOPE_REFS,
        *[
            path
            for sim in SIMS
            for path in (
                f"{ZS_INPUTS}/{sim}/filtered_feature_bc_matrix/matrix.mtx.gz",
                f"{ZS}/numbat/{sim}/prep/{sim}.count_mat.rds",
                f"{ZS}/numbat/{sim}/prep/{sim}_allele_counts.tsv.gz",
                f"{ZS_INPUTS}/{sim}/cnamaste_sample_sheet.tsv",
            )
        ],
        f"{OV}/numbat/{SAMPLE}/prep/{SAMPLE}.count_mat.rds",
        f"{OV}/numbat/{SAMPLE}/prep/{SAMPLE}_allele_counts.tsv.gz",
        f"{OV_INPUTS}/{SAMPLE}/cnamaste_sample_sheet.tsv",
        CELLTYPES,
    ]


def uv_run(env, repo):
    """`uv run` in `env`, the uv environment of `repo`, without syncing it."""
    return f"UV_PROJECT_ENVIRONMENT={ENVS[env]} {BIN}/uv run --no-sync --project {REPOS[repo]}"


UV_CNAMASTE = uv_run("CNAMASTE_ENV", "CNAMASTE_REPO")
UV_ANALYSIS = uv_run("CNAMASTE_ANALYSIS_ENV", "CNAMASTE_ANALYSIS_REPO")

shell.prefix(f"set -euo pipefail; mkdir -p {TMP}; export TMPDIR={TMP}; ")


wildcard_constraints:
    repo="|".join(REPOS),
    env="|".join(ENVS),
    case="|".join(CASES),
    sim=r"[^/]+",
    s=r"[^/]+",
    label=r"[^/]+",
    solver="|".join(SOLVERS),
    seed=r"\d+",
    name="|".join(CALICOST_INPUTS) or "-",
    tag="|".join(CALICOST_RUNS) or "-",


##################################################
# Targets

ANALYSIS_ROOT = REPOS["CNAMASTE_ANALYSIS_REPO"]
RESULTS = f"{OUT}/benchmark_results"
# NB: AnalysisOptions' default names (scripts/analysis/options.py)
SIM_DATASET = "calicost-sim"
REAL_DATASET = "ovarian"
TEXTS = ["table", "metric"]
FIGURES = ["combined.pdf", "diagnosis.pdf", "appendix.pdf"]
ANALYSIS_OUTPUTS = [
    *[f"{RESULTS}/text/{name}.pdf" for name in TEXTS],
    *[f"{RESULTS}/plots/{name}" for name in FIGURES],
]


rule all:
    input:
        ANALYSIS_OUTPUTS,


rule setup:
    input:
        f"{SETUP}/inputs.ok",
        expand(f"{SETUP}/{{key}}.ok", key=[*REPOS, *ENVS]),


rule prepare:
    input:
        prepared_outputs(),


rule run_methods:
    input:
        [path for case in CASES for method in METHODS for path in baseline_outputs(case, method)],
        expand(f"{CN}/qa/qa.seed{{seed}}.done", seed=SEEDS),
        f"{OV}/cnamaste/binned_data.done",
        f"{ZS}/simulation.runtime.tsv",
        f"{OV}/runtime.{SAMPLE}.tsv",


rule analysis:
    input:
        ANALYSIS_OUTPUTS,


##################################################
# Setup


def required_inputs():
    """(config key, path) of every input the steps read."""
    found = [
        ("CALICOST_ZENODO_SIM_ROOT", f"{ZENODO}/{d}")
        for d in ("simulated_data_related", "nomixing_calicost_related", "numbat_related")
    ]
    found += [
        ("OVARIAN_SPACERANGER_ROOT", p)
        for p in (
            f"{B16}/filtered_feature_bc_matrix.h5",
            f"{B16}/filtered_feature_bc_matrix/matrix.mtx.gz",
            f"{B16}/spatial/scalefactors_json.json",
            f"{B16}/spatial/tissue_positions.parquet",
            f"{B16}/spatial/tissue_hires_image.png",
            f"{B16}/analysis/clustering/gene_expression_graphclust/clusters.csv",
            f"{OUTS}/spatial",
            f"{OUTS}/barcode_mappings.parquet",
            BAM16 if BAM16 == BAM16_SPACERANGER else f"{OUTS}/possorted_genome_bam.bam",
        )
    ]
    found += [("PATH_TO_REFDATA", f"{C['PATH_TO_REFDATA']}/{p}") for p in ("fasta/genome.fa", "genes/genes.gtf.gz")]
    found += [
        (key, C[key])
        for key in (
            "PATH_TO_SNP_PANEL",
            "PATH_TO_PHASE_PANEL",
            "PATH_TO_EAGLE_MAP",
            "PATH_TO_ORDERED_GENES",
            "FORMAT_BB_SCRIPT",
        )
    ]
    found += [("CONDA_ROOT", f"{CONDA}/etc/profile.d/conda.sh")]
    return found


rule check_inputs:
    output:
        f"{SETUP}/inputs.ok",
    log:
        f"{LOGS}/setup/check_inputs.log",
    run:
        missing = [f"{key}: {path}" for key, path in required_inputs() if not os.path.exists(path)]
        with open(log[0], "w") as f:
            f.write("".join(f"missing {m}\n" for m in missing) or "every input found\n")
        if missing:
            listed = "\n  ".join(missing)
            raise ValueError(f"inputs missing, set them in configs/cna_maste_reproducibility.yaml:\n  {listed}")
        shell("touch {output}")


rule clone_repo:
    output:
        f"{SETUP}/{{repo}}.ok",
    log:
        f"{LOGS}/setup/clone_repo.{{repo}}.log",
    benchmark:
        f"{LOGS}/setup/clone_repo.{{repo}}.tsv"
    params:
        path=lambda wc: REPOS[wc.repo],
        url=lambda wc: C["REPOS"][wc.repo]["url"],
        rev=lambda wc: C["REPOS"][wc.repo]["rev"],
    shell:
        """
        exec >{log} 2>&1
        if [ ! -d {params.path}/.git ]; then
            git clone {params.url} {params.path}
            git -C {params.path} checkout {params.rev}
        fi
        want=$(git -C {params.path} rev-parse --verify -q {params.rev}^0 || true)
        [ -n "$want" ] || want=$(git -C {params.path} rev-parse --verify -q origin/{params.rev}^0)
        echo "{wildcards.repo}: {params.path} at $(git -C {params.path} describe --always --dirty)"
        echo "{params.rev} is $want"
        if [ "$(git -C {params.path} rev-parse HEAD)" != "$want" ]; then
            echo "{params.path} is not at {params.rev}: check it out, or set {wildcards.repo} to another clone" >&2
            exit 1
        fi
        touch {output}
        """


def env_repo(env):
    """The repository an environment is built from, or None."""
    recipe = C["ENV_RECIPES"][env]
    if recipe.startswith("uv:"):
        return recipe.split(":")[1]
    return {"CALICOST_ENV": "CALICOST_REPO", "CLONALSCOPE_ENV": "CLONALSCOPE_REPO"}.get(
        env, next((key for key in REPOS if "{" + key + "}" in recipe), None)
    )


def env_command(env):
    """The shell command that builds `env` where missing, then installs a repository into it where needed."""
    path, recipe = ENVS[env], C["ENV_RECIPES"][env]
    if recipe.startswith("uv:"):
        _, repo, *group = recipe.split(":")
        groups = " ".join(f"--group {g}" for g in group)
        return (
            f"if [ -x {path}/bin/python ]; then echo '{path} exists'; else "
            f"UV_PROJECT_ENVIRONMENT={path} {BIN}/uv sync --project {REPOS[repo]} --locked --all-extras {groups}; fi"
        )
    command = (
        f"if [ -d {path}/conda-meta ]; then echo '{path} exists'; else "
        f"{MAMBA} env create -y -p {path} -f {under_repo(recipe.format(**REPOS))}; fi"
    )
    if env == "CALICOST_ENV":
        command += f"; {path}/bin/python -c 'import calicost' || {path}/bin/pip install -e {REPOS['CALICOST_REPO']}"
    if env == "CLONALSCOPE_ENV":
        command += (
            f"; {path}/bin/Rscript -e 'library(Clonalscope)' || "
            f"{path}/bin/Rscript -e 'remotes::install_local(\"{REPOS['CLONALSCOPE_REPO']}\", upgrade = \"never\")'"
        )
    return command


rule build_env:
    input:
        lambda wc: [f"{SETUP}/{env_repo(wc.env)}.ok"] if env_repo(wc.env) else [],
    output:
        f"{SETUP}/{{env}}.ok",
    log:
        f"{LOGS}/setup/build_env.{{env}}.log",
    benchmark:
        f"{LOGS}/setup/build_env.{{env}}.tsv"
    params:
        command=lambda wc: env_command(wc.env),
    shell:
        "{{ {params.command}; }} > {log} 2>&1; touch {output}"


##################################################
# Prepare: configs


# NB: the configs CalicoST's scripts read through scripts/baselines/common/config.sh ($CONFIGS)
rule calicost_configs:
    input:
        under_repo(C["CALICOST_CONFIG"]),
    output:
        datasets=f"{CFG}/baselines/datasets.yaml",
        metadata=f"{CFG}/baselines/metadata.yaml",
        sim=f"{CFG}/baselines/calicost.calicost-sim.yaml",
        ovarian=f"{CFG}/baselines/calicost.ovarian.yaml",
    log:
        f"{LOGS}/prepare/calicost_configs.log",
    run:
        write_yaml(
            output.datasets,
            {
                "datasets": {
                    "calicost-sim": {"deposit": ZENODO, "sim_dir": SIM_DIR},
                    "ovarian": {"samples": {SAMPLE: {"ugp_bb": {16: BB16}}}},
                }
            },
        )
        write_yaml(
            output.metadata,
            {
                "references": {"calicost": CALICOST_REFERENCES},
                "software": {"envs": {"calicost": ENVS["CALICOST_ENV"]}, "calicost_repo": REPOS["CALICOST_REPO"]},
            },
        )
        write_yaml(output.sim, {"dataset": "calicost-sim", "out_root": f"{ZS}/calicost"})
        write_yaml(
            output.ovarian,
            {
                "dataset": "ovarian",
                "staged_dir": f"{OV}/calicost/input",
                "outdir": f"{OV}/calicost",
                "inputs": {name: {"sample": SAMPLE, **spec} for name, spec in CALICOST_INPUTS.items()},
                "configuration_cna": CALICOST["ovarian"]["configuration_cna"],
                "runs": CALICOST_RUNS,
            },
        )
        with open(log[0], "w") as f:
            f.write("".join(f"wrote {path}\n" for path in output))


rule cnamaste_configs:
    output:
        expand(f"{CFG}/cnamaste.{{case}}.yaml", case=CASES),
    log:
        f"{LOGS}/prepare/cnamaste_configs.log",
    params:
        cnamaste=json.dumps(CNAMASTE, sort_keys=True),
    run:
        for case in CASES:
            write_yaml(
                f"{CFG}/cnamaste.{case}.yaml", {"references": CALICOST_REFERENCES, **CNAMASTE[case]["config"]}
            )
        with open(log[0], "w") as f:
            f.write("".join(f"wrote {path}\n" for path in output))


##################################################
# Prepare: references


# NB: Clonalscope's 200kb bins, chromosome-arm segments and data-raw/ at its pinned commit, shared by both cases
rule clonalscope_refs:
    input:
        f"{SETUP}/CLONALSCOPE_ENV.ok",
    output:
        CLONALSCOPE_REFS,
    log:
        f"{LOGS}/prepare/clonalscope_refs.log",
    benchmark:
        f"{LOGS}/prepare/clonalscope_refs.tsv"
    shell:
        "ENV={ENVS[CLONALSCOPE_ENV]} bash {BASELINES}/clonalscope/build_refs_clonalscope.sh {CLONALSCOPE_REF} "
        "> {log} 2>&1"


rule celltyping_markers:
    output:
        f"{OUT}/references/omnipath_panglaodb_annotations.tsv",
    log:
        f"{LOGS}/prepare/celltyping_markers.log",
    params:
        url=C["CELLTYPING_MARKERS_URL"],
    shell:
        "curl -fsSL -o {output} '{params.url}' > {log} 2>&1"


##################################################
# Prepare: simulations, each a folder of the deposit's simulated_data_related/


# NB: Ensembl IDs from Cell Ranger GRCh38-2020-A, which Clonalscope's gene BED joins on
rule sim_10x_matrix:
    input:
        f"{SETUP}/inputs.ok",
        f"{SETUP}/CALICOST_ENV.ok",
        h5ad=f"{SIM_DIR}/{{sim}}/filtered_feature_bc_matrix.h5ad",
    output:
        f"{ZS_INPUTS}/{{sim}}/filtered_feature_bc_matrix/matrix.mtx.gz",
    log:
        f"{LOGS}/prepare/sim_10x_matrix/{{sim}}.log",
    params:
        od=lambda wc: f"{ZS_INPUTS}/{wc.sim}/filtered_feature_bc_matrix",
    shell:
        in_env("CALICOST_ENV")
        + "python -I {BASELINES}/common/h5ad_to_10x.py {input.h5ad} {C[PATH_TO_ORDERED_GENES]} {params.od} > {log} 2>&1"


rule sim_numbat_count_mat:
    input:
        f"{SETUP}/NUMBAT_ENV.ok",
        matrix=rules.sim_10x_matrix.output[0],
    output:
        f"{ZS}/numbat/{{sim}}/prep/{{sim}}.count_mat.rds",
    log:
        f"{LOGS}/prepare/sim_numbat_count_mat/{{sim}}.log",
    params:
        mdir=lambda wc: f"{ZS_INPUTS}/{wc.sim}/filtered_feature_bc_matrix",
    shell:
        in_env("NUMBAT_ENV") + "Rscript {BASELINES}/numbat/build_count_mat.R {params.mdir} {output} > {log} 2>&1"


# NB: the deposit's own numbat_allele_counts.tsv.gz, Numbat's allele input as the deposit ran it
rule sim_numbat_allele_counts:
    input:
        f"{SETUP}/inputs.ok",
        counts=f"{SIM_DIR}/{{sim}}/numbat_allele_counts.tsv.gz",
    output:
        f"{ZS}/numbat/{{sim}}/prep/{{sim}}_allele_counts.tsv.gz",
    log:
        f"{LOGS}/prepare/sim_numbat_allele_counts/{{sim}}.log",
    shell:
        "cp -v {input.counts} {output} > {log} 2>&1"


# NB: run_cnamaste's paths.sample_sheet, the simulation folder its spaceranger_dir and snp_dir
rule sim_cnamaste_sample_sheet:
    input:
        f"{SETUP}/inputs.ok",
        h5ad=f"{SIM_DIR}/{{sim}}/filtered_feature_bc_matrix.h5ad",
    output:
        f"{ZS_INPUTS}/{{sim}}/cnamaste_sample_sheet.tsv",
    log:
        f"{LOGS}/prepare/sim_cnamaste_sample_sheet/{{sim}}.log",
    params:
        d=lambda wc: f"{SIM_DIR}/{wc.sim}/",
    shell:
        "printf 'bam\\tsample_id\\tspaceranger_dir\\tsnp_dir\\nNULL\\t%s\\t%s\\t%s\\n' "
        "{wildcards.sim} {params.d} {params.d} > {output} 2> {log}"


##################################################
# Prepare: the ovarian sample

if BAM16 != BAM16_SPACERANGER:

    rule ovarian_bam_016um:
        input:
            f"{SETUP}/inputs.ok",
            f"{SETUP}/UGP_REPO.ok",
        output:
            bam=BAM16,
            bai=f"{BAM16}.bai",
        log:
            f"{LOGS}/prepare/ovarian_bam_016um.log",
        benchmark:
            f"{LOGS}/prepare/ovarian_bam_016um.tsv"
        shell:
            "{PY} {REPOS[UGP_REPO]}/resources/scripts/convert_bam_visium3prime.py {OUTS}/barcode_mappings.parquet "
            "{OUTS}/possorted_genome_bam.bam {output.bam} square_002um square_016um > {log} 2>&1"


rule ugp_configs:
    input:
        f"{SETUP}/inputs.ok",
        f"{SETUP}/UGP_REPO.ok",
        BAM16,
        under_repo(C["UGP_SAMPLE_TEMPLATE"]),
        under_repo(C["UGP_CONFIG_TEMPLATE"]),
    output:
        sample=f"{CFG}/ugp/sample.json",
        config=f"{CFG}/ugp/config.yaml",
    log:
        f"{LOGS}/prepare/ugp_configs.log",
    run:
        record = dict(UGP_SAMPLE_TEMPLATE["samples"][0], sample_id=SAMPLE)
        record["files"] = {
            "alignment": BAM16,
            "alignment_index": f"{BAM16}.bai",
            "barcodes": f"{B16}/filtered_feature_bc_matrix/barcodes.tsv.gz",
            "matrix_h5": f"{B16}/filtered_feature_bc_matrix.h5",
            "tissue_positions": f"{B16}/spatial/tissue_positions.parquet",
            "scalefactors": f"{B16}/spatial/scalefactors_json.json",
        }
        with open(output.sample, "w") as f:
            json.dump({**UGP_SAMPLE_TEMPLATE, "samples": [record]}, f, indent=2)
        ugp, refdata = REPOS["UGP_REPO"], C["PATH_TO_REFDATA"]
        write_yaml(
            output.config,
            {
                **UGP_CONFIG_TEMPLATE,
                "sample_id": SAMPLE,
                "sample_file": output.sample,
                "reference": f"{refdata}/fasta/genome.fa",
                "gtf_file": f"{refdata}/genes/genes.gtf.gz",
                "genome_size": f"{ugp}/resources/data/hg38.chrom.sizes",
                "region_bed": f"{ugp}/resources/data/hg38.regions.bed",
                "blacklist_bed": f"{ugp}/resources/data/hg38-blacklist.v2.bed.gz",
                "gene_blacklist_file": f"{ugp}/resources/data/ig_gene_list.txt",
                "snp_panel": C["PATH_TO_SNP_PANEL"],
                "phasing_panel": C["PATH_TO_PHASE_PANEL"],
                "gmap_path": C["PATH_TO_EAGLE_MAP"],
                "threads": {k: min(int(v), CORES) for k, v in UGP_CONFIG_TEMPLATE["threads"].items()},
            },
        )
        with open(log[0], "w") as f:
            f.write("".join(f"wrote {path}\n" for path in output))


rule ugp_ovarian:
    input:
        config=rules.ugp_configs.output.config,
        env=f"{SETUP}/UNIVERSAL_GENOTYPING_PIPELINE_ENV.ok",
    output:
        **UGP_OUTPUTS,
    log:
        f"{LOGS}/prepare/ugp_ovarian.log",
    benchmark:
        f"{LOGS}/prepare/ugp_ovarian.tsv"
    threads: CORES
    resources:
        mem_mb=MEM_MB,
    params:
        env=ENVS["UNIVERSAL_GENOTYPING_PIPELINE_ENV"],
        ugp=REPOS["UGP_REPO"],
    shell:
        "PATH={params.env}/bin:{CONDA}/condabin:$PATH snakemake --profile {params.ugp}/profile "
        "-s {params.ugp}/workflow/Snakefile --configfile {input.config} --directory {UGP_DIR} --cores {threads} "
        "--conda-prefix {params.env}/conda --rerun-incomplete > {log} 2>&1"


# NB: Numbat on 32um bins: the 16um matrix and UGP's 16um pileup summed (run_numbat's stats::hclust stops above
#     65536 spots)
rule ovarian_numbat_aggregate_bins:
    input:
        f"{SETUP}/NUMBAT_ENV.ok",
        matrix=f"{B16}/filtered_feature_bc_matrix/matrix.mtx.gz",
        pileup=UGP_OUTPUTS["pileup"],
    output:
        count_mat=f"{OV}/numbat/{SAMPLE}/prep/{SAMPLE}.count_mat.rds",
        pileup=f"{NUMBAT_PILEUP}/cellSNP.tag.AD.mtx",
    log:
        f"{LOGS}/prepare/ovarian_numbat_aggregate_bins.log",
    benchmark:
        f"{LOGS}/prepare/ovarian_numbat_aggregate_bins.tsv"
    params:
        bin_um=NUMBAT["ovarian"]["bin_um"],
    shell:
        in_env("NUMBAT_ENV") + "Rscript {BASELINES}/numbat/aggregate_bins.R {B16}/filtered_feature_bc_matrix "
        "{UGP_PILEUP} {params.bin_um} {output.count_mat} {NUMBAT_PILEUP} > {log} 2>&1"


# NB: numbat:::preprocess_allele on the binned pileup and UGP's phasing
rule ovarian_numbat_allele_counts:
    input:
        f"{SETUP}/NUMBAT_ENV.ok",
        pileup=rules.ovarian_numbat_aggregate_bins.output.pileup,
        phase=UGP_OUTPUTS["phase"],
    output:
        f"{OV}/numbat/{SAMPLE}/prep/{SAMPLE}_allele_counts.tsv.gz",
    log:
        f"{LOGS}/prepare/ovarian_numbat_allele_counts.log",
    benchmark:
        f"{LOGS}/prepare/ovarian_numbat_allele_counts.tsv"
    shell:
        in_env("NUMBAT_ENV")
        + "Rscript {BASELINES}/numbat/build_allele_counts.R {NUMBAT_PILEUP} {UGP_DIR}/phase {C[PATH_TO_EAGLE_MAP]} "
        "{output} > {log} 2>&1"


# NB: the slide images and the 16um positions, the spatial/ UGP's bins are formatted with
rule ovarian_spatial_016um:
    input:
        f"{SETUP}/inputs.ok",
    output:
        directory(f"{OV_INPUTS}/{SAMPLE}/spatial_016um"),
    log:
        f"{LOGS}/prepare/ovarian_spatial_016um.log",
    shell:
        "{{ mkdir -p {output} && cp -v {OUTS}/spatial/* {B16}/spatial/* {output}/; }} > {log} 2>&1"


# NB: UGP's 16um bins in CalicoST's input layout, CalicoST's and cnamaste's input
rule ovarian_format_bb:
    input:
        f"{SETUP}/inputs.ok",
        h5ad=UGP_OUTPUTS["h5ad"],
        spatial=rules.ovarian_spatial_016um.output[0],
    output:
        directory(BB16),
    log:
        f"{LOGS}/prepare/ovarian_format_bb.log",
    params:
        python=C.get("FORMAT_BB_PYTHON") or PY,
    shell:
        "{params.python} {C[FORMAT_BB_SCRIPT]} --bb-dir {UGP_DIR}/bb --assay VISIUM3prime --out-dir {output} "
        "--spatial {DATASET_ID}={input.spatial} --barcode-mappings {OUTS}/barcode_mappings.parquet > {log} 2>&1"


rule ovarian_cnamaste_sample_sheet:
    input:
        BB16,
    output:
        f"{OV_INPUTS}/{SAMPLE}/cnamaste_sample_sheet.tsv",
    log:
        f"{LOGS}/prepare/ovarian_cnamaste_sample_sheet.log",
    shell:
        "printf 'bam\\tsample_id\\tspaceranger_dir\\tsnp_dir\\nNULL\\t%s\\t%s/\\t%s/\\n' "
        "{SAMPLE} {input} {input} > {output} 2> {log}"


rule celltyping_ovarian:
    input:
        f"{SETUP}/inputs.ok",
        f"{SETUP}/CELLTYPING_ENV.ok",
        markers=rules.celltyping_markers.output[0],
    output:
        directory(CELLTYPES),
    log:
        f"{LOGS}/prepare/celltyping_ovarian.log",
    benchmark:
        f"{LOGS}/prepare/celltyping_ovarian.tsv"
    shell:
        "{ENVS[CELLTYPING_ENV]}/bin/python {REPO}/scripts/celltyping/annotate_celltypes.py {B16} {input.markers} "
        "{output} --sample-id {SAMPLE} > {log} 2>&1"


##################################################
# Run: Numbat, per case and sample (scripts/baselines/numbat/)


rule numbat_run:
    input:
        f"{SETUP}/NUMBAT_ENV.ok",
        count_mat=f"{OUT}/{{case}}/numbat/{{s}}/prep/{{s}}.count_mat.rds",
        allele=f"{OUT}/{{case}}/numbat/{{s}}/prep/{{s}}_allele_counts.tsv.gz",
    output:
        touch(f"{OUT}/{{case}}/numbat/{{s}}/{{label}}/run_numbat/done"),
    log:
        f"{OUT}/{{case}}/numbat/logs/run_numbat.{{s}}.{{label}}.log",
    benchmark:
        f"{OUT}/{{case}}/numbat/bench/run_numbat.{{s}}.{{label}}.tsv"
    threads: lambda wc: NUMBAT[CASES[wc.case]]["ncores"]
    params:
        od=lambda wc: f"{OUT}/{wc.case}/numbat/{wc.s}/{wc.label}/run_numbat",
        max_entropy=lambda wc: NUMBAT[CASES[wc.case]]["max_entropy"],
    shell:
        in_env("NUMBAT_ENV")
        + "Rscript {BASELINES}/numbat/run.R {input.count_mat} {input.allele} {params.od} {threads} "
        "{params.max_entropy} > {log} 2>&1"


# NB: run_numbat's outputs -> <sample>.RNA.cells.tsv and cluster_states.tsv
rule numbat_postprocess:
    input:
        rules.numbat_run.output[0],
    output:
        f"{OUT}/{{case}}/numbat/postprocess/{{label}}/{{s}}/{{s}}.RNA.cells.tsv",
    log:
        f"{OUT}/{{case}}/numbat/logs/postprocess.{{s}}.{{label}}.log",
    params:
        rd=lambda wc: f"{OUT}/{wc.case}/numbat/{wc.s}/{wc.label}/run_numbat",
        od=lambda wc: f"{OUT}/{wc.case}/numbat/postprocess/{wc.label}/{wc.s}",
    shell:
        in_env("NUMBAT_ENV")
        + "python {BASELINES}/numbat/postprocess.py {wildcards.s} {params.rd} {params.od} > {log} 2>&1"


##################################################
# Run: Clonalscope, per case and sample, the Visium tumor_spots_tutorial (scripts/baselines/clonalscope/)


# NB: a coarser bin's matrix, the 16um matrix summed into it
rule clonalscope_aggregate_matrix:
    input:
        f"{SETUP}/CLONALSCOPE_ENV.ok",
        matrix=f"{B16}/filtered_feature_bc_matrix/matrix.mtx.gz",
    output:
        f"{OV}/clonalscope/{{s}}/input/filtered_feature_bc_matrix/matrix.mtx.gz",
    log:
        f"{OV}/clonalscope/logs/aggregate_matrix.{{s}}.log",
    params:
        bin_um=lambda wc: CLONALSCOPE_AGGREGATE[wc.s],
        od=lambda wc: f"{OV}/clonalscope/{wc.s}/input/filtered_feature_bc_matrix",
    shell:
        in_env("CLONALSCOPE_ENV")
        + "Rscript {BASELINES}/common/aggregate_10x_matrix.R {B16}/filtered_feature_bc_matrix {params.bin_um} "
        "{params.od} > {log} 2>&1"


# NB: one single-threaded job per sample: BayesNonparCluster is a single-threaded Gibbs sampler
rule clonalscope_run:
    input:
        f"{SETUP}/CLONALSCOPE_ENV.ok",
        CLONALSCOPE_REFS,
        matrix=lambda wc: f"{clonalscope_matrix(wc.case, wc.s)}/matrix.mtx.gz",
    output:
        f"{OUT}/{{case}}/clonalscope/{{s}}/{{label}}/raw_output/{{s}}.cells.tsv",
    log:
        f"{OUT}/{{case}}/clonalscope/logs/run_tutorial.{{s}}.{{label}}.log",
    benchmark:
        f"{OUT}/{{case}}/clonalscope/bench/run_tutorial.{{s}}.{{label}}.tsv"
    params:
        mdir=lambda wc: clonalscope_matrix(wc.case, wc.s),
        od=lambda wc: f"{OUT}/{wc.case}/clonalscope/{wc.s}/{wc.label}/raw_output",
    shell:
        in_env("CLONALSCOPE_ENV")
        + "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript {BASELINES}/clonalscope/tutorial.R {wildcards.s} "
        "{params.mdir} {CLONALSCOPE_DATA_RAW} {params.od} > {log} 2>&1"


# NB: RunCovCluster's result_final -> per-cluster copies; a run that left no Cov_obj.rds (Clonalscope stopped) is
#     skipped, as its log says
rule clonalscope_decode_states:
    input:
        cells=rules.clonalscope_run.output[0],
    output:
        touch(f"{OUT}/{{case}}/clonalscope/{{s}}/{{label}}/decode_states.done"),
    log:
        f"{OUT}/{{case}}/clonalscope/logs/decode_states.{{s}}.{{label}}.log",
    params:
        rds=lambda wc: f"{OUT}/{wc.case}/clonalscope/{wc.s}/{wc.label}/raw_output/{wc.s}.Cov_obj.rds",
        prefix=lambda wc: f"{OUT}/{wc.case}/clonalscope/{wc.s}/{wc.label}/postprocess/{wc.s}.scRNA",
        ploidy=lambda wc: CLONALSCOPE[CASES[wc.case]].get("ploidy", 2),
    shell:
        """
        if [ ! -e {params.rds} ]; then
            echo "no {params.rds}: Clonalscope stopped, skipped" >{log}
            exit 0
        fi
        {IN_CLONALSCOPE} OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript {BASELINES}/clonalscope/decode_states.R \
            {wildcards.s} {params.rds} {input.cells} {params.prefix} {params.ploidy} >{log} 2>&1
        """


##################################################
# Run: CalicoST (scripts/baselines/calicost/)


# NB: the deposit's own configfile0..4 per simulation, paths remapped; parse_input.py, then the rects in parallel
rule calicost_sim:
    input:
        f"{SETUP}/inputs.ok",
        f"{SETUP}/CALICOST_ENV.ok",
        f"{SETUP}/CALICOST_REPO.ok",
        rules.calicost_configs.output,
    output:
        touch(f"{ZS}/calicost/{{sim}}/repro.done"),
    log:
        f"{LOGS}/run/calicost_sim/{{sim}}.log",
    benchmark:
        f"{LOGS}/run/calicost_sim/{{sim}}.tsv"
    threads: len(CALICOST["calicost-sim"]["rects"])
    params:
        rects=" ".join(map(str, CALICOST["calicost-sim"]["rects"])),
    shell:
        "{CALICOST_SHELL} bash {BASELINES}/calicost/repro.sh {wildcards.sim} {ZS}/calicost '{params.rects}' "
        "> {log} 2>&1"


# NB: a staged copy of UGP's formatted bins, unique_snp_ids.npy re-saved for CalicoST's numpy 1.24
rule calicost_stage_input:
    input:
        f"{SETUP}/CALICOST_ENV.ok",
        rules.calicost_configs.output,
        source=lambda wc: CALICOST_INPUT_SOURCES[wc.name],
    output:
        directory(f"{OV}/calicost/input/{{name}}"),
    log:
        f"{LOGS}/run/calicost_stage_input.{{name}}.log",
    shell:
        "{CALICOST_SHELL} bash {BASELINES}/calicost/stage_input.sh {input.source} {output} > {log} 2>&1"


# NB: one run of configs/baselines/calicost.yaml: its configuration_cna written, then run under a cgroup MemoryMax
rule calicost_ovarian:
    input:
        f"{SETUP}/CALICOST_REPO.ok",
        rules.calicost_configs.output,
        staged=lambda wc: f"{OV}/calicost/input/{CALICOST_RUNS[wc.tag]['input']}",
    output:
        touch(f"{OV}/calicost/{{tag}}/run.done"),
    log:
        f"{LOGS}/run/calicost_ovarian.{{tag}}.log",
    benchmark:
        f"{LOGS}/run/calicost_ovarian.{{tag}}.tsv"
    resources:
        mem_mb=MEM_MB,
    params:
        configuration=f"{LOGS}/baselines/calicost_{{tag}}.configuration_cna",
        mem_max=CALICOST["ovarian"]["mem_max"] or f"{C['MEM_GB']}G",
    shell:
        "{{ {CALICOST_SHELL} {PY} -I {BASELINES}/calicost/ovarian_config.py write {wildcards.tag} "
        "{params.configuration} && {CALICOST_SHELL} bash {BASELINES}/calicost/run.sh {wildcards.tag} "
        "{params.configuration} {params.mem_max}; }} "
        "> {log} 2>&1"


##################################################
# Run: cnamaste on the simulations, as cnamaste's batch_zenodo_run.sh runs each job


def cnamaste_over_rides(solver, seed):
    """batch_zenodo_run.sh's overrides of one job: its seed, label solver, initializer and decoder."""
    label_solver, initializer, decoder = solver.split("-")
    decoder_key = {"lattice": "hmm_ab", "calicost": "calicost"}[decoder]
    return (
        f"-o hmrf.random_seed={seed} -o hmrf.label_solver={label_solver} -o hmm.initializer={initializer} "
        f"-o ab_decode.decoder={decoder_key}"
    )


# NB: a run writes run_<id>/ under paths.output_dir; its id is read from the one config_<id>.yaml it wrote
RUN_CNAMASTE = (
    "mkdir -p {params.od}; started=$(mktemp {params.od}/.started.XXXXXX); trap 'rm -f \"$started\"' EXIT; "
    "OMP_NUM_THREADS={threads} NUMBA_NUM_THREADS={threads} " + UV_CNAMASTE + " run_cnamaste {input.config} "
    "-o paths.sample_sheet={input.sheet} -o paths.output_dir={params.od} {params.over_rides} > {log} 2>&1; "
    "find {params.od} -mindepth 2 -maxdepth 2 -path '*/run_*/config_*.yaml' -newer \"$started\" "
    "| sed 's|.*/config_||; s|[.]yaml$||' > {output}; "
    'test "$(wc -l < {output})" -eq 1'
)
PLOT_CNAMASTE = (
    "rid=$(cat {input.rid}); " + UV_CNAMASTE + " run_cnamaste_plot {params.od}/run_$rid/config_$rid.yaml > {log} 2>&1"
)


rule cnamaste_sim:
    input:
        f"{SETUP}/CNAMASTE_ENV.ok",
        sheet=rules.sim_cnamaste_sample_sheet.output[0],
        config=f"{CFG}/cnamaste.zenodo_sim.yaml",
    output:
        f"{CN_RUNS}/{{sim}}/{{solver}}/run_id.seed{{seed}}.txt",
    log:
        f"{LOGS}/run/cnamaste_sim/{{sim}}.{{solver}}.seed{{seed}}.log",
    benchmark:
        f"{CN}/bench/run_cnamaste.{{sim}}.{{solver}}.seed{{seed}}.tsv"
    threads: CNAMASTE["zenodo_sim"]["threads"]
    params:
        od=lambda wc: f"{CN_RUNS}/{wc.sim}/{wc.solver}",
        over_rides=lambda wc: cnamaste_over_rides(wc.solver, wc.seed),
    shell:
        RUN_CNAMASTE


rule cnamaste_sim_plot:
    input:
        rid=rules.cnamaste_sim.output,
        env=f"{SETUP}/CNAMASTE_ENV.ok",
    output:
        touch(f"{CN_RUNS}/{{sim}}/{{solver}}/plot.seed{{seed}}.done"),
    log:
        f"{LOGS}/run/cnamaste_sim_plot/{{sim}}.{{solver}}.seed{{seed}}.log",
    benchmark:
        f"{LOGS}/run/cnamaste_sim_plot/{{sim}}.{{solver}}.seed{{seed}}.tsv"
    params:
        od=lambda wc: f"{CN_RUNS}/{wc.sim}/{wc.solver}",
    shell:
        PLOT_CNAMASTE


rule cnamaste_sim_qa:
    input:
        lambda wc: expand(f"{CN_RUNS}/{{sim}}/{{solver}}/plot.seed{wc.seed}.done", sim=SIMS, solver=SOLVERS),
        config=f"{CFG}/cnamaste.zenodo_sim.yaml",
    output:
        touch(f"{CN}/qa/qa.seed{{seed}}.done"),
    log:
        f"{LOGS}/run/cnamaste_sim_qa.seed{{seed}}.log",
    benchmark:
        f"{LOGS}/run/cnamaste_sim_qa.seed{{seed}}.tsv"
    params:
        baselines=" ".join(["calicost", *[s for s in SOLVERS if s != QA_SOLVER]]),
        ploidy=CNAMASTE["zenodo_sim"]["calicost_ploidy"],
    shell:
        """
        ln -sfn {SIM_DIR} {CN}/simulated_data_related
        ln -sfn {ZENODO}/nomixing_calicost_related {CN}/nomixing_calicost_related
        rc=0
        for baseline in {params.baselines}; do
            {UV_CNAMASTE} run_sim_qa --root {CN} --config {input.config} -o hmrf.random_seed={wildcards.seed} \
                --solver {QA_SOLVER} --baseline_solver $baseline --calicost_ploidy {params.ploidy} --one_shot \
                --overwrite || rc=$?
        done >{log} 2>&1
        exit $rc
        """


##################################################
# Run: cnamaste on the ovarian sample, cnamaste's defaults


rule cnamaste_ovarian:
    input:
        f"{SETUP}/CNAMASTE_ENV.ok",
        sheet=rules.ovarian_cnamaste_sample_sheet.output[0],
        config=f"{CFG}/cnamaste.ovarian.yaml",
    output:
        f"{OV}/cnamaste/run_id.txt",
    log:
        f"{LOGS}/run/cnamaste_ovarian.log",
    benchmark:
        f"{OV}/cnamaste/bench/run_cnamaste.{SAMPLE}.tsv"
    threads: CNAMASTE["ovarian"]["threads"]
    resources:
        mem_mb=MEM_MB,
    params:
        od=f"{OV}/cnamaste",
        over_rides="",
    shell:
        RUN_CNAMASTE


rule cnamaste_ovarian_plot:
    input:
        rid=rules.cnamaste_ovarian.output,
        env=f"{SETUP}/CNAMASTE_ENV.ok",
    output:
        touch(f"{OV}/cnamaste/plot.done"),
    log:
        f"{LOGS}/run/cnamaste_ovarian_plot.log",
    benchmark:
        f"{LOGS}/run/cnamaste_ovarian_plot.tsv"
    params:
        od=f"{OV}/cnamaste",
    shell:
        PLOT_CNAMASTE


# NB: binned_data.npz, which the analysis' cnamaste page reads and cnamaste no longer writes
rule cnamaste_ovarian_binned_data:
    input:
        rid=rules.cnamaste_ovarian.output,
        plot=rules.cnamaste_ovarian_plot.output,
    output:
        touch(f"{OV}/cnamaste/binned_data.done"),
    log:
        f"{LOGS}/run/cnamaste_ovarian_binned_data.log",
    benchmark:
        f"{LOGS}/run/cnamaste_ovarian_binned_data.tsv"
    params:
        od=f"{OV}/cnamaste",
    shell:
        "rid=$(cat {input.rid}); " + UV_CNAMASTE + " python {BASELINES}/cnamaste/write_binned_data.py "
        "{params.od}/run_$rid {params.od}/run_$rid/binned_data.npz > {log} 2>&1"


##################################################
# Run: runtime tables (scripts/metrics/runtime.py)


rule sim_runtime:
    input:
        [path for method in ("numbat", "clonalscope") for path in baseline_outputs("zenodo_sim", method)],
        expand(f"{CN_RUNS}/{{sim}}/{{solver}}/run_id.seed{{seed}}.txt", sim=SIMS, solver=SOLVERS, seed=SEEDS),
    output:
        f"{ZS}/simulation.runtime.tsv",
    log:
        f"{LOGS}/run/sim_runtime.log",
    run:
        from scripts.metrics.runtime import write_simulation_runtime

        with open(log[0], "w") as f, contextlib.redirect_stdout(f):
            write_simulation_runtime(SIM_DIR, ZS, f"{ZENODO}/numbat_related", output[0])


rule vishd_runtime:
    input:
        [path for method in METHODS for path in baseline_outputs("ovarian", method)],
    output:
        f"{OV}/runtime.{SAMPLE}.tsv",
    log:
        f"{LOGS}/run/vishd_runtime.log",
    run:
        from scripts.metrics.runtime import write_vishd_runtime

        with open(log[0], "w") as f, contextlib.redirect_stdout(f):
            write_vishd_runtime(OV, f"{LOGS}/baselines", output[0])


##################################################
# Analysis


def run_dir(path):
    """`path`/run_<id>, the id its run_id file names."""
    with open(f"{path}/run_id.txt" if os.path.isdir(path) else path) as f:
        return f"{os.path.dirname(f.name)}/run_{f.read().strip()}"


rule benchmark_yaml:
    input:
        rules.run_methods.input,
        CELLTYPES,
        f"{OV}/cnamaste/plot.done",
    output:
        benchmark=f"{CFG}/benchmark.yaml",
        config=f"{CFG}/config.yaml",
    log:
        f"{LOGS}/analysis/benchmark_yaml.log",
    run:
        numbat, clonalscope = NUMBAT["ovarian"]["run_label"], CLONALSCOPE["ovarian"]["run_label"]
        clonalscope_sample = SAMPLE if 32 not in CLONALSCOPE["ovarian"]["bins_um"] else f"{SAMPLE}_032um"
        calicost_runs = {
            "calicost" if t == CALICOST_ANALYSIS_RUN else f"calicost_{t.removeprefix(SAMPLE + '_')}": {
                "dir": f"{OV}/calicost/{t}"
            }
            for t in CALICOST_RUNS
        }
        vishd = {
            "runtime": f"{OV}/runtime.{SAMPLE}.tsv",
            "ugp": {
                "dir": UGP_DIR,
                "spots": UGP_OUTPUTS["spots"],
                "spaceranger_dir": B16,
                "bb_dir": {16: f"{UGP_DIR}/bb"},
            },
            "numbat": {"dir": f"{OV}/numbat/{SAMPLE}/{numbat}/run_numbat"},
            **calicost_runs,
            "cnamaste": {"dir": run_dir(f"{OV}/cnamaste")},
            "clonalscope": {"dir": f"{OV}/clonalscope/{clonalscope_sample}/{clonalscope}/postprocess"},
            "cell_type": {"labels": f"{CELLTYPES}/spot_annotation.tsv"},
        }
        seed = SEEDS[0]
        cnamaste_sim = {}
        for solver in SOLVERS:
            run = os.path.basename(run_dir(f"{CN_RUNS}/{SIMS[0]}/{solver}/run_id.seed{seed}.txt"))
            key = "cnamaste" if solver == QA_SOLVER else f"cnamaste-{solver}"
            cnamaste_sim[key] = {"dir": f"{CN_RUNS}/{{sample_id}}/{solver}/{run}", "solver": solver}
        calicost_sim = (
            f"{ZENODO}/nomixing_calicost_related" if ANALYSIS["calicost_sim"] == "deposit" else f"{ZS}/calicost"
        )
        write_yaml(
            output.benchmark,
            {
                "vishd": {SAMPLE: vishd},
                "simulation": {
                    "groundtruth_dir": SIM_DIR,
                    "runtime": f"{ZS}/simulation.runtime.tsv",
                    "numbat": {
                        "dir": f"{ZS}/numbat/{{sample_id}}/{NUMBAT['calicost-sim']['run_label']}/run_numbat"
                    },
                    "clonalscope": {
                        "dir": f"{ZS}/clonalscope/{{sample_id}}/{CLONALSCOPE['calicost-sim']['run_label']}/postprocess"
                    },
                    "calicost": {"dir": f"{calicost_sim}/{{sample_id}}"},
                    **cnamaste_sim,
                },
            },
        )
        write_yaml(
            output.config,
            {
                "reference": {
                    "bed_file": f"{REPO}/data/hg38.regions.bed",
                    "refdata": {"gene_list": {"gtf_file": f"{C['PATH_TO_REFDATA']}/genes/genes.gtf.gz"}},
                }
            },
        )
        with open(log[0], "w") as f:
            f.write("".join(f"wrote {path}\n" for path in output))


def figure_options(benchmark_yaml):
    """The figure step's inputs, read from `benchmark_yaml`'s vishd.<OVARIAN_SAMPLE_ID>, as analysis options."""
    entry = read_yaml(benchmark_yaml)["vishd"][SAMPLE]
    path = lambda p: under_repo(str(p))
    bb = entry["ugp"].get("bb_dir", {})
    return {
        "spaceranger_dir": path(entry["ugp"]["spaceranger_dir"]),
        "spots": path(entry["ugp"]["spots"]),
        "bb_dir": [f"{side}={path(d)}" for side, d in sorted(bb.items())],
        "celltype": path(entry["cell_type"]["labels"]),
        "cnamaste_labels": path(f"{entry['cnamaste']['dir']}/clone_labels_int.tsv.gz"),
        "cnamaste_seglevel": path(f"{entry['cnamaste']['dir']}/cnv_diploid_seglevel.tsv"),
        "calicost_dir": f"{CALICOST_LABEL}={path(entry['calicost']['dir'])}",
        "numbat_dir": path(entry["numbat"]["dir"]),
        "clonalscope_dir": path(entry["clonalscope"]["dir"]),
    }


rule analysis_options:
    input:
        rules.benchmark_yaml.output.benchmark,
        rules.benchmark_yaml.output.config,
    output:
        f"{CFG}/analysis_options.yaml",
    log:
        f"{LOGS}/analysis/analysis_options.log",
    run:
        benchmark_yaml, config_yaml = input
        write_yaml(
            output[0],
            {
                "sample": SAMPLE,
                "benchmark_yaml": benchmark_yaml,
                "config_yaml": config_yaml,
                "outdir": RESULTS,
                **figure_options(benchmark_yaml),
            },
        )
        shell('echo "cnamaste: {ANALYSIS_ROOT} ($(git -C {ANALYSIS_ROOT} describe --always --dirty))" > {log}')


def analysis_step(module, function):
    """The shell command of one analysis step: scripts.analysis.<module>.<function> on the analysis options, in the
    cnamaste analysis environment, its outputs' directories made first."""
    return (
        f"mkdir -p {RESULTS}/tsvs/gene_level_cnp {RESULTS}/text {RESULTS}/plots {RESULTS}/pdf-pages; "
        f"CNAMASTE_ROOT={ANALYSIS_ROOT} PYTHONPATH={ANALYSIS_ROOT}/python:{REPO} {UV_ANALYSIS} "
        "python - > {log} 2>&1 <<'EOF'\n"
        "from scripts.analysis.options import get_cnamaste_root, read_options\n"
        f"from scripts.analysis.{module} import {function}\n"
        "get_cnamaste_root()\n"
        f"{function}(read_options('{{input.options}}'))\n"
        "EOF"
    )


ANALYSIS_ENV = [f"{SETUP}/CNAMASTE_ANALYSIS_REPO.ok", f"{SETUP}/CNAMASTE_ANALYSIS_ENV.ok"]


rule analysis_simulation_metrics:
    input:
        ANALYSIS_ENV,
        options=rules.analysis_options.output[0],
    output:
        [f"{RESULTS}/tsvs/{SIM_DATASET}.tsv.gz", f"{RESULTS}/tsvs/{SIM_DATASET}.median.tsv.gz"],
    log:
        f"{LOGS}/analysis/simulation_metrics.log",
    benchmark:
        f"{LOGS}/analysis/simulation_metrics.tsv"
    shell:
        analysis_step("simulation", "write_simulation_metrics")


rule analysis_real_data_metrics:
    input:
        ANALYSIS_ENV,
        options=rules.analysis_options.output[0],
    output:
        [f"{RESULTS}/tsvs/{REAL_DATASET}.tsv.gz", f"{RESULTS}/tsvs/{REAL_DATASET}.json"],
    log:
        f"{LOGS}/analysis/real_data_metrics.log",
    benchmark:
        f"{LOGS}/analysis/real_data_metrics.tsv"
    shell:
        analysis_step("real_data", "write_real_data_metrics")


rule analysis_tables:
    input:
        ANALYSIS_ENV,
        rules.analysis_simulation_metrics.output,
        rules.analysis_real_data_metrics.output,
        options=rules.analysis_options.output[0],
    output:
        [f"{RESULTS}/text/{name}.{ext}" for name in TEXTS for ext in ("tex", "pdf")],
    log:
        f"{LOGS}/analysis/tables.log",
    benchmark:
        f"{LOGS}/analysis/tables.tsv"
    shell:
        analysis_step("tables", "write_tables")


rule analysis_figures:
    input:
        ANALYSIS_ENV,
        rules.analysis_simulation_metrics.output,
        rules.analysis_real_data_metrics.output,
        options=rules.analysis_options.output[0],
    output:
        [f"{RESULTS}/plots/{name}" for name in FIGURES],
    log:
        f"{LOGS}/analysis/figures.log",
    benchmark:
        f"{LOGS}/analysis/figures.tsv"
    shell:
        analysis_step("figures", "write_figures")
