# The CalicoST simulations' inputs for every baseline, or for one: `methods` = all (default), numbat, clonalscope or
# calicost.
#   numbat: each simulation's 10x matrix directory from its h5ad, then <numbat outdir>/<sim>/prep/<sim>.count_mat.rds
#       and <sim>_allele_counts.tsv.gz, the deposit's numbat_allele_counts.tsv.gz as is (read by numbat/run.smk).
#   clonalscope: each simulation's 10x matrix directory, and Clonalscope's reference directory
#       (clonalscope/build_refs_clonalscope.sh), the parent of references.clonalscope_data_raw.
#   calicost: none; CalicoST reads the deposit's simulation folders as they are (calicost/repro.sh).
#
# Config: configs/datasets.yaml (sim_dir, inputs_dir) and configs/metadata.yaml (references.ordered_genes,
#   references.clonalscope_data_raw, the calicost, numbat and clonalscope envs); the numbat outdir from
#   configs/numbat.calicost-sim.yaml. A method config passed with --configfile is not needed and changes nothing.
#
# Run from baselines/ (scripts/run_baseline.sh <method> calicost-sim runs it with methods=<method>):
#   snakemake -s scripts/common/prepare_sim_inputs.smk --configfile configs/datasets.yaml configs/metadata.yaml \
#       --cores 8 [--config methods=numbat]

import yaml

config.setdefault("dataset", "calicost-sim")

include: "samples.smk"

METHODS = ["numbat", "clonalscope", "calicost"]
method = str(config.get("methods", "all"))
if method not in ["all", *METHODS]:
    raise ValueError(f"methods={method}: one of all, {', '.join(METHODS)}")
SELECTED = METHODS if method == "all" else [method]

with open(os.path.join(os.path.dirname(os.path.dirname(COMMON)), "configs", f"numbat.{config['dataset']}.yaml")) as f:
    NUMBAT_OUT = yaml.safe_load(f)["outdir"]
DATA_RAW = REFERENCES["clonalscope_data_raw"]
REF_DIR = os.path.dirname(DATA_RAW)
CLONALSCOPE_REFS = [
    f"{REF_DIR}/hg38.200kb.bins.bed",
    f"{REF_DIR}/seg_table.arm.tsv",
    *[f"{DATA_RAW}/{f}" for f in ["sizes.cellranger-GRCh38-1.0.0.txt", "cyclegenes.rds", "hg38.genes.bed",
                                  "cytoarm_table_hg38.txt"]],
]
TARGETS = {
    "numbat": [
        *[f"{MATRIX_DIR[s]}/matrix.mtx.gz" for s in SAMPLES],
        *[f"{NUMBAT_OUT}/{s}/prep/{s}.count_mat.rds" for s in SAMPLES],
        *[f"{NUMBAT_OUT}/{s}/prep/{s}_allele_counts.tsv.gz" for s in SAMPLES],
    ],
    "clonalscope": [*[f"{MATRIX_DIR[s]}/matrix.mtx.gz" for s in SAMPLES], *CLONALSCOPE_REFS],
    "calicost": [],
}

rule all:
    input: sorted({t for m in SELECTED for t in TARGETS[m]})

rule h5ad_to_10x:
    input: f"{DATASET['sim_dir']}/{{s}}/filtered_feature_bc_matrix.h5ad"
    output: f"{DATASET['inputs_dir']}/{{s}}/filtered_feature_bc_matrix/matrix.mtx.gz"
    params:
        od = f"{DATASET['inputs_dir']}/{{s}}/filtered_feature_bc_matrix",
        py = f"{SOFTWARE['envs']['calicost']}/bin/python",
        genes = REFERENCES["ordered_genes"],
    log: f"{DATASET['inputs_dir']}/{{s}}/h5ad_to_10x.log"
    shell:
        "{params.py} -I {COMMON}/h5ad_to_10x.py {input} {params.genes} {params.od} > {log} 2>&1"

rule numbat_count_mat:
    input: lambda wc: f"{MATRIX_DIR[wc.s]}/matrix.mtx.gz"
    output: f"{NUMBAT_OUT}/{{s}}/prep/{{s}}.count_mat.rds"
    params:
        mdir = lambda wc: MATRIX_DIR[wc.s],
        rscript = f"{SOFTWARE['envs']['numbat']}/bin/Rscript",
        script = os.path.join(os.path.dirname(COMMON), "numbat", "build_count_mat.R"),
    log: f"{NUMBAT_OUT}/logs/count_mat.{{s}}.log"
    shell: "{params.rscript} {params.script} {params.mdir} {output} > {log} 2>&1"

rule numbat_allele_counts:
    input: f"{DATASET['sim_dir']}/{{s}}/numbat_allele_counts.tsv.gz"
    output: f"{NUMBAT_OUT}/{{s}}/prep/{{s}}_allele_counts.tsv.gz"
    shell: "cp {input} {output}"

rule clonalscope_refs:
    output: CLONALSCOPE_REFS
    params:
        script = os.path.join(os.path.dirname(COMMON), "clonalscope", "build_refs_clonalscope.sh"),
        env = SOFTWARE["envs"]["clonalscope"],
        ref_dir = lambda w, output: os.path.dirname(output[0]),
    log: f"{REF_DIR}/build_refs_clonalscope.log"
    shell:
        "ENV={params.env} bash {params.script} {params.ref_dir} > {log} 2>&1"
