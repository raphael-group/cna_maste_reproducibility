# Clonalscope run: Visium tumor_spots_tutorial (no paired DNA) per sample, one single-threaded job each.
# FindNormalReference(pca) -> arm segments -> RunCovCluster -> MalignantAssignment
# (scripts/clonalscope/tutorial.R). Optional config `aggregate` builds a coarse-bin matrix first.
#
# Run from baselines/ (scripts/run_baseline.sh clonalscope <case> runs it):
#   snakemake -s scripts/clonalscope/run.smk --configfile configs/datasets.yaml configs/metadata.yaml \
#       configs/clonalscope.<case>.yaml --cores 8

include: "../common/samples.smk"

OUT = config["outdir"]
LABEL = config["run_label"]
RAW = f"{OUT}/{{s}}/{LABEL}/raw_output"
AGG = config.get("aggregate", {})
MATRIX_DIR.update({s: f"{OUT}/{s}/input/filtered_feature_bc_matrix" for s in AGG})

rule all:
    input: expand(f"{RAW}/{{s}}.cells.tsv", s=SAMPLES)

if AGG:
    rule aggregate_matrix:
        input: lambda wc: f"{DATASET['samples'][AGG[wc.s]['source']]['matrix_016um']}/matrix.mtx.gz"
        output: f"{OUT}/{{s}}/input/filtered_feature_bc_matrix/matrix.mtx.gz"
        wildcard_constraints: s = "|".join(AGG)
        params:
            src = lambda wc: DATASET["samples"][AGG[wc.s]["source"]]["matrix_016um"],
            bin_um = lambda wc: AGG[wc.s]["bin_um"],
            od = f"{OUT}/{{s}}/input/filtered_feature_bc_matrix",
        log: f"{OUT}/logs/aggregate_matrix.{{s}}.log"
        shell: "Rscript {COMMON}/aggregate_10x_matrix.R {params.src} {params.bin_um} {params.od} > {log} 2>&1"

rule run_tutorial:
    input: lambda wc: f"{MATRIX_DIR[wc.s]}/matrix.mtx.gz"
    output: f"{RAW}/{{s}}.cells.tsv"
    params: mdir = lambda wc: MATRIX_DIR[wc.s], od = RAW, raw = REFERENCES["clonalscope_data_raw"]
    threads: 1
    log: f"{OUT}/logs/run_tutorial.{{s}}.{LABEL}.log"
    benchmark: f"{OUT}/bench/run_tutorial.{{s}}.{LABEL}.tsv"
    shell:
        "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "Rscript {SCRIPTS}/tutorial.R {wildcards.s} {params.mdir} "
        "{params.raw} {params.od} > {log} 2>&1"
