# Clonalscope state decoding: RunCovCluster result_final -> per-cluster relative CN (Uest) and
# TCN = round(ploidy * Uest), plus cells / coverage / segments TSVs (scripts/clonalscope/decode_states.R).
# Samples whose run left no Cov_obj.rds (Clonalscope errored) are skipped.
#
# Run from baselines/ (scripts/run_baseline.sh clonalscope <case> runs it):
#   snakemake -s scripts/clonalscope/decode_states.smk --configfile configs/datasets.yaml configs/metadata.yaml \
#       configs/clonalscope.<case>.yaml --cores 8

include: "../common/samples.smk"

OUT = config["outdir"]
LABEL = config["run_label"]
RAW = f"{OUT}/{{s}}/{LABEL}/raw_output"
PP = f"{OUT}/{{s}}/{LABEL}/postprocess/{{s}}.scRNA"
DONE = [s for s in SAMPLES if os.path.exists(f"{OUT}/{s}/{LABEL}/raw_output/{s}.Cov_obj.rds")]

rule all:
    input: expand(f"{PP}.profile.tsv", s=DONE)

rule decode_states:
    input:
        rds = f"{RAW}/{{s}}.Cov_obj.rds",
        cells = f"{RAW}/{{s}}.cells.tsv",
    output: f"{PP}.profile.tsv"
    params: prefix = PP, ploidy = config.get("ploidy", 2)
    log: f"{OUT}/logs/decode_states.{{s}}.{LABEL}.log"
    shell:
        "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "Rscript {SCRIPTS}/decode_states.R {wildcards.s} {input.rds} {input.cells} "
        "{params.prefix} {params.ploidy} > {log} 2>&1"
