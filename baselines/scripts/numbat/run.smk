# Numbat run: run_numbat() with ref_hca, genome hg38, defaults except max_entropy (scripts/numbat/run.R).
# Results go to <outdir>/<sample>/<run_label>/run_numbat; prep/ inputs are shared across labels.
#
# Run from baselines/ (scripts/run_baseline.sh numbat <case> runs it):
#   snakemake -s scripts/numbat/run.smk --configfile configs/datasets.yaml configs/metadata.yaml \
#       configs/numbat.<case>.yaml --cores 12

include: "../common/samples.smk"

OUT = config["outdir"]
LABEL = config["run_label"]
RUN = f"{OUT}/{{s}}/{LABEL}/run_numbat"

rule all:
    input: expand(f"{RUN}/done", s=SAMPLES)

rule run_numbat:
    input:
        cm = f"{OUT}/{{s}}/prep/{{s}}.count_mat.rds",
        allele = f"{OUT}/{{s}}/prep/{{s}}_allele_counts.tsv.gz",
    output: touch(f"{RUN}/done")
    params: od = RUN
    threads: config["ncores"]
    log: f"{OUT}/logs/run_numbat.{{s}}.{LABEL}.log"
    benchmark: f"{OUT}/bench/run_numbat.{{s}}.{LABEL}.tsv"
    shell:
        "Rscript {SCRIPTS}/run.R {input.cm} {input.allele} {params.od} {threads} "
        "{config[max_entropy]} > {log} 2>&1"
