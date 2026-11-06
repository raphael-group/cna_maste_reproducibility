# Numbat postprocess: run_numbat outputs -> <sample>.RNA.cells.tsv + cluster_states.tsv
# under <outdir>/postprocess/<run_label>/<sample>/ (scripts/numbat/postprocess.py).
#
# Run from baselines/ (scripts/run_baseline.sh numbat <case> runs it):
#   snakemake -s scripts/numbat/postprocess.smk --configfile configs/datasets.yaml configs/metadata.yaml \
#       configs/numbat.<case>.yaml --cores 1

include: "../common/samples.smk"

OUT = config["outdir"]
LABEL = config["run_label"]
RUN = f"{OUT}/{{s}}/{LABEL}/run_numbat"
PP = f"{OUT}/postprocess/{LABEL}/{{s}}"

rule all:
    input: expand(f"{PP}/{{s}}.RNA.cells.tsv", s=SAMPLES)

rule postprocess:
    input: f"{RUN}/done"
    output: f"{PP}/{{s}}.RNA.cells.tsv"
    params: rd = RUN, od = PP
    log: f"{OUT}/logs/postprocess.{{s}}.{LABEL}.log"
    shell: "python {SCRIPTS}/postprocess.py {wildcards.s} {params.rd} {params.od} > {log} 2>&1"
