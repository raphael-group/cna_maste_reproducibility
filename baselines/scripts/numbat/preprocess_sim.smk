# Numbat preprocess (CalicoST simulations): count matrix from the converted 10x directory and the
# deposit's simulated numbat_allele_counts.tsv.gz as is. Writes the same prep/ files as numbat/preprocess_ovarian.smk.
#
# Run from baselines/ (scripts/run_baseline.sh numbat calicost-sim runs it):
#   snakemake -s scripts/numbat/preprocess_sim.smk --configfile configs/datasets.yaml configs/metadata.yaml \
#       configs/numbat.calicost-sim.yaml --cores 8

include: "../common/samples.smk"

OUT = config["outdir"]

rule all:
    input:
        expand(f"{OUT}/{{s}}/prep/{{s}}_allele_counts.tsv.gz", s=SAMPLES),
        expand(f"{OUT}/{{s}}/prep/{{s}}.count_mat.rds", s=SAMPLES),

rule count_mat:
    input: lambda wc: f"{MATRIX_DIR[wc.s]}/matrix.mtx.gz"
    output: f"{OUT}/{{s}}/prep/{{s}}.count_mat.rds"
    params: mdir = lambda wc: MATRIX_DIR[wc.s]
    log: f"{OUT}/logs/count_mat.{{s}}.log"
    shell: "Rscript {SCRIPTS}/build_count_mat.R {params.mdir} {output} > {log} 2>&1"

rule allele_counts:
    input: f"{DATASET['sim_dir']}/{{s}}/numbat_allele_counts.tsv.gz"
    output: f"{OUT}/{{s}}/prep/{{s}}_allele_counts.tsv.gz"
    shell: "cp {input} {output}"
