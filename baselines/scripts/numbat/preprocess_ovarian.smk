# Numbat preprocess (ovarian): sum the 16um count matrix and UGP's 16um pileup into <bin_um> bins, then
# the allele-count table from that pileup + UGP's phasing (numbat:::preprocess_allele).
#
# Run from baselines/ (scripts/run_baseline.sh numbat ovarian runs it):
#   snakemake -s scripts/numbat/preprocess_ovarian.smk --configfile configs/datasets.yaml configs/metadata.yaml \
#       configs/numbat.ovarian.yaml --cores 2

include: "../common/samples.smk"

OUT = config["outdir"]
BIN = int(config["bin_um"])
PU = f"{OUT}/{{s}}/prep/pileup_{BIN:03d}um"

rule all:
    input:
        expand(f"{OUT}/{{s}}/prep/{{s}}_allele_counts.tsv.gz", s=SAMPLES),
        expand(f"{OUT}/{{s}}/prep/{{s}}.count_mat.rds", s=SAMPLES),

rule aggregate_bins:
    input:
        mtx = lambda wc: f"{MATRIX_DIR[wc.s]}/matrix.mtx.gz",
        ad = lambda wc: f"{DATASET['samples'][wc.s]['ugp_pileup']}/cellSNP.tag.AD.mtx",
    output:
        cm = f"{OUT}/{{s}}/prep/{{s}}.count_mat.rds",
        ad = f"{PU}/cellSNP.tag.AD.mtx",
    params: mdir = lambda wc: MATRIX_DIR[wc.s], pu = lambda wc: DATASET["samples"][wc.s]["ugp_pileup"], od = PU
    log: f"{OUT}/logs/aggregate_bins.{{s}}.log"
    benchmark: f"{OUT}/bench/aggregate_bins.{{s}}.tsv"
    shell: "Rscript {SCRIPTS}/aggregate_bins.R {params.mdir} {params.pu} {BIN} {output.cm} {params.od} > {log} 2>&1"

rule allele_counts:
    input:
        ad = f"{PU}/cellSNP.tag.AD.mtx",
        ph = lambda wc: f"{DATASET['samples'][wc.s]['ugp_phase']}/chr1.vcf.gz",
    output: f"{OUT}/{{s}}/prep/{{s}}_allele_counts.tsv.gz"
    params: pu = PU, ph = lambda wc: DATASET["samples"][wc.s]["ugp_phase"], gmap = REFERENCES["eagle_gmap"]
    log: f"{OUT}/logs/allele_counts.{{s}}.log"
    benchmark: f"{OUT}/bench/allele_counts.{{s}}.tsv"
    shell: "Rscript {SCRIPTS}/build_allele_counts.R {params.pu} {params.ph} {params.gmap} {output} > {log} 2>&1"
