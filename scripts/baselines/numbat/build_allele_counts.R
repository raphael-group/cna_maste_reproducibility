#!/usr/bin/env Rscript
# Numbat allele-count table from UGP's cellSNP-lite pileup and Eagle2-phased VCFs.
#
# Runpeng Luo (2026-10-24)
#
# Dependencies: numbat (>= 1.5.2), data.table, dplyr, stringr, Matrix.
#
# Usage:
#   Rscript numbat/build_allele_counts.R <pileup_dir> <phase_dir> <gmap> <out_tsv_gz>
#           pileup_dir  # UGP pileup/<assay>_<rep>/ (cellSNP.base.vcf.gz, cellSNP.tag.{AD,DP}.mtx, cellSNP.samples.tsv)
#           phase_dir   # UGP phase/ (chr{1..22}.vcf.gz, phased GT in one sample column)
#           gmap        # Eagle2 genetic map (genetic_map_hg38_withX.txt.gz)
#           out_tsv_gz  # numbat df_allele table
#
# Inputs
#   pileup_dir: per-spot ALT (AD) and total (DP) counts at UGP's genotyped SNPs.
#   phase_dir: UGP Eagle2-phased genotypes; only 0|1 and 1|0 sites are kept by numbat.
# Outputs:
#   out_tsv_gz: cell, snp_id, CHROM, POS, cM, REF, ALT, AD, DP, GT, gene (heterozygous SNPs only).
# Notes/References:
#   Reproduces the final "Generate allele count dataframes" block of numbat inst/bin/pileup_and_phase.R
#   (numbat 1.5.2): same vcf/gmap parsing, numbat:::preprocess_allele with gtf_hg38, same GT filter.
#   Only the pileup and phasing inputs differ: they come from UGP instead of numbat's own
#   cellSNP-lite + numbat:::genotype + Eagle2 run.
#   https://github.com/kharchenkolab/numbat/blob/main/inst/bin/pileup_and_phase.R
suppressPackageStartupMessages({
    library(numbat); library(data.table); library(dplyr); library(stringr); library(Matrix); library(glue)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: numbat/build_allele_counts.R <pileup_dir> <phase_dir> <gmap> <out_tsv_gz>")
pu_dir <- args[1]; phase_dir <- args[2]; gmap <- args[3]; out_file <- args[4]

genetic_map <- fread(gmap) %>%
    setNames(c("CHROM", "POS", "rate", "cM")) %>%
    group_by(CHROM) %>%
    mutate(start = POS, end = c(POS[2:length(POS)], POS[length(POS)])) %>%
    ungroup()

vcf_phased <- lapply(1:22, function(chr) {
        fread(glue("{phase_dir}/chr{chr}.vcf.gz"), skip = "#CHROM") %>%
            rename(CHROM = `#CHROM`) %>%
            mutate(CHROM = str_remove(CHROM, "chr"))
    }) %>%
    Reduce(rbind, .) %>%
    mutate(CHROM = factor(CHROM, unique(CHROM)))
gt_col <- colnames(vcf_phased)[10]
message("phased VCF: ", nrow(vcf_phased), " SNPs; GT column '", gt_col, "'")

vcf_pu <- fread(glue("{pu_dir}/cellSNP.base.vcf.gz"), skip = "#CHROM") %>%
    rename(CHROM = `#CHROM`) %>%
    mutate(CHROM = str_remove(CHROM, "chr"))
AD <- readMM(glue("{pu_dir}/cellSNP.tag.AD.mtx"))
DP <- readMM(glue("{pu_dir}/cellSNP.tag.DP.mtx"))
cell_barcodes <- fread(glue("{pu_dir}/cellSNP.samples.tsv"), header = FALSE) %>% pull(V1)
message("pileup: ", nrow(vcf_pu), " SNPs x ", length(cell_barcodes), " spots")

df <- numbat:::preprocess_allele(
    sample = gt_col, vcf_pu = vcf_pu, vcf_phased = vcf_phased, AD = AD, DP = DP,
    barcodes = cell_barcodes, gtf = gtf_hg38, gmap = genetic_map
) %>%
    filter(GT %in% c("1|0", "0|1"))

dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
fwrite(df, out_file, sep = "\t")
message("wrote ", out_file, ": ", nrow(df), " rows, ", length(unique(df$snp_id)), " het SNPs, ",
        length(unique(df$cell)), " spots")
