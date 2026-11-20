#!/usr/bin/env Rscript
# Run numbat (scRNA mode) on one sample: run_numbat() with the ref_hca reference, defaults except max_entropy if given.
#
# Runpeng Luo (2026-10-24)
#
# Dependencies: numbat (>= 1.5.2), data.table, Matrix.
#
# Usage:
#   Rscript numbat/run.R <count_mat_rds> <allele_counts_tsv_gz> <out_dir> <ncores> [max_entropy]
#           count_mat_rds         # gene x spot dgCMatrix from numbat/build_count_mat.R
#           allele_counts_tsv_gz  # <sample>_allele_counts.tsv.gz from numbat/build_allele_counts.R
#           out_dir               # run_numbat() output directory
#           ncores                # run_numbat() ncores (runtime only)
#           max_entropy           # run_numbat() max_entropy; omitted = package default (0.5)
#
# Inputs
#   count_mat_rds: raw UMI counts, gene-symbol rownames.
#   allele_counts_tsv_gz: phased het-SNP allele counts per spot.
# Parameters
#   lambdas_ref = numbat::ref_hca (the package's built-in normal expression reference),
#   genome = "hg38"; max_entropy only when given; every other run_numbat() argument is the package default.
# Outputs:
#   out_dir: numbat native outputs (clone_post_*.tsv, joint_post_*.tsv, segs_consensus_*.tsv,
#     bulk_clones_final.tsv.gz, plots).
# Notes/References:
#   https://kharchenkolab.github.io/numbat/articles/numbat.html
#   max_entropy 0.8 for Visium: https://kharchenkolab.github.io/numbat/articles/spatial-rna.html
suppressPackageStartupMessages({ library(numbat); library(Matrix); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(4, 5)) stop("usage: numbat/run.R <count_mat_rds> <allele_counts_tsv_gz> <out_dir> <ncores> [max_entropy]")
count_mat <- readRDS(args[1])
df_allele <- fread(args[2], data.table = FALSE)
out_dir <- args[3]; ncores <- as.integer(args[4])
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message("count_mat: ", nrow(count_mat), " genes x ", ncol(count_mat), " spots; df_allele: ",
        nrow(df_allele), " rows over ", length(unique(df_allele$cell)), " spots; ncores = ", ncores)

extra <- if (length(args) == 5) list(max_entropy = as.numeric(args[5])) else list()
message("non-default run_numbat args: ", if (length(extra)) paste(names(extra), extra, sep = "=", collapse = ", ") else "none")
do.call(run_numbat, c(list(count_mat = count_mat, lambdas_ref = numbat::ref_hca, df_allele = df_allele,
                         genome = "hg38", ncores = ncores, out_dir = out_dir), extra))
message("run_numbat done")
