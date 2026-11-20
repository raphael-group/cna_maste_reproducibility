#!/usr/bin/env Rscript
# Gene x spot UMI count matrix (dgCMatrix, gene-symbol rownames) for numbat run_numbat().
#
# Runpeng Luo (2026-10-24)
#
# Dependencies: Matrix.
#
# Usage:
#   Rscript numbat/build_count_mat.R <matrix_dir> <out_rds>
#           matrix_dir  # 10x filtered_feature_bc_matrix/ (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
#           out_rds     # output dgCMatrix
#
# Inputs
#   matrix_dir: Space Ranger binned counts; features column 2 holds gene symbols.
# Outputs:
#   out_rds: dgCMatrix, rows = gene symbols made unique with make.unique (Seurat::Read10X
#     convention, so they match numbat's symbol-keyed ref_hca), columns = spot barcodes.
# Notes/References:
#   https://kharchenkolab.github.io/numbat/articles/numbat.html (count_mat input spec)
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: numbat/build_count_mat.R <matrix_dir> <out_rds>")
matrix_dir <- args[1]; out_rds <- args[2]

mtx <- as(readMM(file.path(matrix_dir, "matrix.mtx.gz")), "CsparseMatrix")
barcodes <- read.table(file.path(matrix_dir, "barcodes.tsv.gz"), sep = "\t", stringsAsFactors = FALSE)[, 1]
features <- read.table(file.path(matrix_dir, "features.tsv.gz"), sep = "\t", stringsAsFactors = FALSE)
dimnames(mtx) <- list(make.unique(features[, 2]), barcodes)
dir.create(dirname(out_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(mtx, out_rds)
message("wrote ", out_rds, ": ", nrow(mtx), " genes x ", ncol(mtx), " spots")
