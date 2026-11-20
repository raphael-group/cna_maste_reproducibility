#!/usr/bin/env Rscript
# Sum a Visium HD 16um 10x count matrix into coarser square bins, written as a 10x matrix directory.
#
# Runpeng Luo (2026-10-25)
#
# Dependencies: Matrix.
#
# Usage:
#   Rscript common/aggregate_10x_matrix.R <matrix_dir> <bin_um> <out_dir>
#           matrix_dir  # Space Ranger square_016um filtered_feature_bc_matrix/
#           bin_um      # target bin size in um, a multiple of 16 (e.g. 32)
#           out_dir     # output 10x directory
#
# Inputs
#   matrix_dir: matrix.mtx.gz, barcodes.tsv.gz (s_016um_<row>_<col>-1), features.tsv.gz.
# Outputs:
#   out_dir: matrix.mtx.gz (genes x coarse bins, summed counts), barcodes.tsv.gz
#     (s_<bin_um:03d>um_<row:05d>_<col:05d>-1, sorted by row then col), features.tsv.gz (copied).
# Notes/References:
#   Same binning rule as numbat/aggregate_bins.R and results/aggregate_bins/README.md:
#   16um bin (r, c) goes to coarse bin (r %/% k, c %/% k), k = bin_um / 16.
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: common/aggregate_10x_matrix.R <matrix_dir> <bin_um> <out_dir>")
matrix_dir <- args[1]; bin_um <- as.integer(args[2]); out_dir <- args[3]
stopifnot(bin_um %% 16 == 0, bin_um > 16)
k <- bin_um %/% 16

bc16 <- read.table(file.path(matrix_dir, "barcodes.tsv.gz"), stringsAsFactors = FALSE)[, 1]
rc <- regmatches(bc16, regexec("^s_016um_(\\d+)_(\\d+)-1$", bc16))
if (any(lengths(rc) != 3)) stop("unexpected 16um barcode format")
r <- as.integer(vapply(rc, `[`, "", 2)) %/% k
c <- as.integer(vapply(rc, `[`, "", 3)) %/% k
coarse <- sprintf("s_%03dum_%05d_%05d-1", bin_um, r, c)
levels <- sort(unique(coarse))
M <- sparseMatrix(i = seq_along(bc16), j = match(coarse, levels), x = 1,
                  dims = c(length(bc16), length(levels)))

cnt <- as(readMM(file.path(matrix_dir, "matrix.mtx.gz")), "CsparseMatrix")
agg <- drop0(cnt %*% M)
stopifnot(sum(agg) == sum(cnt))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
mtx_tmp <- file.path(out_dir, "matrix.mtx")
writeMM(as(agg, "dgCMatrix"), mtx_tmp)
system2("gzip", c("-f", mtx_tmp))
gz <- gzfile(file.path(out_dir, "barcodes.tsv.gz"), "w"); writeLines(levels, gz); close(gz)
file.copy(file.path(matrix_dir, "features.tsv.gz"), file.path(out_dir, "features.tsv.gz"), overwrite = TRUE)
message("aggregated ", ncol(cnt), " 16um -> ", ncol(agg), " ", bin_um, "um bins; wrote ", out_dir)
