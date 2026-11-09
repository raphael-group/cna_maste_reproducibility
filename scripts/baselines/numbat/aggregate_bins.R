#!/usr/bin/env Rscript
# Sum 16um Visium HD bins into coarser square bins: count matrix and cellSNP-lite pileup.
#
# Runpeng Luo (2026-10-24)
#
# Dependencies: Matrix.
#
# Usage:
#   Rscript numbat/aggregate_bins.R <matrix_dir> <pileup_dir> <bin_um> <out_count_rds> <out_pileup_dir>
#           matrix_dir      # Space Ranger square_016um filtered_feature_bc_matrix/
#           pileup_dir      # UGP cellSNP-lite pileup over the 16um barcodes
#           bin_um          # target bin size in um, a multiple of 16 (e.g. 32)
#           out_count_rds   # gene x coarse-bin dgCMatrix (gene-symbol rownames, make.unique)
#           out_pileup_dir  # cellSNP-lite layout over the coarse bins
#
# Inputs
#   matrix_dir: 16um UMI counts; barcodes s_016um_<row>_<col>-1.
#   pileup_dir: cellSNP.base.vcf.gz, cellSNP.tag.{AD,DP,OTH}.mtx, cellSNP.samples.tsv.
# Outputs:
#   out_count_rds: counts summed over the 16um bins of each coarse bin.
#   out_pileup_dir: cellSNP.tag.{AD,DP,OTH}.mtx summed the same way, cellSNP.samples.tsv with the
#     coarse barcodes, cellSNP.base.vcf.gz copied (SNP rows and pseudobulk INFO do not change).
# Notes/References:
#   Binning rule of results/aggregate_bins/README.md (scripts/aggregate_bins.py): 16um bin (r, c)
#   goes to coarse bin (r %/% k, c %/% k), k = bin_um / 16, barcode s_<bin_um:03d>um_<row:05d>_<col:05d>-1.
#   UMIs are counted per 16um barcode, so a UMI split across two 16um bins is counted twice.
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("usage: numbat/aggregate_bins.R <matrix_dir> <pileup_dir> <bin_um> <out_count_rds> <out_pileup_dir>")
matrix_dir <- args[1]; pileup_dir <- args[2]; bin_um <- as.integer(args[3])
out_rds <- args[4]; out_pu <- args[5]
stopifnot(bin_um %% 16 == 0, bin_um > 16)
k <- bin_um %/% 16

coarse_of <- function(bc) {
    rc <- regmatches(bc, regexec("^s_016um_(\\d+)_(\\d+)-1$", bc))
    if (any(lengths(rc) != 3)) stop("unexpected 16um barcode format")
    r <- as.integer(vapply(rc, `[`, "", 2)); c <- as.integer(vapply(rc, `[`, "", 3))
    sprintf("s_%03dum_%05d_%05d-1", bin_um, r %/% k, c %/% k)
}

# 16um x coarse 0/1 indicator; coarse bins sorted by (row, col)
indicator <- function(bc16, levels) {
    j <- match(coarse_of(bc16), levels)
    sparseMatrix(i = seq_along(bc16), j = j, x = 1, dims = c(length(bc16), length(levels)),
                 dimnames = list(bc16, levels))
}

bc16 <- read.table(file.path(matrix_dir, "barcodes.tsv.gz"), stringsAsFactors = FALSE)[, 1]
levels <- sort(unique(coarse_of(bc16)))
features <- read.table(file.path(matrix_dir, "features.tsv.gz"), sep = "\t", stringsAsFactors = FALSE)
cnt <- as(readMM(file.path(matrix_dir, "matrix.mtx.gz")), "CsparseMatrix")
cnt32 <- drop0(cnt %*% indicator(bc16, levels))
dimnames(cnt32) <- list(make.unique(features[, 2]), levels)
cnt32 <- as(cnt32, "dgCMatrix")
stopifnot(sum(cnt32) == sum(cnt))
dir.create(dirname(out_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(cnt32, out_rds)
message("count matrix: ", ncol(cnt), " 16um -> ", ncol(cnt32), " ", bin_um, "um bins; wrote ", out_rds)

pu_bc <- read.table(file.path(pileup_dir, "cellSNP.samples.tsv"), stringsAsFactors = FALSE)[, 1]
M <- indicator(pu_bc, levels)
dir.create(out_pu, showWarnings = FALSE, recursive = TRUE)
for (tag in c("AD", "DP", "OTH")) {
    x <- readMM(file.path(pileup_dir, sprintf("cellSNP.tag.%s.mtx", tag)))
    y <- drop0(as(x, "CsparseMatrix") %*% M)
    stopifnot(sum(y) == sum(x))
    writeMM(as(y, "dgCMatrix"), file.path(out_pu, sprintf("cellSNP.tag.%s.mtx", tag)))
}
writeLines(levels, file.path(out_pu, "cellSNP.samples.tsv"))
file.copy(file.path(pileup_dir, "cellSNP.base.vcf.gz"), file.path(out_pu, "cellSNP.base.vcf.gz"), overwrite = TRUE)
message("pileup: ", length(pu_bc), " 16um -> ", length(levels), " ", bin_um, "um bins; wrote ", out_pu)
