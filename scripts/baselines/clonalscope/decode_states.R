#!/usr/bin/env Rscript
# Decode one Clonalscope run (RunCovCluster result_final): cluster states (Uest), TCN and benchmark TSVs.
#
# Runpeng Luo (2026-10-31)
#
# Dependencies: Clonalscope (07dd0c5).
#
# Usage:
#   Rscript clonalscope/decode_states.R <sample> <cov_obj_rds> <cells_tsv> <out_prefix> [ploidy]
#           sample        # sample id, used in log lines
#           cov_obj_rds   # <run_dir>/<sample>.Cov_obj.rds written by clonalscope/tutorial.R
#           cells_tsv     # <run_dir>/<sample>.cells.tsv (BARCODE, malignant, cluster)
#           out_prefix    # e.g. <pp_dir>/<sample>.scRNA; files are <out_prefix>.<kind>
#           ploidy        # baseline total copy number of the normal reference (default 2)
#
# Inputs
#   cov_obj_rds: list with result_final = (df_obj, result); df_obj$df is the spot x segment
#     coverage-ratio matrix the clustering ran on, df_obj$seg_table_filtered its segments, and
#     result is AssignCluster output (Zest per spot, Uest cluster x segment).
#   cells_tsv: per-spot cluster (Zest) and MalignantAssignment T/N call.
# Parameters
#   ploidy: TCN = round(ploidy * Uest), floored at 0.
# Outputs:
#   out_prefix.cells.tsv: BARCODE, malignant, cluster (copied from cells_tsv).
#   out_prefix.coverage.tsv.gz: BARCODE x REGION coverage ratio (df_obj$df).
#   out_prefix.segments.tsv: #CHR, START, END, LENGTH, REGION; rows match coverage columns.
#   out_prefix.cluster_states.tsv: REGION x cluster<Zest> Uest (only clusters present in Zest).
#   out_prefix.profile.tsv: CLUSTER, MALIGNANT, #CHR, START, END, REGION, REL_CN (Uest), TCN.
# Notes/References:
#   https://github.com/seasoncloud/Clonalscope/blob/identify_normal_cells_noWGS/R/RunCovCluster.R
#   Clonalscope gives no integer copy number for coverage-only data: AssignCluster returns Uest,
#   the per-cluster mean relative coverage (1 = same as the normal reference). RunCovCluster encodes
#   gain/loss priors as 1.5/0.5, i.e. ratio = CN / 2 for a diploid reference; TCN applies that
#   relation in reverse and is a derived value, not a Clonalscope output.
#   Uest rows are indexed by Zest cluster id; rows of unused clusters are NA and are dropped.
#   Spots removed by PrepCovMatrix (ngene / ncell filters) are absent from every output.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(4, 5)) stop("usage: clonalscope/decode_states.R <sample> <cov_obj_rds> <cells_tsv> <out_prefix> [ploidy]")
sample <- args[1]; rds <- args[2]; cells_tsv <- args[3]; prefix <- args[4]
ploidy <- if (length(args) == 5) as.numeric(args[5]) else 2
dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("[", sample, "] wrote ", path)
}

fin <- readRDS(rds)$result_final
cov <- as.matrix(fin$df_obj$df)
seg <- fin$df_obj$seg_table_filtered
res <- fin$result
cells <- read.table(cells_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
regions <- colnames(res$Uest)
if (!identical(regions, colnames(cov))) stop("Uest and coverage matrix columns differ")
if (!all(cells$BARCODE %in% rownames(cov))) stop("cells.tsv barcodes missing from the coverage matrix")

write_tsv(cells, paste0(prefix, ".cells.tsv"))

con <- gzfile(paste0(prefix, ".coverage.tsv.gz"), "w")
write.table(data.frame(BARCODE = rownames(cov), cov, check.names = FALSE), con,
            sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
message("[", sample, "] wrote ", prefix, ".coverage.tsv.gz (", nrow(cov), " spots x ", ncol(cov), " segments)")

# NB: Uest columns are named chr<chr>-<start>-<end> from seg_table_filtered
coord <- do.call(rbind, strsplit(sub("^chr", "", regions), "-", fixed = TRUE))
si <- match(paste(coord[, 1], as.numeric(coord[, 2])), paste(seg$chr, as.numeric(seg$start)))
if (anyNA(si)) stop("segments not found in seg_table_filtered$chrr: ", paste(regions[is.na(si)], collapse = ","))
segs <- data.frame(`#CHR` = seg$chr[si], START = as.numeric(coord[, 2]), END = as.numeric(coord[, 3]),
                   LENGTH = seg$length[si], REGION = regions, check.names = FALSE)
write_tsv(segs, paste0(prefix, ".segments.tsv"))

ids <- sort(unique(as.integer(cells$cluster)))
U <- res$Uest[ids, , drop = FALSE]
states <- data.frame(REGION = regions, t(U), check.names = FALSE)
colnames(states)[-1] <- paste0("cluster", ids)
write_tsv(states, paste0(prefix, ".cluster_states.tsv"))

mal <- tapply(cells$malignant, cells$cluster, function(x) names(which.max(table(x))))
prof <- do.call(rbind, lapply(seq_along(ids), function(i) data.frame(
  CLUSTER = ids[i], MALIGNANT = mal[as.character(ids[i])], segs[, c("#CHR", "START", "END", "REGION")],
  REL_CN = U[i, ], TCN = pmax(0, round(ploidy * U[i, ])), check.names = FALSE)))
write_tsv(prof, paste0(prefix, ".profile.tsv"))
message("[", sample, "] ", length(ids), " clusters, ", nrow(segs), " segments, ",
        sum(prof$TCN != ploidy), " non-", ploidy, " cluster-segment calls")
