#!/usr/bin/env Rscript
# Clonalscope on one Visium sample, following the tumor_spots_tutorial (no paired DNA) call for call.
#
# Runpeng Luo (2026-10-30)
#
# Dependencies: Clonalscope (07dd0c5), Matrix.
#
# Usage:
#   Rscript clonalscope/tutorial.R <sample> <matrix_dir> <data_raw_dir> <out_dir>
#           sample        # sample id, used in file names and log lines
#           matrix_dir    # 10x filtered_feature_bc_matrix/ (features column 1 = Ensembl ID, column 2 = symbol)
#           data_raw_dir  # Clonalscope data-raw/ at the installed commit (sizes, cyclegenes, genes bed, cytoarm table)
#           out_dir       # RunCovCluster dir_path and result files
#
# Inputs
#   matrix_dir: gene x spot UMI counts.
#   data_raw_dir: sizes.cellranger-GRCh38-1.0.0.txt, cyclegenes.rds, hg38.genes.bed, cytoarm_table_hg38.txt.
# Parameters
#   FindNormalReference(method = "pca"); arm-level seg_table_filtered;
#   RunCovCluster(alpha = 5, beta = 5, var_pt = 0.99, var_pt_ctrl = 0.99, include = "all",
#     alpha_source = "all", ctrl_region = NULL, breaks = 50, prep_mode = "intersect", est_cap = 2,
#     clust_mode = "all"), every other argument at its default; MalignantAssignment(cutoff = 0.5).
# Outputs:
#   out_dir/<sample>.cell_types_binary.tsv: barcode, type (FindNormalReference normal spots = "normal").
#   out_dir/<sample>.Cov_obj.rds: RunCovCluster result.
#   out_dir/<sample>.cells.tsv: BARCODE, malignant (T/N from MalignantAssignment), cluster (Zest).
# Notes/References:
#   https://github.com/seasoncloud/Clonalscope/blob/identify_normal_cells_noWGS/samples/V11Y04-378-A1/tumor_spots_tutorial/README.md
#   Only file paths differ from the tutorial; SpatialPlot is not called.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: clonalscope/tutorial.R <sample> <matrix_dir> <data_raw_dir> <out_dir>")
sample <- args[1]; mdir <- args[2]; raw <- args[3]; dir_path <- args[4]
suppressPackageStartupMessages({ library(Clonalscope); library(Matrix) })
dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()

mtx <- readMM(file.path(mdir, "matrix.mtx.gz"))
barcodes <- read.table(file.path(mdir, "barcodes.tsv.gz"), stringsAsFactors = F, sep = "\t", header = F)
features <- read.table(file.path(mdir, "features.tsv.gz"), stringsAsFactors = F, sep = "\t", header = F)
colnames(mtx) <- barcodes[, 1]
rownames(mtx) <- features[, 1]
gene_symbols <- features[, 2]
message("[", sample, "] loaded ", nrow(mtx), " genes x ", ncol(mtx), " spots")

initial_normal_spots <- FindNormalReference(counts = mtx, gene_symbols = gene_symbols, method = "pca")
celltype <- barcodes; rownames(celltype) <- barcodes[, 1]
celltype$annotation <- "tumor"
celltype[initial_normal_spots, "annotation"] <- "normal"
write.table(setNames(celltype, c("barcode", "type")), file.path(dir_path, paste0(sample, ".cell_types_binary.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("[", sample, "] FindNormalReference(pca): ", length(initial_normal_spots), " normal spots")

size <- read.table(file.path(raw, "sizes.cellranger-GRCh38-1.0.0.txt"), stringsAsFactors = F)
cyclegenes <- readRDS(file.path(raw, "cyclegenes.rds"))
bed <- read.table(file.path(raw, "hg38.genes.bed"), sep = "\t", header = T)
chrarm <- read.table(file.path(raw, "cytoarm_table_hg38.txt"), stringsAsFactors = F, sep = "\t", header = T)
chrarm <- chrarm[order(as.numeric(chrarm[, 1]), as.numeric(chrarm[, 3])), ]
bin_bed <- chrarm[, -2]
seg_table_filtered <- data.frame("chr" = bin_bed[, 1], "start" = as.numeric(bin_bed[, 2]),
                                 "end" = as.numeric(bin_bed[, 3]), "states" = 1,
                                 "length" = as.numeric(bin_bed[, 3]) - as.numeric(bin_bed[, 2]),
                                 "mean" = 0, "var" = 0, "Var1" = 1:nrow(bin_bed), "Freq" = 50000,
                                 "chrr" = paste0(bin_bed[, 1], ":", bin_bed[, 2]), stringsAsFactors = F)

Input_filtered <- FilterFeatures(mtx = mtx, barcodes = barcodes, features = features, cyclegenes = cyclegenes)
rm(mtx); rm(barcodes); rm(features)

message("[", sample, "] RunCovCluster start")
Cov_obj <- RunCovCluster(mtx = Input_filtered$mtx, barcodes = Input_filtered$barcodes,
                         features = Input_filtered$features, bed = bed, alpha = 5, beta = 5,
                         celltype0 = celltype, var_pt = 0.99, var_pt_ctrl = 0.99, include = "all",
                         alpha_source = "all", ctrl_region = NULL,
                         seg_table_filtered = seg_table_filtered, size = size,
                         dir_path = dir_path, breaks = 50, prep_mode = "intersect", est_cap = 2,
                         clust_mode = "all")
saveRDS(Cov_obj, file.path(dir_path, paste0(sample, ".Cov_obj.rds")))
message("[", sample, "] RunCovCluster done")

final_tn_res <- MalignantAssignment(Cov_obj, cutoff = 0.5)
zest <- Cov_obj$result_final$result$Zest
tn <- final_tn_res$final_assignment[names(zest)]
cells <- data.frame(BARCODE = names(zest), malignant = ifelse(tn == "Tumor", "T", "N"), cluster = zest)
write.table(cells, file.path(dir_path, paste0(sample, ".cells.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
message("[", sample, "] wrote cells.tsv: ", sum(cells$malignant == "T"), " tumor / ", sum(cells$malignant == "N"),
        " normal, ", length(unique(zest)), " clusters (", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), ")")
