#!/usr/bin/env bash
# Build Clonalscope's ref/ directory: 200kb bins, chromosome-arm seg_table, data-raw at the pinned commit.
#
# Runpeng Luo (2026-11-08)
#
# Dependencies:
#   bedtools, Rscript (base R), curl; defaults to the clonalscope_env binaries.
#
# Usage:
#   bash scripts/baselines/clonalscope/build_refs_clonalscope.sh <ref_dir> [sizes] [regions_bed]
#           ref_dir       # output directory, created if absent
#           sizes         # chrom<TAB>length (default: the repository's data/hg38.chrom.sizes)
#           regions_bed   # non-centromere/non-telomere keep-regions (default: the repository's data/hg38.regions.bed)
#
# Inputs
#   sizes, regions_bed: genome sizes and keep-regions, hg38; the repository's data/ copies equal
#       copytyping-benchmark/reference/hg38.sizes and hg38.regions.bed in their first three columns.
# Parameters
#   BIN_SIZE: 200000 bp bins.
#   CLONALSCOPE_SHA: Clonalscope commit for data-raw/, matching the installed package.
#   ENV: conda env prefix providing bedtools and Rscript (env var; default clonalscope_env).
# Outputs:
#   <ref_dir>/hg38.200kb.bins.bed: 200kb windows clipped to keep-regions, sorted in sizes-file order.
#   <ref_dir>/seg_table.arm.tsv: one row per autosome arm (chr start end states length mean var chrr), states=1.
#   <ref_dir>/clonalscope_data-raw_07dd0c5/: sizes.cellranger-GRCh38-1.0.0.txt, cyclegenes.rds,
#       hg38.genes.bed, cytoarm_table_hg38.txt from Clonalscope data-raw/.
# Notes/References:
#   Bins and arm table reproduce copytyping-benchmark workflows/clonalscope/scripts/build_ref_window_bed.sh
#   and workflows/clonalscope_nobulk/scripts/build_arm_segtable.R.
#   Arm segments per the Clonalscope no-WGS tutorial:
#   https://github.com/seasoncloud/Clonalscope/tree/main/samples/P5847/scRNA
#   The arm table warns "q-arm end != chrom size" for every autosome: the keep-regions stop short of the telomere.
set -euo pipefail

REF_DIR=${1:?usage: build_refs_clonalscope.sh <ref_dir> [sizes] [regions_bed]}
DATA=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../data" && pwd)
SIZES=${2:-$DATA/hg38.chrom.sizes}
REGIONS_BED=${3:-$DATA/hg38.regions.bed}
ENV=${ENV:-/n/fs/ragr-research/users/rl6004/envs/clonalscope_env}
BEDTOOLS=$ENV/bin/bedtools
RSCRIPT=$ENV/bin/Rscript
BIN_SIZE=200000
CLONALSCOPE_SHA=07dd0c576b6754d18dda86ec8905d8a9be5e77b0
DATA_RAW=$REF_DIR/clonalscope_data-raw_${CLONALSCOPE_SHA:0:7}

mkdir -p "$REF_DIR" "$DATA_RAW"

echo "[bins] ${BIN_SIZE} bp from $SIZES, clipped to $REGIONS_BED"
# NB: sizes-file order, not lexicographic; downstream bedtools -sorted needs it.
"$BEDTOOLS" makewindows -g "$SIZES" -w "$BIN_SIZE" \
    | "$BEDTOOLS" intersect -a - -b "$REGIONS_BED" \
    | "$BEDTOOLS" sort -g "$SIZES" -i - \
    > "$REF_DIR/hg38.200kb.bins.bed"

echo "[arm_segtable] from $REGIONS_BED"
"$RSCRIPT" -e '
args <- commandArgs(trailingOnly = TRUE)
regions_bed <- args[1]; sizes_file <- args[2]; out_tsv <- args[3]
regions <- read.table(regions_bed, sep = "\t", stringsAsFactors = FALSE)[, 1:3]
colnames(regions) <- c("chr", "start", "end")
regions <- regions[regions$chr %in% paste0("chr", 1:22), , drop = FALSE]
regions$chr <- as.integer(sub("^chr", "", regions$chr))
regions <- regions[order(regions$chr, regions$start), , drop = FALSE]
sizes <- read.table(sizes_file, stringsAsFactors = FALSE, col.names = c("chr", "length"))
chr_len <- setNames(sizes$length, sub("^chr", "", sizes$chr))
last_end <- tapply(regions$end, regions$chr, max)
mismatch <- which(as.numeric(last_end) != as.numeric(chr_len[names(last_end)]))
if (length(mismatch) > 0)
  message("[arm_segtable] WARNING: q-arm end != chrom size for chr ",
          paste(names(last_end)[mismatch], collapse = ", "))
seg_table <- data.frame(
  chr = regions$chr, start = regions$start, end = regions$end,
  states = 1, length = regions$end - regions$start, mean = 0, var = 0,
  chrr = paste0(regions$chr, ":", regions$start), stringsAsFactors = FALSE)
write.table(seg_table, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("[arm_segtable] wrote ", out_tsv, " (", nrow(seg_table), " arms over ",
        length(unique(seg_table$chr)), " autosomes)")
' "$REGIONS_BED" "$SIZES" "$REF_DIR/seg_table.arm.tsv"

echo "[data-raw] Clonalscope@${CLONALSCOPE_SHA}"
for f in sizes.cellranger-GRCh38-1.0.0.txt cyclegenes.rds hg38.genes.bed cytoarm_table_hg38.txt; do
    curl -fsSL -o "$DATA_RAW/$f" \
        "https://raw.githubusercontent.com/seasoncloud/Clonalscope/${CLONALSCOPE_SHA}/data-raw/$f"
done

echo "[done] $REF_DIR"