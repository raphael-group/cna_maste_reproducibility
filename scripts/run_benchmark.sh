#!/usr/bin/env bash
# Run the benchmark on H1-N3KXRGQ into output/benchmark_results/, the previous run's results there deleted first.
#
# Author: Runpeng Luo (2026-11-04)
#
# Dependencies:
#   a cnamaste repository at $CNAMASTE_ROOT, its environment and `luo` dependency group
#   (`uv sync --locked --all-extras --group luo` there); uv; pdflatex; git, optional, for the commit logged.
#
# Usage:
#   CNAMASTE_ROOT=<cnamaste> bash <benchmark>/scripts/run_benchmark.sh [run_benchmark.py option ...]
#           option    # added to run_benchmark.py's command line, e.g. --cnamaste_normal_fraction 1=0.38
#
# Inputs
#   configs/benchmark.local.yaml, configs/config.local.yaml: run_benchmark.py's runs, references and genes.
#   input_data/real_data/: the sample's UGP bins, Space Ranger outputs, cell types and each method's run.
# Parameters
#   CNAMASTE_ROOT: environment variable, the cnamaste repository's root; `cnamaste` is imported from its python/.
# Outputs:
#   output/benchmark_results/: the run's tsvs/, text/, plots/ and pdf-pages/ (run_benchmark.py), replacing the
#       previous run's.
# Notes/References:
#   README.md, Run. The benchmark is the directory above this script, imported as `scripts.*`.
set -euo pipefail

: "${CNAMASTE_ROOT:?set CNAMASTE_ROOT to the root of a cnamaste repository, README.md Run}"
CNAMASTE_ROOT=$(cd "$CNAMASTE_ROOT" && pwd)
BENCHMARK=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export CNAMASTE_ROOT
export PYTHONPATH="$CNAMASTE_ROOT/python:$BENCHMARK${PYTHONPATH:+:$PYTHONPATH}"
COMMIT=$(git -C "$CNAMASTE_ROOT" describe --always --dirty 2>/dev/null || echo "no git commit")
echo "cnamaste: $CNAMASTE_ROOT ($COMMIT)"
O=$BENCHMARK/output/benchmark_results
B=$BENCHMARK/input_data/real_data
H=$B/output_10092026/H1-N3KXRGQ
CNAMASTE=$B/cnamaste/N3KXRGQ-ac88605f/run_ac88605f

rm -rf "$O"
mkdir -p "$O/tsvs/gene_level_cnp" "$O/plots" "$O/pdf-pages" "$O/text"
echo "Run directory: $O"

uv run --project "$CNAMASTE_ROOT" --group luo python "$BENCHMARK/scripts/run_benchmark.py" --sample H1-N3KXRGQ \
    --spaceranger_dir "$H/bb-032um" \
    --spots "$H/bb/MSR300/VISIUM3prime/barcodes.tsv.gz" \
    --bb_dir "16=$H/bb" --bb_dir "32=$H/bb-032um" \
    --celltype "$B/celltyping/H1-N3KXRGQ/celltyping-16um-all_spots-10x_graphclust_panglaodb_ulm/spot_annotation.tsv" \
    --cnamaste_labels "$CNAMASTE/clone_labels_int.tsv.gz" \
    --cnamaste_seglevel "$CNAMASTE/cnv_diploid_seglevel.tsv" \
    --calicost_dir "msu20=$B/calicost/H1-N3KXRGQ_016um_msu20" \
    --numbat_dir "$B/numbat/H1-N3KXRGQ-numbat-max_entropy_08" \
    --clonalscope_dir "$B/clonalscope/H1-N3KXRGQ_032um-visium_tutorial/postprocess" \
    --outdir "$O" "$@"
echo "Done: $O"
