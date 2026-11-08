#!/usr/bin/env bash
# Run one benchmarked method on one case: its Snakemake stages in the method's conda env, with status and
# per-minute resource logs; CalicoST runs its own runner.
#
# Runpeng Luo (2026-11-03)
#
# Dependencies: conda (scripts/common/config.sh), snakemake in each method env (configs/metadata.yaml software.envs).
#
# Usage:
#   [CORES=<n>] bash scripts/run_baseline.sh <method> <case> [snakemake args, e.g. -n]
#           method  # clonalscope | numbat | calicost
#           case    # ovarian | calicost-sim (configs/<method>.<case>.yaml)
#
# Inputs
#   configs/datasets.yaml, configs/metadata.yaml, configs/<method>.<case>.yaml, merged by --configfile.
# Outputs:
#   logs/<method>_<case>.<stage>.log: snakemake log per stage.
#   logs/<method>_<case>.status: START / END line per stage with exit code.
#   logs/<method>_<case>.resources.tsv: per-minute RSS / CPU of the stage's process tree.
#   Method outputs under the config's outdir (README.md).
# Notes/References:
#   Stages run in order and stop at the first failure; a stage is scripts/{<method>,common}/<stage>.smk, logged by
#   its file name. CORES sets snakemake --cores (default nproc). --rerun-triggers mtime keeps finished
#   outputs when only the shell command text changed. calicost runs scripts/calicost/ovarian.sh or
#   scripts/calicost/sim.sh, which take no snakemake args. The common/prepare_sim_inputs stage prepares the
#   simulations' inputs of <method> alone (--config methods=<method>).
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/common/config.sh

method=$1; run_case=$2; shift 2
tag=${method}_${run_case}
case $tag in
    numbat_ovarian) stages=(numbat/preprocess_ovarian numbat/run numbat/postprocess) ;;
    numbat_calicost-sim) stages=(common/prepare_sim_inputs numbat/run numbat/postprocess) ;;
    clonalscope_ovarian) stages=(clonalscope/run clonalscope/decode_states) ;;
    clonalscope_calicost-sim) stages=(common/prepare_sim_inputs clonalscope/run clonalscope/decode_states) ;;
    calicost_ovarian) exec bash scripts/calicost/ovarian.sh ;;
    calicost_calicost-sim) exec bash scripts/calicost/sim.sh ;;
    *) echo "no case ${run_case} for method ${method}" >&2; exit 2 ;;
esac
configs=("$CONFIGS/datasets.yaml" "$CONFIGS/metadata.yaml" "$CONFIGS/${method}.${run_case}.yaml")

source "$CONDA/etc/profile.d/conda.sh"
conda activate "$(get "c['software']['envs']['${method}']" "$CONFIGS/metadata.yaml")"
mkdir -p logs
for smk in "${stages[@]}"; do
    stage=$(basename "$smk")
    echo "$(date '+%F %T') START ${stage}" >> logs/${tag}.status
    bash scripts/common/resource_sampler.sh $$ "${tag}" "${stage}" >> logs/${tag}.resources.tsv &
    sampler=$!
    only=()
    if [[ $stage == prepare_sim_inputs ]]; then only=(--config "methods=${method}"); fi
    set +e
    snakemake -s "scripts/${smk}.smk" ${only[@]+"${only[@]}"} --configfile "${configs[@]}" --cores "${CORES:-$(nproc)}" -p \
        --rerun-incomplete --rerun-triggers mtime "$@" > logs/${tag}.${stage}.log 2>&1
    rc=$?
    set -e
    kill "$sampler" 2>/dev/null || true
    echo "$(date '+%F %T') END ${stage} rc=${rc}" >> logs/${tag}.status
    if (( rc != 0 )); then exit $rc; fi
done
