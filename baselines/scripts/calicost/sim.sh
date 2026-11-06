#!/usr/bin/env bash
# CalicoST on the CalicoST simulations with the deposit's configfiles, several simulations in parallel.
#
# Runpeng Luo (2026-11-02)
#
# Dependencies: conda base python with pyyaml (scripts/common/config.sh); the calicost env and CalicoST clone
#   (configs/metadata.yaml software); GNU xargs, GNU time.
#
# Usage:
#   bash scripts/calicost/sim.sh    # or: bash scripts/run_baseline.sh calicost calicost-sim
#
# Inputs
#   configs/datasets.yaml (deposit, sim_dir), configs/calicost.calicost-sim.yaml.
# Outputs:
#   <out_root>/<sim>/: configfile{0..4}, parsed_inputs/, clone3_rectangle<r>_w1.0/.
#   <out_root>/logs/<sim>.{prepare,run<r>}.{log,time.txt}, <out_root>/logs/driver.<sim>.log.
# Notes/References:
#   A simulation whose every requested run<r>.time.txt reports "Exit status: 0" is skipped.
#   Deposit Snakefile rules prepare_calicost_data, run_calicost: https://doi.org/10.5281/zenodo.14175627
set -euo pipefail
cd "$(dirname "$0")/../.."
source scripts/common/config.sh
configs=("$CONFIGS/datasets.yaml" "$CONFIGS/calicost.calicost-sim.yaml")

sim_dir=$(get 'c["datasets"][c["dataset"]]["sim_dir"]' "${configs[@]}")
out_root=$(get 'c["out_root"]' "${configs[@]}")
rects=$(get '" ".join(map(str, c["rects"]))' "${configs[@]}")
max_jobs=$(get 'c["max_jobs"]' "${configs[@]}")
sims=$(get '"" if c["sims"] == "all" else " ".join(c["sims"])' "${configs[@]}")
if [[ -z $sims ]]; then
    sims=$(find "$sim_dir" -mindepth 1 -maxdepth 1 -type d -name 'numcnas*' -printf '%f\n' | sort)
fi

todo=()
for s in $sims; do
    ok=1
    for r in $rects; do
        grep -q "Exit status: 0" "$out_root/logs/$s.run$r.time.txt" 2>/dev/null || ok=0
    done
    if (( ok )); then echo "done: $s"; else todo+=("$s"); fi
done
echo "$(date '+%F %T') ${#todo[@]} simulations to run, $max_jobs at a time"
mkdir -p "$out_root/logs"
printf '%s\n' "${todo[@]}" | xargs -r -P "$max_jobs" -I{} \
    bash -c 'bash "$0" {} "$1" "$2" > "$1/logs/driver.{}.log" 2>&1' \
    scripts/calicost/repro.sh "$out_root" "$rects"
