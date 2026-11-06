#!/usr/bin/env bash
# CalicoST on Visium HD 3prime ovarian H1-N3KXRGQ: stage inputs, then each configured run in order.
#
# Runpeng Luo (2026-11-02)
#
# Dependencies: conda base python with pyyaml (scripts/common/config.sh); the calicost env and CalicoST clone
#   (configs/metadata.yaml software); systemd user scopes (scripts/calicost/run.sh).
#
# Usage:
#   bash scripts/calicost/ovarian.sh    # or: bash scripts/run_baseline.sh calicost ovarian
#
# Inputs
#   configs/datasets.yaml, configs/metadata.yaml, configs/calicost.ovarian.yaml.
# Outputs:
#   <staged_dir>/<input>/: staged CalicoST input (skipped if present).
#   <outdir>/<tag>/: CalicoST outputs.
#   logs/calicost_<tag>.configuration_cna: the run's configuration (scripts/calicost/ovarian_config.py).
#   logs/calicost_<tag>.{status,log,resources.tsv,time.txt}.
# Notes/References:
#   A run whose logs/calicost_<tag>.status ends with "END run rc=0" is skipped.
#   https://github.com/raphael-group/CalicoST
set -euo pipefail
cd "$(dirname "$0")/../.."
source scripts/common/config.sh
plan() { "$CONDA/bin/python" -I scripts/calicost/ovarian_config.py "$@"; }

while read -r name src dst; do
    if [[ -d $dst ]]; then echo "staged: $name"; continue; fi
    bash scripts/calicost/stage_input.sh "$src" "$dst"
done < <(plan inputs)

mem_max=$(get 'c.get("mem_max") or ""' "$CONFIGS/calicost.ovarian.yaml")
mkdir -p logs
for tag in $(plan runs); do
    if tail -n 1 "logs/calicost_${tag}.status" 2>/dev/null | grep -q "END run rc=0$"; then
        echo "done: $tag"; continue
    fi
    plan write "$tag" "logs/calicost_${tag}.configuration_cna"
    echo "$(date '+%F %T') run $tag"
    bash scripts/calicost/run.sh "$tag" "logs/calicost_${tag}.configuration_cna" $mem_max
done
