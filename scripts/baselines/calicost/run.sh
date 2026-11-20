#!/usr/bin/env bash
# Run CalicoST calicost_main.py on one configuration under a hard wall-time and memory limit, with resource logging.
#
# Runpeng Luo (2026-10-25)
#
# Dependencies: the calicost env and CalicoST clone ($CONFIGS/metadata.yaml software); common/limit.sh (systemd user
#   scopes with the memory controller delegated, GNU timeout); GNU time.
#
# Usage:
#   bash scripts/baselines/calicost/run.sh <tag> <configuration_cna> <max_walltime> <max_rss>
#           tag                # log name, e.g. H1-N3KXRGQ_016um_msu50
#           configuration_cna  # CalicoST configuration (calicost/ovarian_config.py write)
#           max_walltime       # wall-time limit, a GNU timeout duration, e.g. 24h (common/limit.sh)
#           max_rss            # memory limit, a systemd MemoryMax, e.g. 512G (common/limit.sh)
#
# Inputs
#   configuration_cna: CalicoST configuration_cna with paths filled in.
# Outputs:
#   <output_dir from config>/: CalicoST outputs (created if missing).
#   $LOGS/calicost_<tag>.status: START / END lines with the limits, exit code, wall time and the limit that stopped
#       the run, where one did (common/limit.sh).
#   $LOGS/calicost_<tag>.log: CalicoST stdout and stderr.
#   $LOGS/calicost_<tag>.resources.tsv: per-minute RSS / CPU of the run (common/resource_sampler.sh).
#   $LOGS/calicost_<tag>.time.txt: /usr/bin/time -v (wall time, max RSS).
# Notes/References:
#   Logs go to $LOGS (common/config.sh). Command from the CalicoST README:
#   OMP_NUM_THREADS=1 python <CalicoST directory>/src/calicost/calicost_main.py -c configuration_cna
#   https://github.com/raphael-group/CalicoST
set -euo pipefail
cfg=$(realpath "$2")
cd "$(dirname "$0")/.."
source common/config.sh

tag=$1
max_walltime=$3
max_rss=$4
calicost=$(get 'c["software"]["calicost_repo"]' "$CONFIGS/metadata.yaml")
env=$(get 'c["software"]["envs"]["calicost"]' "$CONFIGS/metadata.yaml")
outdir=$(awk -F' : ' '$1 == "output_dir" {print $2}' "$cfg")
mkdir -p "$outdir"

name=calicost_${tag}
bash common/resource_sampler.sh $$ "${name}" run >> "$LOGS/${name}.resources.tsv" &
sampler=$!
set +e
bash common/limit.sh "$max_walltime" "$max_rss" "$LOGS/${name}.status" \
    /usr/bin/time -v -o "$LOGS/${name}.time.txt" \
    env OMP_NUM_THREADS=1 "${env}/bin/python" "${calicost}/src/calicost/calicost_main.py" -c "$cfg" \
    > "$LOGS/${name}.log" 2>&1
rc=$?
set -e
kill "$sampler" 2>/dev/null || true
exit $rc
