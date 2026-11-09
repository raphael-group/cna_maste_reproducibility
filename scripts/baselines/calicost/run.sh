#!/usr/bin/env bash
# Run CalicoST calicost_main.py on one configuration under a hard memory cap, with resource logging.
#
# Runpeng Luo (2026-10-25)
#
# Dependencies: the calicost env and CalicoST clone ($CONFIGS/metadata.yaml software); systemd user scopes
#   with the memory controller delegated (MemoryMax is enforced by the cgroup OOM killer).
#
# Usage:
#   bash scripts/baselines/calicost/run.sh <tag> <configuration_cna> [mem_max]
#           tag                # log name, e.g. H1-N3KXRGQ_016um_msu50
#           configuration_cna  # CalicoST configuration (calicost/ovarian_config.py write)
#           mem_max            # cgroup MemoryMax; default half of MemTotal, in GB (e.g. 503G)
#
# Inputs
#   configuration_cna: CalicoST configuration_cna with paths filled in.
# Outputs:
#   <output_dir from config>/: CalicoST outputs (created if missing).
#   $LOGS/calicost_<tag>.status: START / END lines with exit code; a SIGKILL (rc 137) under the cap is marked as OOM.
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
mem_max=${3:-$(awk '/MemTotal/ {printf "%dG", $2/1048576/2}' /proc/meminfo)}
calicost=$(get 'c["software"]["calicost_repo"]' "$CONFIGS/metadata.yaml")
env=$(get 'c["software"]["envs"]["calicost"]' "$CONFIGS/metadata.yaml")
outdir=$(awk -F' : ' '$1 == "output_dir" {print $2}' "$cfg")
mkdir -p "$outdir"

name=calicost_${tag}
echo "$(date '+%F %T') START run MemoryMax=${mem_max}" >> "$LOGS/${name}.status"
bash common/resource_sampler.sh $$ "${name}" run >> "$LOGS/${name}.resources.tsv" &
sampler=$!
set +e
systemd-run --user --scope -q -p MemoryMax="${mem_max}" -p MemorySwapMax=0 \
    /usr/bin/time -v -o "$LOGS/${name}.time.txt" \
    env OMP_NUM_THREADS=1 "${env}/bin/python" "${calicost}/src/calicost/calicost_main.py" -c "$cfg" \
    > "$LOGS/${name}.log" 2>&1
rc=$?
set -e
kill "$sampler" 2>/dev/null || true
note=""
if [[ $rc -eq 137 ]] || grep -q "terminated by signal 9" "$LOGS/${name}.time.txt" 2>/dev/null; then
    note=" (killed: OOM at MemoryMax=${mem_max})"
fi
echo "$(date '+%F %T') END run rc=${rc}${note}" >> "$LOGS/${name}.status"
exit $rc
