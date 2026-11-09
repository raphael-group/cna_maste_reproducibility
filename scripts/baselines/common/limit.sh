#!/usr/bin/env bash
# Run a command under a hard wall-time and memory limit: GNU timeout inside a systemd user scope with MemoryMax.
#
# Runpeng Luo (2026-11-09)
#
# Dependencies: systemd user scopes with the memory controller delegated; GNU coreutils timeout.
#
# Usage:
#   bash scripts/baselines/common/limit.sh <max_walltime> <max_rss> <status> <command> [argument ...]
#           max_walltime  # wall-time limit, a GNU timeout duration: <N>[s|m|h|d]
#           max_rss       # memory limit, the scope's MemoryMax with swap off: <N>[K|M|G|T]
#           status        # file the START and END lines are appended to
#           command       # the run, with its arguments
#
# Outputs:
#   status: a START line with the limits; an END line with the exit code, the wall time, and the limit that stopped
#       the run, where one did.
#   Exit code: the command's; 124 where max_walltime stopped it (timeout's SIGTERM), 137 where it was killed (SIGKILL:
#       memory over max_rss, or a run still alive 10 minutes after max_walltime's SIGTERM).
# Notes/References:
#   The limits are BASELINE_LIMITS of configs/cna_maste_reproducibility.yaml; a stopped run fails its rule.
#   MemoryMax, MemorySwapMax: https://www.freedesktop.org/software/systemd/man/latest/systemd.resource-control.html
#   timeout signals the command's process group, its children included:
#   https://www.gnu.org/software/coreutils/manual/html_node/timeout-invocation.html
set -euo pipefail

max_walltime=$1
max_rss=$2
status=$3
shift 3

seconds() {
    local n=${1%[smhd]}
    case $1 in
    *d) echo $((n * 86400)) ;;
    *h) echo $((n * 3600)) ;;
    *m) echo $((n * 60)) ;;
    *) echo "$n" ;;
    esac
}

limit=$(seconds "$max_walltime")
mkdir -p "$(dirname "$status")"
echo "$(date '+%F %T') START max_walltime=${max_walltime} max_rss=${max_rss}" >> "$status"
start=$(date +%s)
set +e
systemd-run --user --scope -q -p MemoryMax="$max_rss" -p MemorySwapMax=0 \
    timeout --kill-after=10m "$max_walltime" "$@"
rc=$?
set -e
elapsed=$(($(date +%s) - start))
note=""
if [[ $rc -ne 0 ]] && ((elapsed >= limit)); then
    note=" (stopped: max_walltime ${max_walltime})"
elif [[ $rc -eq 137 ]]; then
    note=" (killed: max_rss ${max_rss})"
fi
echo "$(date '+%F %T') END rc=${rc} elapsed=${elapsed}s${note}" >> "$status"
exit $rc
