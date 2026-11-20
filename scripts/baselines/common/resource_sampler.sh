#!/usr/bin/env bash
# Print one line per minute: total RSS (GB), %CPU and process count of a process tree.
#
# Runpeng Luo (2026-10-24)
#
# Usage:
#   bash scripts/baselines/common/resource_sampler.sh <root_pid> <method> <stage> [interval_s]
#           root_pid    # parent of the tree to sum over
#           method      # label column
#           stage       # label column
#           interval_s  # sampling interval, default 60
#
# Outputs:
#   stdout: time  method  stage  rss_gb  cpu_pct  nproc  sys_avail_gb  top_cmd
#     top_cmd is the command of the largest-RSS process in the tree.
root=$1; method=$2; stage=$3; dt=${4:-60}
while kill -0 "$root" 2>/dev/null; do
    pids=$(pstree -p "$root" | grep -o '([0-9]\+)' | tr -d '()' | sort -u | paste -sd,)
    read -r rss cpu n top < <(ps -o rss=,pcpu=,args= -p "$pids" 2>/dev/null | awk '
        {r+=$1; c+=$2; n++; if ($1>m) {m=$1; $1=""; $2=""; t=substr($0,3,80)}}
        END {gsub(/[ \t]+/, "_", t); printf "%.2f %.0f %d %s\n", r/1048576, c, n, t}')
    avail=$(awk '/MemAvailable/ {printf "%.0f", $2/1048576}' /proc/meminfo)
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$(date '+%F %T')" "$method" "$stage" "$rss" "$cpu" "$n" "$avail" "$top"
    sleep "$dt"
done
