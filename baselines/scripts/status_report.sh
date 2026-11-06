#!/usr/bin/env bash
# One-screen status of every benchmark run: stage state, elapsed time, memory, progress markers.
#
# Runpeng Luo (2026-10-25)
#
# Usage:
#   bash scripts/status_report.sh
#
# Inputs
#   logs/<tag>.status, logs/<tag>.resources.tsv: written by scripts/run_baseline.sh / scripts/calicost/run.sh.
#   numbat run_numbat/log.txt, CalicoST *_nstates*_*.npz: progress markers in the ovarian output dirs
#   (configs/numbat.ovarian.yaml, configs/calicost.ovarian.yaml).
# Outputs:
#   stdout: one block per run; free system memory at the end.
cd "$(dirname "$0")/.."
source scripts/common/config.sh
NUMBAT_RUN=$(get '"/".join([c["outdir"], c["samples"][0], c["run_label"], "run_numbat"])' "$CONFIGS/numbat.ovarian.yaml")
CALICOST=$(get 'c["outdir"]' "$CONFIGS/calicost.ovarian.yaml")
PY=$(get 'c["software"]["envs"]["calicost"]' "$CONFIGS/metadata.yaml")/bin/python
now=$(date +%s)
echo "status at $(date '+%F %T')"

hm() { local s=$1; printf "%dd%02dh%02dm" $((s/86400)) $((s%86400/3600)) $((s%3600/60)); }

for st in logs/*.status; do
    tag=$(basename "$st" .status)
    last=$(tail -1 "$st")
    first=$(head -1 "$st" | awk '{print $1" "$2}')
    if [[ $last == *END* ]]; then
        state="finished: ${last#* * }"
        el=$(( $(date -d "$(echo "$last" | awk '{print $1" "$2}')" +%s) - $(date -d "$first" +%s) ))
    else
        state="running: ${last#* * }"
        el=$(( now - $(date -d "$first" +%s) ))
    fi
    mem=$(tail -1 logs/${tag}.resources.tsv 2>/dev/null | awk -F'\t' '{print $4" GiB RSS (sampled "$1")"}')
    peak=$(awk -F'\t' '{if ($4 > x) x = $4} END {printf "%.1f", x}' logs/${tag}.resources.tsv 2>/dev/null)
    echo "- ${tag}: ${state} | elapsed $(hm $el) | now ${mem:-n/a} | peak ${peak} GiB"
    case $tag in
        numbat_ovarian)
            l=$(grep INFO $NUMBAT_RUN/log.txt 2>/dev/null | grep -v -E 'Mem used|\] $' | tail -1 | cut -c1-110)
            echo "    last numbat step: ${l}";;
        calicost_*)
            run=${tag#calicost_}
            for f in $(ls -tr $CALICOST/$run/clone3_rectangle0_w1.0/*_nstates*_*.npz 2>/dev/null); do
                n=$($PY -c "import numpy as np,sys; d=np.load(sys.argv[1],allow_pickle=True); print(int(d['num_iterations']) if 'num_iterations' in d else 'final')" "$f" 2>/dev/null)
                echo "    $(basename $f): rounds=${n}"
            done
            ls $CALICOST/$run/clone3_rectangle0_w1.0/clone_labels.tsv >/dev/null 2>&1 && echo "    clone_labels.tsv written";;
        clonalscope*)
            s=H1-N3KXRGQ; [[ $tag == *032um ]] && s=H1-N3KXRGQ_032um
            done_rules=$(grep -h 'Finished jobid' logs/${tag}.*.log 2>/dev/null | sed 's/.*Rule: //; s/)//' | paste -sd, )
            echo "    finished rules: ${done_rules:-none}";;
    esac
done
echo "system MemAvailable: $(awk '/MemAvailable/ {printf "%.0f GiB", $2/1048576}' /proc/meminfo)"
