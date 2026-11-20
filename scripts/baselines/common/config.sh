# Shell access to the baseline configs, sourced by CalicoST's shell scripts under scripts/baselines/calicost/.
#
# Runpeng Luo (2026-11-03)
#
# Dependencies: a python with pyyaml (PY, default conda's base python).
#
# Usage:
#   CONFIGS=<dir> [CONDA=<conda root>] [PY=<python>] [LOGS=<dir>] source common/config.sh    # from scripts/baselines/
#   get '<python expression of c>' <yaml> [yaml ...]
#           expression  # evaluated on c, the yaml files merged by top-level key, and printed
#           yaml        # e.g. "$CONFIGS/metadata.yaml"
#
# Inputs
#   CONFIGS: the baseline configs, datasets.yaml, metadata.yaml and calicost.<case>.yaml.
# Outputs:
#   CONFIGS, CONDA (conda root), PY and LOGS (run logs, created), exported; get: prints the expression.
# Notes/References:
#   run_cna_maste_reproducibility.smk writes CONFIGS (<OUTDIR>/configs/baselines) from
#   configs/cna_maste_reproducibility.yaml and sets CONDA, PY (its root environment's python, with pyyaml and
#   numpy >= 2) and LOGS (<OUTDIR>/logs/baselines). PY reads every other path from CONFIGS.
CONFIGS=$(cd "${CONFIGS:?set CONFIGS to <OUTDIR>/configs/baselines (run_cna_maste_reproducibility.smk)}" && pwd)
CONDA=${CONDA:-$(conda info --base)}
PY=${PY:-$CONDA/bin/python}
LOGS=${LOGS:-$CONFIGS/../../logs/baselines}
mkdir -p "$LOGS"
LOGS=$(cd "$LOGS" && pwd)
export CONFIGS CONDA PY LOGS

get() {
    "$PY" -I -c '
import sys, yaml
c = {}
for f in sys.argv[2:]:
    c.update(yaml.safe_load(open(f)))
print(eval(sys.argv[1]))
' "$@"
}
