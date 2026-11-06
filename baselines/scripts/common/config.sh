# Shell access to configs/, sourced by the shell scripts under scripts/.
#
# Runpeng Luo (2026-11-03)
#
# Dependencies: conda base python with pyyaml.
#
# Usage:
#   source scripts/common/config.sh    # from baselines/
#   get '<python expression of c>' <yaml> [yaml ...]
#           expression  # evaluated on c, the yaml files merged by top-level key, and printed
#           yaml        # e.g. "$CONFIGS/metadata.yaml"
#
# Outputs:
#   CONDA: conda root; CONFIGS: absolute path of configs/; get: prints the expression.
# Notes/References:
#   CONDA is the one path read before YAML can be parsed; its base python reads every other path from
#   configs/. Set CONDA in the environment to override it.
CONDA=${CONDA:-/n/fs/ragr-research/users/rl6004/miniconda3}
CONFIGS=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../configs" && pwd)

get() {
    "$CONDA/bin/python" -I -c '
import sys, yaml
c = {}
for f in sys.argv[2:]:
    c.update(yaml.safe_load(open(f)))
print(eval(sys.argv[1]))
' "$@"
}
