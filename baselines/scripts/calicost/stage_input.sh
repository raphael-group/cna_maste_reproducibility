#!/usr/bin/env bash
# Copy a CalicoST-formatted input dir and re-save unique_snp_ids.npy so numpy 1.x can load it.
#
# Runpeng Luo (2026-10-25)
#
# Dependencies: conda base python (numpy >= 2, scripts/common/config.sh) and the calicost env python (numpy 1.24,
#   CalicoST's pin; configs/metadata.yaml software.envs.calicost).
#
# Usage:
#   bash scripts/calicost/stage_input.sh <src_dir> <dst_dir>
#           src_dir  # formatted input (results/format_bb layout), left unchanged
#           dst_dir  # copy used as CalicoST spaceranger_dir and snp_dir
#
# Outputs:
#   dst_dir: byte copy of src_dir except unique_snp_ids.npy, rewritten as the same object array of
#     the same strings, pickled by numpy 1.24.
# Notes/References:
#   An object .npy is a pickle; one written by numpy 2 references numpy._core, which numpy 1.x lacks
#   (ModuleNotFoundError: No module named 'numpy._core' in CalicoST utils_IO.load_data).
set -euo pipefail
src=$1; dst=$2
source "$(dirname "$0")/../common/config.sh"
base_py=$CONDA/bin/python
cc_py=$(get 'c["software"]["envs"]["calicost"]' "$CONFIGS/metadata.yaml")/bin/python

mkdir -p "$dst"
cp -r "$src"/. "$dst"/
"$base_py" -c "import numpy as np, sys; ids = np.load(sys.argv[1], allow_pickle=True); open(sys.argv[2], 'w').write('\n'.join(map(str, ids)) + '\n')" \
    "$src/unique_snp_ids.npy" "$dst/unique_snp_ids.txt"
"$cc_py" -c "import numpy as np, sys; ids = np.array(open(sys.argv[1]).read().split(), dtype=object); np.save(sys.argv[2], ids, allow_pickle=True)" \
    "$dst/unique_snp_ids.txt" "$dst/unique_snp_ids.npy"
"$base_py" -c "import numpy as np, sys; a = np.load(sys.argv[1], allow_pickle=True); b = np.load(sys.argv[2], allow_pickle=True); assert a.dtype == b.dtype == object and list(a) == list(b); print('unique_snp_ids identical:', len(a))" \
    "$src/unique_snp_ids.npy" "$dst/unique_snp_ids.npy"
rm "$dst/unique_snp_ids.txt"
