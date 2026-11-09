#!/usr/bin/env bash
# Rerun CalicoST on one deposit simulation with the deposit's own configfiles (paths remapped only).
#
# Runpeng Luo (2026-11-02)
#
# Dependencies: PY, a python with pyyaml (common/config.sh); the calicost env with the CalicoST clone
#   installed editable ($CONFIGS/metadata.yaml software); GNU time.
#
# Usage:
#   bash scripts/baselines/calicost/repro.sh <sim> [out_root] [rects]
#           sim       # deposit simulation id, e.g. numcnas6.3_cnasize1e7_ploidy2_random0
#           out_root  # default out_root of $CONFIGS/calicost.calicost-sim.yaml
#           rects     # initializations to run, default "0 1 2 3 4"; a sanity check runs only the
#                     # deposit's chosen one (method_results/calicost/best_init.tsv)
#
# Inputs
#   $CONFIGS/datasets.yaml (deposit, sim_dir), $CONFIGS/metadata.yaml (references.calicost, software).
#   <sim_dir>/<sim>/: simulated counts, used as spaceranger_dir and snp_dir (read only).
#   <deposit>/nomixing_calicost_related/<sim>/configfile{0..4}: the deposit's CalicoST configs.
# Parameters
#   Every config value is kept from the deposit configfile except spaceranger_dir, snp_dir,
#   output_dir and the four resource files (geneticmap, hgtable, filtergenelist, filterregion),
#   which point to $CONFIGS/metadata.yaml references.calicost (CalicoST/GRCh38_resources, same content as
#   the files the deposit used).
# Outputs:
#   <out_root>/<sim>/configfile{0..4}, initial_phase.npz, parsed_inputs/, clone3_rectangle{r}_w1.0/
#   <out_root>/logs/<sim>.{prepare,run<r>}.{log,time.txt}: stdout/stderr and /usr/bin/time -v.
# Notes/References:
#   Deposit: https://doi.org/10.5281/zenodo.14175627 (Snakefile rules prepare_calicost_data, run_calicost).
#   Commands match the deposit Snakefile:
#     OMP_NUM_THREADS=1 python <calicost>/src/calicost/parse_input.py -c configfile0
#     OMP_NUM_THREADS=1 python <calicost>/src/calicost/calicost_main.py -c configfile<r>
#   The selected initializations run in parallel, as independent jobs in the deposit.
set -euo pipefail
source "$(dirname "$0")/../common/config.sh"

sim=$1
out_root=${2:-$(get 'c["out_root"]' "$CONFIGS/calicost.calicost-sim.yaml")}
dep=$(get 'c["datasets"]["calicost-sim"]["deposit"]' "$CONFIGS/datasets.yaml")
sim_dir=$(get 'c["datasets"]["calicost-sim"]["sim_dir"]' "$CONFIGS/datasets.yaml")
calicost=$(get 'c["software"]["calicost_repo"]' "$CONFIGS/metadata.yaml")
py=$(get 'c["software"]["envs"]["calicost"]' "$CONFIGS/metadata.yaml")/bin/python
# NB: key=path per CalicoST resource file
resources=$(get '" ".join(f"{k}={v}" for k, v in c["references"]["calicost"].items())' "$CONFIGS/metadata.yaml")
rects=${3:-0 1 2 3 4}
out=$out_root/$sim
logs=$out_root/logs
mkdir -p "$out" "$logs"

for r in 0 1 2 3 4; do
  OMP_NUM_THREADS=1 "$py" -I - "$dep/nomixing_calicost_related/$sim/configfile$r" "$out/configfile$r" \
    "$sim_dir/$sim" "$out" $resources <<'EOF'
import sys

from calicost.arg_parse import read_configuration_file, write_config_file

src, dst, indir, outdir, *resources = sys.argv[1:]
cfg = read_configuration_file(src)
cfg["spaceranger_dir"] = cfg["snp_dir"] = indir
cfg["output_dir"] = outdir
cfg.update(r.split("=", 1) for r in resources)
write_config_file(dst, cfg)
EOF
done

echo "$(date '+%F %T') $sim prepare"
/usr/bin/time -v -o "$logs/$sim.prepare.time.txt" env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  "$py" "$calicost/src/calicost/parse_input.py" -c "$out/configfile0" > "$logs/$sim.prepare.log" 2>&1

pids=()
for r in $rects; do
  /usr/bin/time -v -o "$logs/$sim.run$r.time.txt" env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    "$py" "$calicost/src/calicost/calicost_main.py" -c "$out/configfile$r" > "$logs/$sim.run$r.log" 2>&1 &
  pids+=($!)
done
rc=0
for p in "${pids[@]}"; do wait "$p" || rc=1; done
echo "$(date '+%F %T') $sim done rc=$rc"
exit $rc
