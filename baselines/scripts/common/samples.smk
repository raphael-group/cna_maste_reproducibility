# Dataset, samples, 10x matrix dir per sample and script dirs, shared by every snakefile.
#
# Config: configs/datasets.yaml, configs/metadata.yaml and configs/<method>.<case>.yaml, merged by
#   --configfile; `dataset` names the entry of datasets.
# Dataset with sim_dir (CalicoST simulations): every numcnas* folder under sim_dir is a sample; its
#   10x matrix is <inputs_dir>/<sample>/filtered_feature_bc_matrix (written by common/prepare_sim_inputs.smk).
# Otherwise (ovarian): `samples` of the method config; a sample of the dataset reads its 16um matrix.
# COMMON is scripts/common/ (shared scripts); SCRIPTS is the running snakefile's dir (scripts/<method>/).

import os

DATASET = config["datasets"][config["dataset"]]
REFERENCES = config["references"]
SOFTWARE = config["software"]
COMMON = os.path.abspath(str(workflow.current_basedir))
SCRIPTS = os.path.abspath(workflow.basedir)

if "sim_dir" in DATASET:
    SAMPLES = sorted(
        d for d in os.listdir(DATASET["sim_dir"])
        if d.startswith("numcnas") and os.path.isdir(os.path.join(DATASET["sim_dir"], d))
    )
    MATRIX_DIR = {s: f"{DATASET['inputs_dir']}/{s}/filtered_feature_bc_matrix" for s in SAMPLES}
else:
    SAMPLES = list(config["samples"])
    MATRIX_DIR = {s: DATASET["samples"][s]["matrix_016um"] for s in SAMPLES if s in DATASET["samples"]}

wildcard_constraints:
    s = r"[^/]+"
