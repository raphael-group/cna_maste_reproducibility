# CalicoST simulation h5ad -> 10x matrix directory, read by numbat and clonalscope.
#
# Config: the dataset's sim_dir and inputs_dir; references.ordered_genes (Cell Ranger GRCh38-2020-A
# genes/ordered_genes.tsv, the Ensembl IDs Clonalscope joins on); the calicost env (has anndata).
#
# Run from baselines/ (scripts/run_baseline.sh <method> calicost-sim runs it first):
#   snakemake -s scripts/common/prepare_sim_inputs.smk --configfile configs/datasets.yaml configs/metadata.yaml \
#       configs/<method>.calicost-sim.yaml --cores 8

include: "samples.smk"

rule all:
    input: [f"{MATRIX_DIR[s]}/matrix.mtx.gz" for s in SAMPLES]

rule h5ad_to_10x:
    input: f"{DATASET['sim_dir']}/{{s}}/filtered_feature_bc_matrix.h5ad"
    output: f"{DATASET['inputs_dir']}/{{s}}/filtered_feature_bc_matrix/matrix.mtx.gz"
    params:
        od = f"{DATASET['inputs_dir']}/{{s}}/filtered_feature_bc_matrix",
        py = f"{SOFTWARE['envs']['calicost']}/bin/python",
        genes = REFERENCES["ordered_genes"],
    log: f"{DATASET['inputs_dir']}/{{s}}/h5ad_to_10x.log"
    shell:
        "{params.py} -I {COMMON}/h5ad_to_10x.py {input} {params.genes} {params.od} > {log} 2>&1"
