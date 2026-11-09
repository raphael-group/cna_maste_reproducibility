# cna_maste_reproducibility

![Status](https://img.shields.io/badge/status-in%20development-yellow)
![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)

Code to reproduce the analysis and figures in the paper [Inference of allele-specific copy number aberrations in high-resolution spatial transcriptomics](todo). The `cna_maste` github repo is [here](https://github.com/raphael-group/cna_maste).

```sh
# NB: requirements: Linux x86_64; conda with mamba; GNU time and systemd user scopes with the memory controller
#     delegated (CalicoST's runs); pdflatex (TeX Live); read access to the private cna_maste repository
uname -sm  # NB: Linux x86_64
/usr/bin/time -v true
grep -w memory /sys/fs/cgroup/user.slice/user-$(id -u).slice/user@$(id -u).service/cgroup.controllers
command -v pdflatex && kpsewhich geometry.sty amsmath.sty amssymb.sty booktabs.sty

# NB: conda and mamba, if missing: Miniforge, https://github.com/conda-forge/miniforge
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh

git clone https://github.com/raphael-group/cna_maste_reproducibility.git
cd cna_maste_reproducibility
mamba env create -p /path/to/envs/cna_maste_reproducibility -f configs/envs/cna_maste_reproducibility.yaml
conda activate /path/to/envs/cna_maste_reproducibility

# NB: first set the /path/to/ values of configs/cna_maste_reproducibility.yaml (docs/reference.md, Configuration)
snakemake -s run_cna_maste_reproducibility.smk --cores 32 --resources mem_mb=262144

# NB: one step and the steps before it: setup, prepare, run_methods or analysis
snakemake -s run_cna_maste_reproducibility.smk --cores 32 --resources mem_mb=262144 prepare
```

## Docs

| Doc | Holds |
| --- | --- |
| [docs/reference.md](docs/reference.md) | Environments, configuration, steps (baselines, cnamaste, analysis), metrics, outputs (tables, figures) and tests |
| [docs/celltyping.md](docs/celltyping.md) | The ovarian sample's cell types: method, inputs and outputs |
