# cna_maste_reproducibility

[![CI](https://github.com/raphael-group/cna_maste_reproducibility/actions/workflows/ci.yml/badge.svg)](https://github.com/raphael-group/cna_maste_reproducibility/actions/workflows/ci.yml)
![Status](https://img.shields.io/badge/status-in%20development-yellow)
![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)

Code to reproduce the analysis and figures in the manuscript [Inference of allele-specific copy number aberrations in high-resolution spatial transcriptomics](under_preparation). The `cna_maste` github repo is [here](https://github.com/raphael-group/cna_maste).

## Environment
```sh
# NB: requirements: Linux x86_64; conda with mamba; GNU time and timeout, and systemd user scopes with the memory
#     controller delegated (the baseline runs' limits); pdflatex (TeX Live); read access to the private cna_maste repository
uname -sm  # NB: Linux x86_64
/usr/bin/time -v true && timeout 1 true
grep -w memory /sys/fs/cgroup/user.slice/user-$(id -u).slice/user@$(id -u).service/cgroup.controllers
command -v pdflatex && kpsewhich geometry.sty amsmath.sty amssymb.sty booktabs.sty

# NB: conda and mamba, if missing: Miniforge, https://github.com/conda-forge/miniforge
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

## Configuration & Run
Copy the "Set" sections of `configs/cna_maste_reproducibility.yaml` into your own config, set its paths, and pass it
with `--configfile`. Its keys override the default, which stays as shipped, so the pinned "Fixed" settings keep coming
from this repository.

```sh
git clone https://github.com/raphael-group/cna_maste_reproducibility.git
cd cna_maste_reproducibility
mamba env create -p /path/to/envs/cna_maste_reproducibility -f configs/envs/cna_maste_reproducibility.yaml
conda activate /path/to/envs/cna_maste_reproducibility

# NB: your config, the default's "Set" sections; then set its /path/to/ values, CORES and MEM_GB (docs/reference.md)
awk '/^# Fixed/{exit} {print}' configs/cna_maste_reproducibility.yaml > /path/to/my_config.yaml

# NB: --cores and --resources mem_mb match your CORES and MEM_GB * 1024
snakemake -s run_cna_maste_reproducibility.smk --configfile /path/to/my_config.yaml --cores 32 --resources mem_mb=262144
```

## Docs

| Doc | Holds |
| --- | --- |
| [docs/reference.md](docs/reference.md) | Environments, configuration, steps (baselines, cnamaste, analysis), metrics, outputs (tables, figures) and tests |
| [docs/celltyping.md](docs/celltyping.md) | The ovarian sample's cell types: method, inputs and outputs |
