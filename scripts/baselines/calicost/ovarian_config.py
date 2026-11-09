#!/usr/bin/env python
"""One CalicoST ovarian run's configuration_cna, from the baseline configs.

Runpeng Luo (2026-11-03)

Dependencies: pyyaml.

Usage:
    python -I calicost/ovarian_config.py write <tag> <out_file>
            tag       # a key of runs in $CONFIGS/calicost.ovarian.yaml
            out_file  # configuration_cna to write

Inputs:
    $CONFIGS/datasets.yaml: UGP bb-<XXX>um-formatted dir of each sample and bin size (ugp_bb).
    $CONFIGS/metadata.yaml: CalicoST GRCh38_resources files (references.calicost).
    $CONFIGS/calicost.ovarian.yaml: staged_dir, outdir, inputs, configuration_cna, runs.
Outputs:
    out_file: "key : value" lines: spaceranger_dir and snp_dir (the run's staged input),
        output_dir (<outdir>/<tag>), the resource files, configuration_cna, then the run's own keys.
Notes/References:
    run_cna_maste_reproducibility.smk rule calicost_ovarian writes each run's configuration_cna with this.
    CalicoST reads configuration_cna line by line over its defaults, so key order does not matter and a
    key left out keeps its default: read_configuration_file in src/calicost/arg_parse.py,
    https://github.com/raphael-group/CalicoST/blob/c1abcae3e3657e01e547ee4529e3b9d039221453/src/calicost/arg_parse.py
"""

import os
import sys

import yaml

CONFIGS = os.environ["CONFIGS"]
CONFIG_FILES = ["datasets.yaml", "metadata.yaml", "calicost.ovarian.yaml"]


def get_config() -> dict:
    """Get the ovarian CalicoST config.

    Returns:
        $CONFIGS/datasets.yaml, $CONFIGS/metadata.yaml and $CONFIGS/calicost.ovarian.yaml, merged by top-level key.
    """
    config: dict = {}
    for name in CONFIG_FILES:
        with open(os.path.join(CONFIGS, name)) as f:
            config.update(yaml.safe_load(f))
    return config


def get_inputs(config: dict) -> dict[str, tuple[str, str]]:
    """Get the inputs to stage.

    Args:
        config: Merged config (get_config).

    Returns:
        Input name -> (UGP bb-<XXX>um-formatted source dir, staged dir).
    """
    samples = config["datasets"][config["dataset"]]["samples"]
    return {
        name: (samples[spec["sample"]]["ugp_bb"][spec["bin_um"]], f"{config['staged_dir']}/{name}")
        for name, spec in config["inputs"].items()
    }


def get_configuration(config: dict, tag: str) -> dict[str, str]:
    """Get one run's configuration_cna.

    Args:
        config: Merged config (get_config).
        tag: Run tag, a key of runs.

    Returns:
        CalicoST key -> value as written, in file order.
    """
    run = dict(config["runs"][tag])
    staged = get_inputs(config)[run.pop("input")][1]
    entries = {"spaceranger_dir": staged, "snp_dir": staged, "output_dir": f"{config['outdir']}/{tag}"}
    entries.update(config["references"]["calicost"])
    entries.update(config["configuration_cna"])
    entries.update(run)
    return {key: str(value) for key, value in entries.items()}


def main(args: list[str]) -> None:
    """Write one run's configuration_cna.

    Args:
        args: ["write", tag, out_file].
    """
    if len(args) != 3 or args[0] != "write":
        sys.exit("usage: calicost/ovarian_config.py write <tag> <out_file>")
    with open(args[2], "w") as f:
        f.writelines(f"{key} : {value}\n" for key, value in get_configuration(get_config(), args[1]).items())


if __name__ == "__main__":
    main(sys.argv[1:])
