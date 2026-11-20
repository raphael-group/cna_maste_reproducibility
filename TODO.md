# cna_maste_reproducibility TODO

## Moving the cnamaste pointers to its public repository

The workflow clones cnamaste from `REPOS.CNAMASTE_REPO` in `configs/cna_maste_reproducibility.yaml`: the private
michaelJwilson/cna_maste, pinned at `ddc060b`, the commit of its `cleanup3` branch. `README.md` already links
raphael-group/cna_maste. `git grep -n "michaelJwilson\|private\|ddc060b\|cleanup3\|42b96ff"` lists every line below.

Always:

- [ ] `configs/cna_maste_reproducibility.yaml`, `REPOS.CNAMASTE_REPO`: `url` the public repository; `rev` the public
      commit or tag, a commit as for the other repositories, not a branch; drop the `# NB: private` line.
- [ ] Drop the private-access notes: `README.md`'s prerequisite "read access to the private cna_maste repository";
      `docs/reference.md`'s IMPORTANT bullet on cloning with GitHub credentials; its `CNAMASTE_REPO` row (repository
      link and revision).
- [ ] `docs/reference.md`'s `CNAMASTE_ENV` row: the new revision.

Where the public history is not the private one (squashed or new), the commits named as sources do not exist there;
restate each as the public commit with the same code. The ported code itself does not change.

- [ ] `42b96ff`, the source of the port: the header of every `scripts/cnamaste_port/` module, `pyproject.toml`,
      `ruff.toml`, `docs/reference.md`.
- [ ] `ddc060b` and `cleanup3`, the source of the shipped cnamaste configs and of the solver page:
      `configs/cna_maste_reproducibility.yaml` comments, `run_cna_maste_reproducibility.smk`, `docs/reference.md`,
      `scripts/plots/simulation.py`, `scripts/baselines/cnamaste/write_binned_data.py`.

Where the new `rev` holds code other than `ddc060b`'s, rerun the setup and cnamaste rules, and check what the workflow
takes from cnamaste still holds:

- [ ] the commands `run_cnamaste`, `run_cnamaste_plot`, `run_sim_qa` in `CNAMASTE_ENV`, built from the clone's
      `uv.lock`;
- [ ] the solver library builds `infra/expansion/build_expansion.sh` and `infra/opengm/build_opengm.sh`;
- [ ] `cnamaste.io.get_labelled_counts`, which `scripts/baselines/cnamaste/write_binned_data.py` calls;
- [ ] the keys of `CNAMASTE.zenodo_sim.config` and `CNAMASTE.ovarian.config` in
      `configs/cna_maste_reproducibility.yaml`, which `run_cnamaste` reads over its
      `python/cnamaste/config/config_defaults.yaml`.

Unaffected: the analysis, the tests and CI import no cnamaste (`tests/test_no_cnamaste_imports.py`) and never clone
it. `clone_repo` checks an existing clone's HEAD against `rev`, not its remote, so a clone of the private repository
passes while it holds that commit.
