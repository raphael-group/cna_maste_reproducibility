# cna_maste_reproducibility TODO

## Moving the cnamaste pointers to its public repository - done

`REPOS.CNAMASTE_REPO` in `configs/cna_maste_reproducibility.yaml` clones raphael-group/cna_maste at `v0.1.0`
(`1ab2ee5`), whose tree is that of `ddc060b`, the private `cleanup3` commit the workflow pinned before.

- [x] `REPOS.CNAMASTE_REPO`: `url` the public repository, `rev` `v0.1.0`'s commit; the `# NB: private` line dropped.
- [x] The private-access notes dropped from `README.md` and `docs/reference.md`; its `CNAMASTE_ENV` and
      `CNAMASTE_REPO` rows name `v0.1.0`.
- [x] `cleanup3` and `ddc060b` restated as `v0.1.0`, which holds the same tree: the configs' comments,
      `run_cna_maste_reproducibility.smk`, `docs/reference.md`, `scripts/plots/simulation.py`,
      `scripts/baselines/cnamaste/write_binned_data.py`.
- [x] `42b96ff` kept as the port's source: `v0.1.0` changes every module the port copies, so no public commit holds
      that code; `docs/reference.md` names it a commit of cnamaste's development repository, before `v0.1.0`.
- [x] The method code is unchanged: with `ddc060b`'s tree, `v0.1.0` has the same `run_cnamaste`, `run_cnamaste_plot`,
      `run_sim_qa`, solver library builds, `cnamaste.io.get_labelled_counts` and config keys.

## A run directory set up before the move

- [ ] Replace the clone at its `CNAMASTE_REPO` with a clone of `v0.1.0`: the private clone fails `clone_repo`'s check.
      The new setup markers are then newer than the cnamaste runs, which Snakemake's default `mtime` trigger reruns;
      `-n` lists the jobs first, and `--rerun-triggers` without `mtime` leaves them.
