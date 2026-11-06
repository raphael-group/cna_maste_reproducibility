#!/usr/bin/env python
"""Per-simulation runtime and memory of the benchmarked numbat and clonalscope runs, plus numbat-zenodo runtime.

Runpeng Luo (2026-11-01)

Dependencies: pandas (calicost_env).

Usage:
    python -I sim_runtime.py <sim_dir> <bench_out> <deposit_numbat_dir> <out_tsv>
            sim_dir    # simulated_data_related/ (one subdirectory per simulation)
            bench_out  # benchmark output root (numbat/, clonalscope/)
            deposit_numbat_dir  # the deposit's numbat_related/ (published numbat runs)
            out_tsv    # one row per (simulation_id, method)

Inputs:
    numbat: <bench_out>/numbat/bench/run_numbat.<sim>.max_entropy_08.tsv
    clonalscope: <bench_out>/clonalscope/bench/run_tutorial.<sim>.visium_tutorial.tsv
    numbat-zenodo: <deposit_numbat_dir>/<sim>/outs/log.txt (numbat run log)
Outputs:
    out_tsv: simulation_id, method, Runtime [mins], Memory [GB].
Notes/References:
    Snakemake benchmark files: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#benchmark-rules
    Runtime is wall-clock seconds (s) / 60. Memory is max_pss / 1024 (MB -> GB); PSS counts
    pages shared by forked workers once, unlike max_rss.
    Only the method's main run job is counted, not shared input preprocessing.
    A simulation without a benchmark file (job failed) gets NA.
    numbat-zenodo runtime = last minus first [YYYY-mm-dd HH:MM:SS] timestamp in log.txt, which
    has one-second resolution and covers run_numbat() only; its memory is NA because the log
    reports R heap usage at checkpoints, not peak process memory.
"""

import os
import re
import sys
from datetime import datetime, timezone

import pandas as pd

BENCH = {
    "numbat": "numbat/bench/run_numbat.{sim}.max_entropy_08.tsv",
    "clonalscope": "clonalscope/bench/run_tutorial.{sim}.visium_tutorial.tsv",
}

TS = re.compile(r"\[(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\]")


def log_minutes(log: str) -> float:
    """Minutes between the first and last timestamp of a numbat log.txt.

    Args:
        log: numbat log.txt path.

    Returns:
        Minutes, or NaN if the log is missing or has fewer than two timestamps.
    """
    if not os.path.exists(log):
        return float("nan")
    with open(log) as fh:
        ts = TS.findall(fh.read())
    if len(ts) < 2:
        return float("nan")
    t0, t1 = (datetime.strptime(t, "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc) for t in (ts[0], ts[-1]))
    return (t1 - t0).total_seconds() / 60


def main(sim_dir: str, out: str, dep: str, out_tsv: str) -> None:
    """Collect one row per simulation and method.

    Args:
        sim_dir: simulated_data_related directory.
        out: Benchmark output root.
        dep: Deposit numbat_related directory.
        out_tsv: Output table.
    """
    sims = sorted(s for s in os.listdir(sim_dir) if os.path.exists(os.path.join(sim_dir, s, "truth_clone_labels.tsv")))
    rows = []
    for sim in sims:
        for method, pat in BENCH.items():
            f = os.path.join(out, pat.format(sim=sim))
            rt = mem = float("nan")
            if os.path.exists(f):
                b = pd.read_csv(f, sep="\t").iloc[-1]
                rt, mem = b["s"] / 60, b["max_pss"] / 1024
            rows.append(
                {
                    "simulation_id": sim,
                    "method": method,
                    "Runtime [mins]": round(rt, 2),
                    "Memory [GB]": round(mem, 3),
                }
            )
        rows.append(
            {
                "simulation_id": sim,
                "method": "numbat-zenodo",
                "Runtime [mins]": round(log_minutes(f"{dep}/{sim}/outs/log.txt"), 2),
                "Memory [GB]": float("nan"),
            }
        )
    res = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    res.to_csv(out_tsv, sep="\t", index=False, na_rep="NA")
    print(res.groupby("method")[["Runtime [mins]", "Memory [GB]"]].describe().T)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit("usage: sim_runtime.py <sim_dir> <bench_out> <deposit_numbat_dir> <out_tsv>")
    main(*sys.argv[1:])
