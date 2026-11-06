#!/usr/bin/env python
"""Runtime and peak memory of every benchmarked method run on Visium HD H1-N3KXRGQ.

Runpeng Luo (2026-11-02)

Dependencies: pandas.

Usage:
    python -I vishd_runtime.py <bench_out> <log_dir> <out_tsv>
            bench_out  # Visium HD benchmark output root (clonalscope/, numbat/)
            log_dir    # benchmark/logs (CalicoST /usr/bin/time -v files)
            out_tsv    # one row per method run

Inputs:
    numbat: <bench_out>/numbat/bench/run_numbat.H1-N3KXRGQ[.<run_label>].tsv
    clonalscope visium_tutorial: <bench_out>/clonalscope/bench/run_tutorial.<s>.visium_tutorial.tsv
    calicost: <log_dir>/calicost_<tag>.time.txt
Outputs:
    out_tsv: method, setting, resolution [um], Runtime [mins], Memory [GB], source.
Notes/References:
    Snakemake benchmark files: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#benchmark-rules
    Snakemake runs: runtime = s / 60, memory = max_pss / 1024 (MB -> GB), as scripts/sim_runtime.py.
    PSS counts pages shared by forked workers once, unlike max_rss.
    CalicoST: runtime = "Elapsed (wall clock) time", memory = "Maximum resident set size" (kB),
    from /usr/bin/time -v; CalicoST runs as one process (OMP_NUM_THREADS=1), so RSS is used.
    A run without a timing record (not finished, or not run by this pipeline) gets NA.
"""

import os
import re
import sys

import pandas as pd

S16, S32 = "H1-N3KXRGQ", "H1-N3KXRGQ_032um"
CALICOST = [
    ("msu10", 16, "H1-N3KXRGQ_016um_msu10"),
    ("msu15", 16, "H1-N3KXRGQ_016um_msu15"),
    ("msu20", 16, "H1-N3KXRGQ_016um_msu20"),
    ("msu25", 16, "H1-N3KXRGQ_016um_msu25"),
    ("msu50", 16, "H1-N3KXRGQ_016um_msu50"),
    ("msu50", 32, "H1-N3KXRGQ_032um"),
]


def snakemake_bench(files: list[str]) -> tuple[float, float]:
    """Summed runtime and largest PSS over Snakemake benchmark files.

    Args:
        files: Benchmark tsv paths; the last row of each is used.

    Returns:
        (minutes, GB), or NaNs if any file is missing.
    """
    if not all(os.path.exists(f) for f in files):
        return float("nan"), float("nan")
    rows = [pd.read_csv(f, sep="\t").iloc[-1] for f in files]
    return sum(r["s"] for r in rows) / 60, max(r["max_pss"] for r in rows) / 1024


def gnu_time(f: str) -> tuple[float, float]:
    """Wall time and max RSS from a /usr/bin/time -v report.

    Args:
        f: time -v output path.

    Returns:
        (minutes, GB), or NaNs if the file is missing.
    """
    if not os.path.exists(f):
        return float("nan"), float("nan")
    with open(f) as fh:
        txt = fh.read()
    wall = re.search(r"Elapsed \(wall clock\) time .*: ([\d:.]+)", txt).group(1)
    secs = 0.0
    for part in wall.split(":"):
        secs = secs * 60 + float(part)
    rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", txt).group(1))
    return secs / 60, rss / 1024**2


def main(out: str, log_dir: str, out_tsv: str) -> None:
    """Collect one row per method run.

    Args:
        out: Visium HD benchmark output root.
        log_dir: benchmark/logs directory.
        out_tsv: Output table.
    """
    cb, nb = f"{out}/clonalscope/bench", f"{out}/numbat/bench"
    runs = [
        ("numbat", "max_entropy_08", 32, [f"{nb}/run_numbat.{S16}.max_entropy_08.tsv"]),
        ("numbat", "default", 32, [f"{nb}/run_numbat.{S16}.tsv"]),
    ]
    for res, s in ((16, S16), (32, S32)):
        runs.append(
            (
                "clonalscope",
                "visium_tutorial",
                res,
                [f"{cb}/run_tutorial.{s}.visium_tutorial.tsv"],
            )
        )
    rows = []
    for method, setting, res, files in runs:
        rt, mem = snakemake_bench(files)
        rows.append((method, setting, res, rt, mem, files[-1]))
    for setting, res, tag in CALICOST:
        f = f"{log_dir}/calicost_{tag}.time.txt"
        rt, mem = gnu_time(f)
        rows.append(("calicost", setting, res, rt, mem, f))
    df = pd.DataFrame(
        rows,
        columns=[
            "method",
            "setting",
            "resolution [um]",
            "Runtime [mins]",
            "Memory [GB]",
            "source",
        ],
    ).round({"Runtime [mins]": 2, "Memory [GB]": 3})
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False, na_rep="NA")
    print(df.drop(columns="source").to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("usage: vishd_runtime.py <bench_out> <log_dir> <out_tsv>")
    main(*sys.argv[1:])
