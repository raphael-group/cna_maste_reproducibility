"""The runtime tables the metrics read: each method run's runtime and peak memory, simulations and H1-N3KXRGQ.

Runpeng Luo (2026-11-01)

Dependencies: pandas.

Functions (run_cna_maste_reproducibility.smk rules sim_runtime and vishd_runtime):
    write_simulation_runtime(sim_dir, bench_out, deposit_numbat_dir, out_tsv): sim_dir, simulated_data_related/;
        bench_out, the simulations' output root (numbat/, clonalscope/, cnamaste/); deposit_numbat_dir, the
        deposit's numbat_related/.
    write_vishd_runtime(bench_out, log_dir, out_tsv): bench_out, the Visium HD sample's output root (clonalscope/,
        numbat/, cnamaste/); log_dir, the baseline logs (CalicoST /usr/bin/time -v files).

Inputs:
    simulation:
        numbat: <bench_out>/numbat/bench/run_numbat.<sim>.max_entropy_08.tsv
        clonalscope: <bench_out>/clonalscope/bench/run_tutorial.<sim>.visium_tutorial.tsv
        cnamaste: <bench_out>/cnamaste/bench/run_cnamaste.<sim>.expansion-lattice-lattice.seed0.tsv
            (run_cna_maste_reproducibility.smk rule cnamaste_sim, the solver and seed run_sim_qa scores)
        numbat-zenodo: <deposit_numbat_dir>/<sim>/outs/log.txt (numbat run log)
    vishd:
        numbat: <bench_out>/numbat/bench/run_numbat.H1-N3KXRGQ[.<run_label>].tsv
        clonalscope visium_tutorial: <bench_out>/clonalscope/bench/run_tutorial.<s>.visium_tutorial.tsv
        calicost: <log_dir>/calicost_<tag>.time.txt
        cnamaste: <bench_out>/cnamaste/bench/run_cnamaste.H1-N3KXRGQ.tsv, 16 um bins, its setting the run id of
            <bench_out>/cnamaste/run_id.txt (real_data_metrics.get_run_runtime finds it in the run's run_<id>)
Outputs:
    simulation out_tsv: simulation_id, method, Runtime [mins], Memory [GB]; simulation_metrics.read_runtime reads it.
    vishd out_tsv: method, setting, resolution [um], Runtime [mins], Memory [GB], source; the benchmark yaml's
        vishd.<sample>.runtime.
Notes/References:
    Snakemake benchmark files: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#benchmark-rules
    Snakemake runs: runtime is wall-clock seconds (s) / 60, memory max_pss / 1024 (MB -> GB); PSS counts pages
    shared by forked workers once, unlike max_rss. Only the method's main run job is counted, not shared input
    preprocessing.
    CalicoST: runtime = "Elapsed (wall clock) time", memory = "Maximum resident set size" (kB), from
    /usr/bin/time -v; CalicoST runs as one process (OMP_NUM_THREADS=1), so RSS is used.
    A run without a timing record (failed, not finished, or not run by this pipeline) gets NA.
    numbat-zenodo runtime = last minus first [YYYY-mm-dd HH:MM:SS] timestamp in log.txt, which has one-second
    resolution and covers run_numbat() only; its memory is NA because the log reports R heap usage at checkpoints,
    not peak process memory.
"""

import os
import re
from datetime import datetime, timezone

import pandas as pd

BENCH = {
    "numbat": "numbat/bench/run_numbat.{sim}.max_entropy_08.tsv",
    "clonalscope": "clonalscope/bench/run_tutorial.{sim}.visium_tutorial.tsv",
    "cnamaste": "cnamaste/bench/run_cnamaste.{sim}.expansion-lattice-lattice.seed0.tsv",
}

TS = re.compile(r"\[(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\]")

S16, S32 = "H1-N3KXRGQ", "H1-N3KXRGQ_032um"
CALICOST = [
    ("msu10", 16, "H1-N3KXRGQ_016um_msu10"),
    ("msu15", 16, "H1-N3KXRGQ_016um_msu15"),
    ("msu20", 16, "H1-N3KXRGQ_016um_msu20"),
    ("msu25", 16, "H1-N3KXRGQ_016um_msu25"),
    ("msu50", 16, "H1-N3KXRGQ_016um_msu50"),
    ("msu50", 32, "H1-N3KXRGQ_032um"),
]


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


def write_simulation_runtime(sim_dir: str, out: str, dep: str, out_tsv: str) -> None:
    """Write one row per simulation and method.

    Args:
        sim_dir: simulated_data_related directory.
        out: The simulations' output root.
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


def write_vishd_runtime(out: str, log_dir: str, out_tsv: str) -> None:
    """Write one row per method run on the Visium HD sample.

    Args:
        out: The Visium HD sample's output root.
        log_dir: The baseline logs directory.
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
    with open(f"{out}/cnamaste/run_id.txt") as f:
        runs.append(("cnamaste", f.read().strip(), 16, [f"{out}/cnamaste/bench/run_cnamaste.{S16}.tsv"]))
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
