#!/usr/bin/env python
"""Per-simulation numbat input sizes, step timings and resources for finished runs.

Runpeng Luo (2026-10-30)

Dependencies: pandas, scipy (calicost_env).

Usage:
    python -I numbat/run_stats.py <numbat_out> <run_label> <out_tsv>
            numbat_out  # benchmark numbat/ output root
            run_label   # config run_label (e.g. max_entropy_08)
            out_tsv     # one row per finished simulation

Inputs:
    <numbat_out>/<sim>/prep/<sim>_allele_counts.tsv.gz
    <numbat_out>/<sim>/<run_label>/run_numbat/{log.txt, segs_consensus_*.tsv}
    <numbat_out>/logs/run_numbat.<sim>.<run_label>.log (genes kept)
    <numbat_out>/bench/run_numbat.<sim>.<run_label>.tsv (snakemake benchmark)
Outputs:
    out_tsv: spots, genes_kept (log "number of genes left"), het_snps, allele_reads,
        segs_iter{1,2} and cnv_segs_iter{1,2}, per-cell step minutes per iteration,
        wall_min, cpu_h, max_rss_gb, max_pss_gb.
Notes/References:
    max_rss sums forked workers' shared pages (mclapply), so max_pss is the closer memory figure.
"""

import os
import re
import sys
from datetime import datetime, timezone

import pandas as pd


def step_minutes(log: str) -> list[float]:
    """Minutes between each "Evaluating CNV per cell" and the next "All cells succeeded".

    Args:
        log: numbat log.txt path.

    Returns:
        One duration per iteration.
    """
    ts = re.compile(r"\[(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\] (.*)")
    out, start = [], None
    with open(log) as f:
        lines = f.readlines()
    for line in lines:
        m = ts.search(line)
        if not m:
            continue
        t = datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
        if "Evaluating CNV per cell" in m.group(2):
            start = t
        elif start and ("cells succeeded" in m.group(2) or "failed" in m.group(2)):
            out.append((t - start).total_seconds() / 60)
            start = None
    return out


def main(root: str, label: str, out_tsv: str) -> None:
    """Collect stats for every finished simulation.

    Args:
        root: numbat output root.
        label: run_label.
        out_tsv: Output table.
    """
    rows = []
    for sim in sorted(os.listdir(root)):
        rd = f"{root}/{sim}/{label}/run_numbat"
        bench = f"{root}/bench/run_numbat.{sim}.{label}.tsv"
        if not (os.path.exists(f"{rd}/done") and os.path.exists(bench)):
            continue
        b = pd.read_csv(bench, sep="\t").iloc[0]
        al = pd.read_csv(
            f"{root}/{sim}/prep/{sim}_allele_counts.tsv.gz",
            sep="\t",
            usecols=["cell", "snp_id", "DP"],
        )
        with open(f"{root}/logs/run_numbat.{sim}.{label}.log") as f:
            rule_log = f.read()
        genes = re.search(r"number of genes left: (\d+)", rule_log)
        row = {
            "sim": sim,
            "spots": al["cell"].nunique(),
            "genes_kept": int(genes.group(1)) if genes else None,
            "het_snps": al["snp_id"].nunique(),
            "allele_reads": int(al["DP"].sum()),
        }
        for i in (1, 2):
            f = f"{rd}/segs_consensus_{i}.tsv"
            if os.path.exists(f):
                s = pd.read_csv(f, sep="\t")
                row[f"segs_iter{i}"] = len(s)
                row[f"cnv_segs_iter{i}"] = int((s["cnv_state"] != "neu").sum())
        for i, m in enumerate(step_minutes(f"{rd}/log.txt"), 1):
            row[f"percell_min_iter{i}"] = round(m, 1)
        row.update(
            wall_min=round(b["s"] / 60, 1),
            cpu_h=round(b["cpu_time"] / 3600, 1),
            max_rss_gb=round(b["max_rss"] / 1024, 1),
            max_pss_gb=round(b["max_pss"] / 1024, 1),
        )
        rows.append(row)
    res = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    res.to_csv(out_tsv, sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(res.drop(columns="sim").assign(sim=res["sim"].str.replace("_ploidy2", "")).to_string(index=False))
    if len(res) >= 3:
        num = res.drop(columns="sim").select_dtypes("number")
        num = num.loc[:, num.nunique() > 1]
        print("\nSpearman correlation with wall_min / max_pss_gb:")
        print(num.corr(method="spearman")[["wall_min", "max_pss_gb"]].round(2).to_string())


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("usage: numbat/run_stats.py <numbat_out> <run_label> <out_tsv>")
    main(*sys.argv[1:])
