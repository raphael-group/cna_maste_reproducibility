#!/usr/bin/env python3
"""Numbat postprocess: run_numbat outputs -> standardized benchmark deliverables.

Copied from copytyping-benchmark workflows/numbat/scripts/postprocess_numbat.py. Change: a cell
without a -RNA/-ATAC suffix (numbat scRNA mode, e.g. Visium spots) is the RNA modality, so scRNA-mode
runs emit <sample>.RNA.{cells,cluster_states}.tsv.

TWO runs, one per modality (RNA, ATAC). Numbat's combined-bin mode scores each nucleus twice (a -RNA
and a -ATAC pseudo-cell) and assigns each an INDEPENDENT clone; numbat prescribes no merge, so we
keep them apart. Both modalities SHARE the clone CN profile -- only the cell -> clone map differs.

There is no bulk-guided variant: run_numbat is now called with `n_cut` (config/numbat/params.json,
run.n_cut = 3 -> K = 4 clones incl. the diploid normal), so the run's OWN clones already sit at the
target K and numbat's clone pseudobulks are pooled over exactly those clones. (K is a target, not a
guarantee -- the collapse stops early if no finite-cost contraction remains, and a clone with no
cells drops out.)

Clone CN = the CLONE PSEUDOBULK profile (`bulk_clones_final.tsv.gz`, `cnv_state_post` per (clone,
segment)), taken verbatim. After convergence numbat re-pools each clone's cells and re-runs the HMM
on that pseudobulk (R/main.R, post-loop), and that is what it reports as the clone CNV profile
(bulk_clones_final.png). NOT the phylogeny genotype (`clone_post.GT_opt`), which can disagree: an
event present in a clone's pseudobulk may not be placed on its branch (numbat FAQ: "filtered out due
to high entropy (weak evidence) in single cells").

Numbat's aneuploid posterior `p_cnv` gives the tumour/normal call (>0.5 -> malignant); the normal
clone is the all-neutral (empty GT) clade (fallback: lowest mean p_cnv). Native numbat events ->
6-state taxonomy: del->del, loh->cnLOH, bamp->WGD, bdel->bdel, amp->amp, else neu.

Stdlib only. Usage:
  python numbat/postprocess.py  sample  run_dir  out_dir
    run_dir   output/<sample>/run_numbat   (holds clone_post_*.tsv, bulk_clones_final.tsv.gz)

Outputs (out_dir), per modality <m> in {RNA, ATAC}:
  <sample>.<m>.cells.tsv           BARCODE, malignant (N/T), clone
  <sample>.<m>.cluster_states.tsv  segment x clone CN states (#CHR, START, END, raw_cn_<clone>, transformed_cn_<clone>)

The normal clone is not written out: it is the clone whose cells carry malignant == "N" in cells.tsv.
"""

import csv
import glob
import gzip
import math
import os
import re
import sys

_CHR_ORDER = {c: i for i, c in enumerate([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])}
# native numbat event -> 6-state (allele-aware)
_NUMBAT_STATE = {"neu": "neu", "del": "del", "amp": "amp", "loh": "cnLOH", "bamp": "WGD", "bdel": "bdel"}


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def canon_bc(bc):
    return re.sub(r"_(R|S)\d+$", "", re.sub(r"-(RNA|ATAC)$", "", str(bc)))


def chrname(c):
    c = str(c)
    return c if c.startswith("chr") else f"chr{c}"


def _fnum(x):
    try:
        v = float(x)
        return -1.0 if math.isnan(v) else v  # NB: NaN sorts lowest
    except (TypeError, ValueError):
        return -1.0


def _read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _highest(run_dir, stem):
    """{iter -> path} for run_dir/<stem>_<n>.tsv."""
    out = {}
    for p in glob.glob(os.path.join(run_dir, f"{stem}_*.tsv")):
        m = re.search(rf"{stem}_(\d+)", p)
        if m:
            out[int(m.group(1))] = p
    return out


# ---------------------------------------------------------------- cells
def _modality(cell):
    """RNA | ATAC -- the pseudo-cell suffix numbat's combined-bin mode carries; RNA if absent (scRNA mode)."""
    m = re.search(r"-(RNA|ATAC)$", str(cell))
    return m.group(1) if m else "RNA"


def _cells_by_modality(rows):
    """clone_post rows -> {modality -> {barcode: {'clone', 'p_cnv', 'malignant'}}}.

    Numbat's combined-bin mode scores each nucleus TWICE -- a -RNA and a -ATAC pseudo-cell
    (inst/bin/input_prep.R binCnt: cbind(rna, atac)) -- and emits an INDEPENDENT clone posterior
    for each. Numbat prescribes no merge, and their posteriors are not on a comparable scale
    (ATAC's saturate at 1.0 far more often than RNA's), so each modality stays its own run.
    """
    out = {}
    for r in rows:
        if r.get("clone_opt") in (None, "", "NA"):
            continue
        mod = _modality(r["cell"])
        if mod is None:
            continue
        pcnv = _fnum(r.get("p_cnv"))
        out.setdefault(mod, {})[canon_bc(r["cell"])] = {
            "clone": f"C{r['clone_opt']}",
            "clone_opt": r["clone_opt"],
            "p_cnv": pcnv,
            "malignant": "T" if pcnv > 0.5 else "N",
        }
    return out


# ---------------------------------------------------------------- cluster states (joint per-clone CN)
def _clone_gt(rows):
    """clone_post rows -> {clone_opt(str) -> GT_opt string} (one per clone)."""
    gt = {}
    for r in rows:
        cl = r.get("clone_opt")
        if cl in (None, "", "NA"):
            continue
        gt.setdefault(str(cl), (r.get("GT_opt") or ""))
    return gt


def _bulk_clones_path(run_dir):
    """bulk_clones_final.tsv.gz, else the highest-iteration bulk_clones_<i>.tsv.gz (None if neither).

    A run that dies in the phylogeny step never reaches numbat's post-loop pseudobulk block, so
    `_final` is absent (HT973N1: "No CNV remains after filtering by entropy in single cells"). The
    per-iteration file is the same object from the last completed iteration.
    """
    final = os.path.join(run_dir, "bulk_clones_final.tsv.gz")
    if os.path.exists(final):
        return final
    it = {}
    for p in glob.glob(os.path.join(run_dir, "bulk_clones_*.tsv.gz")):
        m = re.search(r"bulk_clones_(\d+)\.tsv\.gz$", p)
        if m:
            it[int(m.group(1))] = p
    return it[max(it)] if it else None


def _bulk_clone_rows(path):
    """Clone pseudobulk profiles -> long rows (chrom,start,end,clone,raw,state).

    numbat re-pools each clone's cells and re-runs the HMM on that pseudobulk (main.R, post-loop),
    writing one `cnv_state_post` per (clone, segment); `sample` is the clone id. Taken verbatim --
    no min_LLR filter (numbat applies one only when PLOTTING, in plot_psbulk).
    """
    with gzip.open(path, "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out, seen = [], set()
    for r in rows:
        seg, clone = r.get("seg"), r.get("sample")
        if not seg or seg == "NA" or clone in (None, "", "NA"):
            continue
        key = (clone, seg)
        if key in seen:
            continue
        seen.add(key)
        raw = str(r.get("cnv_state_post") or "neu")
        out.append(
            (
                chrname(r["CHROM"]),
                int(r["seg_start"]),
                int(r["seg_end"]),
                f"C{clone}",
                raw,
                _NUMBAT_STATE.get(raw, "neu"),
            )
        )
    return out


def _write_cluster_states(rows, path):
    """long rows -> wide [#CHR, START, END, raw_cn_<clone>, transformed_cn_<clone>], genome-sorted."""
    clones = list(dict.fromkeys(r[3] for r in rows))
    seg_key, raw, st = [], {}, {}
    for chrom, s, e, clone, rw, state in rows:
        k = (chrom, s, e)
        if k not in raw:
            seg_key.append(k)
            raw[k] = {}
            st[k] = {}
        raw[k][clone] = rw
        st[k][clone] = state
    seg_key.sort(key=lambda k: (_CHR_ORDER.get(k[0], 99), k[1]))
    header = ["#CHR", "START", "END"] + [f"raw_cn_{c}" for c in clones] + [f"transformed_cn_{c}" for c in clones]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for k in seg_key:
            w.writerow(list(k) + [raw[k].get(c, "neu") for c in clones] + [st[k].get(c, "neu") for c in clones])


def _normal_clone(clone_gt, cells):
    """Normal clone = the all-neutral (empty GT) clade; ties/none -> lowest mean p_cnv."""
    empty = [c for c, gt in clone_gt.items() if not gt.strip()]
    if len(empty) == 1:
        return f"C{empty[0]}"
    # per-clone mean p_cnv
    agg = {}
    for d in cells.values():
        agg.setdefault(d["clone"], []).append(d["p_cnv"])
    if empty:
        cand = {f"C{c}" for c in empty}
        agg = {k: v for k, v in agg.items() if k in cand} or agg
    if not agg:
        return None
    return min(agg, key=lambda k: sum(agg[k]) / len(agg[k]))


# ---------------------------------------------------------------- emit one variant
def emit(sample, cells, clone_gt, cn_rows, modality, out_dir):
    """cells.tsv + cluster_states.tsv for one modality. malignant: numbat's own p_cnv > 0.5, EXCEPT
    the normal clone's cells, which are forced to N so the clone -> N/T map stays consistent."""
    normal = _normal_clone(clone_gt, cells)
    K = len({d["clone"] for d in cells.values()})
    pre = os.path.join(out_dir, f"{sample}.{modality}")

    with open(f"{pre}.cells.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["BARCODE", "malignant", "clone"])
        for bc in sorted(cells):
            w.writerow([bc, cells[bc]["malignant"], cells[bc]["clone"]])

    if cn_rows:
        _write_cluster_states(cn_rows, f"{pre}.cluster_states.tsv")
    log(f"[{sample}/{modality}] K={K} normal={normal} -> {pre}.{{cells,cluster_states}}.tsv")


def main():
    args = sys.argv[1:]
    if len(args) < 3:
        sys.exit(__doc__)
    sample, run_dir, out_dir = args[:3]
    os.makedirs(out_dir, exist_ok=True)

    posts = _highest(run_dir, "clone_post")
    if not posts:
        sys.exit(f"[{sample}] no clone_post_*.tsv in {run_dir}")
    rows = _read_tsv(posts[max(posts)])
    gt = _clone_gt(rows)

    bc = _bulk_clones_path(run_dir)
    if bc is None:
        sys.exit(f"[{sample}] no bulk_clones_*.tsv.gz in {run_dir}")
    cn_rows = _bulk_clone_rows(bc)
    log(
        f"[{sample}] CN from {os.path.basename(bc)} "
        f"({len({r[3] for r in cn_rows})} clones x {len({r[:3] for r in cn_rows})} segs)"
    )

    for mod, cells in sorted(_cells_by_modality(rows).items()):
        emit(sample, cells, gt, cn_rows, mod, out_dir)


if __name__ == "__main__":
    main()
