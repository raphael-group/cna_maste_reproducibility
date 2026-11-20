"""UGP's 16 um bb summed to coarser bins (`scripts/utils/aggregate_bb.py`): the bins numbat/aggregate_bins.R and
common/aggregate_10x_matrix.R make, each count kept.

Referee: the coarse bin of 16 um bin (r, c) at side X is (r // k, c // k), k = X / 16 (aggregate_bins.R); the sums
counted by hand on four bins.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

anndata = pytest.importorskip("anndata")

from scripts.utils.aggregate_bb import get_coarse_bins, main

ASSAY = "VISIUM3prime"


def write_bb(root: Path) -> list[str]:
    """A bb dir of four 16 um bins, three in coarse bin (0, 0) and one in (1, 1); the h5ad rows reversed."""
    keys = [f"s_016um_{r:05d}_{c:05d}-1_U1_{ASSAY}" for r, c in [(0, 0), (0, 1), (1, 0), (2, 3)]]
    (root / "unit" / ASSAY).mkdir(parents=True)
    (root / "multi_snp" / ASSAY).mkdir(parents=True)
    with gzip.open(root / "unit" / ASSAY / "barcodes.tsv.gz", "wt") as f:
        f.write("\n".join(keys) + "\n")
    a = sp.csr_matrix(np.array([[1, 2, 3, 4], [0, 1, 0, 5]], dtype=np.int64))
    sp.save_npz(root / "multi_snp" / ASSAY / "bb.Aallele.npz", a)
    sp.save_npz(root / "multi_snp" / ASSAY / "bb.Ballele.npz", a * 2)
    pd.DataFrame({"#CHR": ["chr1", "chr2"], "START": [1, 1], "END": [9, 9]}).to_csv(
        root / "multi_snp" / ASSAY / "bb.tsv.gz", sep="\t", index=False
    )
    x = sp.csr_matrix(np.array([[1.0, 0.0], [2.0, 1.0], [0.0, 3.0], [4.0, 4.0]], dtype=np.float32))
    var = pd.DataFrame({"gene_symbol": ["G1", "G2"]}, index=pd.Index(["g1", "g2"], dtype=object))
    var["gene_symbol"] = var["gene_symbol"].astype(object)
    rows = anndata.AnnData(X=x, obs=pd.DataFrame(index=pd.Index(keys, dtype=object)), var=var)[::-1].copy()
    rows.write_h5ad(root / f"{ASSAY}.h5ad")
    return keys


def test_coarse_bins_follow_the_r_div_k_rule() -> None:
    """A 16 um bin (r, c) goes to (r // 2, c // 2) at 32 um, (r // 4, c // 4) at 64 um; a side not above it, or not
    a multiple of it, is refused."""
    keys = [f"s_016um_00003_00005-1_U1_{ASSAY}", f"s_016um_00004_00007-1_U1_{ASSAY}"]
    assert list(get_coarse_bins(keys, 32)[0]) == [
        f"s_032um_00001_00002-1_U1_{ASSAY}",
        f"s_032um_00002_00003-1_U1_{ASSAY}",
    ]
    assert list(get_coarse_bins(keys, 64)[0]) == [
        f"s_064um_00000_00001-1_U1_{ASSAY}",
        f"s_064um_00001_00001-1_U1_{ASSAY}",
    ]
    with pytest.raises(ValueError, match="not a multiple"):
        get_coarse_bins(keys, 24)
    with pytest.raises(ValueError, match="not a multiple"):
        get_coarse_bins(keys, 16)


def test_bb_is_summed_per_coarse_bin(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Each coarse bin's allele and RNA counts are its bins' sums, its h5ad rows in the barcodes' order; the blocks
    are copied."""
    write_bb(tmp_path / "bb")
    out = tmp_path / "bb-032um"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "aggregate_bb.py",
            "--bb-dir",
            str(tmp_path / "bb"),
            "--assay",
            ASSAY,
            "--bin-um",
            "32",
            "--out-dir",
            str(out),
        ],
    )
    main()

    coarse = [f"s_032um_00000_00000-1_U1_{ASSAY}", f"s_032um_00001_00001-1_U1_{ASSAY}"]
    assert pd.read_csv(out / "unit" / ASSAY / "barcodes.tsv.gz", header=None)[0].tolist() == coarse
    a = sp.load_npz(out / "multi_snp" / ASSAY / "bb.Aallele.npz").toarray()
    np.testing.assert_array_equal(a, [[6, 4], [1, 5]])
    np.testing.assert_array_equal(sp.load_npz(out / "multi_snp" / ASSAY / "bb.Ballele.npz").toarray(), a * 2)
    assert (out / "multi_snp" / ASSAY / "bb.tsv.gz").is_file()
    assert not (out / "multi_snp" / ASSAY / "bb.Tallele.npz").exists()

    summed = anndata.read_h5ad(out / f"{ASSAY}.h5ad")
    assert summed.obs_names.tolist() == coarse
    np.testing.assert_array_equal(summed.X.toarray(), [[3.0, 4.0], [4.0, 4.0]])
    assert summed.obs["n_bins_016um"].tolist() == [3, 1]
    assert summed.var_names.tolist() == ["g1", "g2"]


def test_existing_output_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """An existing output directory is not overwritten."""
    write_bb(tmp_path / "bb")
    (tmp_path / "bb-032um").mkdir()
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "aggregate_bb.py",
            "--bb-dir",
            str(tmp_path / "bb"),
            "--assay",
            ASSAY,
            "--bin-um",
            "32",
            "--out-dir",
            str(tmp_path / "bb-032um"),
        ],
    )
    with pytest.raises(FileExistsError):
        main()
