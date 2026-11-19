"""Ported from cna_maste `python/cnamaste/truth.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim less the methods it never calls, with their imports
pointed at `scripts.cnamaste_port`.

The truth for a simulated sample: its true clones, their spot coordinates and allele-specific copies.
"""

from __future__ import annotations

import io
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scripts.cnamaste_port.io import read_text_table


@dataclass(frozen=True)
class Truth:
    """A CalicoST sample's planted clones and copies.

    Args: `barcodes`, `coords`, `labels`, `clones`, `profile`, each documented at its attribute.
    """

    barcodes: np.ndarray
    """Barcode per spot, as `truth_clone_labels.tsv` orders them."""
    coords: np.ndarray
    """`(n, 2)` float `x, y` per barcode."""
    labels: np.ndarray
    """Planted clone per barcode: `normal` 0, `clone_k` `k + 1`, as `remap_clone_num` numbers them."""
    clones: tuple[str, ...]
    profile: pd.DataFrame
    """`chr start end`, then `<clone>_A_copy`, `<clone>_B_copy` per segment."""

    @property
    def n_clones(self) -> int:
        return len(self.clones)


def _truth_table(sample: Path | Mapping[str, bytes], name: str, **kwargs: Any) -> pd.DataFrame:
    """The truth table `name` (tab-separated) from a sample directory, plain or `.gz`, or from its bytes by name.

    Args: `sample` a directory or name -> bytes; `name` the truth file; `kwargs` to `pandas.read_csv`.
    Returns: the frame.
    """
    if isinstance(sample, Path):
        return read_text_table(sample / name, sep="\t", **kwargs)
    return pd.read_csv(io.BytesIO(sample[name]), sep="\t", **kwargs)


def read_truth(sample: Path | Mapping[str, bytes]) -> Truth:
    """`sample`'s truth files: the planted clone labels, spot coordinates and allele-specific copies.

    Args: `sample` a directory, or the truth files' bytes by name.
    Returns: the `Truth`.
    """
    table = _truth_table(sample, "truth_clone_labels.tsv", index_col=0).reset_index()
    table.columns = ["barcode", "clone", "x", "y", *table.columns[4:]]
    others = sorted((n for n in table["clone"].unique() if n != "normal"), key=lambda n: int(n.removeprefix("clone_")))
    clones = ("normal", *others)
    index = {c: i for i, c in enumerate(clones)}
    return Truth(
        barcodes=table["barcode"].to_numpy().astype(str),
        coords=table[["x", "y"]].to_numpy(dtype=np.float64),
        labels=table["clone"].map(index).to_numpy(dtype=np.int64),
        clones=clones,
        profile=_truth_table(sample, "truth_acn_profile.tsv"),
    )
