"""Ported from cna_maste `python/cnamaste/annotation.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim, with their imports pointed at `scripts.cnamaste_port`.
"""

from pathlib import Path

import pandas as pd
from scripts.cnamaste_port.config.clock import start_time
from scripts.cnamaste_port.logger import get_logger

logger = get_logger(__name__, start_time=start_time)


def load_celltype(celltype_file: str | Path) -> pd.DataFrame:
    """Read the cell type per spot from a tab-separated annotation file with `barcode` and `cell_type` columns.

    Args: `celltype_file` the annotation file, e.g. `spot_annotation.tsv`, plain or `.gz`.
    Returns: frame of barcode, cell_type, one row per spot, in file order.
    Raises: ValueError where either column is missing.
    """
    df_celltype = pd.read_csv(celltype_file, sep="\t", usecols=["barcode", "cell_type"], dtype=str)[
        ["barcode", "cell_type"]
    ]

    spots_per_celltype = df_celltype["cell_type"].value_counts(dropna=False)

    logger.info(
        f"Read {len(df_celltype):_} spots of {len(spots_per_celltype)} cell types from {celltype_file}:\n"
        f"{spots_per_celltype.to_string()}"
    )

    return df_celltype
