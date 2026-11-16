"""Ported from cna_maste `python/cnamaste/io.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim, with their imports pointed at `scripts.cnamaste_port`.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from scripts.cnamaste_port.config.clock import start_time
from scripts.cnamaste_port.logger import get_logger
from scripts.cnamaste_port.utils import find_plain_or_gz

logger = get_logger(__name__, start_time=start_time)


pl.Config.set_tbl_cols(-1)


def read_text_table(path, **kwargs):
    """A delimited text table (`.csv`, `.tsv`, `.txt`) at `path` or `path.gz`, by `pandas.read_csv`.

    pandas infers gzip from the `.gz` suffix, so the compressed file is read as the plain one.

    Args: `path` the plain file's name; `kwargs` to `pandas.read_csv`.
    Returns: the frame.
    Raises: FileNotFoundError where neither exists.
    """
    found = find_plain_or_gz(path)
    if found is None:
        raise FileNotFoundError(f"{path} (or {path}.gz)")
    return pd.read_csv(found, **kwargs)


def get_spatial_positions(spaceranger_dir, filter_in_tissue=True):
    # TODO x,y vs row,col?  sub-pixel position?
    names = ("barcode", "in_tissue", "x", "y", "pixel_row", "pixel_col")

    if find_plain_or_gz(f"{spaceranger_dir}/spatial/tissue_positions.csv"):
        df_this_pos = read_text_table(
            f"{spaceranger_dir}/spatial/tissue_positions.csv",
            sep=",",
            header=0,
            names=names,
        )

        logger.info(f"Reading {spaceranger_dir}/spatial/tissue_positions.csv")

    elif find_plain_or_gz(f"{spaceranger_dir}/spatial/tissue_positions_list.csv"):
        df_this_pos = read_text_table(
            f"{spaceranger_dir}/spatial/tissue_positions_list.csv",
            sep=",",
            header=None,
            names=names,
        )

        logger.info(f"Reading {spaceranger_dir}/spatial/tissue_positions_list.csv")

    elif Path(f"{spaceranger_dir}/spatial/tissue_positions.parquet").exists():
        # NB 11,222,500 rows vs 4,992 rows for visium HD.
        #    see https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/spatial-outputs
        #
        #
        #    native columns :  barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres.
        #
        #    visium HD : x2 6.5mm capture areas, array_row = row coordinate of the spot in the array from 0 to 77.
        #                                        array_col = given orange crate [sic], i.e. hexagonal, for 6.55 this uses even numbers from 0 to 126 for even rows
        #                                                   and odd numbers from 1 to 127 for odd rows.
        #                                        pxl_row_in_fullres: The row pixel coordinate of the center of the spot in the full resolution image
        #                                        pxl_col_in_fullres: The column pixel coordinate of the center of the spot in the full resolution image. This is the x coordinate in pixel space.
        """
        df_this_pos = (
            pl.scan_parquet(f"{spaceranger_dir}/spatial/tissue_positions.parquet")
            .rename({"pxl_row_in_fullres": "y", "pxl_col_in_fullres": "x"})
            .select(["barcode", "in_tissue", "x", "y"])
            # .filter(pl.col("in_tissue") == True)
            .collect()
            .with_columns(pl.col("barcode").alias("square_002um"))
        )

        # NB native columns: square_002um, square_008um, square_016um, cell_id, in_nucleus, in_cell
        # TODO CHECK in_cell
        df_this_pos = (
            pl.scan_parquet(f"{spaceranger_dir}/spatial/barcode_mappings.parquet")
            .filter(pl.col("in_cell") == True)
            .collect()
            .join(
                df_this_pos.select(["square_002um", "x", "y", "in_tissue"]),
                on="square_002um",
                how="left",
            )
            .select(["square_002um", "cell_id", "x", "y", "in_tissue"])
            .group_by("cell_id")
            .agg(
                [
                    pl.col("x").mean().alias("x"),
                    pl.col("y").mean().alias("y"),
                    pl.col("square_002um").n_unique().alias("num_square_002um"),
                    pl.col("in_tissue").cast(pl.Boolean).any().alias("in_tissue"),
                ]
            )
            .with_columns(pl.col("cell_id").alias("barcode"))
        )
        """
        # NB x is pxl_col_in_fullres, increasing left to right, and y pxl_row_in_fullres, increasing top to bottom: the
        #    image's sense, opposite to a plot's y. Kept as read, in pixels, as `coords`; `physical_coords` turns them
        #    to um with y up (`platforms.visium.pixel_physical`), and the kNN and every spatial page read those.
        df_this_pos = (
            pl.scan_parquet(f"{spaceranger_dir}/spatial/tissue_positions.parquet")
            .rename({"pxl_row_in_fullres": "y", "pxl_col_in_fullres": "x"})
            .filter(pl.col("in_tissue") == True)
            .collect()
        )

        logger.info(f"Read {spaceranger_dir}/spatial/tissue_positions.parquet:\n{df_this_pos}")

        df_this_pos = df_this_pos.to_pandas()
    else:
        logger.error(f"No spatial coordinate file @ {spaceranger_dir}.")
        raise RuntimeError()

    # TODO alignment defined for in_tissue == True only?
    if filter_in_tissue:
        logger.warning(
            f"Filtering spatial positions to in_tissue == True (retained {100.0 * np.mean(df_this_pos.in_tissue)}%)."
        )

        result = df_this_pos[df_this_pos.in_tissue == True]
    else:
        result = df_this_pos

    # NB x,y positions for each barcode in this sample.
    return result
