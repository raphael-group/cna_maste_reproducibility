"""Visium HD square bins: barcodes read as grid places, placed on the image, given to finer bins."""

import numpy as np
import pandas as pd

VISIUMHD_GRID_UM = 6_700
"""Visium HD's grid of bins, um a side: 3,350 2 um squares, from the capture area's corner; a bin of side `s` at
`(row, col)` lies on it where `row * s` and `col * s` are less."""

VISIUMHD_BARCODE = r"^s_(\d{3})um_(\d{5})_(\d{5})-1$"
"""A Visium HD square bin's barcode, `s_<side>um_<row>_<col>-1`: its side in um, then its row and column on the grid of
bins that side."""


def get_visiumhd_bins(barcodes):
    """Each Visium HD bin's side, row and column, read from its barcode.

    Args: `barcodes` the bins' barcodes.
    Returns: (n_bins, 3) int, `(side_um, row, col)`.
    Raises: `ValueError` naming a barcode that is not a Visium HD square bin's.
    """
    found = pd.Series(np.asarray(barcodes, dtype=str)).str.extract(VISIUMHD_BARCODE)
    if found.isna().any().any():
        raise ValueError(f"not Visium HD square bins, e.g. {pd.Series(barcodes)[found[0].isna().to_numpy()].iloc[0]!r}")
    return found.astype(np.int64).to_numpy()


def get_visiumhd_bin_pixels(barcodes, positions):
    """Each Visium HD bin's centre in full-resolution pixels, `(pxl_col, pxl_row)`, for bins of any side.

    A bin of side `s` at `(row, col)` is centred `((col + 0.5) s, (row + 0.5) s)` um on the capture area, every side on
    one grid; the affine map from those um to the image's pixels is fitted on `positions`' bins, whose pixels Space
    Ranger gives, and taken to `barcodes`.

    Args: `barcodes` the bins to place; `positions` `io.get_spatial_positions`' frame of one Space Ranger directory,
        `barcode` and its pixels `x` (pxl_col_in_fullres), `y` (pxl_row_in_fullres).
    Returns: (n_bins, 2) pixels, `(pxl_col, pxl_row)`.
    """

    def capture_um(found):
        side, row, col = found.T
        return np.column_stack([(col + 0.5) * side, (row + 0.5) * side, np.ones(len(found))]).astype(np.float64)

    fit = np.linalg.lstsq(
        capture_um(get_visiumhd_bins(positions["barcode"])), positions[["x", "y"]].to_numpy(np.float64), rcond=None
    )[0]
    return capture_um(get_visiumhd_bins(barcodes)) @ fit


def get_visiumhd_finer_labels(labels, side_um):
    """Each label of a Visium HD bin given to every bin of side `side_um` inside it.

    Bins of every side share one grid from the capture area's corner, so a bin of side `s` at `(row, col)` holds the
    `f^2` bins of side `s / f` at rows `row * f ... row * f + f - 1` and columns alike, those on the grid
    (`VISIUMHD_GRID_UM`): a bin at the grid's far edge runs past it, and holds fewer.

    Args: `labels` a label per bin, indexed by barcode, the bins of one side; `side_um` the finer bins' side in um.
    Returns: the labels indexed by the finer bins' barcodes, each bin's in turn.
    Raises: `ValueError` where the bins are of several sides, or `side_um` does not divide theirs.
    """
    side, row, col = get_visiumhd_bins(labels.index).T
    sides = np.unique(side)
    if len(sides) != 1 or sides[0] % side_um:
        raise ValueError(f"bins of one side, a multiple of {side_um} um; found {sides.tolist()} um")
    f = int(sides[0] // side_um)
    i, j = np.divmod(np.arange(f * f), f)
    rows, cols = (row[:, None] * f + i).ravel(), (col[:, None] * f + j).ravel()
    on = (rows * side_um < VISIUMHD_GRID_UM) & (cols * side_um < VISIUMHD_GRID_UM)
    barcodes = [f"s_{side_um:03d}um_{r:05d}_{c:05d}-1" for r, c in zip(rows[on], cols[on])]
    return pd.Series(np.repeat(labels.to_numpy(), f * f)[on], index=barcodes, name=labels.name)
