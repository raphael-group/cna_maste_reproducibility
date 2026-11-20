"""Ported from cna_maste `python/cnamaste/platforms/visium.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim, with their imports pointed at `scripts.cnamaste_port`.

Visium and Visium HD: the array's indices and the bins' pixels as physical positions, and which a Space Ranger
directory holds.
"""

import json
from pathlib import Path

import numpy as np


def pixel_physical(coords, microns_per_pixel):
    """Visium HD's bin pixels `(pxl_col, pxl_row)` as physical positions `(x, y)` in um: x the column, left to right,
    and y the row, up, as the slide is seen; full-resolution pixels are square, `microns_per_pixel` a side.

    Args: `coords` (n_spots, 2) pixels, `(pxl_col_in_fullres, pxl_row_in_fullres)` as `io.get_spatial_positions`
    reads them; `microns_per_pixel` Space Ranger's `scalefactors_json.json` value.
    Returns: (n_spots, 2) `(x, y)` in um.
    """
    pixels = np.asarray(coords, dtype=np.float64)
    return np.column_stack([pixels[:, 0], -pixels[:, 1]]) * float(microns_per_pixel)


def get_visiumhd_microns_per_pixel(spaceranger_dir):
    """Visium HD's full-resolution pixel size: `spatial/scalefactors_json.json`'s `microns_per_pixel`.

    Args: `spaceranger_dir` a slice's Space Ranger directory.
    Returns: um per pixel.
    Raises: `ValueError` where the file or the key is missing: a bin's pixels then have no size in um.
    """
    path = Path(f"{spaceranger_dir}/spatial/scalefactors_json.json")
    found = json.loads(path.read_text()).get("microns_per_pixel") if path.is_file() else None
    if found is None:
        raise ValueError(f"{path} gives no microns_per_pixel: Visium HD's bin pixels have no size in um.")
    return float(found)
