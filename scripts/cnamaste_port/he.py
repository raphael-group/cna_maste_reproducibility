"""Ported from cna_maste `python/cnamaste/he.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim, with their imports pointed at `scripts.cnamaste_port`.
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
from scipy.spatial import cKDTree
from scripts.cnamaste_port.config.clock import start_time
from scripts.cnamaste_port.logger import get_logger

logger = get_logger(__name__, start_time=start_time)


def join_tables_xy(first: pl.DataFrame, second: pl.DataFrame, columns_to_merge: list[str]) -> pl.DataFrame:
    """Join `second` to `first` on `columns_to_merge`, matching each row of `first` to the closest (x, y) in `second`.

    Args: `first`, `second` frames with `x`, `y`; `columns_to_merge` the columns of `second` to carry.
    Returns: `first` with those columns joined.
    """
    first_xy = np.column_stack([first["x"].to_numpy(), first["y"].to_numpy()])
    second_xy = np.column_stack([second["x"].to_numpy(), second["y"].to_numpy()])

    tree = cKDTree(second_xy)

    # NB for each row in first, find closest in second.
    distances, indices = tree.query(first_xy)

    new_columns = [pl.Series(col, second[col].to_numpy()[indices]) for col in columns_to_merge]

    new_columns = new_columns + [pl.Series("dist", distances)]

    return first.with_columns(new_columns)


def get_he_image(spaceranger_dir, res="hires", pos=None, num_labels=4, max_side=None):
    assert res in ("lowres", "hires")
    # scalefactor = target_size / max (original image height, original image width),
    #
    # e.g. {
    #    "spot_diameter_fullres": 58.45684684229273,
    #    "bin_size_um": 16.0,
    #    "microns_per_pixel": 0.2737061758251425,
    #    "tissue_lowres_scalef": 0.0071874363,
    #    "fiducial_diameter_fullres": 1205.6724661222877,
    #    "tissue_hires_scalef": 0.071874365,
    #    "regist_target_img_scalef": 0.071874365
    # }
    sf_path = Path(f"{spaceranger_dir}/spatial/scalefactors_json.json")
    img_path = Path(f"{spaceranger_dir}/spatial/tissue_{res}_image.png")

    if not sf_path.exists() or not img_path.exists():
        logger.warning(f"Could not find H&E image or scalefactors at {spaceranger_dir}/spatial/")
        return pos

    with open(sf_path, "r") as ff:
        scalefactors = json.load(ff)

    scalefactor = scalefactors[f"tissue_{res}_scalef"]

    # NB realizes a (H, W, C) numpy array, i.e. (382, 600, 3) for low and (3818, 6000, 3) for high (HD @ 6.5mm).
    tissue_image = plt.imread(img_path)

    # NB every k-th row and column where the image is larger than a panel draws: each kept pixel keeps its place,
    #    at k times the scale factor's pixel size.
    if max_side is not None:
        k = int(np.ceil(max(tissue_image.shape[:2]) / max_side))
        tissue_image, scalefactor = tissue_image[::k, ::k], scalefactor / k

    # NB 11,222,500 rows in the 2 um version of the file for Visium HD
    # full_width = np.ceil(6000 / 0.071874365) = 83_479.0
    # full_height = np.ceil(3818 / 0.071874365) = 53_121.0

    # NB barcode positions correspond to 5_641.821 < y < 30_099.837, 1_190.55 < x < 5_641.82 for 175_561 spots.
    #    i.e. the aligned portion of the original image.
    H, W, C = tissue_image.shape
    rows, cols = np.meshgrid(np.arange(H), np.arange(W), indexing="ij")

    tissue_frame = pl.DataFrame(
        {
            "red": tissue_image[:, :, 0].flatten(),
            "green": tissue_image[:, :, 1].flatten(),
            "blue": tissue_image[:, :, 2].flatten(),
            "array_row": rows.flatten(),
            "array_col": cols.flatten(),
        }
    )

    # NB 0.0 < x < 83_339.869; 0 < y < 53_009.165
    tissue_frame = tissue_frame.with_columns(
        (pl.col("array_col") / scalefactor).alias("x"),
        (pl.col("array_row") / scalefactor).alias("y"),
    )

    if pos is not None:
        # NB limited to in_tissue=True
        pos = pl.from_pandas(pos)

        columns = ("red", "green", "blue")
        # NB each spot takes its nearest pixel in the image's own frame, full-resolution pixels: a tissue positions
        #    file's pixel_col / pixel_row, not its x / y, which are the array's indices (`io.get_spatial_positions`).
        on_pixels = "pixel_col" in pos.columns and "pixel_row" in pos.columns
        joined = pos.with_columns(pl.col("pixel_col").alias("x"), pl.col("pixel_row").alias("y")) if on_pixels else pos
        tissue_frame = join_tables_xy(joined, tissue_frame, columns)
        if on_pixels:
            tissue_frame = tissue_frame.with_columns(pos["x"], pos["y"])

    def crop_values(x):
        return (x - np.min(x)) / (np.max(x) - np.min(x))

    rgb = tissue_frame.select(["red", "green", "blue"]).to_numpy()

    # TODO standard lib?
    gray = 0.2125 * rgb[:, 0] + 0.7154 * rgb[:, 1] + 0.0721 * rgb[:, 2]
    cropped_gray = crop_values(gray)

    percentiles = np.linspace(0.0, 100.0, 1 + num_labels)
    bins = np.percentile(np.sort(cropped_gray.flatten()), percentiles)

    # NB the inner edges alone: num_labels labels, 0 to num_labels - 1, as plot_he colours them.
    labels = np.digitize(cropped_gray, bins=bins[1:-1])

    # TODO rename he_label.
    # NB 0.0 < x < 83_339.869; 0 < y < 53_009.165
    tissue_frame = tissue_frame.with_columns(
        pl.Series("gray", gray),
        pl.Series("cropped_gray", cropped_gray),
        pl.Series("label", labels),
    )

    logger.info(f"Merged with h&e with result:\n{tissue_frame}")

    # NB column by column: polars' to_pandas needs pyarrow, which cnamaste does not declare.
    return pd.DataFrame({name: tissue_frame[name].to_numpy() for name in tissue_frame.columns})


SLIDE_PIXELS = 6_000
"""Pixels along a slide panel's longer side, at most: Space Ranger's hires image whole (2,000 px a side; 6,000 on
Visium HD, 22.9M pixels), so `combined.pdf` embeds the slide at the resolution the run was given and zooms to it."""


def pixels_to_spots(pos):
    """The affine map from full-resolution pixels to a slice's spot coordinates, fitted on its spots.

    A tissue positions file gives each spot both: `x, y` the array's indices, `pixel_col, pixel_row` the pixels.

    Args: `pos` `io.get_spatial_positions`' frame.
    Returns: the (2, 3) matrix taking `(pixel_col, pixel_row, 1)` to `(x, y)`; None where `pos` has no pixel columns.
    """
    if pos is None or "pixel_col" not in pos.columns or "pixel_row" not in pos.columns:
        return None
    pixels = np.column_stack([pos["pixel_col"], pos["pixel_row"], np.ones(len(pos))]).astype(np.float64)
    spots = np.column_stack([pos["x"], pos["y"]]).astype(np.float64)
    return np.linalg.lstsq(pixels, spots, rcond=None)[0].T


def slide_image_and_meta(spaceranger_dirs, pos=None):
    """The slide drawn beside a run's spots: `get_he_image` of its one slice, at most `SLIDE_PIXELS` a side, one row
    per pixel carrying the image (`red`, `green`, `blue`) and what the drawing reads it by -- `array_row`,
    `array_col`, `(x, y)` on the spots' coordinates, and `gray`, `cropped_gray` and `label`, its H&E class.

    Space Ranger's hires image where present, else its lowres; each pixel's full-resolution place is taken to the spots'
    coordinates by the map the slice's own spots give (`pixels_to_spots`), and left in pixels where `pos` carries none.

    Args: `spaceranger_dirs` the run's slices' Space Ranger directories; `pos` the slice's `io.get_spatial_positions`.
    Returns: the per-pixel frame `layout.draw_slide` draws; None for several slices, or where the slice has no image.
    """
    dirs = list(dict.fromkeys(str(d) for d in spaceranger_dirs))
    if len(dirs) != 1:
        logger.info(f"No slide panel: {len(dirs)} slices, and a panel draws one.")
        return None
    res = "hires" if Path(f"{dirs[0]}/spatial/tissue_hires_image.png").exists() else "lowres"
    frame = get_he_image(dirs[0], res=res, pos=None, max_side=SLIDE_PIXELS)
    if frame is None:
        return None
    affine = pixels_to_spots(pos)
    if affine is not None:
        mapped = affine @ np.vstack([frame["x"], frame["y"], np.ones(len(frame))])
        frame["x"], frame["y"] = mapped[0], mapped[1]
    logger.info(
        f"Slide panel: {dirs[0]} {res}, {frame['array_row'].max() + 1} x {frame['array_col'].max() + 1} pixels."
    )
    return frame
