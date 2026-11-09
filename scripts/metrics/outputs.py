"""Where the benchmark's scripts write: tables under `<outdir>/tsvs/`, gzipped, figures under `<outdir>/plots/`, each
page of a figure alone under `<outdir>/pdf-pages/`, LaTeX tables under `<outdir>/text/`.

- `get_table_dir`, `get_table_path` -- the tables' directory, and a table's path, `<outdir>/tsvs/<name>.tsv.gz`;
- `get_plot_dir`, `get_plot_path` -- the figures' directory, and a figure's path, `<outdir>/plots/<name>`;
- `get_page_path` -- a figure's page alone, `<outdir>/pdf-pages/<stem>-P<page>.pdf` or `<label>.pdf`;
- `get_text_path` -- a LaTeX table's path, `<outdir>/text/<name>.tex`;
- `get_json_path` -- a summary's path, `<outdir>/tsvs/<name>.json`, beside its table;
- `write_table` -- a table, tab-separated, gzipped with its header's time 0, so equal tables give equal bytes;
- `write_json` -- a summary, indented, keys in its order, NaN written `null`.

None creates a directory: `<outdir>/tsvs/`, `<outdir>/plots/`, `<outdir>/pdf-pages/` and `<outdir>/text/` exist
beforehand.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from cnamaste.config import start_time
from cnamaste.logger import get_logger

logger = get_logger("outputs", start_time=start_time)

TABLE_DIR = "tsvs"
"""The tables' directory under an `outdir`."""

PLOT_DIR = "plots"
"""The figures' directory under an `outdir`."""

PAGE_DIR = "pdf-pages"
"""The directory of each figure's pages, a PDF apiece, under an `outdir`."""

TEXT_DIR = "text"
"""The LaTeX tables' directory under an `outdir`."""


def get_table_dir(outdir: Path) -> Path:
    """Get the tables' directory, `<outdir>/tsvs`.

    Args:
        outdir: the analysis' `outdir`.

    Returns:
        The directory.
    """
    return Path(outdir) / TABLE_DIR


def get_table_path(outdir: Path, name: str) -> Path:
    """Get a table's path, `<outdir>/tsvs/<name>.tsv.gz`.

    Args:
        outdir: the analysis' `outdir`.
        name: the table's name, e.g. `simulation_metrics`.

    Returns:
        The path.
    """
    return get_table_dir(outdir) / f"{name}.tsv.gz"


def get_plot_dir(outdir: Path) -> Path:
    """Get the figures' directory, `<outdir>/plots`.

    Args:
        outdir: the analysis' `outdir`.

    Returns:
        The directory.
    """
    return Path(outdir) / PLOT_DIR


def get_plot_path(outdir: Path, name: str) -> Path:
    """Get a figure's path, `<outdir>/plots/<name>`.

    Args:
        outdir: the analysis' `outdir`.
        name: the figure's file name, e.g. `appendix.pdf`.

    Returns:
        The path.
    """
    return get_plot_dir(outdir) / name


def get_page_path(outdir: Path, name: str, page: int, label: str = "") -> Path:
    """Get a figure's page alone, `<outdir>/pdf-pages/<stem>-P<page>.pdf`, or `<label>.pdf` for a labelled page.

    Args:
        outdir: the analysis' `outdir`.
        name: the figure's file name, e.g. `combined.pdf`.
        page: the page, from 1.
        label: the page's own name, e.g. `main-fig-cnp-compare`, or empty for none.

    Returns:
        The path.
    """
    return Path(outdir) / PAGE_DIR / f"{label or f'{Path(name).stem}-P{page}'}.pdf"


def get_text_path(outdir: Path, name: str) -> Path:
    """Get a LaTeX table's path, `<outdir>/text/<name>.tex`.

    Args:
        outdir: the analysis' `outdir`.
        name: the tables' name, e.g. `table`.

    Returns:
        The path.
    """
    return Path(outdir) / TEXT_DIR / f"{name}.tex"


def get_json_path(outdir: Path, name: str) -> Path:
    """Get a summary's path, `<outdir>/tsvs/<name>.json`, beside its table.

    Args:
        outdir: the analysis' `outdir`.
        name: the summary's name, e.g. `ovarian`.

    Returns:
        The path.
    """
    return get_table_dir(outdir) / f"{name}.json"


def get_json_value(value: Any) -> Any:
    """Get a value as JSON writes it: a NumPy scalar as Python's, NaN or NA as None.

    Args:
        value: the value.

    Returns:
        The value.
    """
    if pd.isna(value):
        return None
    return value.item() if hasattr(value, "item") else value


def write_json(summary: dict[str, Any], path: Path) -> Path:
    """Write a summary as JSON, indented, keys in its order, logged.

    Args:
        summary: the summary, its values from `get_json_value`.
        path: a `.json` path, e.g. `get_json_path`'s.

    Returns:
        The path.
    """
    path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    logger.info(f"Wrote {path}.")
    return path


def write_table(table: pd.DataFrame, path: Path, **kwargs: Any) -> Path:
    """Write a table, tab-separated, gzipped with its header's time 0, logged.

    Args:
        table: the table, its index left out.
        path: a `.tsv.gz` path, e.g. `get_table_path`'s.
        kwargs: to `pandas.DataFrame.to_csv`, e.g. `float_format`.

    Returns:
        The path.
    """
    table.to_csv(path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0}, **kwargs)
    logger.info(f"Wrote {path}: {len(table)} rows, {table.shape[1]} columns.")
    return path
