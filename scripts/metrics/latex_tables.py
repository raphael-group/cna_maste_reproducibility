"""The benchmark's LaTeX table: a dataset's runs' spot coverage rate, false-positive rate, runtime and memory, one row per
run, a booktabs tabular in a `table` float, in a minipage of `TABLE_WIDTH`, captioned and labelled (`CAPTION`,
`LABEL`), typeset beside it as a PDF (`write_latex_table`).

Dependencies: `pdflatex` with the `geometry`, `amsmath`, `amssymb` and `booktabs` packages (TeX Live). The document
that inputs the table needs `\\usepackage{booktabs}`.

- `get_table_columns`, `get_table_tex`: the metric columns, each name over its unit or its second line and its decimals,
  and the tabular, `$-$` where not recorded;
- `get_table_float`, `get_fragment_tex`: the tabular in its float, and the float after a comment naming it, for
  `\\input`;
- `get_article_tex`, `write_pdf`: a document typesetting a fragment on letter pages, and its PDF.

A row per run in method order, cnamaste, CalicoST, Numbat, Clonalscope (`plots.style.get_ordered_runs`): `Method`, its
method (`read_method_input.get_method_name`), and `Resolution`, the side in um of the bins it ran on, centred, its unit
under its name (`get_run_cells`).
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd
from scripts.metrics.outputs import get_table_path, get_text_path
from scripts.metrics.qc_metrics import get_false_positive_column
from scripts.metrics.read_method_input import get_method_name
from scripts.metrics.simulation_metrics import get_runtime_columns
from scripts.plots.style import get_ordered_runs

from cnamaste.config import start_time
from cnamaste.logger import get_logger

logger = get_logger("latex_tables", start_time=start_time)

MISSING = "$-$"
"""A cell with no recorded value."""

TEX_NOTE = r"% Written by scripts/analysis/tables.py; needs \usepackage{booktabs}."
"""The table's first line."""

LATEX_ENGINE = "pdflatex"
"""The engine a fragment's PDF is typeset with."""

TABLE_WIDTH = r"0.85\linewidth"
"""The table's minipage width."""

CAPTION = (
    "Spot coverage rate, false-positive rate, runtime and peak memory of each method on the Visium HD {dataset} cancer "
    "sample ({sample})."
)
"""The table's caption, its `{dataset}` and `{sample}` filled (`get_fragment_tex`)."""

LABEL = "tab:data_runtime"
"""The table's label."""


def get_table_columns() -> dict[str, tuple[str, str, int]]:
    """Get the table's metric columns, in order, each by its wide-table column: its name, the line under the name (its
    unit, or the name's second line) and the decimals printed: `spot-coverage-rate` and the false-positive column
    (`qc_metrics.get_false_positive_column`) to 0.001, runtime to 0.1 and memory to 0.01 (`get_runtime_columns`)."""
    runtime, memory = get_runtime_columns()
    return {
        "spot-coverage-rate": ("Spot coverage", "rate", 3),
        get_false_positive_column(): ("False-positive", "rate", 3),
        runtime: ("Runtime", "[mins]", 1),
        memory: ("Memory", "[GB]", 2),
    }


def get_value_cell(value: float, decimals: int) -> str:
    """Get a value's cell: the value to `decimals` places, a negative one's sign the math minus `$-$`, `MISSING` where
    NA.

    Args:
        value: the value.
        decimals: the decimals printed.

    Returns:
        The cell.
    """
    if pd.isna(value):
        return MISSING
    text = f"{value:.{decimals}f}"
    return f"$-${text[1:]}" if text.startswith("-") else text


def get_run_cells(method: str, resolution: str) -> list[str]:
    """Get a run's `Method` and `Resolution` cells: its method's name (`read_method_input.get_method_name`), and the
    side in um of the bins it ran on.

    Args:
        method: the run's method.
        resolution: the side of the bins it ran on, `<side>um`.

    Returns:
        The two cells.
    """
    return [get_method_name(method), str(int(str(resolution).removesuffix("um")))]


def get_table_tex(wide: pd.DataFrame) -> str:
    """Get the table as a booktabs tabular: `Method`, `Resolution`, then the metric columns (`get_table_columns`), each
    name over its unit or second line, then a row per run; every column but `Method` centred.

    Args:
        wide: the wide table, a row per run, `method`, `resolution` and the metric columns.

    Returns:
        The tabular, ending in a newline.
    """
    columns = get_table_columns()
    names, notes, _ = zip(*columns.values())
    lines = [
        r"\begin{tabular}{l" + "c" * (len(columns) + 1) + "}",
        r"\toprule",
        " & ".join(["Method", "Resolution", *names]) + r" \\",
        " & ".join(["", r"[$\mu$m]", *notes]) + r" \\",
        r"\midrule",
    ]
    for values in wide.to_dict("records"):
        cells = [get_value_cell(values[column], decimals) for column, (_, _, decimals) in columns.items()]
        lines.append(" & ".join([*get_run_cells(values["method"], values["resolution"]), *cells]) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def get_table_float(tabular: str, caption: str, label: str) -> str:
    """Get a tabular in its float: a `table` placed `[!h]`, centred, in a minipage of `TABLE_WIDTH`, the caption and
    label under the tabular.

    Args:
        tabular: the tabular, ending in a newline.
        caption: its caption.
        label: its label, e.g. `tab:data_runtime`.

    Returns:
        The float, ending in a newline.
    """
    return (
        "\n".join(
            [
                r"\begin{table}[!h]",
                r"\centering",
                rf"\begin{{minipage}}{{{TABLE_WIDTH}}}",
                tabular.rstrip("\n"),
                rf"\caption{{{caption}}}",
                rf"\label{{{label}}}",
                r"\end{minipage}",
                r"\end{table}",
            ]
        )
        + "\n"
    )


def get_fragment_tex(tabular: str, wide: pd.DataFrame, dataset: str) -> str:
    """Get the table as a fragment for `\\input`: `TEX_NOTE`, a comment naming the sample, then the tabular in its
    float (`get_table_float`), `CAPTION` and `LABEL`.

    Args:
        tabular: the tabular, ending in a newline.
        wide: the wide table, its `sample` in the comment and caption.
        dataset: the wide table's name, e.g. `ovarian`, in the caption.

    Returns:
        The fragment, ending in a newline.
    """
    sample = str(wide["sample"].iloc[0])
    caption = CAPTION.format(dataset=dataset, sample=sample)
    return "\n".join([TEX_NOTE, f"% {sample}", get_table_float(tabular, caption, LABEL)])


def get_article_tex(fragment: str) -> str:
    """Get a document typesetting a fragment on letter pages, with `amsmath`, `amssymb` and `booktabs`, lines allowed
    3 em of stretch so long words such as DOIs break without overfull lines. It writes no
    `.aux` (`\\nofiles`), and its PDF holds no date, trailer ID or pdfTeX banner, so equal text gives equal bytes.

    Args:
        fragment: the fragment's path, as `\\input` reads it.

    Returns:
        The document, one line.
    """
    return (
        r"\documentclass[11pt]{article}\usepackage[letterpaper,margin=1in]{geometry}"
        r"\usepackage{amsmath,amssymb,booktabs}\emergencystretch=3em\nofiles\pdfinfoomitdate=1\pdftrailerid{}"
        r"\pdfsuppressptexinfo=-1"
        rf"\begin{{document}}\input{{{fragment}}}\end{{document}}"
    )


def write_pdf(tex_path: Path, document: str) -> Path:
    """Typeset a fragment into a PDF beside it, `<stem>.pdf`, by `document`, with `LATEX_ENGINE` run in its directory;
    the engine's log is `<stem>.log` there.

    Args:
        tex_path: the fragment, `<stem>.tex`.
        document: the document inputting it, one line, e.g. `get_article_tex`'s.

    Returns:
        The PDF's path.

    Raises:
        RuntimeError: where the engine fails.
    """
    command = [
        LATEX_ENGINE,
        "-interaction=nonstopmode",
        "-halt-on-error",
        f"-jobname={tex_path.stem}",
        document,
    ]
    done = subprocess.run(command, cwd=tex_path.parent, capture_output=True, text=True, check=False)
    if done.returncode:
        raise RuntimeError(f"{LATEX_ENGINE} failed on {tex_path}: see {tex_path.with_suffix('.log')}")
    pdf = tex_path.with_suffix(".pdf")
    logger.info(f"Wrote {pdf}.")
    return pdf


def write_latex_table(outdir: Path, dataset: str, name: str) -> Path:
    """Read a dataset's wide table and write its table, `text/<name>.tex`, the tabular in its float, a row per run in
    method order (`plots.style.get_ordered_runs`), and its PDF, `text/<name>.pdf`, with pdflatex's log,
    `text/<name>.log` (`write_pdf`).

    Args:
        outdir: a benchmark run's directory, the wide table in its `tsvs/` (`outputs.get_table_path`).
        dataset: the wide table's name, e.g. `ovarian`.
        name: the table's name under `text/`, e.g. `table`.

    Returns:
        The fragment's path.
    """
    wide = get_ordered_runs(pd.read_csv(get_table_path(outdir, dataset), sep="\t"))
    path = get_text_path(outdir, name)
    path.write_text(get_fragment_tex(get_table_tex(wide), wide, dataset))
    logger.info(f"Wrote {path}: {len(wide)} runs.")
    write_pdf(path, get_article_tex(path.name))
    return path
