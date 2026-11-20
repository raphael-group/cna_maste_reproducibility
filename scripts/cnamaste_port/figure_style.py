"""Ported from cna_maste `python/cnamaste/figure_style.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim, with their imports pointed at `scripts.cnamaste_port`.

One style for every figure cnamaste writes: the face, the text size, the page's width, the palettes.

`combined.png` (`layout.plot_combined_page`) is where these were settled, and this module is them taken out of it so
a figure written on its own is the same figure: the stated serif face at `[tool.cnamaste.figures]`'s `font_size`,
every text at that one size, a page the text column wide, at `DPI`, and clones in `CLONE_PALETTE` with H&E classes
in `HE_PALETTE`; a spatial scatter above `RASTERIZE_ABOVE` spots drawn as one raster (`rasterize`).

Style only. Nothing here decides what is drawn -- which bins, which segments, which states -- and a plotter that
adopts it draws exactly what it drew before, at the page's size and in the page's face.
"""

from __future__ import annotations

import contextlib
from typing import Any, Iterator

import matplotlib as mpl
from matplotlib.text import Text

DEFAULT = {"family": "serif", "font": "STIXGeneral", "mathtext": "stix", "font_size": 8.0}
"""`[tool.cnamaste.figures]` as shipped; used where no `pyproject.toml` is found. cna_maste's `pyproject.toml` at
42b96ff states these same values."""


PAPER_WIDTH = 8.5 - 2 * 1.01
"""Inches: a text column, US letter less 1.01 in margins, as port's `figure_style.PAPER_WIDTH`."""


TEXT_HEIGHT = 11.0 - 2 * 1.01 - 58.0 / 72.27
"""Inches: the text block less its running head and folio, as port's `figure_style.TEXT_HEIGHT`."""


CAPTION_ROOM = 1.0
"""Inches of the text block a page leaves for its caption."""


PAGES = {"third": 1 / 3, "half": 1 / 2, "full": 1.0}
"""The supported page heights, each a share of `TEXT_HEIGHT - CAPTION_ROOM`."""


DPI = 300


CLONE_PALETTE = "rocket"
"""The clones' palette, wherever clones are coloured."""


def style() -> dict[str, Any]:
    """`DEFAULT`, the `[tool.cnamaste.figures]` cna_maste ships.

    Ported change: cnamaste reads `[tool.cnamaste.figures]` from its own `pyproject.toml`; this repository has no
    such table, and cna_maste's at 42b96ff equals `DEFAULT`, so the port returns it.

    Returns: the figure style mapping.
    """
    return dict(DEFAULT)


@contextlib.contextmanager
def page_style() -> Iterator[float]:
    """Set the stated face for the block, over matplotlib's DejaVu Serif.

    Yields: the font size.
    """
    stated = style()
    family = stated["family"]
    with mpl.rc_context(
        {"font.family": family, f"font.{family}": [stated["font"]], "mathtext.fontset": stated["mathtext"]}
    ):
        yield float(stated["font_size"])


def page_size(page: str) -> tuple[float, float]:
    """Return a page's size in inches, width by height.

    Args: `page` "third", "half" or "full".
    Returns: `(width, height)`.
    """
    return PAPER_WIDTH, PAGES[page] * (TEXT_HEIGHT - CAPTION_ROOM)


def column_size(height: float) -> tuple[float, float]:
    """Inches, a figure the text column wide and `height` tall: what a standalone figure is drawn at.

    Args: `height` the figure's own height in inches.
    Returns: `(PAPER_WIDTH, height)`.

    A figure sized for a slide, 15 or 20 inches across, is a figure whose text is unreadable at the width it is
    printed. Its height is the plotter's own, so the rows it tuned keep their inches and only the width moves:
    scaling the height by the same factor as the width squeezes a row below what its own labels need, which
    matplotlib reports as `tight_layout cannot make Axes height small enough`.

    Panels side by side are the exception a caller states: there the height is what keeps their aspect, so the
    caller passes `PAPER_WIDTH / <the old ratio>` rather than the old height.
    """
    return PAPER_WIDTH, float(height)


def stagger(middles: Any, names: Any, size: float, gap: float = 1.0, lines: int = 2) -> list[int]:
    """Which line each label takes so a name does not meet its neighbour: the first line it fits on.

    `combined.png`'s rule (`layout.GenomicAxis.frame`), here so every genomic axis staggers the same way, and with
    the second line tracked as well as the first: the rule as it was let two names that had both been pushed down
    meet each other, which is the collision it exists to prevent. A chromosome name is about `0.6 * size` points
    per character, so a name is half that each side of its middle.

    Args: `middles` each label's centre in points along the axis; `names` the labels; `size` the text size in
        points; `gap` the points to leave between two names on one line; `lines` the lines available.
    Returns: one line index per label, `0 .. lines - 1`.
    """
    right = [-float("inf")] * lines
    taken = []
    for middle, name in zip(middles, names):
        half = 0.3 * size * len(str(name))
        fits = [k for k in range(lines) if middle - half >= right[k] + gap]
        # NB the first line it fits on; where it fits on none, the emptiest, which is the least it can overlap.
        line = fits[0] if fits else min(range(lines), key=lambda k: right[k])
        right[line] = middle + half
        taken.append(line)
    return taken


def uniform(figure: Any, size: float) -> None:
    """Set every text on the figure to `size` points.

    Args: `figure` the page; `size` in points.
    """
    for text in figure.findobj(Text):
        text.set_fontsize(size)


def restyle(figure: Any) -> None:
    """Set every text on the figure to the stated size, in the stated family; `write_fig` calls it on the way out.

    Args: `figure` the figure being written.

    The family is set per text, and generically ("serif"), rather than left to `rcParams`: a `Text` resolves its
    font when it is made, so a figure built while `rcParams` said DejaVu Serif keeps DejaVu Serif however the
    figure is saved. Generic, because the face itself comes from the `page_style` the save runs under, which is
    also what makes the mathtext match.
    """
    stated = style()
    for text in figure.findobj(Text):
        text.set_fontsize(stated["font_size"])
        text.set_fontfamily(stated["family"])
        # NB mathtext picks its fontset per text as well, and the same way: a $\log \mu$ laid out under the
        #    default set stays in DejaVu Sans beside STIX prose, which is the pair that reads as a mistake.
        text.set_math_fontfamily(stated["mathtext"])
