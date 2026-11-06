---
last-modified-date: 2025-10 (main guide: May 2025; Matplotlib vectorization addendum: October 2025)
author: Ben Raphael
contributors: Matthew Reyna, Richard Zhang, and other group members
converted-to-markdown: 2026-11-03
scope: Figure, axis, legend, caption and table rules extracted from "Tips for Writing Manuscripts"
---

# Figure, Axis and Table Rules

❌ = do not do this. ✅ = do this instead.
Rules marked **[added]** are not in the original guide. They follow its principles, and the group should confirm them.

---

## 1. Panels and Layout (FP)

- **FP-1.** Split multi-panel figures into **labeled panels**: a, b, c or A, B, C. Use the journal's default labeling style.
- **FP-2.** Do **not** put panels together unlabeled and refer to them as "Left", "Middle", "Right" in the caption.
- **FP-3.** Assemble multi-panel figures and add panel labels **manually** (PowerPoint, Illustrator, etc.). LaTeX `subfloat`/`subfigure` generally produce poor results.
- **FP-4.** Avoid **unnecessary arrows or callouts** for labeling figure elements.
- **FP-5. [added]** Lay a page out **tight**: size it to its content, the text column wide and a fixed height per
  row of panels, not a fixed page height the panels stretch or shrink to fill; let the layout engine place the
  panels with a fixed gap between them, wide enough that no panel's labels, ticks or letter reach another's; save it
  cropped to its content with a small margin.

```python
fig, axes = plt.subplots(rows, 2, figsize=(text_width, row_height * rows), layout="constrained")
fig.get_layout_engine().set(w_pad=0.1, h_pad=0.16)  # inches: gap left/right, and room for a panel letter
fig.savefig("figure.pdf", bbox_inches="tight", pad_inches=0.02)
```

---

## 2. Axes (FA)

- **FA-1. LABEL AXES IN PLOTS.** (Worth repeating: label axes in plots.)
- **FA-2.** Use **comically large fonts** for axis labels, tick labels and all other figure text, so they stay readable when the figure is shrunk to manuscript size. The default font size of most plotting programs is too small.
- **FA-3.** When showing a *y* = *x* relationship, use the **same scale on both axes**. Otherwise the diagonal becomes *y* = *mx* for an arbitrary slope *m* chosen by the plotting program.

---

## 3. Legends (FL)

- **FL-1.** Place legends so they are **not the first item a reader sees**. Readers go top to bottom and left to right, so legends do not go at the **top or left** of a panel.
- **FL-2. [added]** Legend text is **no larger than the axis labels**. An oversized legend draws the eye first and breaks the same principle as FL-1. The data and axes come first; the legend supports them.

---

## 4. Color and Style (FC)

- **FC-1.** Use a **consistent color scheme across all figures**.
  - Example: when several figures compare methods, each method keeps the **same color and line/marker style** in every figure.
- **FC-2. [added]** Order the methods **cnamaste, CalicoST, Numbat, Clonalscope**, left to right and top to bottom, in every figure: panels, bars, boxes, heatmap and confusion-map rows and columns, and keys.
  - A figure without one of the methods keeps the others in this order.
  - The runs of one method, e.g. a parameter sweep, stay together in that method's place.

---

## 5. Heatmaps (FH)

- **FH-1. Use heatmaps judiciously.** They are wildly overused. They give a **qualitative** view of structure in tabular data, not a quantitative one. Do not ask readers to draw quantitative conclusions from a heatmap.
- **FH-2. Sort rows/columns** so the heatmap shows structure, unless there is a biological or logical order that must be preserved.
- **FH-3. Use "light mode"** (zero = white). Most plotting programs default to "dark mode" (zero = dark), which suits a monitor but not a white manuscript page.

```python
import seaborn as sns

# Light-mode heatmap: zero maps to white
sns.heatmap(data, cmap="Blues", vmin=0)  # sequential data
sns.heatmap(data, cmap="RdBu_r", center=0)  # signed data: zero = white
```

---

## 6. Figure Captions and Legends (FT)

- **FT-1.** Describe **what the figure contains**, referring to specific parts, labels, colors, etc. as needed.
- **FT-2.** For multi-panel figures, write a description for **each subpanel**, matching its label ((a), (b), (c)) or location ((Left), (Right)).
- **FT-3.** Describe **only what is shown**: no interpretation and no speculation on causes or implications.
  - ❌ "Our method outperforms the other method."
- **FT-4.** Rarely or never use the word "figure" in a caption.
  - ❌ "This figure shows…"
  - ❌ "The left panel shows…"
- **FT-5.** Do not use the plot type as the description; describe what the plot shows.
  - ❌ "A boxplot of …", "A violin plot of …"
  - ✅ "Distribution of the …"

---

## 7. Referring to Figures in the Text (FR)

- **FR-1.** Write **"(Figure X)"**, not "(see Figure X)".
- **FR-2.** Keep the **story consistent** between the main text and supporting figures. At minimum, the **order of figure panels matches the order of sentences** in the main text.
- **FR-3.** In a reviewer response, **include the changed figures** in the response itself, so the reviewer does not need to open the manuscript.

---

## 8. File Format and Export (FX)

- **FX-1.** Figures must be **high-quality vector graphics (PDF/SVG)**. Do not create bitmaps by pasting graphics or converting Windows formats to PDF.
  - Check: zoom in on the figure. If text looks fuzzy, vectorization was lost.
- **FX-2.** Some older programs export SVG, EPS or PS (also vector formats). `pdflatex` prefers PDF, so convert these with command-line tools.
- **FX-3. Exception for very many points** (more than ~10K, e.g., scatter plots with millions of points): vector files become huge and slow to render. Use a high-resolution raster image (PNG, lossless) at **≥ 300 DPI**, preferably **≥ 600 DPI** unless you know the final print size.
- **FX-4.** Many journals require **fully vectorized camera-ready figures, including text as editable vector elements**. Matplotlib does not do this by default, even when saving to PDF. Always save as in MP-3 to avoid remaking figures.
- **FX-5.** Exported PDFs can be imported into **Adobe Illustrator or PowerPoint** with vectorization preserved.

---

## 9. Matplotlib / Seaborn Code (MP)

- **MP-1.** Match the LaTeX default font (Computer Modern) instead of Matplotlib's default DejaVu Sans:

```python
import matplotlib.pyplot as plt

plt.rc("font", family="serif", serif=["Computer Modern Roman"])
```

- **MP-2.** To write LaTeX formulas in titles, labels or legends (e.g., `r'$R^2$'`), enable LaTeX text rendering:

```python
plt.rc("text", usetex=True)

# With usetex=True, a literal dollar sign must be escaped: r'\$' or '\\$'
# (the backslash must escape itself in non-raw strings). Other LaTeX quirks also apply.
```

- **MP-3.** Save **fully vectorized** output, shapes **and** text, with at least 300 DPI:

```python
### Header of script or notebook
plt.rcParams["pdf.fonttype"] = 42  # editable text in PDF
plt.rcParams["ps.fonttype"] = 42  # same for PS/EPS
plt.rcParams["svg.fonttype"] = "none"  # keep text as text in SVG

### When saving the figure
plt.savefig(f"{filename}.pdf", dpi=300, bbox_inches="tight", transparent=True)
```

- **MP-4.** Raster fallback for more than ~10K points (FX-3):

```python
plt.savefig(f"{filename}.png", dpi=600, bbox_inches="tight")
```

- **MP-5. [added]** Set font sizes once, globally, so that figure text is large (FA-2) and legend text is no larger than axis labels (FL-2). Choose the values for your figure's final printed size.

```python
plt.rcParams.update(
    {
        "axes.labelsize": 20,  # axis labels: the largest text
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 16,  # no larger than axes.labelsize
        "legend.title_fontsize": 16,
        "axes.titlesize": 20,
    }
)
```

- **MP-6. [added]** Use the same axis scale for *y* = *x* comparisons (FA-3), and keep legends away from the top-left (FL-1):

```python
fig, ax = plt.subplots()
ax.scatter(x, y)
lims = [min(x.min(), y.min()), max(x.max(), y.max())]
ax.plot(lims, lims, color="gray", linestyle="--")  # y = x
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_aspect("equal")
ax.set_xlabel("True value")  # FA-1: always label axes
ax.set_ylabel("Inferred value")
ax.legend(loc="lower right")
```

- **MP-7. [added]** Keep method colors consistent across figures (FC-1) by defining them once and importing them everywhere:

```python
# colors.py — import in every plotting script
METHOD_COLORS = {
    "OurMethod": "#1f77b4",
    "MethodA": "#ff7f0e",
    "MethodB": "#2ca02c",
}
```

---

## 10. Figures in LaTeX (LF)

- **LF-1.** Do **not** use the `h` placement specifier; use `tbp`.
- **LF-2.** Do **not** specify the `.pdf` extension in `\includegraphics`.

```latex
\begin{figure}[tbp]
  \centering
  \includegraphics[width=0.5\textwidth]{filename}
  \caption{Distribution of …}
  \label{fig:example}
\end{figure}
```

- **LF-3.** LaTeX is sensitive to PDF **bounding boxes**. Crop PDFs before including them, using `pdfcrop`, Preview (Mac) or Illustrator.
- **LF-4.** Useful command-line tools (Linux/Mac; installed on department machines):

```bash
pdfcrop input.pdf output.pdf                                # auto-crop margins / fix bounding box
pdftk input_a.pdf input_b.pdf input_c.pdf cat output.pdf    # merge PDFs into one multipage PDF
convert input.pdf output.png                                # ImageMagick conversion (does NOT preserve vectors)
```

---

## 11. Tables (TB) [added]

The original guide has no table-specific rules. These apply its figure, caption and results rules to tables.

- **TB-1.** **Label every column and row**, with units, just as FA-1 requires labeled axes.
- **TB-2.** The table **caption describes what the table contains**, without interpretation (as in FT-3).
  - ❌ "This table shows…" (as in FT-4)
  - ❌ "Our method achieves the best accuracy."
  - ✅ "Adjusted Rand index (ARI) of each method on simulated data and breast cancer genomes."
- **TB-3.** Use the **same method names and order** as in the figures (as in FC-1 and FC-2).
- **TB-4.** Report **quantities, not adjectives**, with consistent precision within a column (as in the results rule "quantify your claims").
- **TB-5.** **Define every acronym** used in the table, in the table or its caption (as in the rule to write out acronyms at first use).
- **TB-6.** Prefer a table over a heatmap when readers need exact values (FH-1).
- **TB-7.** In LaTeX, use `booktabs` rules and no vertical lines:

```latex
\usepackage{booktabs}

\begin{table}[tbp]
  \centering
  \caption{Adjusted Rand index (ARI) of each method on simulated data.}
  \label{tab:ari}
  \begin{tabular}{lrr}
    \toprule
    Method    & ARI  & Runtime (min) \\
    \midrule
    OurMethod & 0.90 & 12 \\
    MethodA   & 0.60 & 45 \\
    \bottomrule
  \end{tabular}
\end{table}
```

---

## Quick Checklist

**Layout and axes**
- [ ] Panels labeled (a), (b), (c) in journal style; assembled outside LaTeX
- [ ] Every axis labeled; all text comically large
- [ ] Same scale on both axes for *y* = *x* comparisons
- [ ] No unnecessary arrows or callouts
- [ ] Page tight: text column wide, a fixed height per row, fixed gaps, cropped to content [added]

**Legends and color**
- [ ] Legends not at the top or left of a panel
- [ ] Legend text no larger than axis labels [added]
- [ ] Each method keeps the same color/style across all figures
- [ ] Methods ordered cnamaste, CalicoST, Numbat, Clonalscope, left to right and top to bottom [added]

**Heatmaps**
- [ ] Used only for qualitative structure
- [ ] Rows/columns sorted unless a natural order exists
- [ ] Light mode (zero = white)

**Captions and text**
- [ ] Each panel described; no interpretation
- [ ] No "This figure shows", no "The left panel shows", no "A boxplot of"
- [ ] "(Figure X)", not "(see Figure X)"; panel order matches text order

**Files**
- [ ] Vector PDF with editable text (`pdf.fonttype = 42`), ≥ 300 DPI
- [ ] PNG at ≥ 300–600 DPI only for > ~10K points
- [ ] PDFs cropped; `[tbp]`; no `.pdf` extension in `\includegraphics`

**Tables [added]**
- [ ] Columns labeled with units; acronyms defined
- [ ] Caption descriptive, not interpretive
- [ ] Method names and order match the figures; `booktabs`, no vertical lines
