"""The benchmark's metrics as LaTeX text, `text/<name>.tex` and its PDF (`write_metric_text`): a paragraph per metric,
defined as `simulation_metrics` and `real_data_metrics` compute it, each followed by a `% shown in:` line naming the
pages and panels that draw it, a paper figure's page by its name (`plots.pages.PAGE_NAMES`), another by its PDF and
page, `diagnosis-P<k>`, as `analysis.figures` writes them, or `text/table.tex`.

- `get_metric_texts` -- each metric's paragraph: its name, its dataset, its table columns, its definition and the
  pages drawing it;
- `get_missing_columns` -- the metric columns the tables score that no paragraph defines;
- `get_metric_tex` -- the paragraphs as one fragment for `\\input`;
- `write_metric_text` -- the fragment and its PDF, refused where a column scored has no paragraph.
"""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

from scripts.cnamaste_port.config import start_time
from scripts.cnamaste_port.logger import get_logger
from scripts.metrics.latex_tables import get_article_tex, write_pdf
from scripts.metrics.outputs import get_text_path
from scripts.metrics.real_data_metrics import get_summary_metrics
from scripts.metrics.simulation_metrics import get_level_metrics, get_runtime_columns
from scripts.plots.pages import get_page_name

logger = get_logger("metric_text", start_time=start_time)

TEX_NOTE = r"% Written by scripts/analysis/tables.py; needs \usepackage{amsmath,amssymb}."
"""The fragment's first line."""

DATASETS = {"simulation": "Simulations", "real_data": "Real data"}
"""Each dataset's section, in order."""

SIMULATION_INTRO = r"""We score each run of a method on a simulation against the true clones and copy numbers of the
simulation.
Let $\mathcal{S}$ be the set of simulated spots, $\hat{\mathcal{S}} \subseteq \mathcal{S}$ the spots the run assigns to
a clone, normal or tumor, and $\ell_i$ and $\hat{\ell}_i$ the true clone of spot $i$ and the clone the run assigns to
it, from \texttt{clone\_labels\_int.tsv.gz} for cnamaste and \texttt{clone\_labels.tsv} for CalicoST.
Let $\mathcal{G}$ be the set of genes inside the autosomal chromosome arms, and
$\hat{\mathcal{G}} \subseteq \mathcal{G}$ the genes at which the run calls a CNA state in every clone.
A cell $(i,g)$ is a spot $i \in \hat{\mathcal{S}}$ at a gene $g \in \hat{\mathcal{G}}$ that a true segment overlaps.
Let $k_{i,g}$ be the true allele-specific copy number $(A,B)$ of clone $\ell_i$ on the true segment overlapping gene $g$
most, and $\hat{k}_{i,g}$ the one the run calls for clone $\hat{\ell}_i$ on its segment overlapping gene $g$ most, with
$(A,B)$ and $(B,A)$ equal.
Let $s_{i,g}$ and $\hat{s}_{i,g}$ be their CNA states: \texttt{neu}, \texttt{bamp} or \texttt{bdel} where $A = B$ and
$A + B$ is $2$, above $2$ or below $2$, and \texttt{cnloh}, \texttt{amp} or \texttt{del} where $A \neq B$ and $A + B$ is
$2$, above $2$ or below $2$.
Numbat calls $\hat{s}_{i,g}$ itself, and Clonalscope calls a total copy number alone, read as \texttt{del},
\texttt{neu} or \texttt{amp} where it is below, at or above $2$.
Let $P = \{(i,g) : s_{i,g} \neq \texttt{neu}\}$ be the truly altered cells and $Q = \{(i,g) : s_{i,g} = \texttt{neu}\}$
the truly diploid cells, \texttt{neu} being $(1,1)$ alone.
Each rate is NaN where its denominator is $0$.
The figures draw the metrics at the genes; the table also holds them on the segments of the run, a segment in place of a
gene, \texttt{<metric> (segment)} and \texttt{<metric> (gene)}."""
"""The simulations' section's opening paragraph: the spots, genes, cells and CNA states the metrics count."""

REAL_DATA_INTRO = (
    r"A sample's runs are compared without a truth, on UGP's in-tissue spots at 16~$\mu$m: a run on 32~$\mu$m bins "
    r"gives each bin's clone to the 16~$\mu$m spots inside it, so every run is counted on the same spots. A run's "
    r"called spots are the in-tissue spots it labels, its tumor spots those of a clone other than the normal clone 0, "
    r"Clonalscope's normal clusters merged into clone 0. A cell is a called spot at a gene, at its clone's CNA state, "
    r"the state of the clone's segment overlapping the gene most, \texttt{uncalled} where none does. Two runs are "
    r"compared with the reference run, cnamaste, at the finest field both give: $(A,B)$ phase-free where both give "
    r"pairs (cnamaste, CalicoST), the state of total copies where either gives total copies alone (Clonalscope), "
    r"else the CNA state (Numbat)."
)
"""The real-data section's opening paragraph."""

MSE_DEFINITION = r"""We quantify the copy-number error of a method by the mean squared error (MSE) between the
read-depth ratio (RDR) and B-allele frequency (BAF) observed in its clones and those expected from its copy numbers, and
we report each MSE as a skill score (SS) against the diploid prediction.
We score the methods that report copy numbers: Clonalscope reports total copy numbers alone and is scored on RDR only,
and Numbat reports no copy numbers and is not scored.
We define the bins by the adaptive segmentation of the Universal Genotyping Pipeline (UGP) with
\texttt{min\_snp\_reads}~$= 300$, which merges adjacent genomic intervals until each bin holds at least $300$
SNP-covering UMIs over all spots of the slice, and we retain the bins inside the autosomal chromosome arms.
We score each method on the in-tissue $16\,\mu$m spots it assigns to a clone; a method run on $32\,\mu$m bins assigns
the clone of each bin to the $16\,\mu$m spots inside it.

Let $\hat{\mathcal{C}}_m$ be the set of spots a method assigns to clone $m$, $N = \sum_m |\hat{\mathcal{C}}_m|$, and
$w_m = |\hat{\mathcal{C}}_m| / N$ the weight of clone $m$.
For spot $i$ and bin $g$, let $x_{g,i}$ be the UMI count, $b_{g,i}$ the phased B-allele count and $n_{g,i}$ the total
allele count, and let $T_i = \sum_g x_{g,i}$ be the UMI count of spot $i$ over the retained bins.
Let $R$ be the set of spots whose broad cell type is neither epithelial nor low confidence, and
$\hat{\lambda}_g = \sum_{i \in R} x_{g,i} / \sum_{i \in R} T_i$ the baseline proportion of bin $g$, so that
$\sum_g \hat{\lambda}_g = 1$.
The observed RDR and BAF of clone $m$ at bin $g$ are those of the pseudobulk of its spots,
\begin{equation}
\mathrm{RDR}_{g,m} = \frac{\sum_{i \in \hat{\mathcal{C}}_m} x_{g,i}}
{\hat{\lambda}_g \sum_{i \in \hat{\mathcal{C}}_m} T_i},
\qquad
\mathrm{BAF}_{g,m} = \frac{\sum_{i \in \hat{\mathcal{C}}_m} b_{g,i}}{\sum_{i \in \hat{\mathcal{C}}_m} n_{g,i}}.
\end{equation}
Let $\hat{A}_{g,m}$ and $\hat{B}_{g,m}$ be the allele-specific copy numbers a method assigns to clone $m$ at bin $g$,
and $\hat{\mu}_{g,m} = \hat{A}_{g,m} + \hat{B}_{g,m}$ the total copy number; bin $g$ takes the copy numbers of the
segment that overlaps it most in the common refinement of the segments of all scored methods.
The expected RDR and BAF of clone $m$ at bin $g$ are
\begin{equation}
\widehat{\mathrm{RDR}}_{g,m} = \frac{\hat{\mu}_{g,m}}{\sum_{g'} \hat{\lambda}_{g'} \hat{\mu}_{g',m}},
\qquad
\widehat{\mathrm{BAF}}_{g,m} = \frac{\hat{B}_{g,m}}{\hat{\mu}_{g,m}},
\end{equation}
where the sum runs over the bins at which the method assigns copy numbers to clone $m$.
As the assignment of the two haplotypes to the alleles $A$ and $B$ by a method need not match the phasing of the
counts, we orient the expected BAF of all clones jointly at each bin $g$: $\widetilde{\mathrm{BAF}}_{g,m}$ is
$\widehat{\mathrm{BAF}}_{g,m}$ for every clone $m$, or $1 - \widehat{\mathrm{BAF}}_{g,m}$ for every clone $m$,
whichever gives the smaller $\sum_m w_m (\mathrm{BAF}_{g,m} - \widetilde{\mathrm{BAF}}_{g,m})^2$.

A bin is scored on RDR (BAF) where every method scored on RDR (BAF) assigns copy numbers to every clone with finite
observed and expected RDR (BAF), and where at least one scored method assigns a CNA to at least one clone: an
allele-specific copy number other than $(1,1)$ or, for Clonalscope, a total copy number other than $2$.
Let $G_{\mathrm{RDR}}$ and $G_{\mathrm{BAF}}$ be the numbers of bins scored on RDR and BAF.
A bin at which every method assigns $(1,1)$ to every clone is not scored, as every method makes the diploid prediction
there.
We quantify the errors of the total and the allele-specific copy numbers by
\begin{equation}
\mathrm{MSE}(\mathrm{RDR}) = \sum_{g,m} \frac{w_m}{G_{\mathrm{RDR}}}
\left( \log_2 \frac{\mathrm{RDR}_{g,m}}{\widehat{\mathrm{RDR}}_{g,m}} \right)^2,
\quad
\mathrm{MSE}(\mathrm{BAF}) = \sum_{g,m} \frac{w_m}{G_{\mathrm{BAF}}}
(\mathrm{BAF}_{g,m} - \widetilde{\mathrm{BAF}}_{g,m})^2,
\end{equation}
where each sum runs over the bins scored on the channel and the clones of the method."""
"""The mean squared error's paragraphs: the spots, the bins, the observed and expected RDR and BAF, and
`cn_error.score_rdr` and `cn_error.score_baf`."""

SKILL_SCORE_DEFINITION = r"""Let $\mathrm{MSE}_{(1,1)}(\mathrm{RDR})$ and $\mathrm{MSE}_{(1,1)}(\mathrm{BAF})$ be the
errors, on the same bins and clones, of the diploid prediction, which assigns $(1,1)$ to every clone at every bin, so
that $\widehat{\mathrm{RDR}}_{g,m} = 1$ and $\widetilde{\mathrm{BAF}}_{g,m} = 1/2$.
We report the skill scores (Murphy 1988,
doi:\allowbreak 10.1175/\allowbreak 1520-0493(1988)116$<$2417:SSBOTM$>$2.0.CO;2)
\begin{equation}
\mathrm{SS}(\mathrm{RDR}) = 1 - \frac{\mathrm{MSE}(\mathrm{RDR})}{\mathrm{MSE}_{(1,1)}(\mathrm{RDR})},
\qquad
\mathrm{SS}(\mathrm{BAF}) = 1 - \frac{\mathrm{MSE}(\mathrm{BAF})}{\mathrm{MSE}_{(1,1)}(\mathrm{BAF})}.
\end{equation}
A skill score of $1$ is an exact reconstruction, $0$ the error of the diploid prediction, and a negative score an error
larger than that of the diploid prediction."""
"""The skill score's paragraph, `cn_error.get_skill_score`."""


class MetricText(NamedTuple):
    """A metric's paragraph: its `name`, its `dataset`, `simulation` or `real_data`, the table `columns` it is, its
    LaTeX `definition`, and the pages `shown` drawing it, `<page> (<panels>)`."""

    name: str
    dataset: str
    columns: tuple[str, ...]
    definition: str
    shown: tuple[str, ...]


def get_metric_texts() -> tuple[MetricText, ...]:
    """Get each metric's paragraph, in each dataset's order: the simulations' metrics as their page draws them, then
    the real data's, from the spots to the copies.

    Returns:
        The paragraphs.
    """
    simulation, ablation = get_page_name("simulation"), get_page_name("ablation")
    # NB diagnosis.pdf's simulation page in full, every metric and simulated class.
    diagnosis = "diagnosis-P6"
    profiles = get_page_name("profiles")
    return (
        MetricText(
            "Normal specificity",
            "simulation",
            ("specificity",),
            r"The share of truly diploid cells the run calls diploid, "
            r"$|\{(i,g) \in Q : \hat{s}_{i,g} = \texttt{neu}\}| / |Q|$.",
            (f"{diagnosis} (a)", f"{ablation} (b)"),
        ),
        MetricText(
            "Spot coverage rate",
            "simulation",
            ("spot-coverage-rate",),
            r"The share of simulated spots the run assigns to a clone, $|\hat{\mathcal{S}}| / |\mathcal{S}|$, "
            r"\#called / \#spots: the spots its clone ARI scores.",
            (f"{diagnosis} (b)",),
        ),
        MetricText(
            "CNA sensitivity",
            "simulation",
            ("sensitivity",),
            r"The share of truly altered cells the run calls altered, "
            r"$|\{(i,g) \in P : \hat{s}_{i,g} \neq \texttt{neu}\}| / |P|$.",
            (f"{simulation} (a)", f"{diagnosis} (c)", f"{ablation} (a)"),
        ),
        MetricText(
            "Clone ARI",
            "simulation",
            ("clone_ari",),
            r"The adjusted Rand index (Hubert and Arabie 1985, doi:\allowbreak 10.1007/\allowbreak BF01908075) of the "
            r"true clones $\ell_i$ and the clones $\hat{\ell}_i$ of the run, over the spots $i \in \hat{\mathcal{S}}$.",
            (f"{simulation} (b)", f"{diagnosis} (d)", f"{ablation} (d)"),
        ),
        MetricText(
            "CNA recall",
            "simulation",
            ("state_altered",),
            r"The share of truly altered cells the run calls in their true CNA state, "
            r"$|\{(i,g) \in P : \hat{s}_{i,g} = s_{i,g}\}| / |P|$.",
            (f"{simulation} (c)", f"{diagnosis} (e)"),
        ),
        MetricText(
            r"$\mathbb{N}$-CNA recall",
            "simulation",
            ("exact_altered",),
            r"The share of truly altered cells the run calls at their true allele-specific copy number, "
            r"$|\{(i,g) \in P : \hat{k}_{i,g} = k_{i,g}\}| / |P|$; NaN for Numbat and Clonalscope, which call no "
            r"allele-specific copy numbers.",
            (f"{simulation} (d)", f"{diagnosis} (f)", f"{ablation} (c)"),
        ),
        MetricText(
            "Runtime and memory",
            "simulation",
            get_runtime_columns(),
            r"A run's wall-clock minutes and peak memory in GB, from \texttt{simulation.runtime.tsv}: Numbat's and "
            r"Clonalscope's run rule's Snakemake benchmark, its \texttt{s} / 60 and \texttt{max\_pss} / 1024; "
            r"cnamaste's its log's \texttt{Done in}, its memory not recorded.",
            (f"{get_page_name('runtime')} (a), (b)",),
        ),
        MetricText(
            "Spot coverage rate",
            "real_data",
            ("spot-coverage-rate",),
            r"The run's called spots over the in-tissue spots, \#called / \#in-tissue.",
            (
                "text/table.tex",
                (
                    f"{get_page_name('labels')} (c)--(f), each clone's share of the in-tissue spots in its key, the "
                    "shares summing to it"
                ),
                f"{profiles} (a), each profile row's share of the in-tissue spots",
            ),
        ),
        MetricText(
            "Genome composition",
            "real_data",
            ("neu", "del", "amp", "cnloh", "bamp", "bdel", "uncalled"),
            r"The share of the run's cells in each CNA state and \texttt{uncalled}, summing to 1: "
            r"$\sum_k n_k G_{k,s} / (G \sum_k n_k)$, $n_k$ the in-tissue spots of clone $k$, $G_{k,s}$ the genes "
            r"clone $k$ calls in state $s$, $G$ the genes.",
            (f"{get_page_name('composition')} (a), del, amp, cnloh and bamp",),
        ),
        MetricText(
            "False-positive rate",
            "real_data",
            ("false-positive-rate",),
            r"Let $R$ be the set of in-tissue spots whose broad cell type is neither epithelial nor low confidence, "
            r"$\hat{\mathcal{S}}$ the in-tissue spots the run assigns to a clone, and "
            r"$\hat{\mathcal{T}} \subseteq \hat{\mathcal{S}}$ those it assigns to a tumor clone, any clone but the "
            r"normal clone 0. The cell types give normal truth only, so the false-positive tumor spots are "
            r"$\hat{\mathcal{T}} \cap R$. The false-positive rate is the share of the called spots of $R$ that the "
            r"run assigns to a tumor clone, $|\hat{\mathcal{T}} \cap R| / |\hat{\mathcal{S}} \cap R|$: the false "
            r"positives over the negatives, $\mathrm{FP} / (\mathrm{FP} + \mathrm{TN})$ (James et al., An "
            r"Introduction to Statistical Learning, 1st ed., sec. 4.4.3, Table 4.7, doi:\allowbreak "
            r"10.1007/\allowbreak 978-1-4614-7138-7); NaN where $\hat{\mathcal{S}} \cap R$ is empty. An epithelial, "
            r"low-confidence or untyped spot is in neither count.",
            (
                "text/table.tex",
                "diagnosis-P1 (a)--(c), the false-positive spots of Numbat, a CalicoST sweep run and cnamaste",
                "diagnosis-P2 (a)--(c), their total and SNP-covering UMI against the other called spots'",
            ),
        ),
        MetricText(
            "Clone ARI",
            "real_data",
            ("clone_ari",),
            r"The adjusted Rand index of two runs' clones over the in-tissue spots both label, Clonalscope's clusters "
            r"as called; the table's is each run's against the reference run's.",
            (f"{profiles} (c), every two runs",),
        ),
        MetricText(
            "CNA ARI",
            "real_data",
            ("cna_ari",),
            r"The adjusted Rand index of two runs' states over the cells both call a state on, at their field: that of "
            r"their confusion, the cells by (state of the reference, state of the run), as a contingency table, "
            r"\texttt{uncalled} left out.",
            ("no figure; tsvs/<dataset>.tsv.gz and .json",),
        ),
        MetricText(
            "State Jaccard index",
            "real_data",
            (),
            r"Per state $s$ of the reference and $t$ of the run, at their field, over the cells both call: the cells "
            r"where the reference calls $s$ and the run $t$, over those where the reference calls $s$ or the run $t$ "
            r"(Jaccard 1912, doi:\allowbreak 10.1111/\allowbreak j.1469-8137.1912.tb05611.x); NaN where neither state holds such a cell.",
            (f"{get_page_name('state_jaccard')} (a)--(c), a panel per run, named by its field",),
        ),
        MetricText(
            "CNP similarity",
            "real_data",
            (),
            r"Of two runs' tumor clones $i$ and $j$, the altered genes $\mathcal{A}_{ij}$ are those both call a state "
            r"on and at least one calls other than \texttt{neu}; $J_{ij}$ is the share of $\mathcal{A}_{ij}$ where "
            r"their CNA states match, and the runs' CNP similarity is $S(a,b) = \max_{i,j} J_{ij}$, NaN where no pair "
            r"has an altered gene.",
            (f"{profiles} (b), every two runs",),
        ),
        MetricText(
            "Mean squared error",
            "real_data",
            ("mse_rdr", "mse_baf"),
            MSE_DEFINITION,
            ("diagnosis-P3 (a)", "diagnosis-P4 (a), the CalicoST sweep"),
        ),
        MetricText(
            "Skill score",
            "real_data",
            ("ss_rdr", "ss_baf"),
            SKILL_SCORE_DEFINITION,
            (f"{profiles} (d)", "diagnosis-P3 (b)", "diagnosis-P4 (b), the CalicoST sweep"),
        ),
        MetricText(
            "Runtime and memory",
            "real_data",
            get_runtime_columns(),
            r"A run's wall-clock minutes and peak memory in GB, its row of the sample's runtime table, by its method, "
            r"the side of its bins and its setting.",
            ("text/table.tex",),
        ),
    )


def get_missing_columns(texts: tuple[MetricText, ...]) -> list[str]:
    """Get the metric columns the tables score that no paragraph of their dataset defines: the simulations' spot
    coverage rate, per-level metrics and runtime (`simulation_metrics.get_level_metrics`), the real data's summary
    metrics (`real_data_metrics.get_summary_metrics`).

    Args:
        texts: the paragraphs.

    Returns:
        Each such column, `<dataset>: <column>`.
    """
    scored = {
        "simulation": ("spot-coverage-rate", *get_level_metrics(), *get_runtime_columns()),
        "real_data": get_summary_metrics(),
    }
    defined = {(text.dataset, column) for text in texts for column in text.columns}
    return [
        f"{dataset}: {column}"
        for dataset, columns in scored.items()
        for column in columns
        if (dataset, column) not in defined
    ]


def get_metric_tex(texts: tuple[MetricText, ...]) -> str:
    """Get the paragraphs as one fragment for `\\input`: `TEX_NOTE`, then a section per dataset, its opening paragraph
    and a paragraph per metric, its columns after its name, then a `% shown in:` line per page drawing it.

    Args:
        texts: the paragraphs.

    Returns:
        The fragment, ending in a newline.
    """
    intros = {"simulation": SIMULATION_INTRO, "real_data": REAL_DATA_INTRO}
    lines = [TEX_NOTE]
    for dataset, section in DATASETS.items():
        lines += ["", rf"\subsection*{{{section}}}", intros[dataset]]
        for text in (text for text in texts if text.dataset == dataset):
            columns = ", ".join(
                rf"\texttt{{{column.replace('_', chr(92) + '_').replace('#', chr(92) + '#')}}}"
                for column in text.columns
            )
            lines += ["", rf"\paragraph{{{text.name}.}}" + (f" ({columns})" if columns else ""), text.definition]
            lines += [f"% shown in: {page}" for page in text.shown]
    return "\n".join(lines) + "\n"


def write_metric_text(outdir: Path, name: str) -> Path:
    """Write the metrics' paragraphs, `text/<name>.tex` (`get_metric_tex`), and their PDF, `text/<name>.pdf`, with
    pdflatex's log, `text/<name>.log` (`latex_tables.write_pdf`).

    Args:
        outdir: a benchmark run's directory.
        name: the text's name under `text/`, e.g. `metric`.

    Returns:
        The fragment's path.

    Raises:
        ValueError: where a metric column the tables score has no paragraph (`get_missing_columns`).
    """
    texts = get_metric_texts()
    missing = get_missing_columns(texts)
    if missing:
        raise ValueError(f"no metric text defines {missing}: add it to metric_text.get_metric_texts")
    path = get_text_path(outdir, name)
    path.write_text(get_metric_tex(texts))
    logger.info(f"Wrote {path}: {len(texts)} metrics.")
    write_pdf(path, get_article_tex(path.name))
    return path
