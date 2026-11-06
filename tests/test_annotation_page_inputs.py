"""The inputs `scripts/run_benchmark.py` converts: Visium HD bins given to their finer
bins and placed on the image (`platforms.visium`), a Numbat run's clones and clone profile
(`scripts/metrics/read_method_input.py`), spots drawn as lattice cells (`layout.draw_spot_image`), and
segments on cnplot's arms as `(A,B)` states in cnamaste's colours (`layout.get_cnplot_segments`, `layout.get_copy_states`).

Referees: a 32 um bin's 16 um bins are the 2 x 2 block at twice its indices, their centres a quarter of its side
either way of its own; the grid ends at 3,350 2 um squares, so a bin on its far edge holds one; Numbat's documented
outputs (https://kharchenkolab.github.io/numbat/articles/descriptions.html), its clone 1 the normal cells, and
its clone-level pseudobulk profiles.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from scripts.metrics.read_method_input import (
    get_calicost_ploidy,
    read_numbat_labels,
    read_numbat_profile,
)
from scripts.platforms.visium import get_visiumhd_bin_pixels, get_visiumhd_finer_labels


def test_finer_bins_tile_their_bin_and_end_with_the_grid() -> None:
    """Three 32 um bins to 16 um: each label on its own block, each 16 um centre a quarter side from its bin's, under a
    turned, scaled and shifted image as Space Ranger fits one; the bin at (209, 209), past which the 16 um grid stops
    at index 418, keeps one."""
    labels = pd.Series([1, 2, 0], index=["s_032um_00000_00000-1", "s_032um_00005_00007-1", "s_032um_00209_00209-1"])
    fine = get_visiumhd_finer_labels(labels, 16)

    assert fine.groupby(fine).size().to_dict() == {0: 1, 1: 4, 2: 4}
    assert set(fine[fine == 2].index) == {f"s_016um_{r:05d}_{c:05d}-1" for r in (10, 11) for c in (14, 15)}

    turn = np.radians(0.4)
    image = 3.65 * np.array([[np.cos(turn), -np.sin(turn)], [np.sin(turn), np.cos(turn)]])
    grid = [f"s_032um_{r:05d}_{c:05d}-1" for r in range(0, 200, 37) for c in range(0, 200, 41)]
    um = np.array([[(int(b[14:19]) + 0.5) * 32, (int(b[8:13]) + 0.5) * 32] for b in grid])
    xy = um @ image.T + [1_200.0, 5_600.0]
    positions = pd.DataFrame({"barcode": grid, "x": xy[:, 0], "y": xy[:, 1]})

    coarse = get_visiumhd_bin_pixels(labels.index[:2], positions)
    within = get_visiumhd_bin_pixels(fine[fine != 0].index, positions).reshape(2, 4, 2)
    assert np.allclose(within.mean(axis=1), coarse)
    assert np.allclose(np.abs(np.linalg.solve(image, (within - coarse[:, None]).reshape(-1, 2).T)), 8.0)


def test_numbat_clone_profile_is_each_clones_pseudobulk_state(tmp_path: Path) -> None:
    """A run of two iterations: labels read at the last, clones shifted to cnamaste's numbering, the normal clone 0; the
    profile each clone's `cnv_state` per consensus segment of `bulk_clones_final.tsv.gz`, one row per SNP in the file,
    one per segment out, in genome order."""
    pd.DataFrame({"cell": ["a"], "clone_opt": [1]}).to_csv(tmp_path / "clone_post_1.tsv", sep="\t")
    pd.DataFrame({"cell": ["a", "b", "c", "d"], "clone_opt": [1, 2, 3, 3]}).to_csv(
        tmp_path / "clone_post_2.tsv", sep="\t", index=False
    )
    states = {1: ["neu", "neu"], 2: ["neu", "amp"], 3: ["loh", "amp"]}
    pd.DataFrame(
        [
            {"sample": k, "CHROM": chrom, "seg_cons": seg, "seg_start": 100, "seg_end": 200, "cnv_state": s}
            for k, row in states.items()
            for (chrom, seg), s in zip(((2, "2a"), (1, "1a")), row)
            for _ in range(3)
        ]
    ).to_csv(tmp_path / "bulk_clones_final.tsv.gz", sep="\t", index=False)

    assert read_numbat_labels(tmp_path).to_dict() == {"a": 0, "b": 1, "c": 2, "d": 2}
    profile = read_numbat_profile(tmp_path)
    assert profile.to_dict("list") == {
        "CHR": [1, 2], "START": [100, 100], "END": [200, 200],
        "clone0": ["neu", "neu"], "clone1": ["amp", "neu"], "clone2": ["amp", "loh"],
    }  # fmt: skip


def test_segments_are_cut_at_the_arms(tmp_path: Path) -> None:
    """A segment through chr1's centromere, on two arm regions, is drawn either side of it and nowhere in it; clone 0,
    (1, 1) throughout, is not drawn; each state is in cnamaste's profile colour, a B > A state its mirror's, lighter,
    (1,1) white; the key groups the states by BAF, ascending, a state with its phase flip, (0,0) first."""
    pytest.importorskip("cnplot", reason="cnplot is in the luo dependency group: uv sync --group luo")
    import matplotlib.colors as mcolors
    from scripts.plots.layout import (
        get_cnplot_segments,
        get_copy_state_colours,
        get_copy_state_groups,
        get_copy_states,
        get_genome_axis,
        get_profile_clones,
    )

    from cnamaste.palette import COPY_STATES

    (tmp_path / "regions.bed").write_text("chr1\t0\t100\tchr1p\nchr1\t150\t300\tchr1q\n")
    (tmp_path / "sizes").write_text("chr1\t300\n")
    axis = get_genome_axis(str(tmp_path / "regions.bed"), str(tmp_path / "sizes"))
    seglevel = pd.DataFrame(
        {"CHR": [1, 1], "START": [20, 200], "END": [180, 280], "clone0 A": [1, 1], "clone0 B": [1, 1],
         "clone1 A": [2, 1], "clone1 B": [1, 0], "clone2 A": [1, 1], "clone2 B": [2, 1]}
    )  # fmt: skip

    segments = get_cnplot_segments(seglevel, axis)
    states = get_copy_states(segments)
    clones = get_profile_clones(states, integer=True)
    colours = get_copy_state_colours([segments], [clones])

    assert states[["START", "END"]].to_numpy().tolist() == [[20, 100], [150, 180], [200, 280]]
    assert states["clone1"].tolist() == ["(2,1)", "(2,1)", "(1,0)"]
    assert clones == [1, 2]
    assert list(colours) == ["(1,1)", "(1,0)", "(2,1)", "(1,2)"]
    drawn = [
        "(1,1)",
        "(0,0)",
        "(1,0)",
        "(2,0)",
        "(0,2)",
        "(3,0)",
        "(4,1)",
        "(1,3)",
        "(2,1)",
        "(1,2)",
        "(4,2)",
        "(3,2)",
        "(2,3)",
        "(2,2)",
    ]
    assert [list(group) for group in get_copy_state_groups(dict.fromkeys(drawn, "")).values()] == [
        ["(0,0)"],
        ["(1,0)", "(2,0)", "(0,2)", "(3,0)"],
        ["(4,1)"],
        ["(1,3)"],
        ["(2,1)", "(1,2)", "(4,2)"],
        ["(3,2)", "(2,3)"],
        ["(1,1)", "(2,2)"],
    ]
    assert colours["(2,1)"] == COPY_STATES[(2, 1)].lower()
    assert colours["(1,1)"] == "#ffffff"
    flip = mcolors.to_rgb(colours["(1,2)"])
    assert flip == pytest.approx(tuple(0.5 * v + 0.5 for v in mcolors.to_rgb(COPY_STATES[(2, 1)])), abs=1 / 255)


def test_spots_are_one_cell_each_whatever_their_order() -> None:
    """A 3 x 4 lattice turned 0.4 degrees, 16 um pitch, two spots missing: each spot its own cell in its own colour,
    stored as a square of equal pixels, the missing cells clear, and the image the same for the spots in reverse
    order."""
    pytest.importorskip("cnplot", reason="cnplot is in the luo dependency group: uv sync --group luo")
    from scripts.plots.layout import draw_spot_image

    turn = np.radians(0.4)
    rotation = np.array([[np.cos(turn), -np.sin(turn)], [np.sin(turn), np.cos(turn)]])
    grid = np.array([(c, r) for r in range(3) for c in range(4) if (c, r) not in ((1, 1), (3, 2))])
    xy = (16.0 * grid) @ rotation.T + [500.0, -2_000.0]
    rgba = np.random.default_rng(0).uniform(size=(len(grid), 4))

    images = []
    for order in (np.arange(len(grid)), np.arange(len(grid))[::-1]):
        figure, ax = plt.subplots()
        draw_spot_image(ax, xy[order], rgba[order])
        stored = np.asarray(ax.images[0].get_array())
        plt.close(figure)
        repeat = stored.shape[0] // 3
        assert stored.shape == (3 * repeat, 4 * repeat, 4)
        cells = stored[::repeat, ::repeat]
        assert np.array_equal(stored, np.repeat(np.repeat(cells, repeat, axis=0), repeat, axis=1))
        images.append(cells)

    assert images[0].shape == (3, 4, 4)
    assert np.array_equal(images[0], images[1])
    assert np.array_equal(images[0][grid[:, 1], grid[:, 0]], rgba)
    assert np.all(images[0][[1, 2], [1, 3], 3] == 0.0)


def test_calicost_table_is_the_first_holding_every_labelled_clone(tmp_path: Path) -> None:
    """Clones 0-2 labelled. The unconstrained table, holding every clone, is read; leaving out clone 2, the first of
    diploid then tetraploid that holds it is read, diploid where both do; where neither does, the unconstrained."""

    def write(ploidy: str, clones: list[int]) -> None:
        suffix = "" if ploidy == "unconstrained" else f"_{ploidy}"
        columns = {f"clone{k} {a}": [1] for k in clones for a in "AB"}
        pd.DataFrame({"CHR": [1], "START": [0], "END": [10], **columns}).to_csv(
            tmp_path / f"cnv{suffix}_seglevel.tsv", sep="\t", index=False
        )

    pd.DataFrame({"BARCODES": ["a", "b", "c"], "clone_label": [0, 1, 2]}).to_csv(
        tmp_path / "clone_labels.tsv", sep="\t", index=False
    )
    write("unconstrained", [0, 1, 2])
    assert get_calicost_ploidy(tmp_path) == "unconstrained"
    write("unconstrained", [0, 1])
    write("diploid", [0, 1])
    write("tetraploid", [0, 1, 2])
    assert get_calicost_ploidy(tmp_path) == "tetraploid"
    write("diploid", [0, 1, 2])
    assert get_calicost_ploidy(tmp_path) == "diploid"
    write("diploid", [0, 1])
    write("tetraploid", [0, 2])
    assert get_calicost_ploidy(tmp_path) == "unconstrained"
    assert get_calicost_ploidy(tmp_path, fallbacks=()) == "unconstrained"


def test_chromosomes_without_a_cna_are_left_out(tmp_path: Path) -> None:
    """A chromosome is drawn whole where a drawn clone of any profile holds a CNA on either arm: chr1's only drawn CNA
    is an integer (2,1) on chr1p, so chr1q is drawn too; chr2 holds a state table's AMP and an integer (2,0); chr3's
    only CNA is a state table's clone 0, the normal clone, not drawn, so chr3 is left out."""
    pytest.importorskip("cnplot", reason="cnplot is in the luo dependency group: uv sync --group luo")
    from scripts.plots.layout import get_cna_regions, get_cnplot_segments, get_genome_axis

    (tmp_path / "regions.bed").write_text(
        "chr1\t0\t100\tchr1p\nchr1\t150\t300\tchr1q\nchr2\t0\t100\tchr2p\nchr2\t150\t300\tchr2q\n"
        "chr3\t0\t100\tchr3p\nchr3\t150\t300\tchr3q\n"
    )
    (tmp_path / "sizes").write_text("chr1\t300\nchr2\t300\nchr3\t300\n")
    bed, sizes = str(tmp_path / "regions.bed"), str(tmp_path / "sizes")
    states = pd.DataFrame(
        {
            "CHR": [1, 2, 3],
            "START": [150, 20, 20],
            "END": [280, 80, 80],
            "clone0": ["DEL", "NEU", "DEL"],
            "clone1": ["NEU", "AMP", "NEU"],
        }
    )
    ones = [1, 1, 1]
    copies = pd.DataFrame(
        {"CHR": [1, 1, 2], "START": [20, 150, 150], "END": [100, 300, 300], "clone0 A": ones, "clone0 B": ones}
    ).assign(**{"clone1 A": [2, 1, 2], "clone1 B": [1, 1, 0]})

    regions = get_cna_regions({"states": states, "copies": copies, "none": None}, get_genome_axis(bed, sizes), "NEU")
    axis = get_genome_axis(bed, sizes, regions)

    assert regions["NAME"].tolist() == ["chr1p", "chr1q", "chr2p", "chr2q"]
    assert axis.xlab_chrs == ["chr1", "chr2"]
    assert get_cnplot_segments(states, axis)[["#CHR", "START", "END"]].to_numpy().tolist() == [
        ["chr1", 150, 280],
        ["chr2", 20, 80],
    ]


def test_figures_name_a_run_by_its_method_side_and_setting_alone() -> None:
    """A figure names a run by its method, its bins' side and its `msu<k>`: a cnamaste run's id in its key is left
    out, and a CalicoST sweep run keeps its min-snp-umi."""
    pytest.importorskip("cnplot", reason="cnplot is in the luo dependency group: uv sync --group luo")
    from scripts.plots.style import get_run_label

    assert get_run_label("cnamaste-ac88605f", "cnamaste", 16) == get_run_label("cnamaste", "cnamaste", 16)
    assert get_run_label("calicost_016um_msu10", "calicost", 16) == r"$\mathtt{CalicoST}$ (16 $\mu$m, min-snp-umi: 10)"


def test_cnamaste_tracks_pool_each_integer_clones_spots() -> None:
    """Fitted clones 0 and 2 merge into integer clone 0, 1 and 3 into 1: each spot's track is its integer clone's, and
    each integer clone takes its first clone's library shift; a fitted clone split over two integer clones is refused."""
    pytest.importorskip("cnplot", reason="cnplot is in the luo dependency group: uv sync --group luo")
    from scripts.plots.cnamaste_page import get_integer_fit

    barcodes = np.array(["s0", "s1", "s2", "s3", "s4"])
    fit = {"new_assignment": np.array([0, 1, 2, 3, 3]), "new_log_mu_shift": np.array([0.0, 0.1, 0.2, 0.3])}
    labels = pd.Series([1, 1, 0, 0, 1], index=["s4", "s3", "s2", "s0", "s1"])

    integer = get_integer_fit(fit, barcodes, labels)

    assert integer["new_assignment"].tolist() == [0, 1, 0, 1, 1]
    assert integer["new_log_mu_shift"].tolist() == [0.0, 0.1]
    with pytest.raises(ValueError, match="not a merge"):
        get_integer_fit(fit, barcodes, pd.Series([0, 1, 0, 1, 3], index=barcodes))


def test_cnamaste_inputs_are_one_runs_integer_clones() -> None:
    """The labels must be the run's integer clones and the copies the same run's; anything else is refused."""
    pytest.importorskip("cnplot", reason="cnplot is in the luo dependency group: uv sync --group luo")
    from scripts.plots.cnamaste_page import get_cnamaste_input_mismatch

    run, other = Path("run_a"), Path("run_b")
    assert get_cnamaste_input_mismatch(run / "clone_labels_int.tsv.gz", run / "cnv_diploid_seglevel.tsv") == ""
    assert "integer clones" in get_cnamaste_input_mismatch(run / "clone_labels.tsv.gz", run / "cnv_seglevel.tsv")
    assert "run directory" in get_cnamaste_input_mismatch(run / "clone_labels_int.tsv", other / "cnv_seglevel.tsv")


def test_cnamaste_normal_fractions_are_its_logs_unless_given(tmp_path: Path) -> None:
    """Each integer clone takes the logged fraction of the fitted clone naming it, clone 0 held at 1; fractions given
    on the command line replace the log's; a clone the log leaves out is refused."""
    pytest.importorskip("cnplot", reason="cnplot is in the luo dependency group: uv sync --group luo")
    from scripts.plots.cnamaste_page import get_normal_fractions

    line = (
        "INFO - run_integer_copy_decoding:1258 - hmm_ab normal fractions per clone: [1.0, 0.45, 0.0, 0.36], status {}."
    )
    (tmp_path / "cnamaste_x.log").write_text(f"start\n{line}\nDone\n")

    assert get_normal_fractions(tmp_path, [0, 1], {}) == {0: 1.0, 1: 0.45}
    assert get_normal_fractions(tmp_path, [0, 1], {1: 0.3}) == {0: 1.0, 1: 0.3}
    with pytest.raises(ValueError, match="none of clone"):
        get_normal_fractions(tmp_path, [0, 4], {})
