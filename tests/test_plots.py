"""Tests for shiba2corr.plots."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from shiba2corr import plots


def test_apply_font_none_is_noop():
	plots._apply_font(None)  # should not raise


def test_apply_font_sets_rcparams():
	plots._apply_font("DejaVu Sans")
	assert matplotlib.rcParams["font.family"] == ["DejaVu Sans"] or \
		matplotlib.rcParams["font.family"] == "DejaVu Sans"


def test_ensure_plot_dirs_creates(tmp_path):
	png_dir, pdf_dir = plots.ensure_plot_dirs(str(tmp_path))
	assert png_dir.is_dir() and pdf_dir.is_dir()
	assert png_dir.name == "png" and pdf_dir.name == "pdf"


# ---------------------------------------------------------------------------
# DSE Venn drawing
# ---------------------------------------------------------------------------

def test_draw_dse_direction_venn_zero_counts():
	fig, ax = plt.subplots()
	counts = {"target_only": 0, "reference_only": 0, "intersection_count": 0,
			  "target_count": 0, "reference_count": 0}
	plots.draw_dse_direction_venn(ax, counts, "SE", "up")
	plt.close(fig)


def test_draw_dse_direction_venn_happy_up():
	fig, ax = plt.subplots()
	counts = {"target_only": 3, "reference_only": 2, "intersection_count": 1,
			  "target_count": 4, "reference_count": 3}
	plots.draw_dse_direction_venn(ax, counts, "SE", "up",
								   target_color="#111111", reference_color="#222222")
	plt.close(fig)


def test_draw_dse_direction_venn_happy_down():
	fig, ax = plt.subplots()
	counts = {"target_only": 2, "reference_only": 4, "intersection_count": 1,
			  "target_count": 3, "reference_count": 5}
	plots.draw_dse_direction_venn(ax, counts, "MXE", "down")
	plt.close(fig)


def test_draw_dse_direction_venn_no_venn_fallback(monkeypatch):
	monkeypatch.setattr(plots, "VENN_AVAILABLE", False)
	fig, ax = plt.subplots()
	counts = {"target_only": 1, "reference_only": 1, "intersection_count": 1,
			  "target_count": 2, "reference_count": 2}
	plots.draw_dse_direction_venn(ax, counts, "SE", "up")
	plt.close(fig)


def test_draw_dse_direction_venn_grid_empty_returns(tmp_path):
	png = tmp_path / "g.png"
	pdf = tmp_path / "g.pdf"
	# overlap_data with only empty direction dicts -> no valid events -> early return.
	plots.draw_dse_direction_venn_grid({}, "up", png, pdf)
	assert not png.exists() and not pdf.exists()


def test_draw_dse_direction_venn_grid_single_event(tmp_path):
	# n_events == 1 -> rows == cols == 1 branch
	overlap = {"SE": {"up": {"target_only": 1, "reference_only": 1, "intersection_count": 1,
								 "target_count": 2, "reference_count": 2}}}
	png = tmp_path / "g.png"
	pdf = tmp_path / "g.pdf"
	plots.draw_dse_direction_venn_grid(overlap, "up", png, pdf)
	assert png.exists() and pdf.exists()


def test_draw_dse_direction_venn_grid_pads_axes(tmp_path):
	# n_events == 4 -> 2x3 grid with 2 empty axes (pad branch)
	overlap = {f"E{i}": {"up": {"target_only": 1, "reference_only": 1, "intersection_count": 1,
								   "target_count": 2, "reference_count": 2}}
			   for i in range(4)}
	png = tmp_path / "g.png"
	pdf = tmp_path / "g.pdf"
	plots.draw_dse_direction_venn_grid(overlap, "up", png, pdf)
	assert png.exists()


def test_draw_dse_direction_venn_grid_row_or_col_only(tmp_path):
	# n_events == 2 -> rows == 1, cols == 2 -> elif branch
	overlap = {f"E{i}": {"down": {"target_only": 1, "reference_only": 1, "intersection_count": 1,
									 "target_count": 2, "reference_count": 2}}
			   for i in range(2)}
	png = tmp_path / "g.png"
	pdf = tmp_path / "g.pdf"
	plots.draw_dse_direction_venn_grid(overlap, "down", png, pdf)
	assert png.exists()


def test_plot_dse_overlap_venns_skips_empty_direction(tmp_path):
	# Empty direction dict triggers `if not counts: continue` in plot_dse_overlap_venns.
	overlap = {"SE": {"up": {}, "down": {"target_only": 1, "reference_only": 1,
										 "intersection_count": 1, "target_count": 2,
										 "reference_count": 2}}}
	plots.plot_dse_overlap_venns(overlap, str(tmp_path))
	assert not (tmp_path / "plots" / "png" / "dse_venn_SE_up.png").exists()
	assert (tmp_path / "plots" / "png" / "dse_venn_SE_down.png").exists()


def test_plot_dse_overlap_venns_writes_files(tmp_path, synthetic_overlap_data):
	plots.plot_dse_overlap_venns(synthetic_overlap_data, str(tmp_path))
	png_dir = tmp_path / "plots" / "png"
	pdf_dir = tmp_path / "plots" / "pdf"
	# SE has both directions; "all" has only up (down counts all zero -> still drawn).
	assert (png_dir / "dse_venn_SE_up.png").exists()
	assert (pdf_dir / "dse_venn_SE_up.pdf").exists()
	assert (png_dir / "dse_venn_grid_up.png").exists()
	assert (png_dir / "dse_venn_grid_down.png").exists()


# ---------------------------------------------------------------------------
# Weighted correlation scatter plots
# ---------------------------------------------------------------------------

def _scatter_df(n=20, seed=0):
	rng = np.random.default_rng(seed)
	x = rng.uniform(-0.5, 0.5, n)
	y = x + rng.normal(0, 0.05, n)
	return pd.DataFrame({
		"pos_id": [f"p{i}" for i in range(n)],
		"target_dpsi": y,
		"reference_dpsi": x,
		"weight": np.ones(n),
		"weight_normalized": np.ones(n) / n,
	})


def test_plot_weighted_correlation_scatter_writes_files(tmp_path):
	scatter = {"SE": _scatter_df()}
	results = {"SE": {"n_used": 20, "r": 0.9, "p_value_approx": 1e-5}}
	plots.plot_weighted_correlation_scatter(scatter, results, str(tmp_path))
	assert (tmp_path / "plots" / "png" / "weighted_correlation_scatter_SE.png").exists()
	assert (tmp_path / "plots" / "pdf" / "weighted_correlation_scatter_SE.pdf").exists()


def test_plot_weighted_correlation_scatter_skips_empty(tmp_path):
	scatter = {
		"SE": pd.DataFrame(columns=["pos_id", "target_dpsi", "reference_dpsi",
									 "weight", "weight_normalized"]),
		"FIVE": _scatter_df(n=1),  # n_used < 2 branch
	}
	results = {
		"SE": {"n_used": 0, "r": np.nan, "p_value_approx": np.nan},
		"FIVE": {"n_used": 1, "r": np.nan, "p_value_approx": np.nan},
	}
	plots.plot_weighted_correlation_scatter(scatter, results, str(tmp_path))
	# Neither event should have a file
	assert not (tmp_path / "plots" / "png" / "weighted_correlation_scatter_SE.png").exists()
	assert not (tmp_path / "plots" / "png" / "weighted_correlation_scatter_FIVE.png").exists()


def test_plot_weighted_correlation_scatter_singular_kde(tmp_path):
	df = pd.DataFrame({
		"pos_id": ["a", "b", "c"],
		"target_dpsi": [0.1, 0.1, 0.1],
		"reference_dpsi": [0.2, 0.2, 0.2],
		"weight": [1, 1, 1],
		"weight_normalized": [1 / 3, 1 / 3, 1 / 3],
	})
	scatter = {"SE": df}
	results = {"SE": {"n_used": 3, "r": np.nan, "p_value_approx": np.nan}}
	plots.plot_weighted_correlation_scatter(scatter, results, str(tmp_path))
	assert (tmp_path / "plots" / "png" / "weighted_correlation_scatter_SE.png").exists()


def test_plot_weighted_correlation_scatter_p_value_formats(tmp_path):
	# exp == 0 path (p=0.5) and small-p path (p=1e-5) handled within one frame each.
	scatter = {"SE": _scatter_df()}
	results_big_p = {"SE": {"n_used": 20, "r": 0.1, "p_value_approx": 0.5}}
	plots.plot_weighted_correlation_scatter(scatter, results_big_p, str(tmp_path / "a"))
	assert (tmp_path / "a" / "plots" / "png" / "weighted_correlation_scatter_SE.png").exists()

	results_small_p = {"SE": {"n_used": 20, "r": 0.9, "p_value_approx": 1e-5}}
	plots.plot_weighted_correlation_scatter(scatter, results_small_p, str(tmp_path / "b"))
	assert (tmp_path / "b" / "plots" / "png" / "weighted_correlation_scatter_SE.png").exists()

	# NaN p_value -> "NA" branch
	results_nan = {"SE": {"n_used": 20, "r": np.nan, "p_value_approx": np.nan}}
	plots.plot_weighted_correlation_scatter(scatter, results_nan, str(tmp_path / "c"))
	assert (tmp_path / "c" / "plots" / "png" / "weighted_correlation_scatter_SE.png").exists()
