"""Tests for shiba2corr.correlation."""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from shiba2corr import correlation
from shiba2corr.correlation import (
	WEIGHT_SCHEMES,
	compute_weighted_correlation_dpsi_union,
	save_scatter_plot_data,
	save_weighted_correlation_results,
)


@pytest.mark.parametrize("scheme", list(WEIGHT_SCHEMES))
def test_compute_weighted_correlation_all_schemes(
	example_target_dir, example_reference_dir, scheme,
):
	results, scatter = compute_weighted_correlation_dpsi_union(
		example_target_dir, example_reference_dir, weight_scheme=scheme,
	)
	assert "SE" in results and "all" in results
	row = results["all"]
	assert {"n_union", "n_used", "weight_scheme", "r", "p_value_approx", "neff"} <= set(row)
	assert row["weight_scheme"] == scheme
	se_df = scatter["SE"]
	assert {"pos_id", "target_dpsi", "reference_dpsi", "weight", "weight_normalized"} <= set(se_df.columns)


def test_invalid_weight_scheme_raises(example_target_dir, example_reference_dir):
	with pytest.raises(ValueError):
		compute_weighted_correlation_dpsi_union(
			example_target_dir, example_reference_dir, weight_scheme="bogus",
		)


def test_event_list_filter_shrinks_union(
	example_target_dir, example_reference_dir, tmp_path,
):
	# Pick a pos_id known to exist in the target SE table.
	psi = pd.read_csv(
		Path(example_target_dir) / "results" / "splicing" / "PSI_SE.txt", sep="\t"
	)
	keep_id = str(psi["pos_id"].iloc[0])
	event_list = tmp_path / "events.txt"
	event_list.write_text(keep_id + "\n")
	results, _ = compute_weighted_correlation_dpsi_union(
		example_target_dir, example_reference_dir, event_list=str(event_list),
	)
	for et, row in results.items():
		assert row["n_union"] <= 1


def test_min_events_large_forces_nan_pvalue(example_target_dir, example_reference_dir):
	results, _ = compute_weighted_correlation_dpsi_union(
		example_target_dir, example_reference_dir, min_events=10_000,
	)
	for row in results.values():
		assert np.isnan(row["p_value_approx"])


def test_empty_event_list_yields_nan(
	example_target_dir, example_reference_dir, tmp_path,
):
	event_list = tmp_path / "empty.txt"
	event_list.write_text("")
	results, scatter = compute_weighted_correlation_dpsi_union(
		example_target_dir, example_reference_dir, event_list=str(event_list),
	)
	for row in results.values():
		assert row["n_used"] == 0
		assert np.isnan(row["r"])
	# Scatter frames are empty but have the expected columns.
	for df in scatter.values():
		assert df.empty
		assert {"pos_id", "target_dpsi", "reference_dpsi", "weight", "weight_normalized"} <= set(df.columns)


def test_geom_mean_skips_zero_coverage(monkeypatch):
	"""Cover the geom_mean and coverage_mean zero-coverage continue branches."""
	from shiba2corr import load

	def fake_load(*args, **kwargs):
		return {
			"target": {et: {"dpsi": {}, "coverage": {}, "variance": {}} for et in list(load.EVENT_TYPES) + ['all']},
			"reference": {et: {"dpsi": {}, "coverage": {}, "variance": {}} for et in list(load.EVENT_TYPES) + ['all']},
			"union": {et: set() for et in list(load.EVENT_TYPES) + ['all']},
		}

	# Inject one event with zero coverage so geom_mean / coverage_mean skip it.
	def fake_load_nonzero(*args, **kwargs):
		base = fake_load()
		base["target"]["SE"] = {
			"dpsi": {"p1": 0.1, "p2": 0.2},
			"coverage": {"p1": 0, "p2": 0},
			"variance": {"p1": 0.0, "p2": 0.0},
		}
		base["reference"]["SE"] = {
			"dpsi": {"p1": 0.15, "p2": 0.25},
			"coverage": {"p1": 0, "p2": 0},
			"variance": {"p1": 0.0, "p2": 0.0},
		}
		base["union"]["SE"] = {"p1", "p2"}
		return base

	monkeypatch.setattr(correlation.load, "load_union_event_data", fake_load_nonzero)
	results_g, _ = compute_weighted_correlation_dpsi_union(
		"x", "y", weight_scheme="geom_mean",
	)
	assert results_g["SE"]["n_used"] == 0
	results_cm, _ = compute_weighted_correlation_dpsi_union(
		"x", "y", weight_scheme="coverage_mean",
	)
	assert results_cm["SE"]["n_used"] == 0


def test_inverse_variance_skips_nan_variance(monkeypatch):
	from shiba2corr import load

	def fake_load(*args, **kwargs):
		base = {
			"target": {et: {"dpsi": {}, "coverage": {}, "variance": {}} for et in list(load.EVENT_TYPES) + ['all']},
			"reference": {et: {"dpsi": {}, "coverage": {}, "variance": {}} for et in list(load.EVENT_TYPES) + ['all']},
			"union": {et: set() for et in list(load.EVENT_TYPES) + ['all']},
		}
		base["target"]["SE"] = {
			"dpsi": {"p1": 0.1}, "coverage": {"p1": 10}, "variance": {"p1": np.nan},
		}
		base["reference"]["SE"] = {
			"dpsi": {"p1": 0.2}, "coverage": {"p1": 10}, "variance": {"p1": 0.01},
		}
		base["union"]["SE"] = {"p1"}
		return base

	monkeypatch.setattr(correlation.load, "load_union_event_data", fake_load)
	results, _ = compute_weighted_correlation_dpsi_union("x", "y", weight_scheme="inverse_variance")
	assert results["SE"]["n_used"] == 0


def test_skips_nan_dpsi(monkeypatch):
	"""Cover the `if pd.isna(t_dpsi) or pd.isna(r_dpsi): continue` branch."""
	from shiba2corr import load

	def fake_load(*args, **kwargs):
		base = {
			"target": {et: {"dpsi": {}, "coverage": {}, "variance": {}} for et in list(load.EVENT_TYPES) + ['all']},
			"reference": {et: {"dpsi": {}, "coverage": {}, "variance": {}} for et in list(load.EVENT_TYPES) + ['all']},
			"union": {et: set() for et in list(load.EVENT_TYPES) + ['all']},
		}
		base["target"]["SE"] = {
			"dpsi": {"p1": np.nan}, "coverage": {"p1": 10}, "variance": {"p1": 0.01},
		}
		base["reference"]["SE"] = {
			"dpsi": {"p1": 0.2}, "coverage": {"p1": 10}, "variance": {"p1": 0.01},
		}
		base["union"]["SE"] = {"p1"}
		return base

	monkeypatch.setattr(correlation.load, "load_union_event_data", fake_load)
	results, _ = compute_weighted_correlation_dpsi_union("x", "y", weight_scheme="uniform")
	assert results["SE"]["n_used"] == 0


def test_save_weighted_correlation_results(tmp_path):
	results = {
		"SE": {"n_union": 5, "n_used": 4, "r": 0.5, "p_value_approx": 0.01, "neff": 3.5,
			   "weight_scheme": "uniform"},
	}
	out = tmp_path / "out.tsv"
	save_weighted_correlation_results(results, str(out))
	df = pd.read_csv(out, sep="\t")
	assert list(df.columns) == ["event_type", "n_union", "n_used", "r", "p_value_approx", "neff"]
	assert df["event_type"].iloc[0] == "SE"


def test_save_scatter_plot_data_skips_empty(tmp_path):
	scatter = {
		"SE": pd.DataFrame({
			"pos_id": ["a"], "target_dpsi": [0.1], "reference_dpsi": [0.2],
			"weight": [1.0], "weight_normalized": [1.0],
		}),
		"FIVE": pd.DataFrame(columns=["pos_id", "target_dpsi", "reference_dpsi", "weight", "weight_normalized"]),
	}
	save_scatter_plot_data(scatter, str(tmp_path))
	assert (tmp_path / "weighted_correlation_scatter_SE.tsv").exists()
	assert not (tmp_path / "weighted_correlation_scatter_FIVE.tsv").exists()
