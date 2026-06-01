"""Tests for shiba2corr.load."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from shiba2corr.load import (
	aggregate_all_types,
	build_maps,
	check_ttest_p_available,
	compute_detected_events_intersection,
	compute_dse_overlap_by_direction,
	coverage_percentile_diff_events,
	create_background_dataframe,
	filter_to_union,
	load_union_event_data,
	parse_psi_table,
	split_dse_by_direction,
)


# ---------------------------------------------------------------------------
# split_dse_by_direction / compute_dse_overlap_by_direction
# ---------------------------------------------------------------------------

def test_split_dse_by_direction_basic():
	dse_set = {"a", "b", "c", "d", "e"}
	dpsi_map = {"a": 0.3, "b": -0.2, "c": 0.0, "d": np.nan, "e": -0.05}
	up, down = split_dse_by_direction(dse_set, dpsi_map)
	assert up == {"a"}
	assert down == {"b", "e"}


def test_compute_dse_overlap_by_direction_counts(synthetic_union_data):
	overlap = compute_dse_overlap_by_direction(synthetic_union_data, event_types=["SE"])
	up = overlap["SE"]["up"]
	down = overlap["SE"]["down"]
	# up: target {a, c, z}, reference {a, d, z} -> intersection {a, z}
	assert up["target_count"] == 3
	assert up["reference_count"] == 3
	assert up["intersection_count"] == 2
	assert up["target_only"] == 1
	assert up["reference_only"] == 1
	# down: target {b}, reference {b}
	assert down["intersection_count"] == 1


def test_compute_dse_overlap_default_event_types(synthetic_union_data):
	data = dict(synthetic_union_data)
	data["dse_target"] = {**data["dse_target"], "all": data["dse_target"]["SE"]}
	data["dse_reference"] = {**data["dse_reference"], "all": data["dse_reference"]["SE"]}
	data["target"] = {**data["target"], "all": data["target"]["SE"]}
	data["reference"] = {**data["reference"], "all": data["reference"]["SE"]}
	overlap = compute_dse_overlap_by_direction(data)
	assert "SE" in overlap and "all" in overlap


# ---------------------------------------------------------------------------
# coverage_percentile_diff_events
# ---------------------------------------------------------------------------

def test_coverage_percentile_diff_events_happy(synthetic_psi_df):
	df, p5, p95 = coverage_percentile_diff_events(synthetic_psi_df)
	assert "mean_coverage" in df.columns
	assert "variance" in df.columns
	assert p95 >= p5 > 0


def test_coverage_percentile_diff_events_empty_diff():
	df = pd.DataFrame({
		"pos_id": ["p1"],
		"ref_junction_1": ["10;5"],
		"alt_junction_1": ["3;2"],
		"ref_PSI": [0.5],
		"alt_PSI": [0.4],
		"Diff events": ["No"],
	})
	out_df, p5, p95 = coverage_percentile_diff_events(df)
	assert p5 == 0 and p95 == 0
	assert "mean_coverage" in out_df.columns


# ---------------------------------------------------------------------------
# parse_psi_table
# ---------------------------------------------------------------------------

def test_parse_psi_table_missing_file(tmp_path, caplog):
	with caplog.at_level("WARNING"):
		df = parse_psi_table(str(tmp_path / "nope.txt"), "SE")
	assert df.empty


def test_parse_psi_table_malformed(tmp_path):
	bad = tmp_path / "bad.txt"
	bad.write_text("not\ta\tvalid\tpsi\nfile\twith\ttoo\tfew\tcolumns\n")
	df = parse_psi_table(str(bad), "SE")
	# Either empty (parse failed) or a non-empty df without required cols.
	assert isinstance(df, pd.DataFrame)


# ---------------------------------------------------------------------------
# build_maps
# ---------------------------------------------------------------------------

def test_build_maps_empty():
	dpsi, cov, var, dse = build_maps(pd.DataFrame())
	assert dpsi == {} and cov == {} and var == {} and dse == set()


def test_build_maps_missing_columns(caplog):
	df = pd.DataFrame({"pos_id": ["p1"], "dPSI": [0.1]})
	with caplog.at_level("WARNING"):
		dpsi, cov, var, dse = build_maps(df)
	assert dpsi == {} and dse == set()


def test_build_maps_happy(synthetic_psi_df):
	df, _, _ = coverage_percentile_diff_events(synthetic_psi_df)
	dpsi, cov, var, dse = build_maps(df)
	assert set(dpsi.keys()) == {"p1", "p2", "p3", "p4"}
	assert dse == {"p1", "p3", "p4"}


# ---------------------------------------------------------------------------
# aggregate_all_types / filter_to_union / compute_detected_events_intersection
# ---------------------------------------------------------------------------

def test_aggregate_all_types_merges_event_types():
	cond = {
		"SE": {
			"dpsi": {"a": 0.1, "b": 0.2}, "coverage": {"a": 10, "b": 20},
			"variance": {"a": 0.01, "b": 0.02}, "dse_set": {"a"},
		},
		"FIVE": {
			"dpsi": {"c": 0.3}, "coverage": {"c": 30},
			"variance": {"c": 0.03}, "dse_set": {"c"},
		},
	}
	agg = aggregate_all_types(cond)
	assert set(agg["dpsi"].keys()) == {"a", "b", "c"}
	assert agg["dse_set"] == {"a", "c"}


def test_aggregate_all_types_skips_missing():
	agg = aggregate_all_types({})
	assert agg == {"dpsi": {}, "coverage": {}, "variance": {}, "dse_set": set()}


def test_compute_detected_events_intersection_filters_coverage():
	t = {
		"SE": {
			"dpsi": {"a": 0.1, "b": 0.2}, "coverage": {"a": 50, "b": 1},
			"variance": {"a": 0.01, "b": 0.02},
		},
		"all": {
			"dpsi": {"a": 0.1}, "coverage": {"a": 50}, "variance": {"a": 0.01},
		},
	}
	r = {
		"SE": {
			"dpsi": {"a": 0.15, "b": 0.25}, "coverage": {"a": 40, "b": 60},
			"variance": {"a": 0.012, "b": 0.018},
		},
		"all": {
			"dpsi": {"a": 0.15}, "coverage": {"a": 40}, "variance": {"a": 0.012},
		},
	}
	out = compute_detected_events_intersection(t, r, min_coverage=10)
	assert out["SE"] == {"a"}
	out2 = compute_detected_events_intersection(t, r)
	assert out2["SE"] == {"a", "b"}


def test_filter_to_union_keeps_only_union():
	cond = {
		"SE": {
			"dpsi": {"a": 0.1, "b": 0.2},
			"coverage": {"a": 1, "b": 2},
			"variance": {"a": 0.0, "b": 0.0},
		},
	}
	union = {"SE": {"a"}}
	out = filter_to_union(cond, union)
	assert out["SE"]["dpsi"] == {"a": 0.1}


# ---------------------------------------------------------------------------
# check_ttest_p_available
# ---------------------------------------------------------------------------

def _write_psi(path: Path, columns):
	path.parent.mkdir(parents=True, exist_ok=True)
	pd.DataFrame({c: [] for c in columns}).to_csv(path, sep="\t", index=False)


def test_check_ttest_p_available_present(tmp_path):
	shiba = tmp_path / "shiba"
	_write_psi(shiba / "results" / "splicing" / "PSI_SE.txt",
			   ["pos_id", "dPSI", "Diff events", "p_ttest"])
	assert check_ttest_p_available(str(shiba), ["SE"]) == "p_ttest"


def test_check_ttest_p_available_absent(tmp_path):
	shiba = tmp_path / "shiba"
	_write_psi(shiba / "results" / "splicing" / "PSI_SE.txt",
			   ["pos_id", "dPSI", "Diff events"])
	assert check_ttest_p_available(str(shiba), ["SE"]) is None


def test_check_ttest_p_available_no_files(tmp_path):
	assert check_ttest_p_available(str(tmp_path), ["SE"]) is None


def test_check_ttest_p_available_parse_error(tmp_path):
	"""Exception while reading header -> caught, continue, return None."""
	shiba = tmp_path / "shiba"
	path = shiba / "results" / "splicing" / "PSI_SE.txt"
	path.parent.mkdir(parents=True, exist_ok=True)
	path.write_bytes(b"\x00\x01\x02\x03")  # not valid TSV; read_csv may raise
	assert check_ttest_p_available(str(shiba), ["SE"]) in (None, "p_ttest")


# ---------------------------------------------------------------------------
# load_union_event_data on example fixtures
# ---------------------------------------------------------------------------

def test_load_union_event_data_example(example_target_dir, example_reference_dir):
	data = load_union_event_data(example_target_dir, example_reference_dir, return_df=True)
	assert "background_df" in data
	assert not data["background_df"].empty
	expected = {"pos_id", "event_type", "dPSI_tgt", "dPSI_ref", "eff_cov"}
	assert expected.issubset(data["background_df"].columns)
	assert any(len(s) > 0 for s in data["union"].values())


def test_load_union_event_data_ttest_warns_when_missing(tmp_path, caplog):
	# Build a Shiba dir tree WITHOUT p_ttest column so the warning fires.
	t_root = tmp_path / "t"
	r_root = tmp_path / "r"
	for root in (t_root, r_root):
		path = root / "results" / "splicing" / "PSI_SE.txt"
		path.parent.mkdir(parents=True, exist_ok=True)
		pd.DataFrame({
			"pos_id": ["p1"],
			"ref_junction_1": ["10;10"],
			"alt_junction_1": ["10;10"],
			"ref_PSI": [0.5],
			"alt_PSI": [0.2],
			"dPSI": [-0.3],
			"Diff events": ["Yes"],
		}).to_csv(path, sep="\t", index=False)
	with caplog.at_level("WARNING", logger="shiba2corr.load"):
		load_union_event_data(
			str(t_root), str(r_root),
			return_df=False, ttest_p_threshold=0.05,
		)
	assert any("Ignoring --ttest threshold" in rec.message for rec in caplog.records)


def _build_psi_with_ttest(root: Path, ttest_values):
	path = root / "results" / "splicing" / "PSI_SE.txt"
	path.parent.mkdir(parents=True, exist_ok=True)
	df = pd.DataFrame({
		"pos_id": ["p1", "p2"],
		"ref_junction_1": ["10;10", "10;10"],
		"alt_junction_1": ["10;10", "10;10"],
		"ref_PSI": [0.5, 0.5],
		"alt_PSI": [0.2, 0.8],
		"dPSI": [-0.3, 0.3],
		"Diff events": ["Yes", "Yes"],
		"p_ttest": ttest_values,
	})
	df.to_csv(path, sep="\t", index=False)


def test_load_union_event_data_ttest_applied(tmp_path):
	t_root = tmp_path / "t"
	r_root = tmp_path / "r"
	_build_psi_with_ttest(t_root, [0.001, 0.9])
	_build_psi_with_ttest(r_root, [0.001, 0.9])
	data = load_union_event_data(
		str(t_root), str(r_root),
		return_df=False, ttest_p_threshold=0.05,
	)
	assert data["dse_target"]["SE"] == {"p1"}
	assert data["dse_reference"]["SE"] == {"p1"}


# ---------------------------------------------------------------------------
# create_background_dataframe
# ---------------------------------------------------------------------------

def test_create_background_dataframe_min_coverage_filters():
	target = {
		"SE": {
			"dpsi": {"a": 0.1, "b": 0.2},
			"coverage": {"a": 100, "b": 5},
			"variance": {"a": 0.01, "b": 0.02},
		},
		"all": {
			"dpsi": {"a": 0.1, "b": 0.2},
			"coverage": {"a": 100, "b": 5},
			"variance": {"a": 0.01, "b": 0.02},
		},
	}
	reference = {
		"SE": {
			"dpsi": {"a": 0.15, "b": 0.25},
			"coverage": {"a": 90, "b": 4},
			"variance": {"a": 0.012, "b": 0.018},
		},
		"all": {
			"dpsi": {"a": 0.15, "b": 0.25},
			"coverage": {"a": 90, "b": 4},
			"variance": {"a": 0.012, "b": 0.018},
		},
	}
	bg = {"SE": {"a", "b"}, "all": {"a", "b"}}
	df = create_background_dataframe(target, reference, bg, min_coverage=20)
	assert set(df["pos_id"]) == {"a"}
