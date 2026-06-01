"""Shared fixtures for shiba2corr tests."""

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # noqa: E402  must be set before any pyplot import

import numpy as np
import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="session")
def example_target_dir():
	return str(REPO_ROOT / "example" / "target")


@pytest.fixture(scope="session")
def example_reference_dir():
	return str(REPO_ROOT / "example" / "reference")


@pytest.fixture
def tmp_output_dir(tmp_path):
	out = tmp_path / "output"
	out.mkdir()
	return str(out)


@pytest.fixture
def synthetic_union_data():
	"""Tiny union_data dict accepted by compute_dse_overlap_by_direction."""
	return {
		"target": {
			"SE": {
				"dpsi": {"a": 0.2, "b": -0.3, "c": 0.1, "z": 0.4},
				"coverage": {"a": 50, "b": 80, "c": 30, "z": 90},
				"variance": {"a": 0.01, "b": 0.02, "c": 0.015, "z": 0.005},
			},
		},
		"reference": {
			"SE": {
				"dpsi": {"a": 0.25, "b": -0.4, "d": 0.15, "z": 0.35},
				"coverage": {"a": 60, "b": 70, "d": 40, "z": 95},
				"variance": {"a": 0.012, "b": 0.018, "d": 0.02, "z": 0.006},
			},
		},
		"dse_target": {"SE": {"a", "b", "c", "z"}},
		"dse_reference": {"SE": {"a", "b", "d", "z"}},
		"target_unfiltered": None,
		"reference_unfiltered": None,
	}


@pytest.fixture
def synthetic_overlap_data():
	"""Minimal overlap_data dict accepted by plot_dse_overlap_venns."""
	return {
		"SE": {
			"up": {
				"target_count": 5,
				"reference_count": 4,
				"intersection_count": 2,
				"target_only": 3,
				"reference_only": 2,
			},
			"down": {
				"target_count": 3,
				"reference_count": 3,
				"intersection_count": 1,
				"target_only": 2,
				"reference_only": 2,
			},
		},
		"all": {
			"up": {
				"target_count": 5,
				"reference_count": 4,
				"intersection_count": 2,
				"target_only": 3,
				"reference_only": 2,
			},
			"down": {
				"target_count": 0,
				"reference_count": 0,
				"intersection_count": 0,
				"target_only": 0,
				"reference_only": 0,
			},
		},
	}


@pytest.fixture
def synthetic_psi_df():
	"""DataFrame mimicking a Shiba PSI_<EVENT>.txt minimal payload."""
	return pd.DataFrame({
		"pos_id": ["p1", "p2", "p3", "p4"],
		"ref_junction_1": ["10;5", "0;0", "12;8", "20;15"],
		"ref_junction_2": ["8;6", "1;1", "10;10", "15;10"],
		"alt_junction_1": ["3;2", "0;0", "20;18", "30;25"],
		"alt_junction_2": ["4;3", "1;0", "22;19", "28;22"],
		"ref_PSI": [0.6, np.nan, 0.4, 0.55],
		"alt_PSI": [0.3, np.nan, 0.7, 0.25],
		"dPSI": [-0.3, np.nan, 0.3, -0.3],
		"Diff events": ["Yes", "No", "Yes", "Yes"],
	})
