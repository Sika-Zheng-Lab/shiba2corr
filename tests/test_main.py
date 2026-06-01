"""Tests for shiba2corr.main."""

import json
import sys
from pathlib import Path

import pytest

from shiba2corr import main as cli


def test_version_exits(monkeypatch, capsys):
	monkeypatch.setattr(sys, "argv", ["shiba2corr", "--version"])
	with pytest.raises(SystemExit) as exc:
		cli.parse_args()
	assert exc.value.code == 0


def test_missing_required_args_exits(monkeypatch):
	monkeypatch.setattr(sys, "argv", ["shiba2corr"])
	with pytest.raises(SystemExit):
		cli.parse_args()


def test_main_end_to_end(monkeypatch, tmp_path, example_target_dir, example_reference_dir):
	out = tmp_path / "out"
	monkeypatch.setattr(sys, "argv", [
		"shiba2corr",
		"-t", example_target_dir,
		"-r", example_reference_dir,
		"-o", str(out),
		"--verbose",
	])
	cli.main()
	assert (out / "results" / "event_dpsi_weighted_correlation.tsv").exists()
	assert (out / "results" / "weighted_correlation_scatter_SE.tsv").exists()
	assert (out / "report.json").exists()
	# Venn plots produced (no --no-venn)
	png_dir = out / "plots" / "png"
	assert any(png_dir.glob("dse_venn_*_up.png"))
	report = json.loads((out / "report.json").read_text())
	assert report["tool"]["name"] == "shiba2corr"


def test_main_no_venn(monkeypatch, tmp_path, example_target_dir, example_reference_dir):
	out = tmp_path / "out"
	monkeypatch.setattr(sys, "argv", [
		"shiba2corr",
		"-t", example_target_dir,
		"-r", example_reference_dir,
		"-o", str(out),
		"--no-venn",
		"--weight-scheme", "uniform",
	])
	cli.main()
	png_dir = out / "plots" / "png"
	# Scatter plots present, but no Venn plots.
	assert any(png_dir.glob("weighted_correlation_scatter_*.png"))
	assert not list(png_dir.glob("dse_venn_*"))
