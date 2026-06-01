"""Tests for shiba2corr.general."""

import json
import sys
from pathlib import Path

import pytest

from shiba2corr import general


def test_format_command_line(monkeypatch):
	monkeypatch.setattr(sys, "argv", ["shiba2corr", "-t", "a", "-r", "b"])
	assert general.format_command_line() == "shiba2corr -t a -r b"


def test_create_simple_report_shape():
	class Args:
		pass
	report = general.create_simple_report(Args(), 1_700_000_000.0, 1_700_000_005.5, "v0.1.0")
	assert report["tool"]["name"] == "shiba2corr"
	assert report["tool"]["version"] == "v0.1.0"
	assert report["run"]["duration_seconds"] == 5.5
	assert "T" in report["run"]["start_time"]
	assert "T" in report["run"]["end_time"]


def test_write_report_writes_json(tmp_path):
	report = {"tool": {"name": "shiba2corr"}, "run": {}}
	general.write_report(report, str(tmp_path))
	out = Path(tmp_path) / "report.json"
	assert out.exists()
	loaded = json.loads(out.read_text())
	assert loaded["tool"]["name"] == "shiba2corr"


def test_write_report_handles_failure(tmp_path, caplog):
	# Use a path where the parent is a file, not a directory, to trigger the except branch.
	blocker = tmp_path / "blocker"
	blocker.write_text("not a dir")
	with caplog.at_level("WARNING"):
		general.write_report({"x": 1}, str(blocker))
	assert any("Could not write report" in rec.message for rec in caplog.records)
