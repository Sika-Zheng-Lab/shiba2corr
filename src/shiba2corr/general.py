import json
import sys
from datetime import datetime
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


def format_command_line():
	return " ".join(sys.argv)


def create_simple_report(args, start_time, end_time, version):
	report = {
		"tool": {
			"name": "shiba2corr",
			"version": version,
		},
		"run": {
			"command": format_command_line(),
			"start_time": datetime.fromtimestamp(start_time).isoformat(),
			"end_time": datetime.fromtimestamp(end_time).isoformat(),
			"duration_seconds": round(end_time - start_time, 2),
		},
	}
	return report


def write_report(report_dict, output_dir):
	output_path = Path(output_dir) / "report.json"
	try:
		output_path.parent.mkdir(parents=True, exist_ok=True)
		with open(output_path, "w", encoding="utf-8") as f:
			json.dump(report_dict, f, indent=2, ensure_ascii=False)
		logger.info(f"Report saved to {output_path}")
	except Exception as e:
		logger.warning(f"Could not write report: {e}")
