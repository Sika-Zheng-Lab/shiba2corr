"""shiba2corr command-line entry point."""

import argparse
import logging
import os
import time

from . import __version__, correlation, general, load, plots, publication, sensitivity

logger = logging.getLogger(__name__)


def parse_args():
	parser = argparse.ArgumentParser(
		prog="shiba2corr",
		description=(
			f"shiba2corr {__version__} — Weighted dPSI correlation analysis "
			f"for alternative splicing events quantified by Shiba"
		),
	)
	parser.add_argument("-v", "--version", action="version", version=f"shiba2corr {__version__}")
	parser.add_argument("-t", "--target", required=True, help="Target Shiba working directory")
	parser.add_argument("-r", "--reference", required=True, help="Reference Shiba working directory")
	parser.add_argument("-o", "--output", required=True, help="Output directory")
	parser.add_argument("-l", "--event-list", default=None, help="Optional event list file (one pos_id per line)")
	parser.add_argument("--weight-scheme", default="inverse_variance",
						choices=list(correlation.WEIGHT_SCHEMES),
						help="Weighting scheme for the correlation (default: inverse_variance)")
	parser.add_argument("--min-events", type=int, default=3,
						help="Minimum effective sample size required for the p-value (default: 3)")
	parser.add_argument("--ttest", type=float, default=None,
						help="Filter DSEs by t-test P-value threshold (e.g. 0.05). "
							 "Requires 'p_ttest' column in PSI files.")
	parser.add_argument("--target-color", default=None,
						help="Color for Target DSEs in Venn diagrams (e.g. '#E49414FF')")
	parser.add_argument("--reference-color", default=None,
						help="Color for Reference DSEs in Venn diagrams (e.g. '#8087AAFF')")
	parser.add_argument("--font-family", default=None, help="Font family for plots (e.g. 'Arial')")
	parser.add_argument("--no-venn", action="store_true",
						help="Skip the DSE overlap Venn diagrams")
	parser.add_argument("--verbose", action="store_true", help="Increase verbosity")
	parser.add_argument("--sensitivity-analysis", action="store_true",
						help="Repeat correlations on random event subsets and export publication figures")
	parser.add_argument("--subsample-fractions", type=float, nargs="+",
						default=list(sensitivity.DEFAULT_FRACTIONS), help="Retention fractions in (0, 1); default: 0.25 0.5 0.75 0.9")
	parser.add_argument("--subsample-repeats", type=int, default=200, help="Repeats per fraction (default: 200)")
	parser.add_argument("--subsample-seed", type=int, default=0, help="Nonnegative random seed (default: 0)")
	args = parser.parse_args()
	try:
		sensitivity.validate_settings(args.subsample_fractions, args.subsample_repeats, args.subsample_seed)
	except ValueError as exc:
		parser.error(str(exc))
	return args


def main():
	start_time = time.time()
	args = parse_args()

	logging.basicConfig(
		format="[%(asctime)s] %(levelname)7s %(message)s",
		level=logging.DEBUG if args.verbose else logging.INFO,
	)
	logger.info(f"Running shiba2corr ({__version__})")
	logger.debug(f"Arguments: {args}")

	if args.sensitivity_analysis:
		publication.publication_style()  # Fail early if required Arial fonts are missing.
	os.makedirs(os.path.join(args.output, "results"), exist_ok=True)

	logger.info("===== Weighted dPSI correlation =====")
	wcorr_results, scatter_data_dict = correlation.compute_weighted_correlation_dpsi_union(
		args.target, args.reference,
		event_list=args.event_list,
		weight_scheme=args.weight_scheme,
		min_events=args.min_events,
		ttest_p_threshold=args.ttest,
	)

	wcorr_file = os.path.join(args.output, "results", "event_dpsi_weighted_correlation.tsv")
	correlation.save_weighted_correlation_results(wcorr_results, wcorr_file)
	correlation.save_scatter_plot_data(scatter_data_dict, os.path.join(args.output, "results"))

	if args.sensitivity_analysis:
		logger.info("Running random event-selection sensitivity analysis")
		raw, summary = sensitivity.analyze_sensitivity(
			scatter_data_dict, fractions=args.subsample_fractions,
			repeats=args.subsample_repeats, seed=args.subsample_seed, min_events=args.min_events)
		raw.to_csv(os.path.join(args.output, "results", "event_selection_sensitivity.tsv"), sep="\t", index=False)
		summary.to_csv(os.path.join(args.output, "results", "event_selection_sensitivity_summary.tsv"), sep="\t", index=False)
		publication.plot_sensitivity(raw, summary, args.output)

	logger.info("Creating weighted correlation scatter plots")
	plots.plot_weighted_correlation_scatter(
		scatter_data_dict, wcorr_results, args.output,
		font_family=args.font_family,
	)

	if not args.no_venn:
		logger.info("===== DSE overlap Venn diagrams =====")
		union_data = load.load_union_event_data(args.target, args.reference,
												ttest_p_threshold=args.ttest)
		overlap_data = load.compute_dse_overlap_by_direction(union_data)
		plots.plot_dse_overlap_venns(
			overlap_data, args.output,
			target_color=args.target_color, reference_color=args.reference_color,
			font_family=args.font_family,
		)
	else:
		logger.info("Skipping DSE overlap Venn diagrams (--no-venn)")

	end_time = time.time()
	logger.info("shiba2corr finished successfully")
	logger.info(f"Output directory: {args.output}")

	report = general.create_simple_report(args, start_time, end_time, __version__)
	report["parameters"] = vars(args)
	general.write_report(report, args.output)


if __name__ == "__main__":  # pragma: no cover
	main()
