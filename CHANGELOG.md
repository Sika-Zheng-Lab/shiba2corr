# Change log

All notable changes to this shiba2corr project will be documented in this file.

## [v1.1.0] - 2026-10-27

### Added

- Output unweighted Pearson and Spearman coefficients and P values alongside weighted Pearson results by default, using the same usable events.
- Add optional event-selection sensitivity analysis with `--sensitivity-analysis`, sampling events without replacement for each event type and the pooled `all` set.
- Support configurable retention fractions, repeat counts, and random seeds through `--subsample-fractions`, `--subsample-repeats`, and `--subsample-seed` (defaults: 0.25/0.5/0.75/0.9, 200 repeats, and seed 0).
- Export raw and summary sensitivity TSVs with full-data baselines, medians, 2.5–97.5 percentile ranges, and valid repeat counts. These ranges describe subsampling variability, not confidence intervals.
- Generate three-panel sensitivity figures in PDF and 600-dpi PNG formats. Sensitivity figures require Arial regular and bold fonts and report an error if they are unavailable.
- Record resolved run parameters in `report.json`.

### Changed

- Sort events by position ID to make event ordering and seeded subsampling reproducible with the same inputs and software environment.

### Fixed

- Exclude nonfinite dPSI pairs, including infinite values, from correlation analysis.

## [v1.0.0] - 2026-09-14

Initial release.
