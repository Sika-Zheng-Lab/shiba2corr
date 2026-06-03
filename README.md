# shiba2corr

[![GitHub License](https://img.shields.io/github/license/Sika-Zheng-Lab/shiba2corr)](https://github.com/Sika-Zheng-Lab/shiba2corr/blob/main/LICENSE)
[![GitHub Release](https://img.shields.io/github/v/release/Sika-Zheng-Lab/shiba2corr?style=flat)](https://github.com/Sika-Zheng-Lab/shiba2corr/releases)
[![GitHub Release Date](https://img.shields.io/github/release-date/Sika-Zheng-Lab/shiba2corr)](https://github.com/Sika-Zheng-Lab/shiba2corr/releases)
[![Tests](https://github.com/Sika-Zheng-Lab/shiba2corr/actions/workflows/test.yml/badge.svg)](https://github.com/Sika-Zheng-Lab/shiba2corr/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/Sika-Zheng-Lab/shiba2corr/branch/main/graph/badge.svg)](https://codecov.io/gh/Sika-Zheng-Lab/shiba2corr)
[![Create Release](https://github.com/Sika-Zheng-Lab/shiba2corr/actions/workflows/release.yaml/badge.svg)](https://github.com/Sika-Zheng-Lab/shiba2corr/actions/workflows/release.yaml)
[![Publish to PyPI](https://github.com/Sika-Zheng-Lab/shiba2corr/actions/workflows/publish.yaml/badge.svg)](https://github.com/Sika-Zheng-Lab/shiba2corr/actions/workflows/publish.yaml)
[![Python](https://img.shields.io/pypi/pyversions/shiba2corr.svg?label=Python&color=blue)](https://pypi.org/project/shiba2corr/)
[![PyPI](https://img.shields.io/pypi/v/shiba2corr.svg?label=PyPI&color=orange)](https://pypi.org/project/shiba2corr/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/shiba2corr.svg?label=PyPI%20-%20Downloads&color=orange)](https://pypi.org/project/shiba2corr/)
[![Docker](https://img.shields.io/docker/v/naotokubota/shiba2corr?color=blue&label=Docker)](https://hub.docker.com/r/naotokubota/shiba2corr)
[![Docker Pulls](https://img.shields.io/docker/pulls/naotokubota/shiba2corr)](https://hub.docker.com/r/naotokubota/shiba2corr)
[![Docker Image Size](https://img.shields.io/docker/image-size/naotokubota/shiba2corr)](https://hub.docker.com/r/naotokubota/shiba2corr)

Quantify the *reproducibility* of differentially spliced events between two
[Shiba](https://github.com/Sika-Zheng-Lab/Shiba) experiments via a
**weighted Pearson correlation of dPSI values**, computed on the union of
differentially spliced events (DSEs) and weighted by the reliability of each
PSI estimate.

## Overview

Two Shiba runs (a **Target** condition and a **Reference** condition) produce
per-event PSI tables. A naive Pearson correlation of their dPSI values
overweights noisy events with low junction coverage. `shiba2corr` instead
computes a **weighted Pearson correlation** on the union of DSEs from both
runs, where the weight of each event reflects how precisely its PSI was
estimated in each condition.

**Weighted correlation:**

$$r_w = \frac{\mathrm{Cov}_w(x, y)}{\sqrt{\mathrm{Var}_w(x)\ \mathrm{Var}_w(y)}}$$

with unbiased weighted estimators

$$\bar{x}_w = \frac{\sum w_i x_i}{\sum w_i},\quad
\mathrm{Cov}_w(x,y) = \frac{\sum w_i (x_i - \bar{x}_w)(y_i - \bar{y}_w)}{\sum w_i - \frac{\sum w_i^2}{\sum w_i}}$$

and significance assessed via a t-test on the **Kish effective sample size**

$$n_{\text{eff}} = \frac{(\sum w_i)^2}{\sum w_i^2}.$$

<table>
  <thead>
    <tr><th>Weight scheme</th><th>$w_i$</th><th>Use when</th></tr>
  </thead>
  <tbody>
    <tr>
      <td><code>inverse_variance</code> <em>(default)</em></td>
      <td>$1 / (\mathrm{Var}(\mathrm{PSI}_{\mathrm{tgt}}) + \mathrm{Var}(\mathrm{PSI}_{\mathrm{ref}}) + \varepsilon)$</td>
      <td>Beta-derived PSI variances are reliable</td>
    </tr>
    <tr>
      <td><code>geom_mean</code></td>
      <td>$\sqrt{\mathrm{cov}_{\mathrm{tgt}} \cdot \mathrm{cov}_{\mathrm{ref}}}$</td>
      <td>Coverage is the only confidence proxy</td>
    </tr>
    <tr>
      <td><code>coverage_mean</code></td>
      <td>$\tfrac{1}{2}(\mathrm{cov}_{\mathrm{tgt}} + \mathrm{cov}_{\mathrm{ref}})$</td>
      <td>One condition may have low coverage</td>
    </tr>
    <tr>
      <td><code>uniform</code></td>
      <td>$1$</td>
      <td>Sanity-check vs. unweighted Pearson</td>
    </tr>
  </tbody>
</table>

In addition to the correlation table, `shiba2corr` produces:

- per-event-type **KDE-coloured scatter plots** of Target vs Reference dPSI;
- **DSE overlap Venn diagrams** (Target vs Reference, split by dPSI direction)
  for each event type and a combined grid.

## Installation

```bash
pip install shiba2corr
```

For development:

```bash
git clone https://github.com/Sika-Zheng-Lab/shiba2corr.git
cd shiba2corr
pip install -e ".[dev]"
```

## Usage

### Quick start with example data

A small synthetic Shiba result pair is included in [`example/`](example/):

```bash
shiba2corr \
    -t example/target \
    -r example/reference \
    -o example/output
```

### Basic usage

```bash
shiba2corr -t target_shiba_dir -r reference_shiba_dir -o output/
```

### Filter DSEs by t-test P-value

If Shiba was run with t-tests enabled (column `p_ttest` present in the PSI
files), restrict the analysis to events that also pass a P-value threshold:

```bash
shiba2corr -t target/ -r reference/ -o output/ --ttest 0.05
```

### Choose a different weighting scheme

```bash
shiba2corr -t target/ -r reference/ -o output/ --weight-scheme geom_mean
```

### Custom colors and font

```bash
shiba2corr -t target/ -r reference/ -o output/ \
    --target-color "#FF0000FF" \
    --reference-color "#00008BFF" \
    --font-family Arial
```

### Skip Venn diagrams

```bash
shiba2corr -t target/ -r reference/ -o output/ --no-venn
```

## Input File Format

`-t` / `-r` accept the path to a **Shiba working directory**. Within each
directory, `shiba2corr` reads:

```
<shiba_dir>/results/splicing/PSI_SE.txt
<shiba_dir>/results/splicing/PSI_FIVE.txt
<shiba_dir>/results/splicing/PSI_THREE.txt
<shiba_dir>/results/splicing/PSI_MXE.txt
<shiba_dir>/results/splicing/PSI_RI.txt
<shiba_dir>/results/splicing/PSI_MSE.txt
<shiba_dir>/results/splicing/PSI_AFE.txt
<shiba_dir>/results/splicing/PSI_ALE.txt
```

Each PSI table is a tab-separated file with at least the following columns:

| Column | Description |
|---|---|
| `pos_id` | Unique event ID (matched across conditions) |
| `dPSI` | Delta PSI vs. control samples |
| `ref_PSI`, `alt_PSI` | PSI values used to compute `dPSI` |
| `ref_junction*`, `alt_junction*` | Per-sample junction read counts (semicolon-delimited; one column per replicate) |
| `Diff events` | `Yes` / `No`, marking DSEs |
| `p_ttest` *(optional)* | T-test P-value; required only with `--ttest` |

Missing event-type files are skipped with a warning.

### `-l/--event-list` (optional)

A plain-text file with one `pos_id` per line. When provided, the weighted
correlation is restricted to the intersection of the union DSE set and this
list.

## Output Files

Given `-o output/`:

| File | Description |
|---|---|
| `results/event_dpsi_weighted_correlation.tsv` | Per-event-type summary (`n_union`, `n_used`, `r`, `p_value_approx`, `neff`) |
| `results/weighted_correlation_scatter_{TYPE}.tsv` | Per-event-type per-event scatter data (`pos_id`, dPSI values, weights) |
| `plots/png/weighted_correlation_scatter_{TYPE}.png` | KDE-coloured scatter plot (also `.pdf` under `plots/pdf/`) |
| `plots/png/dse_venn_{TYPE}_{up,down}.png` | DSE overlap Venn diagram per event type and direction |
| `plots/png/dse_venn_grid_{up,down}.png` | Combined Venn diagram grid per direction |
| `report.json` | Machine-readable run summary (version, timestamp, command line) |

`TYPE` is one of `SE`, `FIVE`, `THREE`, `MXE`, `RI`, `MSE`, `AFE`, `ALE`,
or `all` (aggregated across event types). `--no-venn` skips all DSE overlap
Venn output.

## CLI Options

| Option | Description |
|---|---|
| `-v`, `--version` | Show version and exit |
| `-t`, `--target` | Target Shiba working directory (required) |
| `-r`, `--reference` | Reference Shiba working directory (required) |
| `-o`, `--output` | Output directory (required) |
| `-l`, `--event-list` | Restrict analysis to events listed in this file (optional) |
| `--weight-scheme` | `inverse_variance` *(default)*, `geom_mean`, `coverage_mean`, `uniform` |
| `--min-events` | Minimum effective sample size for the t-test p-value (default: 3) |
| `--ttest` | Filter DSEs by t-test P-value threshold (requires `p_ttest` column) |
| `--target-color` | Hex color for Target DSEs in Venn diagrams |
| `--reference-color` | Hex color for Reference DSEs in Venn diagrams |
| `--font-family` | Font family for plot text (e.g. `Arial`) |
| `--no-venn` | Skip the DSE overlap Venn diagrams |
| `--verbose` | Enable DEBUG-level logging |

## License

MIT License

## Citation

If you use `shiba2corr` in your research, please cite this repository.

## Contributing

Thank you for wanting to improve shiba2corr! If you have any bugs or
questions, feel free to [open an issue](https://github.com/Sika-Zheng-Lab/shiba2corr/issues)
or pull request. See [CONTRIBUTING.md](CONTRIBUTING.md) for the development
workflow and release process.

## Authors

- Naoto Kubota ([0000-0003-0612-2300](https://orcid.org/0000-0003-0612-2300))
- Sika Zheng ([0000-0002-0573-4981](https://orcid.org/0000-0002-0573-4981))
