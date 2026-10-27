"""Numerical, reproducibility and edge-case checks for event subsampling."""
import sys

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from shiba2corr.correlation import calculate_correlations
from shiba2corr.sensitivity import analyze_sensitivity
from shiba2corr import main as cli, publication


def frame(n=20):
	x = np.arange(n, dtype=float)
	return pd.DataFrame(dict(pos_id=[f'p{i}' for i in range(n)], target_dpsi=x,
		reference_dpsi=x + np.sin(x) * 5, weight=np.ones(n)))


def test_correlations_match_scipy_with_ties():
	x, y = [1, 1, 2, 3, 8], [4, 2, 2, 7, 5]
	result = calculate_correlations(x, y, np.ones(5))
	assert result['r'] == pytest.approx(stats.pearsonr(x, y)[0])
	assert (result['pearson_r'], result['pearson_p_value']) == pytest.approx(stats.pearsonr(x, y))
	assert (result['spearman_rho'], result['spearman_p_value']) == pytest.approx(stats.spearmanr(x, y))


def test_nonuniform_weighted_result():
	x, y, w = np.array([1, 2, 5]), np.array([4, 2, 6]), np.array([1, 2, 3])
	cov = np.cov(x, y, aweights=w)
	result = calculate_correlations(x, y, w)
	assert result['r'] == pytest.approx(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1]))
	assert result['neff'] == pytest.approx(w.sum() ** 2 / (w ** 2).sum())


@pytest.mark.parametrize('x,y', [([], []), ([1], [2]), ([1, 1, 1], [1, 2, 3])])
def test_undefined(x, y):
	r = calculate_correlations(x, y, np.ones(len(x)))
	assert all(np.isnan(r[k]) for k in ('r', 'pearson_r', 'spearman_rho'))


def test_reproducible_and_full_baseline():
	df = frame()
	raw, summary = analyze_sensitivity({'SE': df}, fractions=[0.5], repeats=12, seed=4)
	again, _ = analyze_sensitivity({'SE': df.iloc[::-1]}, fractions=[0.5], repeats=12, seed=4)
	pd.testing.assert_frame_equal(raw, again)
	other, _ = analyze_sensitivity({'SE': df}, fractions=[0.5], repeats=12, seed=5)
	assert not raw.correlation.equals(other.correlation)
	assert len(raw) == 13 * 3
	assert set(raw[raw.fraction == 0.5].n_sample) == {10}
	baseline = raw[(raw.fraction == 1) & (raw.method == 'Pearson')].correlation.iloc[0]
	assert baseline == pytest.approx(stats.pearsonr(df.target_dpsi, df.reference_dpsi)[0])
	for row in summary[summary.fraction < 1].itertuples():
		v = raw[(raw.fraction == row.fraction) & (raw.method == row.method)].correlation
		assert [row.lower, row.median, row.upper] == pytest.approx(v.quantile([.025, .5, .975]))


def test_sampling_without_replacement(monkeypatch):
	from shiba2corr import sensitivity
	original = sensitivity.calculate_correlations
	def checked(x, y, w, min_events):
		assert len(np.unique(x)) == len(x)
		return original(x, y, w, min_events)
	monkeypatch.setattr(sensitivity, 'calculate_correlations', checked)
	analyze_sensitivity({'SE': frame()}, repeats=5)


def test_small_empty_and_constant_sets():
	df = frame(3)
	df['target_dpsi'] = 1.
	raw, summary = analyze_sensitivity({'SE': df, 'RI': frame(0)}, fractions=[.25], repeats=3)
	assert raw.correlation.isna().all()
	assert (summary.n_valid == 0).all()
	assert (raw[raw.fraction == .25].n_sample == 0).all()


@pytest.mark.parametrize('option,value', [('--subsample-fractions', 'nan'),
	('--subsample-fractions', '1'), ('--subsample-repeats', '0'), ('--subsample-seed', '-1')])
def test_invalid_cli(monkeypatch, option, value):
	monkeypatch.setattr(sys, 'argv', ['shiba2corr', '-t', 'x', '-r', 'y', '-o', 'z', option, value])
	with pytest.raises(SystemExit):
		cli.parse_args()


def test_missing_arial_is_explicit(monkeypatch):
	def missing(*args, **kwargs):
		raise ValueError('missing')
	monkeypatch.setattr(publication.font_manager, 'findfont', missing)
	with pytest.raises(RuntimeError, match='Arial'):
		publication.publication_style()


def test_sensitivity_cli_outputs(monkeypatch, tmp_path, example_target_dir, example_reference_dir):
	# Numerical pipeline is portable; real Arial export is checked separately.
	monkeypatch.setattr(publication, 'publication_style', lambda: {})
	monkeypatch.setattr(cli.plots, 'plot_weighted_correlation_scatter', lambda *a, **k: None)
	monkeypatch.setattr(sys, 'argv', ['shiba2corr', '-t', example_target_dir, '-r', example_reference_dir,
		'-o', str(tmp_path), '--sensitivity-analysis', '--subsample-repeats', '3', '--no-venn'])
	cli.main()
	raw = pd.read_csv(tmp_path / 'results/event_selection_sensitivity.tsv', sep='\t')
	assert {'Weighted Pearson', 'Pearson', 'Spearman'} == set(raw.method)
	assert (tmp_path / 'plots/pdf/event_selection_sensitivity_all.pdf').exists()
	assert (tmp_path / 'plots/png/event_selection_sensitivity_all.png').exists()
