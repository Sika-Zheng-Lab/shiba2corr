"""Random event subsampling, conditional on the filtered, usable DSE union."""

import numpy as np
import pandas as pd

from .correlation import calculate_correlations

METHODS = {'r': 'Weighted Pearson', 'pearson_r': 'Pearson', 'spearman_rho': 'Spearman'}
DEFAULT_FRACTIONS = (0.25, 0.5, 0.75, 0.9)


def validate_settings(fractions, repeats, seed):
	if not fractions or any(not np.isfinite(f) or not 0 < f < 1 for f in fractions):
		raise ValueError('Subsample fractions must be finite and strictly between 0 and 1.')
	if len(set(fractions)) != len(fractions):
		raise ValueError('Subsample fractions must be unique.')
	if not isinstance(repeats, int) or repeats < 1:
		raise ValueError('Subsample repeats must be a positive integer.')
	if not isinstance(seed, int) or seed < 0:
		raise ValueError('Subsample seed must be a nonnegative integer.')


def analyze_sensitivity(scatter_data, fractions=DEFAULT_FRACTIONS, repeats=200,
						seed=0, min_events=3):
	"""Sample rows uniformly without replacement; retain original event weights.

	Each type (including pooled `all`) is sampled independently. Sorting makes
	results invariant to input row/dictionary order. The full-data baseline is
	stored once as fraction=1, repeat=0. Small/constant subsets remain NaN.
	"""
	fractions = tuple(fractions)
	validate_settings(fractions, repeats, seed)
	rows = []
	for event_type in sorted(scatter_data):
		df = scatter_data[event_type].sort_values('pos_id').reset_index(drop=True)
		x = df['target_dpsi'].to_numpy(dtype=float)
		y = df['reference_dpsi'].to_numpy(dtype=float)
		w = df['weight'].to_numpy(dtype=float)
		if not (np.isfinite(x).all() and np.isfinite(y).all()
				and np.isfinite(w).all() and (w > 0).all()):
			raise ValueError('Sensitivity input must contain finite pairs and positive finite weights.')
		# Stable independent stream for each event type, unaffected by missing types.
		rng = np.random.default_rng(np.random.SeedSequence([seed, *event_type.encode()]))
		for fraction in (*sorted(fractions), 1.0):
			n = int(np.floor(len(df) * fraction))
			for repeat in (range(1, repeats + 1) if fraction < 1 else [0]):
				idx = rng.choice(len(df), size=n, replace=False) if fraction < 1 else np.arange(n)
				values = calculate_correlations(x[idx], y[idx], w[idx], min_events)
				for key, method in METHODS.items():
					rows.append(dict(event_type=event_type, fraction=fraction,
						repeat=repeat, n_total=len(df), n_sample=n, method=method,
						correlation=values[key], seed=seed))
	columns = ['event_type', 'fraction', 'repeat', 'n_total', 'n_sample', 'method', 'correlation', 'seed']
	raw = pd.DataFrame(rows, columns=columns)
	summaries = []
	for keys, group in raw.groupby(['event_type', 'fraction', 'method'], sort=False):
		values = group['correlation'].dropna()
		q = values.quantile([0.025, 0.5, 0.975]).to_numpy() if len(values) else [np.nan] * 3
		summaries.append(dict(event_type=keys[0], fraction=keys[1], method=keys[2],
			n_total=group.n_total.iloc[0], n_sample=group.n_sample.iloc[0],
			n_repeats=len(group), n_valid=len(values), lower=q[0], median=q[1], upper=q[2]))
	summary = pd.DataFrame(summaries, columns=['event_type', 'fraction', 'method',
		'n_total', 'n_sample', 'n_repeats', 'n_valid', 'lower', 'median', 'upper'])
	return raw, summary
