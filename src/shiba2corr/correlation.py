"""Weighted dPSI correlation analysis."""

import os
import logging
import numpy as np
import pandas as pd
import scipy.stats as stats

from . import load

logger = logging.getLogger(__name__)

WEIGHT_SCHEMES = ('inverse_variance', 'geom_mean', 'coverage_mean', 'uniform')


def compute_weighted_correlation_dpsi_union(target_shiba_dir, reference_shiba_dir,
											event_list=None,
											weight_scheme='inverse_variance',
											min_events=3,
											epsilon=1e-8,
											ttest_p_threshold=None):
	"""Compute the weighted dPSI correlation between Target and Reference over the DSE union.

	Returns (correlation_results, scatter_data_dict).
	"""
	if weight_scheme not in WEIGHT_SCHEMES:
		raise ValueError(f"Unknown weight_scheme '{weight_scheme}'. Choose from {WEIGHT_SCHEMES}.")

	corr_data = load.load_union_event_data(target_shiba_dir, reference_shiba_dir,
										   return_df=False, ttest_p_threshold=ttest_p_threshold)
	target_data = corr_data['target']
	reference_data = corr_data['reference']
	union_dse = corr_data['union']

	if event_list is not None:
		logger.info(f"Filtering union DSEs using provided event list: {event_list}")
		with open(event_list, 'r') as f:
			specified = set(line.strip() for line in f if line.strip())
		for et in union_dse:
			before = len(union_dse[et])
			union_dse[et] = {p for p in union_dse[et] if p in specified}
			logger.info(f"Event type {et}: filtered union DSEs from {before} to {len(union_dse[et])}")

	results = {}
	scatter_data_dict = {}
	event_type_list = list(load.EVENT_TYPES) + ['all']

	for et in event_type_list:
		logger.info(f"Computing weighted correlation for event type: {et}")
		union_pos_ids = union_dse.get(et, set())
		n_union = len(union_pos_ids)

		xs, ys, ws = [], [], []
		scatter_rows = []
		for pos_id in sorted(union_pos_ids):
			t_dpsi = target_data[et]['dpsi'].get(pos_id, np.nan)
			r_dpsi = reference_data[et]['dpsi'].get(pos_id, np.nan)
			if not np.isfinite(t_dpsi) or not np.isfinite(r_dpsi):
				continue

			if weight_scheme == 'inverse_variance':
				t_var = target_data[et]['variance'].get(pos_id, np.nan)
				r_var = reference_data[et]['variance'].get(pos_id, np.nan)
				if pd.isna(t_var) or pd.isna(r_var):
					continue
				weight = 1.0 / (t_var + r_var + epsilon)
			elif weight_scheme == 'geom_mean':
				t_cov = target_data[et]['coverage'].get(pos_id, 0.0)
				r_cov = reference_data[et]['coverage'].get(pos_id, 0.0)
				if t_cov <= 0 or r_cov <= 0:
					continue
				weight = float(np.sqrt(t_cov * r_cov))
			elif weight_scheme == 'coverage_mean':
				t_cov = target_data[et]['coverage'].get(pos_id, 0.0)
				r_cov = reference_data[et]['coverage'].get(pos_id, 0.0)
				if t_cov <= 0 and r_cov <= 0:
					continue
				weight = 0.5 * (t_cov + r_cov)
			else:  # uniform
				weight = 1.0

			if np.isfinite(weight) and weight > 0:
				xs.append(float(t_dpsi))
				ys.append(float(r_dpsi))
				ws.append(float(weight))
				scatter_rows.append({
					'pos_id': pos_id,
					'target_dpsi': float(t_dpsi),
					'reference_dpsi': float(r_dpsi),
					'weight': float(weight),
				})

		n_used = len(xs)
		if scatter_rows:
			df_tmp = pd.DataFrame(scatter_rows)
			wsum = df_tmp['weight'].sum()
			df_tmp['weight_normalized'] = df_tmp['weight'] / wsum if wsum > 0 else 0.0
			scatter_data_dict[et] = df_tmp
		else:
			scatter_data_dict[et] = pd.DataFrame(
				columns=['pos_id', 'target_dpsi', 'reference_dpsi', 'weight', 'weight_normalized'])

		results[et] = {
			'n_union': n_union, 'n_used': n_used, 'weight_scheme': weight_scheme,
			**calculate_correlations(xs, ys, ws, min_events=min_events),
		}

	logger.info(f"Completed weighted dPSI correlation analysis using {weight_scheme} weighting")
	return results, scatter_data_dict


def calculate_correlations(x, y, weights, min_events=3):
	"""Calculate all three correlations on the same finite, positive-weight pairs.

	Undefined correlations (fewer than two pairs or constant values) are NaN.
	min_events applies only to the approximate weighted-correlation P value.
	"""
	x, y, w = (np.asarray(a, dtype=float) for a in (x, y, weights))
	valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
	x, y, w = x[valid], y[valid], w[valid]
	result = dict.fromkeys(('r', 'p_value_approx', 'neff', 'pearson_r',
		'pearson_p_value', 'spearman_rho', 'spearman_p_value'), np.nan)
	if len(x) < 2:
		return result
	w = w / w.max()
	w = w / w.sum()
	neff = 1.0 / np.sum(w ** 2)
	result['neff'] = neff
	if np.ptp(x) == 0 or np.ptp(y) == 0:
		return result
	dx, dy = x - np.sum(w * x), y - np.sum(w * y)
	vx, vy = np.sum(w * dx ** 2), np.sum(w * dy ** 2)
	r = np.clip(np.sum(w * dx * dy) / np.sqrt(vx * vy), -1., 1.)
	result['r'] = r
	if neff >= max(3, min_events):
		tval = r * np.sqrt((neff - 2) / max(1e-12, 1 - r ** 2))
		result['p_value_approx'] = 2 * stats.t.sf(abs(tval), df=neff - 2)
	result['pearson_r'], result['pearson_p_value'] = stats.pearsonr(x, y)
	result['spearman_rho'], result['spearman_p_value'] = stats.spearmanr(x, y)
	return result


def save_weighted_correlation_results(wcorr_results, output_file):
	rows = []
	for et, result in wcorr_results.items():
		rows.append({
			'event_type': et,
			'n_union': result.get('n_union', np.nan),
			'n_used': result.get('n_used', np.nan),
			'r': result.get('r', np.nan),
			'p_value_approx': result.get('p_value_approx', np.nan),
			'neff': result.get('neff', np.nan),
			**{key: result.get(key, np.nan) for key in (
				'pearson_r', 'pearson_p_value', 'spearman_rho', 'spearman_p_value')},
		})
	pd.DataFrame(rows).to_csv(output_file, sep='\t', index=False)
	logger.info(f"Weighted dPSI correlation results saved to {output_file}")


def save_scatter_plot_data(scatter_data_dict, output_dir):
	for et, df in scatter_data_dict.items():
		if df.empty:
			continue
		out = os.path.join(output_dir, f"weighted_correlation_scatter_{et}.tsv")
		df.to_csv(out, sep='\t', index=False)
	logger.info(f"Scatter data saved to {output_dir}/weighted_correlation_scatter_*.tsv")
