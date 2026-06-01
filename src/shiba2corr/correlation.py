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
		for pos_id in union_pos_ids:
			t_dpsi = target_data[et]['dpsi'].get(pos_id, np.nan)
			r_dpsi = reference_data[et]['dpsi'].get(pos_id, np.nan)
			if pd.isna(t_dpsi) or pd.isna(r_dpsi):
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

		if n_used < 2:
			results[et] = {'n_union': n_union, 'n_used': n_used, 'weight_scheme': weight_scheme,
						   'r': np.nan, 'p_value_approx': np.nan, 'neff': np.nan}
			continue

		x = np.asarray(xs)
		y = np.asarray(ys)
		w = np.asarray(ws)
		try:
			w_sum = w.sum()
			w2_sum = (w ** 2).sum()
			denom = max(1e-12, w_sum - (w2_sum / w_sum))
			mx = (w * x).sum() / w_sum
			my = (w * y).sum() / w_sum
			cov = (w * (x - mx) * (y - my)).sum() / denom
			vx = (w * (x - mx) ** 2).sum() / denom
			vy = (w * (y - my) ** 2).sum() / denom
			r = cov / np.sqrt(vx * vy) if vx > 0 and vy > 0 else np.nan
			neff = (w_sum ** 2) / w2_sum if w2_sum > 0 else np.nan
			if not np.isfinite(neff) or neff < max(3, min_events):
				p_value = np.nan
			elif np.isfinite(r):
				denom_t = max(1e-12, 1.0 - r ** 2)
				tval = r * np.sqrt((neff - 2.0) / denom_t)
				p_value = 2.0 * stats.t.sf(np.abs(tval), df=max(1.0, neff - 2.0))
			else:
				p_value = np.nan
		except (ZeroDivisionError, ValueError) as e:
			logger.warning(f"Error calculating correlation for {et}: {e}")
			r, p_value, neff = np.nan, np.nan, np.nan

		results[et] = {
			'n_union': n_union, 'n_used': n_used, 'weight_scheme': weight_scheme,
			'r': r, 'p_value_approx': p_value, 'neff': neff,
		}

	logger.info(f"Completed weighted dPSI correlation analysis using {weight_scheme} weighting")
	return results, scatter_data_dict


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
