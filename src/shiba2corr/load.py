"""Loading and parsing utilities for Shiba PSI tables."""

import os
import pandas as pd
import numpy as np
import logging
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)

EVENT_TYPES = ('SE', 'FIVE', 'THREE', 'MXE', 'RI', 'MSE', 'AFE', 'ALE')
TTEST_P_COLUMNS = ('p_ttest',)


def coverage_percentile_diff_events(psi_file_df: pd.DataFrame) -> Tuple[pd.DataFrame, float, float]:
	"""Add mean_coverage and Beta-derived dPSI variance columns; return df + 5/95 coverage percentiles of DSEs."""
	ref_junction_cols = [c for c in psi_file_df.columns if c.startswith('ref_junction')]
	alt_junction_cols = [c for c in psi_file_df.columns if c.startswith('alt_junction')]
	ref_psi_col = 'ref_PSI'
	alt_psi_col = 'alt_PSI'

	pos_id_values = psi_file_df['pos_id'].values
	ref_junction_col_values = psi_file_df[ref_junction_cols].values
	alt_junction_col_values = psi_file_df[alt_junction_cols].values
	ref_psi_col_values = psi_file_df[ref_psi_col].values
	alt_psi_col_values = psi_file_df[alt_psi_col].values

	mean_coverage_list = []
	variance_list = []

	def beta_variance(psi, coverage):
		if coverage <= 1 or pd.isna(psi) or psi <= 0 or psi >= 1:
			return np.nan
		alpha = psi * coverage
		beta_param = (1 - psi) * coverage
		numerator = alpha * beta_param
		denominator = (alpha + beta_param) ** 2 * (alpha + beta_param + 1)
		return numerator / denominator

	for i in range(len(pos_id_values)):
		ref_coverage = sum(sum(int(v) for v in str(x).split(';')) for x in ref_junction_col_values[i] if pd.notna(x))
		alt_coverage = sum(sum(int(v) for v in str(x).split(';')) for x in alt_junction_col_values[i] if pd.notna(x))
		mean_coverage = (ref_coverage + alt_coverage) / 2 if (ref_coverage + alt_coverage) > 0 else 0
		mean_coverage_list.append(mean_coverage)

		ref_psi = float(ref_psi_col_values[i]) if pd.notna(ref_psi_col_values[i]) else np.nan
		alt_psi = float(alt_psi_col_values[i]) if pd.notna(alt_psi_col_values[i]) else np.nan
		ref_var = beta_variance(ref_psi, ref_coverage)
		alt_var = beta_variance(alt_psi, alt_coverage)
		variance = ref_var + alt_var if (pd.notna(ref_var) and pd.notna(alt_var)) else np.nan
		variance_list.append(variance)

	psi_file_df = psi_file_df.copy()
	psi_file_df['mean_coverage'] = mean_coverage_list
	psi_file_df['variance'] = variance_list
	psi_file_diff_df = psi_file_df[psi_file_df['Diff events'] == 'Yes'].reset_index(drop=True)
	if psi_file_diff_df.empty:
		return (psi_file_df, 0, 0)
	mean_coverage_diff_list = psi_file_diff_df['mean_coverage'].values
	return (psi_file_df, np.percentile(mean_coverage_diff_list, 5), np.percentile(mean_coverage_diff_list, 95))


def parse_psi_table(psi_file_path: str, event_type: str) -> pd.DataFrame:
	if not os.path.exists(psi_file_path):
		logger.warning(f"PSI file not found for {event_type}: {psi_file_path}")
		return pd.DataFrame()
	try:
		df = pd.read_csv(psi_file_path, sep="\t")
		df, _, _ = coverage_percentile_diff_events(df)
		return df
	except Exception as e:
		logger.error(f"Error parsing PSI table for {event_type}: {e}")
		return pd.DataFrame()


def build_maps(df: pd.DataFrame) -> Tuple[Dict[str, float], Dict[str, float], Dict[str, float], Set[str]]:
	dpsi_map, coverage_map, variance_map, dse_set = {}, {}, {}, set()
	if df.empty:
		return dpsi_map, coverage_map, variance_map, dse_set
	required_cols = {'pos_id', 'dPSI', 'mean_coverage', 'variance', 'Diff events'}
	if not required_cols.issubset(df.columns):
		logger.warning(f"Missing required columns: {required_cols - set(df.columns)}")
		return dpsi_map, coverage_map, variance_map, dse_set
	for _, row in df.iterrows():
		pos_id = str(row['pos_id'])
		dpsi_map[pos_id] = float(row['dPSI']) if pd.notna(row['dPSI']) else np.nan
		coverage_map[pos_id] = float(row['mean_coverage']) if pd.notna(row['mean_coverage']) else 0.0
		variance_map[pos_id] = float(row['variance']) if pd.notna(row['variance']) else np.nan
		if row.get('Diff events') == 'Yes':
			dse_set.add(pos_id)
	return dpsi_map, coverage_map, variance_map, dse_set


def split_dse_by_direction(dse_set: Set[str], dpsi_map: Dict[str, float]) -> Tuple[Set[str], Set[str]]:
	up_set, down_set = set(), set()
	for pos_id in dse_set:
		dpsi = dpsi_map.get(pos_id, np.nan)
		if not np.isfinite(dpsi):
			continue
		if dpsi > 0:
			up_set.add(pos_id)
		elif dpsi < 0:
			down_set.add(pos_id)
	logger.debug(f"Split {len(dse_set)} DSEs: {len(up_set)} up, {len(down_set)} down")
	return up_set, down_set


def aggregate_all_types(condition_results: Dict[str, Dict]) -> Dict[str, any]:
	all_dpsi, all_coverage, all_variance = {}, {}, {}
	all_dse_set = set()
	for event_type in EVENT_TYPES:
		if event_type not in condition_results:
			continue
		ed = condition_results[event_type]
		for k, v in ed['dpsi'].items():
			all_dpsi.setdefault(k, v)
		for k, v in ed['coverage'].items():
			all_coverage.setdefault(k, v)
		for k, v in ed['variance'].items():
			all_variance.setdefault(k, v)
		all_dse_set.update(ed['dse_set'])
	return {
		'dpsi': all_dpsi,
		'coverage': all_coverage,
		'variance': all_variance,
		'dse_set': all_dse_set,
	}


def compute_detected_events_intersection(target_results: Dict, reference_results: Dict,
										min_coverage: Optional[float] = None) -> Dict[str, Set[str]]:
	"""Background = events with finite stats in BOTH conditions."""
	detected_intersection = {}
	min_cov = min_coverage or 0.0
	for event_type in list(EVENT_TYPES) + ['all']:
		t = target_results.get(event_type, {})
		r = reference_results.get(event_type, {})
		t_det = {p for p in t.get('dpsi', {})
				 if t.get('coverage', {}).get(p, 0.0) >= min_cov
				 and np.isfinite(t.get('dpsi', {}).get(p, np.nan))
				 and np.isfinite(t.get('variance', {}).get(p, np.nan))}
		r_det = {p for p in r.get('dpsi', {})
				 if r.get('coverage', {}).get(p, 0.0) >= min_cov
				 and np.isfinite(r.get('dpsi', {}).get(p, np.nan))
				 and np.isfinite(r.get('variance', {}).get(p, np.nan))}
		detected_intersection[event_type] = t_det & r_det
	return detected_intersection


def filter_to_union(condition_results: Dict, union_sets: Dict[str, Set[str]]) -> Dict:
	filtered = {}
	for event_type in list(EVENT_TYPES) + ['all']:
		if event_type not in condition_results:
			continue
		keep = union_sets.get(event_type, set())
		ed = condition_results[event_type]
		filtered[event_type] = {
			'dpsi': {p: v for p, v in ed['dpsi'].items() if p in keep},
			'coverage': {p: v for p, v in ed['coverage'].items() if p in keep},
			'variance': {p: v for p, v in ed['variance'].items() if p in keep},
		}
	return filtered


def check_ttest_p_available(shiba_dir: str, event_types: List[str]) -> Optional[str]:
	for event_type in event_types:
		psi_file_path = os.path.join(shiba_dir, "results", "splicing", f"PSI_{event_type}.txt")
		if os.path.exists(psi_file_path):
			try:
				header = pd.read_csv(psi_file_path, sep='\t', nrows=0).columns.tolist()
				for col_name in TTEST_P_COLUMNS:
					if col_name in header:
						return col_name
				return None
			except Exception:
				continue
	return None


def load_union_event_data(target_dir: str,
						 reference_dir: str,
						 event_types: Optional[List[str]] = None,
						 min_coverage: Optional[float] = None,
						 return_df: bool = True,
						 ttest_p_threshold: Optional[float] = None) -> Dict:
	"""Load pos_id->(dPSI, coverage, variance) maps and DSE sets per condition.

	The 'union' key holds the union of Target/Reference DSE sets per event type;
	'background' holds the intersection of detected (finite-stats) events used to
	filter the maps and to compose the optional background_df.
	"""
	if event_types is None:
		event_types = list(EVENT_TYPES)

	logger.info(f"Loading event data from target: {target_dir}, reference: {reference_dir}")

	ttest_col_name = None
	if ttest_p_threshold is not None:
		t_col = check_ttest_p_available(target_dir, event_types)
		r_col = check_ttest_p_available(reference_dir, event_types)
		if t_col is None or r_col is None:
			missing = [n for n, c in (("target", t_col), ("reference", r_col)) if c is None]
			logger.warning(
				f"ttest P-value column not found in {' and '.join(missing)} PSI files. "
				f"Ignoring --ttest threshold ({ttest_p_threshold}).")
			ttest_p_threshold = None
		else:
			ttest_col_name = t_col
			logger.info(f"Applying t-test P-value filter: {ttest_col_name} <= {ttest_p_threshold}")

	def load_single_condition(shiba_dir: str, condition_name: str) -> Dict:
		results = {}
		for event_type in event_types:
			psi_file_path = os.path.join(shiba_dir, "results", "splicing", f"PSI_{event_type}.txt")
			df = parse_psi_table(psi_file_path, event_type)
			dpsi_map, coverage_map, variance_map, dse_set = build_maps(df)
			if ttest_p_threshold is not None and ttest_col_name is not None and ttest_col_name in df.columns:
				ttest_map = dict(zip(df['pos_id'].astype(str), df[ttest_col_name]))
				before = len(dse_set)
				dse_set = {p for p in dse_set
						   if p in ttest_map and pd.notna(ttest_map[p]) and ttest_map[p] <= ttest_p_threshold}
				if before - len(dse_set) > 0:
					logger.debug(f"Filtered {before - len(dse_set)}/{before} DSEs by {ttest_col_name} > {ttest_p_threshold} "
								f"for {event_type} in {condition_name}")
			results[event_type] = {
				'dpsi': dpsi_map,
				'coverage': coverage_map,
				'variance': variance_map,
				'dse_set': dse_set,
			}
		results['all'] = aggregate_all_types(results)
		return results

	target_results = load_single_condition(target_dir, "target")
	reference_results = load_single_condition(reference_dir, "reference")

	background_sets = compute_detected_events_intersection(target_results, reference_results, min_coverage)
	union_sets = {et: target_results.get(et, {}).get('dse_set', set()) | reference_results.get(et, {}).get('dse_set', set())
				  for et in list(EVENT_TYPES) + ['all']}

	target_filtered = filter_to_union(target_results, background_sets)
	reference_filtered = filter_to_union(reference_results, background_sets)

	result = {
		'target': target_filtered,
		'reference': reference_filtered,
		'union': union_sets,
		'background': background_sets,
		'dse_target': {et: target_results[et]['dse_set'] for et in list(EVENT_TYPES) + ['all']},
		'dse_reference': {et: reference_results[et]['dse_set'] for et in list(EVENT_TYPES) + ['all']},
		'target_unfiltered': target_results,
		'reference_unfiltered': reference_results,
	}

	if return_df:
		result['background_df'] = create_background_dataframe(
			target_filtered, reference_filtered, background_sets, min_coverage)

	logger.info(f"Loaded event data: {len(result['target']['all']['dpsi'])} target events, "
				f"{len(background_sets['all'])} background events (intersection of detected events)")
	return result


def create_background_dataframe(target_data: Dict, reference_data: Dict,
								background_sets: Dict[str, Set[str]],
								min_coverage: Optional[float] = None) -> pd.DataFrame:
	"""Create a unified background DataFrame using the detected-event intersection."""
	rows = []
	for event_type in list(EVENT_TYPES) + ['all']:
		if event_type not in background_sets:
			continue
		for pos_id in background_sets[event_type]:
			t_dpsi = target_data.get(event_type, {}).get('dpsi', {}).get(pos_id, np.nan)
			r_dpsi = reference_data.get(event_type, {}).get('dpsi', {}).get(pos_id, np.nan)
			t_cov = target_data.get(event_type, {}).get('coverage', {}).get(pos_id, 0.0)
			r_cov = reference_data.get(event_type, {}).get('coverage', {}).get(pos_id, 0.0)
			t_var = target_data.get(event_type, {}).get('variance', {}).get(pos_id, np.nan)
			r_var = reference_data.get(event_type, {}).get('variance', {}).get(pos_id, np.nan)
			eff_cov = (np.sqrt(t_cov * r_cov) if t_cov > 0 and r_cov > 0 else 0.0)
			if min_coverage is not None and eff_cov < min_coverage:
				continue
			rows.append({
				'pos_id': pos_id, 'event_type': event_type,
				'dPSI_tgt': t_dpsi, 'dPSI_ref': r_dpsi,
				'Var_tgt': t_var, 'Var_ref': r_var,
				'cov_tgt': t_cov, 'cov_ref': r_cov,
				'eff_cov': eff_cov,
			})
	return pd.DataFrame(rows)


def compute_dse_overlap_by_direction(union_data: Dict,
									event_types: Optional[List[str]] = None) -> Dict[str, Dict[str, Dict]]:
	"""Compute DSE overlap counts (Target vs Reference) split by dPSI sign."""
	if event_types is None:
		event_types = list(EVENT_TYPES) + ['all']

	target_data = union_data['target']
	reference_data = union_data['reference']
	dse_target = union_data['dse_target']
	dse_reference = union_data['dse_reference']
	target_unfiltered = union_data.get('target_unfiltered')
	reference_unfiltered = union_data.get('reference_unfiltered')

	overlap_data = {}
	for et in event_types:
		tgt_dse = dse_target.get(et, set())
		ref_dse = dse_reference.get(et, set())
		tgt_dpsi = (target_unfiltered.get(et, {}).get('dpsi', {})
					if target_unfiltered is not None and et in target_unfiltered
					else target_data.get(et, {}).get('dpsi', {}))
		ref_dpsi = (reference_unfiltered.get(et, {}).get('dpsi', {})
					if reference_unfiltered is not None and et in reference_unfiltered
					else reference_data.get(et, {}).get('dpsi', {}))
		tgt_up, tgt_down = split_dse_by_direction(tgt_dse, tgt_dpsi)
		ref_up, ref_down = split_dse_by_direction(ref_dse, ref_dpsi)

		overlap_data[et] = {}
		for direction, tgt_set, ref_set in [('up', tgt_up, ref_up), ('down', tgt_down, ref_down)]:
			intersection = tgt_set & ref_set
			overlap_data[et][direction] = {
				'target_count': len(tgt_set),
				'reference_count': len(ref_set),
				'intersection_count': len(intersection),
				'target_only': len(tgt_set - ref_set),
				'reference_only': len(ref_set - tgt_set),
			}
	return overlap_data
