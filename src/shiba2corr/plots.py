"""Plotting helpers for shiba2corr."""

import math
import logging
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import gaussian_kde
from pathlib import Path

logger = logging.getLogger(__name__)

try:
	from matplotlib_venn import venn2
	VENN_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised only when matplotlib_venn is missing
	VENN_AVAILABLE = False
	logger.warning("matplotlib_venn is not available. Venn diagram features will be disabled.")
	logger.warning("Install with: pip install matplotlib-venn")

TARGET_COLOR = '#E49414FF'
REFERENCE_COLOR = '#8087AAFF'


def _apply_font(font_family):
	"""Apply font_family AFTER plt.style.use('default')."""
	if font_family is None:
		return
	import matplotlib
	matplotlib.rcParams['font.family'] = font_family
	matplotlib.rcParams['mathtext.fontset'] = 'custom'
	matplotlib.rcParams['mathtext.it'] = f'{font_family}:italic'
	matplotlib.rcParams['mathtext.rm'] = font_family


def ensure_plot_dirs(output_dir):
	plot_dir = Path(output_dir) / "plots"
	png_dir = plot_dir / "png"
	pdf_dir = plot_dir / "pdf"
	png_dir.mkdir(parents=True, exist_ok=True)
	pdf_dir.mkdir(parents=True, exist_ok=True)
	return png_dir, pdf_dir


# ---------------------------------------------------------------------------
# DSE overlap Venn diagrams (Target vs Reference, by dPSI direction)
# ---------------------------------------------------------------------------

def draw_dse_direction_venn(ax, counts, event_type, direction,
							target_color=None, reference_color=None):
	target_color = target_color or TARGET_COLOR
	reference_color = reference_color or REFERENCE_COLOR
	if not VENN_AVAILABLE:
		ax.text(0.5, 0.5, 'Venn diagrams not available\n(matplotlib-venn not installed)',
				ha='center', va='center', transform=ax.transAxes, fontsize=12)
		ax.set_title(event_type, fontsize=16, color='#2E2E2E', weight='normal', pad=0.5)
		return

	target_only = counts['target_only']
	reference_only = counts['reference_only']
	intersection = counts['intersection_count']

	if target_only == 0 and reference_only == 0 and intersection == 0:
		ax.text(0.5, 0.5, 'No DSEs detected',
				ha='center', va='center', transform=ax.transAxes, fontsize=12, color='#888888')
		direction_label = 'Up-regulated' if direction == 'up' else 'Down-regulated'
		ax.set_title(f"{event_type} — {direction_label} DSEs",
					 fontsize=16, color='#2E2E2E', weight='normal', pad=0.5)
		ax.set_facecolor('white')
		ax.grid(False)
		return

	venn = venn2(subsets=(target_only, reference_only, intersection),
				 set_labels=None, ax=ax)
	if venn.get_patch_by_id('10'):
		venn.get_patch_by_id('10').set_facecolor(target_color)
		venn.get_patch_by_id('10').set_alpha(0.5)
		venn.get_patch_by_id('10').set_edgecolor(target_color)
		venn.get_patch_by_id('10').set_linewidth(1)
	if venn.get_patch_by_id('01'):
		venn.get_patch_by_id('01').set_facecolor(reference_color)
		venn.get_patch_by_id('01').set_alpha(0.5)
		venn.get_patch_by_id('01').set_edgecolor(reference_color)
		venn.get_patch_by_id('01').set_linewidth(1)
	if venn.get_patch_by_id('11'):
		venn.get_patch_by_id('11').set_facecolor('#B49A5D')
		venn.get_patch_by_id('11').set_alpha(0.6)
		venn.get_patch_by_id('11').set_edgecolor('none')
		venn.get_patch_by_id('11').set_linewidth(1)

	for text in venn.subset_labels:
		if text:
			text.set_fontsize(14)
			text.set_color('black')
			text.set_weight('normal')

	if any([target_only, reference_only, intersection]):
		ax.text(-0.8, -0.6, 'Target DSEs', ha='center', va='center',
				fontsize=14, color=target_color, weight='normal')
		ax.text(0.8, -0.6, 'Reference DSEs', ha='center', va='center',
				fontsize=14, color=reference_color, weight='normal')

	direction_label = 'Up-regulated' if direction == 'up' else 'Down-regulated'
	ax.set_title(f"{event_type} — {direction_label} DSEs",
				 fontsize=16, color='#2E2E2E', weight='normal', pad=0.5)
	ax.set_facecolor('white')
	ax.grid(False)


def draw_single_dse_direction_venn(event_type, direction, counts, png_path, pdf_path,
								   target_color=None, reference_color=None,
								   font_family=None):
	plt.style.use('default')
	_apply_font(font_family)
	fig, ax = plt.subplots(figsize=(6, 5))
	fig.patch.set_facecolor('white')
	draw_dse_direction_venn(ax, counts, event_type, direction,
							target_color=target_color, reference_color=reference_color)
	plt.tight_layout(pad=0.5)
	plt.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
	plt.savefig(pdf_path, bbox_inches='tight', facecolor='white', edgecolor='none')
	plt.close()
	logger.info(f"Saved DSE Venn diagram for {event_type} ({direction}): {png_path}")


def draw_dse_direction_venn_grid(overlap_data, direction, png_path, pdf_path,
								 target_color=None, reference_color=None,
								 font_family=None):
	direction_label = 'Up-regulated' if direction == 'up' else 'Down-regulated'
	valid_events = []
	for event_type, direction_data in overlap_data.items():
		counts = direction_data.get(direction, {})
		if not counts:
			continue
		display_name = "All Events" if event_type == 'all' else event_type
		valid_events.append((display_name, counts))
	if not valid_events:
		return

	n_events = len(valid_events)
	cols = min(3, n_events)
	rows = math.ceil(n_events / cols)

	plt.style.use('default')
	_apply_font(font_family)
	fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4.5))
	fig.patch.set_facecolor('white')
	fig.suptitle(f"{direction_label} DSEs — Target vs Reference",
				 fontsize=18, color='#2E2E2E', weight='normal', y=1.02)

	if rows == 1 and cols == 1:
		axes = [axes]
	elif rows == 1 or cols == 1:
		axes = axes.flatten()
	else:
		axes = axes.flatten()

	for i, (display_name, counts) in enumerate(valid_events):
		draw_dse_direction_venn(axes[i], counts, display_name, direction,
								target_color=target_color, reference_color=reference_color)
	for i in range(n_events, len(axes)):
		axes[i].axis('off')
		axes[i].set_facecolor('white')

	plt.tight_layout(pad=2.0, h_pad=3.0, w_pad=2.0)
	plt.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
	plt.savefig(pdf_path, bbox_inches='tight', facecolor='white', edgecolor='none')
	plt.close()
	logger.info(f"Saved DSE direction Venn grid ({direction}): {png_path}")


def plot_dse_overlap_venns(overlap_data, output_dir,
						   target_color=None, reference_color=None,
						   font_family=None):
	png_dir, pdf_dir = ensure_plot_dirs(output_dir)
	for event_type, direction_data in overlap_data.items():
		display_name = "All Events" if event_type == 'all' else event_type
		file_base = "dse_venn_all" if event_type == 'all' else f"dse_venn_{event_type}"
		for direction in ('up', 'down'):
			counts = direction_data.get(direction, {})
			if not counts:
				continue
			png_path = png_dir / f"{file_base}_{direction}.png"
			pdf_path = pdf_dir / f"{file_base}_{direction}.pdf"
			draw_single_dse_direction_venn(display_name, direction, counts, png_path, pdf_path,
										   target_color=target_color, reference_color=reference_color,
										   font_family=font_family)

	for direction in ('up', 'down'):
		grid_png = png_dir / f"dse_venn_grid_{direction}.png"
		grid_pdf = pdf_dir / f"dse_venn_grid_{direction}.pdf"
		draw_dse_direction_venn_grid(overlap_data, direction, grid_png, grid_pdf,
									 target_color=target_color, reference_color=reference_color,
									 font_family=font_family)
	logger.info("DSE overlap Venn diagram creation completed")


# ---------------------------------------------------------------------------
# Weighted correlation scatter plots
# ---------------------------------------------------------------------------

def plot_weighted_correlation_scatter(scatter_data_dict, correlation_results, output_dir,
									  font_family=None):
	png_dir, pdf_dir = ensure_plot_dirs(output_dir)
	for event_type, corr_df in scatter_data_dict.items():
		if corr_df.empty:
			continue
		stats_d = correlation_results.get(event_type, {})
		n_used = stats_d.get('n_used', 0)
		r = stats_d.get('r', np.nan)
		p_value = stats_d.get('p_value_approx', np.nan)
		if n_used < 2:
			continue

		x = corr_df['reference_dpsi'].values * 100
		y = corr_df['target_dpsi'].values * 100
		try:
			xy = np.vstack([x, y])
			kde = gaussian_kde(xy)
			z = kde(xy)
		except np.linalg.LinAlgError:
			logger.warning(f"KDE failed for {event_type} (singular covariance); using uniform colour")
			z = np.ones_like(x)

		plt.style.use('default')
		_apply_font(font_family)
		fig, ax = plt.subplots(figsize=(5, 4))
		fig.patch.set_facecolor('white')

		sc = ax.scatter(x, y, c=z, cmap='viridis', s=10, edgecolor='none', alpha=0.7)
		sns.regplot(x=x, y=y, scatter=False, color='grey',
					line_kws={'linewidth': 2, 'alpha': 0.7}, ax=ax)

		cbar = plt.colorbar(sc, ax=ax, shrink=0.6, aspect=20)
		cbar.set_label('Kernel density')

		ax.set_xlim(-100, 100)
		ax.set_ylim(-100, 100)
		ax.set_xlabel('Reference dPSI (%)', fontsize=14, color='#2E2E2E')
		ax.set_ylabel('Target dPSI (%)', fontsize=14, color='#2E2E2E')
		ax.axline((0, 0), slope=1, color='lightgray', linestyle='--', alpha=0.5, zorder=0)

		ax.set_title('Weighted dPSI correlation', fontsize=14, color='#2E2E2E',
					  weight='normal', pad=24)
		r_text = f"{r:.3f}" if np.isfinite(r) else "NA"
		if np.isfinite(p_value):
			exp = int(np.floor(np.log10(abs(p_value)))) if p_value != 0 else 0
			mantissa = p_value / (10 ** exp)
			p_text = f"{p_value:.3f}" if exp == 0 else f"{mantissa:.2f}×10$^{{{exp}}}$"
		else:
			p_text = "NA"
		subtitle = f"$n$ = {n_used}, $r$ = {r_text}, $P$ = {p_text}"
		ax.text(0.5, 1.01, subtitle, transform=ax.transAxes, ha='center', va='bottom',
				fontsize=11, color='#2E2E2E')
		ax.grid(True, alpha=0.5)
		ax.set_facecolor('white')

		png_path = png_dir / f"weighted_correlation_scatter_{event_type}.png"
		pdf_path = pdf_dir / f"weighted_correlation_scatter_{event_type}.pdf"
		plt.savefig(png_path, dpi=800, bbox_inches='tight', facecolor='white', edgecolor='none')
		plt.savefig(pdf_path, bbox_inches='tight', facecolor='white', edgecolor='none')
		plt.close()
		logger.info(f"Saved weighted correlation scatter for {event_type}: {png_path}")
	logger.info("Weighted correlation scatter plot creation completed")
