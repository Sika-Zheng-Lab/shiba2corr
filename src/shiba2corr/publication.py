"""Reusable publication styling and sensitivity figure export."""

import logging

import matplotlib as mpl
from matplotlib import font_manager
import matplotlib.pyplot as plt
import numpy as np

from .plots import ensure_plot_dirs
from .sensitivity import METHODS

logger = logging.getLogger(__name__)


def publication_style():
	"""Require real Arial regular/bold fonts; never silently substitute fonts."""
	for weight in ('normal', 'bold'):
		try:
			path = font_manager.findfont(
				font_manager.FontProperties(family='Arial', weight=weight),
				fallback_to_default=False)
		except ValueError as exc:
			raise RuntimeError('Sensitivity figures require installed Arial fonts (regular and bold).') from exc
		font = font_manager.get_font(path)
		if font.family_name != 'Arial' or (weight == 'bold' and 'bold' not in font.style_name.lower()):
			raise RuntimeError(f'Expected genuine Arial {weight}; found {path}.')
		logger.info('Publication font (%s): %s', weight, path)
	return {
		'font.family': 'Arial', 'font.size': 7, 'axes.labelsize': 8,
		'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
		'axes.linewidth': 0.7, 'axes.spines.top': False, 'axes.spines.right': False,
		'axes.grid': False, 'xtick.direction': 'out', 'ytick.direction': 'out',
		'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
		'xtick.major.size': 3, 'ytick.major.size': 3,
		'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
		'text.usetex': False, 'figure.facecolor': 'white', 'axes.facecolor': 'white',
	}


def plot_sensitivity(raw, summary, output_dir):
	"""One 180-mm-wide three-panel figure per event type; no hidden outliers."""
	png_dir, pdf_dir = ensure_plot_dirs(output_dir)
	colors = ('#606060', '#477C9B', '#B87950')
	markers = ('o', 's', '^')
	with mpl.rc_context(publication_style()):
		for event_type, data in raw.groupby('event_type', sort=False):
			fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.75), sharey=True)
			fig.subplots_adjust(left=0.09, right=0.99, bottom=0.31, top=0.83, wspace=0.22)
			for i, (method, ax, color, marker) in enumerate(zip(METHODS.values(), axes, colors, markers)):
				stats = summary[(summary.event_type == event_type) & (summary.method == method)]
				stats = stats[stats.fraction < 1].sort_values('fraction')
				x = np.arange(len(stats))
				for j, row in enumerate(stats.itertuples()):
					values = data[(data.method == method) & (data.fraction == row.fraction)].correlation.dropna()
					# Deterministic horizontal jitter only; every defined replicate shown.
					jitter = np.random.default_rng(0).uniform(-0.16, 0.16, len(values))
					ax.scatter(j + jitter, values, s=3, color=color, alpha=0.22, linewidths=0)
					if row.n_valid:
						ax.errorbar(j, row.median, yerr=[[row.median - row.lower], [row.upper - row.median]],
							fmt=marker, color=color, markersize=3, elinewidth=0.9, capsize=2, zorder=3)
					if row.n_valid < row.n_repeats:
						ax.text(j, -0.97, f'{row.n_valid}/{row.n_repeats}', ha='center', fontsize=6.5)
				baseline = data[(data.method == method) & (data.fraction == 1)].correlation.iloc[0]
				if np.isfinite(baseline):
					ax.axhline(baseline, color='#555555', linewidth=0.8, linestyle='--', zorder=0)
				ax.set_xticks(x)
				ax.set_xticklabels([f'{r.fraction:.0%}\n(n={r.n_sample})' for r in stats.itertuples()])
				ax.set_xlim(-0.5, len(stats) - 0.5)
				ax.set_ylim(-1.05, 1.05)
				ax.set_yticks([-1, -0.5, 0, 0.5, 1])
				ax.set_xlabel('Events retained')
				ax.set_title(method, fontsize=8, pad=7)
				ax.text(-0.12, 1.09, chr(65 + i), transform=ax.transAxes, fontsize=9, fontweight='bold')
			axes[0].set_ylabel('Correlation coefficient')
			fig.text(0.09, 0.96, f'{event_type} | {data.n_total.iloc[0]} eligible events', fontsize=8)
			fig.text(0.09, 0.055, 'Dots: subsamples; markers: median; whiskers: 2.5–97.5 percentiles. Dashed: full event set.', fontsize=6.5)
			fig.text(0.09, 0.015, 'Counts inside panels, when present: defined / attempted correlations. Percentile ranges are not confidence intervals.', fontsize=6.5)
			try:
				fig.savefig(png_dir / f'event_selection_sensitivity_{event_type}.png', dpi=600, bbox_inches='tight')
				fig.savefig(pdf_dir / f'event_selection_sensitivity_{event_type}.pdf', bbox_inches='tight')
			finally:
				plt.close(fig)
