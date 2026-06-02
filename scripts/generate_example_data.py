"""Generate small synthetic Shiba PSI tables under example/{target,reference}/results/splicing/.

The schema mirrors the real Shiba output (see
https://sika-zheng-lab.github.io/Shiba/output/shiba/) including:
- column order (strand/gene_id/gene_name/label after intron coordinates)
- junction-count columns as SCALAR integers (already summed across replicates)
- pos_id of the form ``EVENT@chr@coord1@coord2(...)``
- per-sample ``<sample>_PSI`` columns and a final ``p_ttest`` column

Run once when the example/ data needs to be regenerated. The output files are
committed to the repository so end users don't need to run this script.
"""

import os
from typing import Dict, List

import numpy as np
import pandas as pd

EVENT_TYPES = ("SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE")
N_REPS = 3
N_EVENTS = {"SE": 40, "FIVE": 15, "THREE": 15, "MXE": 10, "RI": 12, "MSE": 8, "AFE": 10, "ALE": 10}
SAMPLE_NAMES = [f"ref{i + 1}" for i in range(N_REPS)] + [f"alt{i + 1}" for i in range(N_REPS)]


def _coords(chrom_idx: int, start: int, length: int) -> str:
	return f"chr{chrom_idx}:{start}-{start + length}"


def _coords_no_chr(start: int, length: int) -> str:
	return f"{start}-{start + length}"


def _scalar_count(rng, mean_cov: float) -> int:
	return int(max(rng.poisson(lam=max(mean_cov, 0.1) * N_REPS), 0))


def _semicolon_counts(rng, mean_covs):
	return ";".join(str(int(max(rng.poisson(lam=max(c, 0.1) * N_REPS), 0))) for c in mean_covs)


def _semicolon_floats(rng, n: int, mean: float, std: float) -> str:
	return ";".join(f"{float(abs(rng.normal(mean, std)) + 0.01):.6g}" for _ in range(n))


def _sample_psi(rng, group_psi: float):
	"""Per-replicate PSI jittered around the group mean."""
	out = []
	for _ in range(N_REPS):
		out.append(float(np.clip(group_psi + rng.normal(0, 0.05), 0.01, 0.99)))
	return out


def _common(event_type, i, rng, dpsi):
	ref_psi = float(np.clip(0.5 - dpsi / 2, 0.02, 0.98))
	alt_psi = float(np.clip(ref_psi + dpsi, 0.02, 0.98))
	dpsi = alt_psi - ref_psi
	diff = "Yes" if abs(dpsi) >= 0.1 else "No"
	q = float(min(1.0, abs(rng.normal(0.04, 0.05)) + (0.1 if diff == "No" else 0.0)))
	p_ttest = float(min(1.0, abs(rng.normal(0.04, 0.05)) + (0.1 if diff == "No" else 0.0)))
	strand = "+" if rng.random() < 0.5 else "-"
	label = "annotated" if rng.random() < 0.8 else "unannotated"
	ref_samples = _sample_psi(rng, ref_psi)
	alt_samples = _sample_psi(rng, alt_psi)
	return {
		"event_id": f"{event_type}_{i:04d}",
		"strand": strand,
		"gene_id": f"GENE{i:04d}",
		"gene_name": f"Gene{i:04d}",
		"label": label,
		"ref_PSI": ref_psi,
		"alt_PSI": alt_psi,
		"dPSI": dpsi,
		"q": q,
		"Diff events": diff,
		"p_ttest": p_ttest,
		"_ref_samples": ref_samples,
		"_alt_samples": alt_samples,
	}


def _trailing(common):
	"""p_maximum/q/Diff events placeholder is filled per-row; here only sample PSIs + p_ttest."""
	cols = {}
	for name, val in zip(SAMPLE_NAMES[: N_REPS], common["_ref_samples"]):
		cols[f"{name}_PSI"] = val
	for name, val in zip(SAMPLE_NAMES[N_REPS:], common["_alt_samples"]):
		cols[f"{name}_PSI"] = val
	cols["p_ttest"] = common["p_ttest"]
	return cols


def _ids_strand_gene(common):
	return {
		"event_id": common["event_id"],
		"strand": common["strand"],
		"gene_id": common["gene_id"],
		"gene_name": common["gene_name"],
		"label": common["label"],
	}


def _row_SE(rng, i, common):
	chrom = (i % 22) + 1
	base = 100000 + i * 5000
	mean_cov = float(rng.integers(20, 120))
	exon = _coords(chrom, base + 200, 100)
	intron_a = _coords(chrom, base, 200)
	intron_b = _coords(chrom, base + 300, 200)
	intron_c = _coords(chrom, base, 500)
	pos_id = (f"SE@chr{chrom}@{base + 200}-{base + 300}"
			  f"@{base}-{base + 500}")
	row = {
		"event_id": common["event_id"],
		"pos_id": pos_id,
		"exon": exon,
		"intron_a": intron_a,
		"intron_b": intron_b,
		"intron_c": intron_c,
		"strand": common["strand"],
		"gene_id": common["gene_id"],
		"gene_name": common["gene_name"],
		"label": common["label"],
		"ref_junction_a": _scalar_count(rng, mean_cov * common["ref_PSI"]),
		"ref_junction_b": _scalar_count(rng, mean_cov * common["ref_PSI"]),
		"ref_junction_c": _scalar_count(rng, mean_cov * (1 - common["ref_PSI"])),
		"ref_PSI": common["ref_PSI"],
		"alt_junction_a": _scalar_count(rng, mean_cov * common["alt_PSI"]),
		"alt_junction_b": _scalar_count(rng, mean_cov * common["alt_PSI"]),
		"alt_junction_c": _scalar_count(rng, mean_cov * (1 - common["alt_PSI"])),
		"alt_PSI": common["alt_PSI"],
		"dPSI": common["dPSI"],
		"OR_junction_a": float(abs(rng.normal(1, 0.5)) + 0.01),
		"p_junction_a": float(min(1.0, abs(rng.normal(0.05, 0.05)))),
		"OR_junction_b": float(abs(rng.normal(1, 0.5)) + 0.01),
		"p_junction_b": float(min(1.0, abs(rng.normal(0.05, 0.05)))),
		"p_maximum": float(min(1.0, abs(rng.normal(0.06, 0.05)))),
		"q": common["q"],
		"Diff events": common["Diff events"],
	}
	row.update(_trailing(common))
	return row


def _row_FIVE_THREE(rng, i, common, event_type):
	chrom = (i % 22) + 1
	base = 100000 + i * 5000
	mean_cov = float(rng.integers(20, 120))
	exon_a = _coords(chrom, base, 200)
	exon_b = _coords(chrom, base + 50, 150)
	intron_a = _coords(chrom, base + 200, 300)
	intron_b = _coords(chrom, base + 200, 250)
	pos_id = (f"{event_type}@chr{chrom}@{base}-{base + 200}"
			  f"@{base + 50}-{base + 200}")
	row = {
		"event_id": common["event_id"],
		"pos_id": pos_id,
		"exon_a": exon_a,
		"exon_b": exon_b,
		"intron_a": intron_a,
		"intron_b": intron_b,
		"strand": common["strand"],
		"gene_id": common["gene_id"],
		"gene_name": common["gene_name"],
		"label": common["label"],
		"ref_junction_a": _scalar_count(rng, mean_cov * common["ref_PSI"]),
		"ref_junction_b": _scalar_count(rng, mean_cov * (1 - common["ref_PSI"])),
		"ref_PSI": common["ref_PSI"],
		"alt_junction_a": _scalar_count(rng, mean_cov * common["alt_PSI"]),
		"alt_junction_b": _scalar_count(rng, mean_cov * (1 - common["alt_PSI"])),
		"alt_PSI": common["alt_PSI"],
		"dPSI": common["dPSI"],
		"OR": float(abs(rng.normal(1, 0.5)) + 0.01),
		"p": float(min(1.0, abs(rng.normal(0.05, 0.05)))),
		"q": common["q"],
		"Diff events": common["Diff events"],
	}
	row.update(_trailing(common))
	return row


def _row_MXE(rng, i, common):
	chrom = (i % 22) + 1
	base = 100000 + i * 5000
	mean_cov = float(rng.integers(20, 120))
	exon_a = _coords(chrom, base + 200, 80)
	exon_b = _coords(chrom, base + 400, 80)
	intron_a1 = _coords(chrom, base, 200)
	intron_a2 = _coords(chrom, base + 280, 200)
	intron_b1 = _coords(chrom, base, 400)
	intron_b2 = _coords(chrom, base + 480, 100)
	pos_id = (f"MXE@chr{chrom}@{base}@{base + 200}-{base + 280}"
			  f"@{base + 400}-{base + 480}@{base + 580}")
	row = {
		"event_id": common["event_id"],
		"pos_id": pos_id,
		"exon_a": exon_a,
		"exon_b": exon_b,
		"intron_a1": intron_a1,
		"intron_a2": intron_a2,
		"intron_b1": intron_b1,
		"intron_b2": intron_b2,
		"strand": common["strand"],
		"gene_id": common["gene_id"],
		"gene_name": common["gene_name"],
		"label": common["label"],
		"ref_junction_a1": _scalar_count(rng, mean_cov * common["ref_PSI"]),
		"ref_junction_a2": _scalar_count(rng, mean_cov * common["ref_PSI"]),
		"ref_junction_b1": _scalar_count(rng, mean_cov * (1 - common["ref_PSI"])),
		"ref_junction_b2": _scalar_count(rng, mean_cov * (1 - common["ref_PSI"])),
		"ref_PSI": common["ref_PSI"],
		"alt_junction_a1": _scalar_count(rng, mean_cov * common["alt_PSI"]),
		"alt_junction_a2": _scalar_count(rng, mean_cov * common["alt_PSI"]),
		"alt_junction_b1": _scalar_count(rng, mean_cov * (1 - common["alt_PSI"])),
		"alt_junction_b2": _scalar_count(rng, mean_cov * (1 - common["alt_PSI"])),
		"alt_PSI": common["alt_PSI"],
		"dPSI": common["dPSI"],
		"OR_junction_a1b1": float(abs(rng.normal(1, 0.5)) + 0.01),
		"p_junction_a1b1": float(min(1.0, abs(rng.normal(0.05, 0.05)))),
		"OR_junction_a1b2": float(abs(rng.normal(1, 0.5)) + 0.01),
		"p_junction_a1b2": float(min(1.0, abs(rng.normal(0.05, 0.05)))),
		"OR_junction_a2b1": float(abs(rng.normal(1, 0.5)) + 0.01),
		"p_junction_a2b1": float(min(1.0, abs(rng.normal(0.05, 0.05)))),
		"OR_junction_a2b2": float(abs(rng.normal(1, 0.5)) + 0.01),
		"p_junction_a2b2": float(min(1.0, abs(rng.normal(0.05, 0.05)))),
		"p_maximum": float(min(1.0, abs(rng.normal(0.06, 0.05)))),
		"q": common["q"],
		"Diff events": common["Diff events"],
	}
	row.update(_trailing(common))
	return row


def _row_RI(rng, i, common):
	chrom = (i % 22) + 1
	base = 100000 + i * 5000
	mean_cov = float(rng.integers(20, 120))
	exon_a = _coords(chrom, base, 100)
	exon_b = _coords(chrom, base + 300, 100)
	exon_c = _coords(chrom, base, 400)
	intron_a = _coords(chrom, base + 100, 200)
	pos_id = f"RI@chr{chrom}@{base + 100}-{base + 300}"
	row = {
		"event_id": common["event_id"],
		"pos_id": pos_id,
		"exon_a": exon_a,
		"exon_b": exon_b,
		"exon_c": exon_c,
		"intron_a": intron_a,
		"strand": common["strand"],
		"gene_id": common["gene_id"],
		"gene_name": common["gene_name"],
		"label": common["label"],
		"ref_junction_a": _scalar_count(rng, mean_cov * common["ref_PSI"]),
		"ref_junction_a_start": _scalar_count(rng, mean_cov * (1 - common["ref_PSI"])),
		"ref_junction_a_end": _scalar_count(rng, mean_cov * (1 - common["ref_PSI"])),
		"ref_PSI": common["ref_PSI"],
		"alt_junction_a": _scalar_count(rng, mean_cov * common["alt_PSI"]),
		"alt_junction_a_start": _scalar_count(rng, mean_cov * (1 - common["alt_PSI"])),
		"alt_junction_a_end": _scalar_count(rng, mean_cov * (1 - common["alt_PSI"])),
		"alt_PSI": common["alt_PSI"],
		"dPSI": common["dPSI"],
		"OR_junction_a_start": float(abs(rng.normal(1, 0.5)) + 0.01),
		"p_junction_a_start": float(min(1.0, abs(rng.normal(0.05, 0.05)))),
		"OR_junction_a_end": float(abs(rng.normal(1, 0.5)) + 0.01),
		"p_junction_a_end": float(min(1.0, abs(rng.normal(0.05, 0.05)))),
		"p_maximum": float(min(1.0, abs(rng.normal(0.06, 0.05)))),
		"q": common["q"],
		"Diff events": common["Diff events"],
	}
	row.update(_trailing(common))
	return row


def _row_MSE(rng, i, common):
	chrom = (i % 22) + 1
	base = 100000 + i * 5000
	mse_n = int(rng.integers(2, 4))
	mean_cov = float(rng.integers(20, 120))
	exon_starts = [base + 200 + k * 300 for k in range(mse_n)]
	exon = ";".join(_coords(chrom, s, 100) for s in exon_starts)
	# inclusive introns: mse_n + 1, plus exclusive intron last entry
	intron_starts = [base + k * 300 for k in range(mse_n + 1)]
	intron = ";".join(_coords(chrom, s, 200) for s in intron_starts)
	intron = intron + ";" + _coords(chrom, base, 500 + 300 * mse_n)
	n_inclusive = mse_n + 1
	ref_means = [mean_cov * common["ref_PSI"]] * n_inclusive + [mean_cov * (1 - common["ref_PSI"])]
	alt_means = [mean_cov * common["alt_PSI"]] * n_inclusive + [mean_cov * (1 - common["alt_PSI"])]
	pos_id = f"MSE@chr{chrom}@{exon.replace(f'chr{chrom}:', '')}@{base}-{base + 500 + 300 * mse_n}"
	row = {
		"event_id": common["event_id"],
		"pos_id": pos_id,
		"mse_n": mse_n,
		"exon": exon,
		"intron": intron,
		"strand": common["strand"],
		"gene_id": common["gene_id"],
		"gene_name": common["gene_name"],
		"label": common["label"],
		"ref_junction": _semicolon_counts(rng, ref_means),
		"ref_PSI": common["ref_PSI"],
		"alt_junction": _semicolon_counts(rng, alt_means),
		"alt_PSI": common["alt_PSI"],
		"dPSI": common["dPSI"],
		"OR_junction": _semicolon_floats(rng, n_inclusive, 1.0, 0.5),
		"p_junction": _semicolon_floats(rng, n_inclusive, 0.05, 0.05),
		"p_maximum": float(min(1.0, abs(rng.normal(0.06, 0.05)))),
		"q": common["q"],
		"Diff events": common["Diff events"],
	}
	row.update(_trailing(common))
	return row


def _row_AFE_ALE(rng, i, common, event_type):
	chrom = (i % 22) + 1
	base = 100000 + i * 5000
	mean_cov = float(rng.integers(20, 120))
	exon_a = _coords(chrom, base + 200, 150)
	exon_b = _coords(chrom, base + 600, 150)
	intron_a = _coords(chrom, base, 300)
	intron_b = _coords(chrom, base + 400, 300)
	pos_id = (f"{event_type}@chr{chrom}@{base + 200}-{base + 350}"
			  f"@{base + 600}-{base + 750}")
	row = {
		"event_id": common["event_id"],
		"pos_id": pos_id,
		"exon_a": exon_a,
		"exon_b": exon_b,
		"intron_a": intron_a,
		"intron_b": intron_b,
		"strand": common["strand"],
		"gene_id": common["gene_id"],
		"gene_name": common["gene_name"],
		"label": common["label"],
		"ref_junction_a": _scalar_count(rng, mean_cov * common["ref_PSI"]),
		"ref_junction_b": _scalar_count(rng, mean_cov * (1 - common["ref_PSI"])),
		"ref_PSI": common["ref_PSI"],
		"alt_junction_a": _scalar_count(rng, mean_cov * common["alt_PSI"]),
		"alt_junction_b": _scalar_count(rng, mean_cov * (1 - common["alt_PSI"])),
		"alt_PSI": common["alt_PSI"],
		"dPSI": common["dPSI"],
		"OR_junction": float(abs(rng.normal(1, 0.5)) + 0.01),
		"p_junction": float(min(1.0, abs(rng.normal(0.05, 0.05)))),
		"p_maximum": float(min(1.0, abs(rng.normal(0.06, 0.05)))),
		"q": common["q"],
		"Diff events": common["Diff events"],
	}
	row.update(_trailing(common))
	return row


def make_psi_table(event_type, n_events, seed, correlated_dpsi=None):
	rng = np.random.default_rng(seed)
	rows = []
	for i in range(n_events):
		if correlated_dpsi is not None:
			dpsi = float(np.clip(correlated_dpsi[i] + rng.normal(0, 0.05), -0.9, 0.9))
		else:
			dpsi = float(np.clip(rng.normal(0, 0.25), -0.9, 0.9))
		common = _common(event_type, i, rng, dpsi)
		if event_type == "SE":
			rows.append(_row_SE(rng, i, common))
		elif event_type in ("FIVE", "THREE"):
			rows.append(_row_FIVE_THREE(rng, i, common, event_type))
		elif event_type == "MXE":
			rows.append(_row_MXE(rng, i, common))
		elif event_type == "RI":
			rows.append(_row_RI(rng, i, common))
		elif event_type == "MSE":
			rows.append(_row_MSE(rng, i, common))
		elif event_type in ("AFE", "ALE"):
			rows.append(_row_AFE_ALE(rng, i, common, event_type))
	return pd.DataFrame(rows)


def main():
	repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	rng = np.random.default_rng(42)

	for condition, _seed_offset in [("target", 0), ("reference", 1000)]:
		out_dir = os.path.join(repo_root, "example", condition, "results", "splicing")
		os.makedirs(out_dir, exist_ok=True)
		for ev in EVENT_TYPES:
			n = N_EVENTS[ev]
			shared = rng.normal(0, 0.25, size=n)
			shared_dpsi = np.clip(shared + rng.normal(0, 0.05, size=n), -0.9, 0.9)
			df = make_psi_table(ev, n, seed=hash((ev, condition)) % (2**32),
								correlated_dpsi=shared_dpsi)
			df.to_csv(os.path.join(out_dir, f"PSI_{ev}.txt"), sep="\t", index=False)
			print(f"Wrote {condition}/{ev}: {len(df)} events")


if __name__ == "__main__":
	main()
