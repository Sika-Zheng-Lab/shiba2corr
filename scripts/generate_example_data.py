"""Generate small synthetic Shiba PSI tables under example/{target,reference}/results/splicing/.

Run once when the example/ data needs to be regenerated. The output files are
committed to the repository so end users don't need to run this script.
"""

import os
import numpy as np
import pandas as pd

EVENT_TYPES = ("SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE")
N_REPS = 3
N_EVENTS = {"SE": 40, "FIVE": 15, "THREE": 15, "MXE": 10, "RI": 12, "MSE": 8, "AFE": 10, "ALE": 10}


def _coverage_string(rng, n_reps, mean_cov):
    counts = rng.poisson(lam=mean_cov, size=n_reps)
    return ";".join(str(int(max(c, 0))) for c in counts)


def make_psi_table(event_type: str, n_events: int, seed: int, correlated_dpsi=None):
    """Create one PSI table. correlated_dpsi (target only) injects shared signal."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_events):
        pos_id = f"{event_type}_{i:04d}"
        mean_cov = float(rng.integers(20, 120))
        if correlated_dpsi is not None:
            dpsi = float(np.clip(correlated_dpsi[i] + rng.normal(0, 0.05), -0.9, 0.9))
        else:
            # base dPSI sampled around 0 with some real signal at the tails
            dpsi = float(np.clip(rng.normal(0, 0.25), -0.9, 0.9))
        ref_psi = float(np.clip(0.5 - dpsi / 2, 0.02, 0.98))
        alt_psi = float(np.clip(ref_psi + dpsi, 0.02, 0.98))
        diff = "Yes" if abs(dpsi) >= 0.1 else "No"
        rows.append({
            "pos_id": pos_id,
            "label": f"gene_{i:04d}",
            "Diff events": diff,
            "dPSI": dpsi,
            "ref_PSI": ref_psi,
            "alt_PSI": alt_psi,
            "ref_junction_1": _coverage_string(rng, N_REPS, mean_cov * (1 - ref_psi)),
            "alt_junction_1": _coverage_string(rng, N_REPS, mean_cov * alt_psi),
            "p_ttest": float(min(1.0, abs(rng.normal(0.04, 0.05)) + (0.1 if diff == "No" else 0.0))),
        })
    return pd.DataFrame(rows)


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    rng = np.random.default_rng(42)

    for condition, seed_offset in [("target", 0), ("reference", 1000)]:
        out_dir = os.path.join(repo_root, "example", condition, "results", "splicing")
        os.makedirs(out_dir, exist_ok=True)
        for ev in EVENT_TYPES:
            n = N_EVENTS[ev]
            # shared signal across target/reference -> positive correlation
            shared = rng.normal(0, 0.25, size=n)
            shared_dpsi = np.clip(shared + rng.normal(0, 0.05, size=n), -0.9, 0.9)
            df = make_psi_table(ev, n, seed=hash((ev, condition)) % (2**32),
                                correlated_dpsi=shared_dpsi)
            df.to_csv(os.path.join(out_dir, f"PSI_{ev}.txt"), sep="\t", index=False)
            print(f"Wrote {condition}/{ev}: {len(df)} events")


if __name__ == "__main__":
    main()
