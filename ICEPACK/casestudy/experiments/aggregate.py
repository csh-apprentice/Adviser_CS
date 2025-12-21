"""
Aggregate multiple run directories into a CSV and generate simple plots.

Assumes each run directory contains:
  - stats.json
  - metrics.json
  - config.json

Example:
  python -m experiments.aggregate --runs-dir runs --out summary
"""
from __future__ import annotations

import argparse
from pathlib import Path
import json

import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-dir", required=True, help="Directory containing run subfolders")
    ap.add_argument("--out", required=True, help="Output directory for summary CSV/plots")
    ap.add_argument("--x", default="fluidity_scale", help="X-axis field for plots (from config)")
    return ap.parse_args()


def load_json(p: Path):
    return json.loads(p.read_text())


def main():
    args = parse_args()
    runs_dir = Path(args.runs_dir)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for run in sorted([p for p in runs_dir.iterdir() if p.is_dir()]):
        cfg_p = run / "config.json"
        stats_p = run / "stats.json"
        met_p = run / "metrics.json"
        if not (cfg_p.exists() and stats_p.exists() and met_p.exists()):
            continue
        cfg = load_json(cfg_p)
        stats = load_json(stats_p)
        metrics = load_json(met_p)
        row = {}
        row.update({f"cfg_{k}": v for k, v in cfg.items() if isinstance(v, (int, float, str))})
        row.update({f"stats_{k}": v for k, v in stats.items() if isinstance(v, (int, float, str))})
        row.update({f"met_{k}": v for k, v in metrics.items() if isinstance(v, (int, float, str))})
        row["run"] = run.name
        rows.append(row)

    if not rows:
        raise SystemExit(f"No valid runs found in {runs_dir}")

    df = pd.DataFrame(rows)
    csv_path = out_dir / "summary.csv"
    df.to_csv(csv_path, index=False)
    print(f"[OK] Wrote: {csv_path}")

    # Simple plots: runtime and flux vs chosen x
    xcol = f"cfg_{args.x}"
    if xcol not in df.columns:
        print(f"[WARN] Column {xcol} not found; skipping plots.")
        return

    # Wall time plot
    if "stats_wall_s" in df.columns:
        plt.figure()
        plt.plot(df[xcol], df["stats_wall_s"], marker="o")
        plt.xlabel(args.x)
        plt.ylabel("wall time (s)")
        plt.title("Runtime vs parameter")
        plt.grid(True)
        plt.savefig(out_dir / "runtime_vs_param.png", dpi=200, bbox_inches="tight")
        plt.close()

    # Flux plot
    if "met_flux_gate_m3_per_yr" in df.columns:
        plt.figure()
        plt.plot(df[xcol], df["met_flux_gate_m3_per_yr"], marker="o")
        plt.xlabel(args.x)
        plt.ylabel("flux across gate (m^3/yr)")
        plt.title("Flux vs parameter")
        plt.grid(True)
        plt.savefig(out_dir / "flux_vs_param.png", dpi=200, bbox_inches="tight")
        plt.close()

    # Max speed plot
    if "met_max_speed_m_per_yr" in df.columns:
        plt.figure()
        plt.plot(df[xcol], df["met_max_speed_m_per_yr"], marker="o")
        plt.xlabel(args.x)
        plt.ylabel("max speed (m/yr)")
        plt.title("Max speed vs parameter")
        plt.grid(True)
        plt.savefig(out_dir / "maxspeed_vs_param.png", dpi=200, bbox_inches="tight")
        plt.close()

    print(f"[OK] Wrote plots to: {out_dir}")


if __name__ == "__main__":
    main()
