"""
Run a single synthetic ice-shelf experiment (forward simulation) and save artifacts.

Example:
  python -m experiments.run_forward --out runs/run_A1 --fluidity-scale 1.0 --dx 5000

This script is intended to be launched by ADVISER as an atomic task.

NOTE on PETSc options:
  Firedrake/PETSc will parse sys.argv for command-line options. If we leave our
  argparse flags (e.g., --dx) in sys.argv, PETSc will warn about "unused database options".
  So we:
    1) parse args early
    2) strip sys.argv down to argv[0]
    3) only then import Firedrake/Icepack code.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
import time


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="Output directory for this run")
    ap.add_argument("--save-fields", action="store_true", help="Write Firedrake field output (fields.pvd + fields/). Default off for benchmarking.")
    ap.add_argument("--R", type=float, default=200e3, help="Domain size parameter (meters)")
    ap.add_argument("--dx", type=float, default=5e3, help="Target mesh size (meters)")
    ap.add_argument("--temperature-K", type=float, default=255.15, help="Temperature in Kelvin")
    ap.add_argument("--fluidity-scale", type=float, default=1.0, help="Multiplier on A(T)")
    ap.add_argument("--num-timesteps", type=int, default=0, help="Thickness evolution timesteps (0 disables)")
    ap.add_argument("--final-time", type=float, default=0.0, help="Final time (years) for thickness evolution")
    ap.add_argument("--gate-id", type=int, default=1, help="Boundary id used as flux gate (default 1)")
    ap.add_argument("--dirichlet-ids", type=str, default="1", help="Comma-separated boundary ids for Dirichlet BCs")
    return ap.parse_args()


def main():
    args = parse_args()

    # Strip argv so PETSc doesn't see our argparse flags.
    sys.argv = sys.argv[:1]

    # Import after argv cleanup
    from .common_shelf_bench import solve_shelf_forward, compute_metrics, save_run

    out_dir = Path(args.out)
    dirichlet_ids = tuple(int(x) for x in args.dirichlet_ids.split(",") if x.strip())

    cfg = {
        "R": args.R,
        "dx": args.dx,
        "temperature_K": args.temperature_K,
        "fluidity_scale": args.fluidity_scale,
        "num_timesteps": args.num_timesteps,
        "final_time": args.final_time,
        "gate_id": args.gate_id,
        "dirichlet_ids": list(dirichlet_ids),
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }

    u, h, stats = solve_shelf_forward(
        out_dir=out_dir,
        R=args.R,
        dx=args.dx,
        temperature_K=args.temperature_K,
        fluidity_scale=args.fluidity_scale,
        num_timesteps=args.num_timesteps,
        final_time=args.final_time,
        dirichlet_ids=dirichlet_ids,
        verbose=True,
    )

    metrics = compute_metrics(u=u, h=h, gate_id=args.gate_id)
    save_run(out_dir=out_dir, cfg=cfg, u=u, h=h, stats=stats, metrics=metrics, save_fields=args.save_fields)

    print(f"[OK] Wrote run artifacts to: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
