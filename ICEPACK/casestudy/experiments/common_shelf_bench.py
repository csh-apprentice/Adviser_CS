"""
ADVISER Icepack case-study scaffold (Synthetic Ice Shelf)

Goal: a small, reproducible Firedrake/Icepack workload that is:
  - container-friendly (no external datasets)
  - param-sweep friendly (embarrassingly parallel)
  - produces paper-ready metrics (runtime, DOFs, flux gate, speed)

Designed to be run *inside* the firedrakeproject/firedrake-vanilla Docker image
with Icepack installed/available.

Based on the Icepack synthetic ice shelf tutorial notebook structure.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Any, Tuple
import time
import json

import gmsh
import firedrake
import icepack
import os


@dataclass
class RunStats:
    wall_s: float
    ndofs_u: int
    ndofs_h: int
    notes: Dict[str, Any]



def make_mesh(out_dir: Path, R: float, dx: float) -> firedrake.Mesh:
    """
    MPI-robust mesh generator:
      - Uses a filesystem lock so only one process runs gmsh + writes the .msh
      - Other processes wait until the final .msh exists, then read it
      - Avoids relying on COMM_WORLD.rank (important if MPI launcher is mismatched)
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    msh_path = out_dir / "ice-shelf.msh"
    tmp_path = out_dir / "ice-shelf.tmp.msh"   # must end with .msh
    lock_path = out_dir / "ice-shelf.msh.lock"

    # Fast path: mesh already exists (e.g., rerun)
    if msh_path.exists() and msh_path.stat().st_size > 0:
        return firedrake.Mesh(str(msh_path))

    # Try to acquire exclusive lock (atomic)
    have_lock = False
    try:
        fd = os.open(str(lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.close(fd)
        have_lock = True
    except FileExistsError:
        have_lock = False

    if have_lock:
        try:
            # Re-check in case someone created it between our checks (rare)
            if not (msh_path.exists() and msh_path.stat().st_size > 0):
                gmsh.initialize()
                gmsh.model.add("ice_shelf")
                try:
                    geo = gmsh.model.geo

                    x1 = geo.add_point(-R, 0, 0, dx)
                    x2 = geo.add_point(+R, 0, 0, dx)
                    center1 = geo.add_point(0, 0, 0, dx)
                    center2 = geo.add_point(0, -4 * R, 0, dx)

                    arcs = [
                        geo.add_circle_arc(x1, center1, x2),
                        geo.add_circle_arc(x2, center2, x1),
                    ]

                    loop = geo.add_curve_loop(arcs)
                    surf = geo.add_plane_surface([loop])

                    # IMPORTANT: synchronize before physical groups
                    geo.synchronize()

                    gmsh.model.add_physical_group(1, [arcs[0]], tag=1)
                    gmsh.model.add_physical_group(1, [arcs[1]], tag=2)
                    gmsh.model.add_physical_group(2, [surf], tag=1)

                    gmsh.model.mesh.generate(2)

                    gmsh.write(str(tmp_path))
                finally:
                    gmsh.finalize()

                # Atomic-ish publish step
                tmp_path.replace(msh_path)
        finally:
            # Release lock
            try:
                lock_path.unlink()
            except FileNotFoundError:
                pass
    else:
        # Wait for the writer to publish the mesh
        # (don’t hang forever; fail loudly if something is wrong)
        deadline = time.time() + 300  # 5 minutes
        while time.time() < deadline:
            if msh_path.exists() and msh_path.stat().st_size > 0:
                break
            time.sleep(0.1)
        else:
            raise RuntimeError(f"Timed out waiting for mesh file: {msh_path}")

    return firedrake.Mesh(str(msh_path))


def build_initial_fields(mesh: firedrake.Mesh, R: float) -> Tuple[firedrake.Function, firedrake.Function]:
    """
    Construct initial thickness h0 and a rough initial velocity guess u0,
    in the spirit of the tutorial.

    NOTE: This is intentionally simple and self-contained; you can swap in
    more realistic fields later (e.g., Pine Island) without changing the
    experiment harness.
    """
    Q = firedrake.FunctionSpace(mesh, "CG", 2)
    V = firedrake.VectorFunctionSpace(mesh, "CG", 2)

    x, y = firedrake.SpatialCoordinate(mesh)

    # A smooth thickness "cap" with a minimum thickness floor.
    hb = 100.0  # m
    h_peak = 600.0  # m
    h_expr = firedrake.max_value(
        firedrake.Constant(hb),
        hb + (h_peak - hb) * firedrake.exp(-((x / R) ** 2 + (y / (2 * R)) ** 2)),
    )

    # A simple inflow-like velocity guess (mostly in y-direction), tapered towards sides.
    u_in = 300.0  # m/yr (scale for demo)
    taper = firedrake.exp(-4.0 * (x / R) ** 2)
    u_expr = firedrake.as_vector((0.0, u_in * taper))

    h0 = firedrake.Function(Q, name="thickness").interpolate(h_expr)
    u0 = firedrake.Function(V, name="velocity").interpolate(u_expr)
    return h0, u0


def solve_shelf_forward(
    out_dir: Path,
    R: float,
    dx: float,
    temperature_K: float,
    fluidity_scale: float,
    num_timesteps: int,
    final_time: float,
    dirichlet_ids=(1,),
    verbose: bool = False,
) -> Tuple[firedrake.Function, firedrake.Function, RunStats]:
    """
    Run a synthetic ice-shelf forward simulation:
      - build mesh
      - set up IceShelf model + FlowSolver
      - diagnostic solve to get velocity
      - optional prognostic thickness evolution loop

    Returns final (u, h, stats).
    """
    out_dir = Path(out_dir)
    t0 = time.time()

    mesh = make_mesh(out_dir, R=R, dx=dx)
    h0, u0 = build_initial_fields(mesh, R=R)

    model = icepack.models.IceShelf()
    solver = icepack.solvers.FlowSolver(model, dirichlet_ids=list(dirichlet_ids))

    # IMPORTANT: icepack FlowSolver expects fluidity to be a Firedrake Constant or Function.
    # If we build it as a UFL product (rate_factor(Constant(T)) * scale), icepack will reject it.
    # So we evaluate the scalar rate factor in Python and wrap it in a Firedrake Constant.
    A_val = float(icepack.rate_factor(float(temperature_K))) * float(fluidity_scale)
    A = firedrake.Constant(A_val)

    # Diagnostic solve: compute velocity for current thickness.
    h = h0.copy(deepcopy=True)
    u = solver.diagnostic_solve(velocity=u0, thickness=h, fluidity=A)

    # Prognostic loop: evolve thickness (accumulation=0 by default) and re-solve.
    if num_timesteps > 0 and final_time > 0:
        dt = float(final_time) / int(num_timesteps)
        a = firedrake.Constant(0.0)
        for _ in range(int(num_timesteps)):
            h = solver.prognostic_solve(
                dt,
                thickness=h,
                velocity=u,
                accumulation=a,
                thickness_inflow=h0,
            )
            u = solver.diagnostic_solve(velocity=u, thickness=h, fluidity=A)

    t1 = time.time()

    V = u.function_space()
    Q = h.function_space()
    stats = RunStats(
        wall_s=t1 - t0,
        ndofs_u=V.dim(),
        ndofs_h=Q.dim(),
        notes={
            "dirichlet_ids": list(dirichlet_ids),
            "R": R,
            "dx": dx,
            "temperature_K": float(temperature_K),
            "fluidity_scale": float(fluidity_scale),
            "num_timesteps": int(num_timesteps),
            "final_time": float(final_time),
        },
    )
    if verbose:
        print("RunStats:", stats)
    return u, h, stats


def compute_metrics(u: firedrake.Function, h: firedrake.Function, gate_id: int = 1) -> Dict[str, float]:
    """
    Compute simple scalar metrics:
      - max speed (m/yr)
      - mean speed (m/yr)
      - flux through boundary gate_id: ∫ (u·n) h ds
    """
    mesh = u.function_space().mesh()
    n = firedrake.FacetNormal(mesh)

    # Speed is a UFL expression; interpolate into a CG space to access .dat
    Q = h.function_space()  # scalar CG space (same as thickness)
    speed_expr = firedrake.sqrt(firedrake.inner(u, u))
    speed_fn = firedrake.Function(Q, name="speed").interpolate(speed_expr)

    max_speed = float(speed_fn.dat.data_ro.max())

    # Mean speed (area average)
    area = firedrake.assemble(1.0 * firedrake.dx(domain=mesh))
    mean_speed = float(firedrake.assemble(speed_fn * firedrake.dx(domain=mesh)) / area)

    # Flux across boundary id (units per 2D idealization)
    ds = firedrake.ds(gate_id, domain=mesh)
    flux = float(firedrake.assemble(firedrake.inner(u, n) * h * ds))

    return {
        "max_speed_m_per_yr": max_speed,
        "mean_speed_m_per_yr": mean_speed,
        "flux_gate_m3_per_yr": flux,
        "gate_id": float(gate_id),
    }



def save_run(out_dir: Path, cfg: Dict[str, Any], u: firedrake.Function, h: firedrake.Function, stats: RunStats, metrics: Dict[str, Any], save_fields: bool = False):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    (out_dir / "config.json").write_text(json.dumps(cfg, indent=2))
    (out_dir / "stats.json").write_text(json.dumps(asdict(stats), indent=2))
    (out_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))

    # Save fields for visualization (ParaView-friendly)
    # You can open these in ParaView or postprocess to make paper figs.
    if save_fields:
        try:
            firedrake.File(str(out_dir / "fields.pvd")).write(u, h)
        except Exception as e:
            # Some environments may not have VTK writers; keep going.
            (out_dir / "write_warning.txt").write_text(f"VTK write failed: {e}\n")
