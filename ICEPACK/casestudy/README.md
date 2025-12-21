# ADVISER + Icepack Case Study (Synthetic Ice Shelf)

This folder contains a small, headless (non-notebook) Icepack/Firedrake workload that is designed
for PEARC-style platform evaluation:

- embarrassingly-parallel parameter sweeps
- clean run artifacts per task (config/stats/metrics + optional .pvd fields)
- aggregation into paper-friendly tables and plots

## Run inside the Firedrake Docker image

From your Mac host:

```bash
docker run --platform=linux/amd64 -it --rm \
  -v "/Users/shihanc/Documents/ADVISER/Local":/workspace \
  -w /workspace \
  firedrakeproject/firedrake-vanilla:2025-01 \
  bash
```

Inside the container:

```bash
source ~/firedrake/bin/activate

# (If needed) install icepack into the venv, or ensure it's already available:
# pip install icepack

cd /workspace/adviser_icepack_case_study

# One run:
python -m experiments.run_forward --out runs/run_A1 --fluidity-scale 0.5 --dx 5000

# Another run:
python -m experiments.run_forward --out runs/run_A2 --fluidity-scale 1.0 --dx 5000

# Aggregate results:
python -m experiments.aggregate --runs-dir runs --out summary --x fluidity_scale
```

Outputs:
- `runs/<runname>/config.json`, `stats.json`, `metrics.json`, and `fields.pvd`
- `summary/summary.csv` + plots `.png`

## Notes for ADVISER integration

- Treat each `python -m experiments.run_forward ...` call as one atomic task.
- Use ADVISER to launch sweeps over `--fluidity-scale`, `--dx`, and optionally `--num-timesteps/--final-time`.
- Collect `summary.csv` and plots as the "paper artifact bundle".
