#!/usr/bin/env bash
set -euo pipefail

# -------------------------
# Usage check
# -------------------------
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <repeat_times> <dx>"
  echo "Example: $0 100 8000"
  exit 1
fi

REPEAT_TIMES="$1"
DX="$2"
TOTAL_RUNS=$((REPEAT_TIMES + 1))

echo "[Running] Case study benchmark"
echo "  repeat_times = ${REPEAT_TIMES}  (warm-up + ${REPEAT_TIMES} measured)"
echo "  dx           = ${DX}"
echo "  total runs   = ${TOTAL_RUNS}"

# -------------------------
# Helper: choose python
# -------------------------
if command -v python >/dev/null 2>&1; then
  PY=python
elif command -v python3 >/dev/null 2>&1; then
  PY=python3
else
  echo "[env] ERROR: neither python nor python3 found in PATH"
  echo "[env] PATH=$PATH"
  exit 127
fi

echo "[env] Using Python: $PY ($(command -v "$PY"))"
$PY --version || true
echo "[env] uname -m: $(uname -m)"

# -------------------------
# Helper: sudo (only if available)
# -------------------------
SUDO=""
if command -v sudo >/dev/null 2>&1; then
  SUDO="sudo"
fi

# -------------------------
# Environment inspection + optional activation
# -------------------------
echo "[env] Inspecting /home/firedrake (if it exists)..."
if [[ -d /home/firedrake ]]; then
  ls -lah /home/firedrake
else
  echo "[env] /home/firedrake does not exist"
fi

echo "[env] Detecting Firedrake venv activate script (optional)..."
ACTIVATED=0
for CANDIDATE in \
  /home/firedrake/firedrake/bin/activate \
  /home/firedrake/firedrake-venv/bin/activate \
  /opt/firedrake/bin/activate \
  /firedrake/bin/activate
do
  if [[ -f "$CANDIDATE" ]]; then
    echo "[env] Activating Firedrake via: $CANDIDATE"
    # shellcheck disable=SC1090
    source "$CANDIDATE"
    ACTIVATED=1
    break
  fi
done
if [[ "$ACTIVATED" -eq 0 ]]; then
  echo "[env] No activate script found (OK for many Firedrake images)"
fi

# Re-evaluate python after activation (some images change PATH)
if command -v python >/dev/null 2>&1; then
  PY=python
elif command -v python3 >/dev/null 2>&1; then
  PY=python3
fi
echo "[env] Python after activation: $PY ($(command -v "$PY"))"
$PY --version || true

echo "[env] Verifying Firedrake import with base Python..."
$PY - <<'EOF'
import sys, platform
print("[env] sys.executable:", sys.executable)
print("[env] platform.machine():", platform.machine())
import firedrake
print("[env] firedrake OK:", firedrake.__file__)
EOF

# -------------------------
# Create an isolated venv for pip installs (avoid apt-managed pip issues)
# Keep Firedrake available via --system-site-packages
# -------------------------
echo "[env] Creating local venv (with system-site-packages so Firedrake remains visible)..."
$PY -m venv --system-site-packages .venv_bench
# shellcheck disable=SC1091
source .venv_bench/bin/activate
PY=python

echo "[env] venv Python: $PY ($(command -v "$PY"))"
$PY --version
$PY -m pip --version

echo "[env] Installing system libs needed by gmsh..."
$SUDO apt-get update
$SUDO apt-get install -y libglu1-mesa libgl1

echo "[env] Installing Python deps into venv..."
$PY -m pip install -U pip
$PY -m pip install gmsh

echo "[env] Installing local editable packages into venv..."
$PY -m pip install -e ./icepack
$PY -m pip install -e ./modelfunc

echo "[env] Installing OpenMPI if needed..."
if ! command -v mpirun >/dev/null 2>&1; then
  $SUDO apt-get update
  $SUDO apt-get install -y openmpi-bin
else
  echo "[env] mpirun already present: $(command -v mpirun)"
fi

# Pin threading (important for reproducibility)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# -------------------------
# Run loop
# -------------------------
cd casestudy

OUT_BASE="../../../adviser_output"
BENCH_DIR="${OUT_BASE}/dx_${DX}"
mkdir -p "${BENCH_DIR}"

echo "[env] Output base: ${BENCH_DIR}"
echo "[env] Starting runs..."

for ((i=0; i<${TOTAL_RUNS}; i++)); do
  if [[ "$i" -eq 0 ]]; then
    echo "[Warm-up] run ${i}"
  else
    echo "[Measured] run ${i}/${REPEAT_TIMES}"
  fi

  $PY -m experiments.run_forward \
    --out "${BENCH_DIR}/trial_$(printf "%03d" "$i")" \
    --dx "${DX}"
done

echo "[Done] Benchmark finished."
echo "Warm-up result: ${BENCH_DIR}/trial_000 (ignore)"
echo "Measured runs:  ${BENCH_DIR}/trial_001 ... trial_$(printf "%03d" "${REPEAT_TIMES}")"
