#!/usr/bin/env bash
set -euo pipefail

##############################################
# MPI Parallel Scaling Test for Icepack
# - Environment setup similar to run_benchmark.sh
# - Warm-up + measured runs per NP
# - Pure MPI (threads pinned to 1)
##############################################

# -------------------------
# Config (edit if you want)
# -------------------------
DX="${1:-2000}"                     # default dx=2000, allow override via arg
NP_LIST="${NP_LIST:-"1 2 4 8"}"     # override via env: NP_LIST="1 2 4 8 16"
MEASURE_RUNS="${MEASURE_RUNS:-1}"   # how many measured runs per NP (default 1)

OUT_BASE="../../../adviser_output/mpi_scaling_dx_${DX}"

echo "=============================================="
echo "  MPI Parallel Scaling Test"
echo "  dx          = ${DX}"
echo "  NP list     = ${NP_LIST}"
echo "  meas. runs  = ${MEASURE_RUNS} per NP"
echo "  Output base = ${OUT_BASE}"
echo "=============================================="

mkdir -p "${OUT_BASE}"

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
# Environment inspection + optional Firedrake activation
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
# Create an isolated venv for pip installs
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

echo "[env] Installing system libs needed by gmsh (OpenGL + X11/font deps)..."
$SUDO apt-get update
$SUDO apt-get install -y \
  libglu1-mesa libgl1 \
  libxft2 libxrender1 libxext6 libsm6 libice6 \
  libfontconfig1

echo "[env] Installing Python deps into venv..."
$PY -m pip install -U pip
$PY -m pip install gmsh

echo "[env] Installing local editable packages into venv..."
$PY -m pip install -e ./icepack
$PY -m pip install -e ./modelfunc

# -------------------------
# Ensure mpiexec exists (OpenMPI)
# -------------------------
echo "[env] Installing (or refreshing) OpenMPI..."
$SUDO apt-get update
$SUDO apt-get install -y openmpi-bin

MPIEXEC=${MPIEXEC:-mpiexec}
echo "[env] MPI launcher: ${MPIEXEC} ($(command -v "${MPIEXEC}"))"

# Allow oversubscription (safe if NP > nproc)
export OMPI_MCA_rmaps_base_oversubscribe=1


# -------------------------
# Avoid thread oversubscription (pure MPI mode)
# -------------------------
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# -------------------------
# Go to casestudy
# -------------------------
cd casestudy
MAX_NP="$(nproc)"
echo "[env] nproc (max MPI ranks without oversubscribe): ${MAX_NP}"


# -------------------------
# Loop over NP values
# -------------------------
for NP in ${NP_LIST}; do
  OUT_DIR="${OUT_BASE}/np_${NP}"
  mkdir -p "${OUT_DIR}"

  echo ""
  echo "----------------------------------------------"
  echo "NP = ${NP}"
  echo "Output dir: ${OUT_DIR}"
  echo "----------------------------------------------"

  ###############################
  # 1) Warm-up run (ignore timing)
  ###############################
  echo "[Warm-up] NP=${NP}"

  "${MPIEXEC}" -n "${NP}" \
    "${PY}" -m experiments.run_forward \
      --out "${OUT_DIR}/warmup" \
      --dx "${DX}" \
      --fluidity-scale 1.0

  #######################################
  # 2) Measured runs (for actual scaling)
  #######################################
  for ((k=1; k<=MEASURE_RUNS; k++)); do
    RUN_NAME="run_$(printf "%02d" "${k}")"
    echo "[Measured] NP=${NP}, ${RUN_NAME}"

    # /usr/bin/time just prints to stderr; stats.json still has wall_s, etc.
      "${MPIEXEC}" -n "${NP}" \
        "${PY}" -m experiments.run_forward \
          --out "${OUT_DIR}/${RUN_NAME}" \
          --dx "${DX}" \
          --fluidity-scale 1.0
  done

  echo "[NP=${NP}] Done. Check ${OUT_DIR}/run_*/stats.json for wall_s."
done

echo ""
echo "=============================================="
echo " All NP values done."
echo " Results under: ${OUT_BASE}"
echo " Each np_<N>/run_*/stats.json has wall_s, etc."
echo " (warmup/ can be ignored for timing analysis)"
echo "=============================================="
