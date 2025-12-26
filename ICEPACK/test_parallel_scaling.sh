#!/usr/bin/env bash
set -euo pipefail

##############################################
# MPI Parallel Scaling Test for Icepack
# - Warm-up + measured run per NP
# - Pure MPI (threads pinned to 1)
##############################################

# -------------------------
# Config (edit if you want)
# -------------------------
DX="${1:-2000}"                     # default dx=2000, allow override via arg
NP_LIST="${NP_LIST:-"1 2 4 8"}"     # override via env: NP_LIST="1 2 4 8 16"
MEASURE_RUNS="${MEASURE_RUNS:-1}"   # how many measured runs per NP (default 1)


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
# Ensure mpiexec exists
# -------------------------
if ! command -v mpiexec >/dev/null 2>&1; then
  echo "[env] mpiexec not found — installing OpenMPI..."
  if command -v sudo >/dev/null 2>&1; then
    sudo apt-get update
    sudo apt-get install -y openmpi-bin
  else
    apt-get update
    apt-get install -y openmpi-bin
  fi
fi

MPIEXEC=${MPIEXEC:-mpiexec}
echo "[env] MPI launcher: ${MPIEXEC} ($(command -v "${MPIEXEC}"))"

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
    python -m experiments.run_forward \
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
    /usr/bin/time -f "[time] NP=${NP} ${RUN_NAME} elapsed=%e sec" \
      "${MPIEXEC}" -n "${NP}" \
        python -m experiments.run_forward \
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
