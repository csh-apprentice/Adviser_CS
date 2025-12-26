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

...

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
else
  echo "[env] ERROR: neither python nor python3 found in PATH after activation"
  echo "[env] PATH=$PATH"
  exit 127
fi

echo "[env] Python after activation: $PY ($(command -v "$PY"))"
$PY --version || true

# -------------------------
# Create + activate a local venv for gmsh / extra deps
# -------------------------
if [[ ! -d ".venv_bench" ]]; then
  echo "[env] Creating local venv .venv_bench..."
  $PY -m venv .venv_bench
fi

# shellcheck disable=SC1091
source .venv_bench/bin/activate

echo "[env] Using benchmark venv Python: $(command -v python)"
PY=python
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
# Force OpenMPI (match run_casestudy.sh)
# -------------------------
echo "[env] Installing (or refreshing) OpenMPI..."
$SUDO apt-get update
$SUDO apt-get install -y openmpi-bin

MPIEXEC=mpiexec
echo "[env] Using OpenMPI launcher: ${MPIEXEC} ($(command -v "${MPIEXEC}"))"

# Allow oversubscribe safely if NP > nproc (harmless if not needed)
export OMPI_MCA_rmaps_base_oversubscribe=1

# Pin threading (important for reproducibility)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# -------------------------
# Detect number of physical cores and set NP accordingly
# -------------------------
echo "[env] Detecting CPU topology for MPI ranks..."

# 1) Get the list of allowed logical CPUs (respects cgroups/affinity)
ALLOWED_LIST=$(grep Cpus_allowed_list /proc/self/status | awk '{print $2}')
echo "[env] Cpus_allowed_list: ${ALLOWED_LIST:-<none>}"

count_logical_from_list() {
  local list="$1"
  local total=0
  IFS=',' read -ra parts <<< "$list"
  for part in "${parts[@]}"; do
    if [[ "$part" == *-* ]]; then
      IFS='-' read -r start end <<< "$part"
      total=$(( total + end - start + 1 ))
    elif [[ -n "$part" ]]; then
      total=$(( total + 1 ))
    fi
  done
  echo "$total"
}

if [[ -n "$ALLOWED_LIST" && "$ALLOWED_LIST" != "0" ]]; then
  LOGICAL_CPUS=$(count_logical_from_list "$ALLOWED_LIST")
else
  # Fallback: count processors in /proc/cpuinfo
  LOGICAL_CPUS=$(grep -c '^processor' /proc/cpuinfo || echo 1)
fi

echo "[env] Logical CPUs visible: ${LOGICAL_CPUS}"

# 2) Determine threads per core (if lscpu is available), else assume 2
THREADS_PER_CORE=2
if command -v lscpu >/dev/null 2>&1; then
  tpc=$(lscpu | awk '/^Thread\(s\) per core:/ {print $4; exit}')
  if [[ -n "$tpc" ]]; then
    THREADS_PER_CORE="$tpc"
  fi
fi
echo "[env] Threads per core (assumed/detected): ${THREADS_PER_CORE}"

# 3) Physical cores = logical / threads_per_core (rounded up)
PHYSICAL_CORES=$(( (LOGICAL_CPUS + THREADS_PER_CORE - 1) / THREADS_PER_CORE ))
if (( PHYSICAL_CORES < 1 )); then
  PHYSICAL_CORES=1
fi

NP="${PHYSICAL_CORES}"
echo "[env] Using NP (MPI ranks) = ${NP} (one rank per physical core)"

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

  ${MPIEXEC} -n "${NP}" \
    "$PY" -m experiments.run_forward \
      --out "${BENCH_DIR}/trial_$(printf "%03d" "$i")" \
      --dx "${DX}"
done

echo "[Done] Benchmark finished."
echo "Warm-up result: ${BENCH_DIR}/trial_000 (ignore)"
echo "Measured runs:  ${BENCH_DIR}/trial_001 ... trial_$(printf "%03d" "${REPEAT_TIMES}")"
