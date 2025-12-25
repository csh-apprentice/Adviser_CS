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
# Environment setup
# -------------------------
source /home/firedrake/firedrake/bin/activate

python -m pip install --upgrade pip
python -m pip install gmsh

pip install -e ./icepack
pip install -e ./modelfunc

sudo apt-get update
sudo apt-get install -y openmpi-bin

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

for ((i=0; i<${TOTAL_RUNS}; i++)); do
  if [[ "$i" -eq 0 ]]; then
    echo "[Warm-up] run ${i}"
  else
    echo "[Measured] run ${i}/${REPEAT_TIMES}"
  fi

  python -m experiments.run_forward \
    --out "${BENCH_DIR}/trial_$(printf "%03d" "$i")" \
    --dx "${DX}"
done

echo "[Done] Benchmark finished."
echo "Warm-up result: ${BENCH_DIR}/trial_000 (ignore)"
echo "Measured runs:  ${BENCH_DIR}/trial_001 ... trial_$(printf "%03d" "${REPEAT_TIMES}")"