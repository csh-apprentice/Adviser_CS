#!/usr/bin/env bash
set -euo pipefail

## Weak scaling version:
## - Reference case: np_ref=1, dx_ref=dy_ref=10km
## - For each np in NP_LIST, dx,dy are set so that work per rank is ~constant:
##   dx = dx_ref / sqrt(np / np_ref)

cd "$HOME/pism-stable/examples/std-greenland"
export PATH="$HOME/pism/bin:$PATH"

./preprocess.sh

# Only the head node runs mpirun tests
if [[ "${ADVISER_NODE_RANK:-}" != "0" ]]; then
  echo "[run] worker node rank=${ADVISER_NODE_RANK:-unknown} idle (waiting for head mpirun)"
  exit 0
fi

echo "[run] head node (weak scaling)"

# --- MPI / host setup ---

if [[ -z "${ADVISER_NODE_IPS:-}" ]]; then
  echo "[error] ADVISER_NODE_IPS is not set; are you running under adviser?"
  exit 1
fi

# Count nodes (handle spaces or newlines)
num_nodes=$(
  echo "${ADVISER_NODE_IPS}" \
    | tr ' ' '\n' \
    | awk 'NF' \
    | wc -l \
    | tr -d ' '
)

slots_per_node="$(nproc || echo 1)"

echo "[mpi] num_nodes=${num_nodes} slots_per_node=${slots_per_node}"

# --- Weak-scaling parameters ---

# Reference case: np_ref=1 uses dx_ref=10 km
np_ref=${NP_REF:-1}
dx_ref_km=${DX_REF_KM:-10}

# NP_LIST is a space-separated list of TOTAL MPI ranks to test in weak scaling.
# Default to classic 1,4,16,64 (can override via env).
if [[ "${NP_LIST:-}" == "" ]]; then
  NP_LIST="1 4 16 64"
fi
echo "[run] NP_LIST (weak scaling) = ${NP_LIST}"
echo "[run] Reference: np_ref=${np_ref}, dx_ref=${dx_ref_km} km"

# Where to copy results so adviser syncs them back
OUT_ROOT="/home/ubuntu/sky_workdir/adviser_output"
mkdir -p "${OUT_ROOT}"

# Helpful MPI env
export OMPI_MCA_rmaps_base_oversubscribe=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

for np in ${NP_LIST}; do
  echo "==============================="
  echo "[run] WEAK scaling run: np=${np}"

  # Require that np is divisible by num_nodes so ppr:ranks_per_node:node makes sense
  if (( np % num_nodes != 0 )); then
    echo "[warn] np=${np} is not divisible by num_nodes=${num_nodes}; skipping this np."
    continue
  fi

  ranks_per_node=$(( np / num_nodes ))
  echo "[mpi] np=${np} ranks_per_node=${ranks_per_node}"

  # --- Compute dx,dy for weak scaling ---
  # dx_km = dx_ref_km / sqrt(np / np_ref)
  # Use bc -l for floating-point + sqrt
  ratio=$(echo "${np} / ${np_ref}" | bc -l)
  dx_km=$(echo "${dx_ref_km} / sqrt(${ratio})" | bc -l)

  # This will look like "10.000000", "5.000000", "2.500000", ...
  echo "[weak] np=${np}: dx=dy=${dx_km} km (target)"

  dx_arg="${dx_km}km"
  dy_arg="${dx_km}km"

  # Build hosts for THIS np, with exactly ranks_per_node slots per node
  hosts="$(
    echo "${ADVISER_NODE_IPS}" \
      | tr ' ' '\n' \
      | awk -v rpn="${ranks_per_node}" 'NF{print $0 ":" rpn}' \
      | paste -sd, -
  )"
  echo "[mpi] hosts=${hosts}"

  # Make a dx tag for filenames, e.g., 10.000000 -> 10p000000
  dx_tag=$(echo "${dx_km}" | tr '.' 'p')

  # Unique file names per np/dx to avoid clashes
  out_file="weak_dx${dx_tag}km_np${np}.nc"
  log_file="weak_dx${dx_tag}km_np${np}.log"

  # --- LOW-I/O PISM RUN with weak scaling dx/dy ---
  mpirun -np "${np}" -H "${hosts}" \
    --map-by "ppr:${ranks_per_node}:node" --bind-to core \
    "$HOME/pism/bin/pism" \
      -i pism_Greenland_5km_v1.1.nc -bootstrap -grid.registration corner \
      -dx "${dx_arg}" -dy "${dy_arg}" \
      -Mz 201 -Mbz 21 -z_spacing equal -Lz 4000 -Lbz 2000 \
      -skip -skip_max 20 -grid.recompute_longitude_and_latitude false \
      -ys -200 -ye 0 \
      -surface given -surface_given_file pism_Greenland_5km_v1.1.nc \
      -front_retreat_file pism_Greenland_5km_v1.1.nc \
      -sia_e 3.0 \
      -o "${out_file}" \
      &> "${log_file}"

  echo "[run] Finished WEAK scaling run np=${np}, out=${out_file}, log=${log_file}"

  # Copy outputs to adviser_output so they sync back
  echo "[run] Copying outputs for np=${np} to ${OUT_ROOT}"
  cp "${log_file}" "${OUT_ROOT}/"

  # tar the NetCDF output (binary) so Adviser will sync it
  tar -czf "${OUT_ROOT}/${out_file%.nc}.tar.gz" "${out_file}"

done

echo "[run] All WEAK scaling runs finished. Outputs in ${OUT_ROOT}"
