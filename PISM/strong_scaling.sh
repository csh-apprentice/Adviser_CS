#!/usr/bin/env bash
set -euo pipefail

## First Run (low-I/O scaling version)

cd "$HOME/pism-stable/examples/std-greenland"
export PATH="$HOME/pism/bin:$PATH"

# Only the head node runs mpirun tests
if [[ "${ADVISER_NODE_RANK:-}" != "0" ]]; then
  echo "[run] worker node rank=${ADVISER_NODE_RANK:-unknown} idle (waiting for head mpirun)"
  exit 0
fi


./preprocess.sh


echo "[run] head node"

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

# This is only used for the default NP_LIST; hosts are now built per-np.
slots_per_node="$(nproc || echo 1)"

echo "[mpi] num_nodes=${num_nodes} slots_per_node=${slots_per_node}"

# NP_LIST is a space-separated list of TOTAL MPI ranks to test, e.g. "1 2 4 8 16"
# If not set, default to using all cores on all nodes once.
if [[ "${NP_LIST:-}" == "" ]]; then
  NP_LIST="$(( num_nodes * slots_per_node ))"
fi
echo "[run] NP_LIST=${NP_LIST}"

# Where to copy results so adviser syncs them back
OUT_ROOT="/home/ubuntu/sky_workdir/adviser_output"
mkdir -p "${OUT_ROOT}"

# Helpful MPI env (optional but nice)
export OMPI_MCA_rmaps_base_oversubscribe=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

for np in ${NP_LIST}; do
  echo "==============================="
  echo "[run] Starting PISM run with np=${np}"

  # Require that np is divisible by num_nodes so ppr:ranks_per_node:node makes sense
  if (( np % num_nodes != 0 )); then
    echo "[warn] np=${np} is not divisible by num_nodes=${num_nodes}; skipping this np."
    continue
  fi

  ranks_per_node=$(( np / num_nodes ))
  echo "[mpi] np=${np} ranks_per_node=${ranks_per_node}"

  # Build hosts for THIS np, with exactly ranks_per_node slots per node
  hosts="$(
    echo "${ADVISER_NODE_IPS}" \
      | tr ' ' '\n' \
      | awk -v rpn="${ranks_per_node}" 'NF{print $0 ":" rpn}' \
      | paste -sd, -
  )"
  echo "[mpi] hosts=${hosts}"

  # Unique file names per np to avoid clashes
  out_file="g10km_10ka_np${np}.nc"
  log_file="out.g10km_10ka_np${np}.log"

  # --- LOW-I/O PISM RUN: only -o, no ts/extra files ---
  mpirun -np "${np}" -H "${hosts}" \
    --map-by "ppr:${ranks_per_node}:node" --bind-to core \
    "$HOME/pism/bin/pism" \
      -i pism_Greenland_5km_v1.1.nc -bootstrap -grid.registration corner \
      -dx 10km -dy 10km \
      -Mz 201 -Mbz 21 -z_spacing equal -Lz 4000 -Lbz 2000 \
      -skip -skip_max 20 -grid.recompute_longitude_and_latitude false \
      -ys -10000 -ye 0 \
      -surface given -surface_given_file pism_Greenland_5km_v1.1.nc \
      -front_retreat_file pism_Greenland_5km_v1.1.nc \
      -sia_e 3.0 \
      -o "${out_file}" \
      &> "${log_file}"

  echo "[run] Finished PISM run with np=${np}, out=${out_file}, log=${log_file}"

  # Copy outputs to adviser_output so they sync back
  echo "[run] Copying outputs for np=${np} to ${OUT_ROOT}"
  cp "${log_file}" "${OUT_ROOT}/"

  # tar the NetCDF output (binary) so Adviser will sync it
  tar -czf "${OUT_ROOT}/${out_file%.nc}.tar.gz" "${out_file}"

done

echo "[run] All NP_LIST runs finished. Outputs in ${OUT_ROOT}"
