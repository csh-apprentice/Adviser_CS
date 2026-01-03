#!/usr/bin/env bash
set -euo pipefail

## First Run (low-I/O scaling version)

cd "$HOME/pism-stable/examples/std-greenland"
export PATH="$HOME/pism/bin:$PATH"

./preprocess.sh

# Only the head node runs mpirun tests
if [[ "${ADVISER_NODE_RANK:-}" != "0" ]]; then
  echo "[run] worker node rank=${ADVISER_NODE_RANK:-unknown} idle (waiting for head mpirun)"
  exit 0
fi

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

slots_per_node="$(nproc)"   # maximum ranks we allow per node
hosts="$(
  echo "${ADVISER_NODE_IPS}" \
    | tr ' ' '\n' \
    | awk -v slots="${slots_per_node}" 'NF{print $0 ":" slots}' \
    | paste -sd, -
)"

echo "[mpi] num_nodes=${num_nodes} slots_per_node=${slots_per_node}"
echo "[mpi] hosts=${hosts}"

# NP_LIST is a space-separated list of TOTAL MPI ranks to test, e.g. "1 2 4 8 16"
# If not set, default to using all cores on all nodes once.
if [[ "${NP_LIST:-}" == "" ]]; then
  NP_LIST="$(( num_nodes * slots_per_node ))"
fi
echo "[run] NP_LIST=${NP_LIST}"

# Where to copy results so adviser syncs them back
OUT_ROOT="/home/ubuntu/sky_workdir/adviser_output/strong_scale"
mkdir -p "${OUT_ROOT}"

for np in ${NP_LIST}; do
  echo "==============================="
  echo "[run] Starting PISM run with np=${np}"

  # Require that np is divisible by num_nodes so ppr:ranks_per_node:node makes sense
  if (( np % num_nodes != 0 )); then
    echo "[warn] np=${np} is not divisible by num_nodes=${num_nodes}; skipping this np."
    continue
  fi

  ranks_per_node=$(( np / num_nodes ))

#   if (( ranks_per_node > slots_per_node )); then
#     echo "[warn] np=${np} implies ranks_per_node=${ranks_per_node} > slots_per_node=${slots_per_node}; skipping to avoid oversubscribe."
#     continue
#   fi

  echo "[mpi] np=${np} ranks_per_node=${ranks_per_node}"

  # Unique file names per np to avoid clashes
  out_file="g20km_10ka_np${np}.nc"
  log_file="out.g20km_10ka_np${np}.log"

  # --- LOW-I/O PISM RUN: only -o, no ts/extra files ---
  mpirun -np "${np}" -H "${hosts}" \
    --map-by "ppr:${ranks_per_node}:node" --bind-to core \
    "$HOME/pism/bin/pism" \
      -i pism_Greenland_5km_v1.1.nc -bootstrap -grid.registration corner \
      -dx 10km -dy 10km \
      -Mz 201 -Mbz 21 -z_spacing equal -Lz 4000 -Lbz 2000 \
      -skip -skip_max 20 -grid.recompute_longitude_and_latitude false \
      -ys -200 -ye 0 \
      -surface given -surface_given_file pism_Greenland_5km_v1.1.nc \
      -front_retreat_file pism_Greenland_5km_v1.1.nc \
      -sia_e 3.0 \
      -o "${out_file}" \
      &> "${log_file}"

  echo "[run] Finished PISM run with np=${np}, out=${out_file}, log=${log_file}"

  # Copy outputs to adviser_output so they sync back
  echo "[run] Copying outputs for np=${np} to ${OUT_ROOT}"
  cp "${out_file}" "${log_file}" "${OUT_ROOT}/"
done

echo "[run] All NP_LIST runs finished. Outputs in ${OUT_ROOT}"
