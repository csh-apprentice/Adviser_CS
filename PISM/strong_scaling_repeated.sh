#!/usr/bin/env bash
set -euo pipefail

## PISM strong scaling (low-I/O) with per-np repeats via NP_LIST + REPEAT_LIST

cd "$HOME/pism-stable/examples/std-greenland"
export PATH="$HOME/pism/bin:$PATH"

# Only the head node runs mpirun tests
if [[ "${ADVISER_NODE_RANK:-}" != "0" ]]; then
  echo "[run] worker node rank=${ADVISER_NODE_RANK:-unknown} idle (waiting for head mpirun)"
  exit 0
fi

# if [[ "${ADVISER_NODE_RANK:-}" != "0" ]]; then
#   echo "[run] worker rank=${ADVISER_NODE_RANK:-unknown} holding (will be released by head)"
#   # Name the process so head can pkill it later
#   exec -a adviser_worker_hold sleep 365d
# fi


./preprocess.sh

echo "[run] head node"

# --- MPI / host setup ---
if [[ -z "${ADVISER_NODE_IPS:-}" ]]; then
  echo "[error] ADVISER_NODE_IPS is not set; are you running under adviser?"
  exit 1
fi


# release_workers() {
#   echo "[run] releasing workers..."
#   # ADVISER_NODE_IPS might be space-separated; normalize to one IP per line
#   while read -r ip; do
#     [[ -z "$ip" ]] && continue
#     # Skip self (head) — harmless if not skipped, but cleaner
#     if [[ "$ip" == "$(hostname -I | awk '{print $1}')" ]]; then
#       continue
#     fi
#     echo "[run] ssh $ip pkill adviser_worker_hold"
#     ssh -o StrictHostKeyChecking=no -o ConnectTimeout=5 "$ip" "pkill -f '^adviser_worker_hold' || true" || true
#   done < <(echo "${ADVISER_NODE_IPS:?ADVISER_NODE_IPS not set}" | tr ' ' '\n')
# }

# # Ensure we release even if the script errors out
# trap release_workers EXIT

# Count nodes (handle spaces or newlines)
num_nodes=$(
  echo "${ADVISER_NODE_IPS}" \
    | tr ' ' '\n' \
    | awk 'NF' \
    | wc -l \
    | tr -d ' '
)

# Only used for default NP_LIST; hosts are built per-np.
slots_per_node="$(nproc || echo 1)"
echo "[mpi] num_nodes=${num_nodes} slots_per_node=${slots_per_node}"

# NP_LIST: space-separated TOTAL MPI ranks to test
if [[ -z "${NP_LIST:-}" ]]; then
  NP_LIST="$(( num_nodes * slots_per_node ))"
fi

# REPEAT_LIST: space-separated repeats aligned with NP_LIST (same length)
# If not provided, defaults to 1 repeat for each np.
if [[ -z "${REPEAT_LIST:-}" ]]; then
  # Build "1 1 1 ..." matching NP_LIST length
  REPEAT_LIST="$(echo "${NP_LIST}" | awk '{for (i=1;i<=NF;i++) printf "1%s", (i==NF?ORS:OFS)}')"
fi

echo "[run] NP_LIST=${NP_LIST}"
echo "[run] REPEAT_LIST=${REPEAT_LIST}"

# Split NP_LIST and REPEAT_LIST into arrays
read -ra NP_ARR <<< "${NP_LIST}"
read -ra REP_ARR <<< "${REPEAT_LIST}"

# Sanity check
if [[ ${#NP_ARR[@]} -ne ${#REP_ARR[@]} ]]; then
  echo "[error] NP_LIST and REPEAT_LIST must have the same number of entries."
  echo "        NP_LIST=${NP_LIST}"
  echo "        REPEAT_LIST=${REPEAT_LIST}"
  echo "        (#NP=${#NP_ARR[@]} #REP=${#REP_ARR[@]})"
  exit 1
fi

# Where to copy results so adviser syncs them back
# NOTE: keep it flat (no subdirs) to avoid FS limitations
OUT_ROOT="/home/ubuntu/sky_workdir/adviser_output"
mkdir -p "${OUT_ROOT}"

# Helpful MPI env
export OMPI_MCA_rmaps_base_oversubscribe=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Helper: 2-digit repeat numbering
pad2() { printf "%02d" "$1"; }

for i in "${!NP_ARR[@]}"; do
  np="${NP_ARR[$i]}"
  nrep="${REP_ARR[$i]}"

  echo "==============================="
  echo "[run] Starting PISM runs with np=${np} (repeats=${nrep})"

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

  # ---- Measured repeats ----
  for r in $(seq 1 "${nrep}"); do
    rep="$(pad2 "${r}")"
    echo "[run] Measured repeat ${r}/${nrep} for np=${np} (rep=${rep})"

    # Unique file names per np and repeat
    out_file="g10km_10ka_np${np}_rep${rep}.nc"
    log_file="out.g10km_10ka_np${np}_rep${rep}.log"

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

    echo "[run] Finished np=${np} rep=${rep}, out=${out_file}, log=${log_file}"

    # Copy outputs to adviser_output so they sync back
    echo "[run] Copying outputs for np=${np} rep=${rep} to ${OUT_ROOT}"
    cp "${log_file}" "${OUT_ROOT}/"

    # Tar the NetCDF output so Adviser will sync it reliably
    tarball="${OUT_ROOT}/${out_file%.nc}.tar.gz"
    tar -czf "${tarball}" "${out_file}"
    echo "[run] Wrote ${tarball}"
  done

done

echo "[run] All NP_LIST runs finished. Outputs in ${OUT_ROOT}"


