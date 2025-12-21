#!/usr/bin/env bash
set -euo pipefail


## First Run


cd $HOME/pism-stable/examples/std-greenland
export PATH="$HOME/pism/bin:$PATH"
./preprocess.sh


# Only head runs mpirun tests
if [[ "${ADVISER_NODE_RANK:-}" != "0" ]]; then
  echo "[run] worker node rank=${ADVISER_NODE_RANK:-unknown} idle (waiting for head mpirun)"
  exit 0
fi

echo "[run] head node"



num_nodes=$(echo "${ADVISER_NODE_IPS:?ADVISER_NODE_IPS not set}" | wc -l | tr -d ' ')
ranks_per_node="${RANKS_PER_NODE:-2}"
total_ranks=$(( num_nodes * ranks_per_node ))
hosts="$(echo "$ADVISER_NODE_IPS" | tr ' ' '\n' | awk 'NF{print $0 ":'"$ranks_per_node"'"}' | paste -sd, -)"
echo "[mpi] num_nodes=${num_nodes} ranks_per_node=${ranks_per_node} total_ranks=${total_ranks}"
echo "[mpi] hosts=${hosts}"
# Pick the default-route IPv4 interface (usually eth0 on AWS)
iface=$(ip -o -4 route show to default | awk '{print $5; exit}')
echo "[run] default-route iface=$iface"
ip -o -4 addr show dev "$iface" || true

# --- Key part: pin OpenMPI to the right NIC & disable IPv6 (NO --mca flags) ---
export OMPI_MCA_btl="tcp,self"
export OMPI_MCA_btl_tcp_if_include="$iface"
export OMPI_MCA_oob_tcp_if_include="$iface"
export OMPI_MCA_oob_tcp_disable_family="ipv6"
# optional: get all help messages (useful while debugging)
export OMPI_MCA_orte_base_help_aggregate="0"

mpirun -np "$total_ranks" -H "$hosts" \
    --map-by "ppr:${ranks_per_node}:node" --bind-to core \
    $HOME/pism/bin/pism \
        -i pism_Greenland_5km_v1.1.nc -bootstrap -grid.registration corner \
        -dx 20km -dy 20km -Mz 101 -Mbz 11 -z_spacing equal -Lz 4000 -Lbz 2000 \
        -skip -skip_max 10 -grid.recompute_longitude_and_latitude false -ys -10000 -ye 0 \
        -surface given -surface_given_file pism_Greenland_5km_v1.1.nc \
        -front_retreat_file pism_Greenland_5km_v1.1.nc -sia_e 3.0 \
        -ts_file ts_g20km_10ka.nc -ts_times -10000:yearly:0 \
        -extra_file ex_g20km_10ka.nc -extra_times -10000:100:0 \
        -extra_vars diffusivity,temppabase,tempicethk_basal,bmelt,tillwat,velsurf_mag,mask,thk,topg,usurf \
        -o g20km_10ka.nc &> out.g20km_10ka &
else
  echo "[run] worker node rank=${ADVISER_NODE_RANK:-unknown} (idle, waiting for mpirun)"
fi