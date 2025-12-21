#!/usr/bin/env bash
set -euo pipefail

echo "[setup] Updating apt and installing prerequisites..."

sudo apt-get update
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  cmake pkg-config make g++ git \
  libfftw3-dev libgsl-dev libnetcdf-dev libudunits2-dev netcdf-bin \
  petsc-dev \
  cdo cmake-curses-gui libpnetcdf-dev libproj-dev libx11-dev nco ncview \
  python3-dev python3-netcdf4 python3-nose python3-numpy python3-pyproj python3-scipy \
  swig

echo "[setup] Setting PETSC_DIR..."
export PETSC_DIR=/usr/lib/petsc
echo 'export PETSC_DIR=/usr/lib/petsc' >> "$HOME/.bashrc"

echo "[setup] Cloning and building PISM..."
cd "$HOME"
if [ ! -d pism-stable ]; then
  git clone https://github.com/pism/pism.git pism-stable
fi

mkdir -p pism-stable/build
cd pism-stable/build


export CC=mpicc
export CXX=mpicxx
cmake -DCMAKE_INSTALL_PREFIX="$HOME/pism" ..

make -j"$(nproc)" install

echo "[setup] Adding PISM to PATH..."
echo 'export PATH=$HOME/pism/bin:$PATH' >> "$HOME/.bashrc"

export PATH="$HOME/pism/bin:$PATH"

echo "[setup] Done. PISM location:"
which pism || echo "pism not on PATH"


## Step 3: Quick “did it work?” tests (local)
# --- Step 3: Quick “did it work?” tests (Adviser/SkyPilot multi-node) ---

# if [ "${ADVISER_NODE_RANK:-}" = "0" ]; then
#   echo "[run] head node"

#   ranks_per_node=${RANKS_PER_NODE:-2}

#   num_nodes=$(echo "${ADVISER_NODE_IPS:?ADVISER_NODE_IPS missing}" | wc -l | tr -d ' ')
#   total_ranks=$((num_nodes * ranks_per_node))

#   hosts=$(echo "$ADVISER_NODE_IPS" | awk -v s=$ranks_per_node '{printf "%s:%d,", $1, s}' | sed 's/,$//')


#   echo "[run] num_nodes=$num_nodes hosts=$hosts total_ranks=$total_ranks"

#   # sanity check placement
#   mpirun -np "$total_ranks" -H "$hosts" \
#     --map-by "ppr:${ranks_per_node}:node" --bind-to core \
#     hostname

# #   # actual PISM run
# #   mpirun -np "$total_ranks" -H "$hosts" \
# #     --map-by "ppr:${ranks_per_node}:node" --bind-to core \
# #     $HOME/pism/bin/pism -test G -y 200
 
#   iface=$(ip -o -4 route show to default | awk '{print $5; exit}')

#   mpirun -np "$total_ranks" -H "$hosts" \
#     --mca btl tcp,self \
#     --mca btl_tcp_if_include "$iface" \
#     --mca oob_tcp_if_include "$iface" \
#     --mca oob_tcp_disable_family ipv6 \
#     --map-by "ppr:${ranks_per_node}:node" --bind-to core \
#     $HOME/pism/bin/pism -test G -y 200
# else
#   echo "[run] worker node rank=${ADVISER_NODE_RANK:-unknown} (idle, waiting for mpirun)"
# fi
# mpiexec -n 4 /home/shihan/pism/bin/pism -test G -y 200

## First Run


# ranks_per_node=${RANKS_PER_NODE:-2}

# num_nodes=$(echo "${ADVISER_NODE_IPS:?ADVISER_NODE_IPS missing}" | wc -l | tr -d ' ')
# total_ranks=$((num_nodes * ranks_per_node))

# hosts=$(echo "$ADVISER_NODE_IPS" | awk -v s=$ranks_per_node '{printf "%s:%d,", $1, s}' | sed 's/,$//')


# echo "[run] num_nodes=$num_nodes hosts=$hosts total_ranks=$total_ranks"


# cd $HOME/pism-stable/examples/std-greenland

# ./spinup.sh