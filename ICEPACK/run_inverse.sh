#!/usr/bin/env bash
set -euo pipefail

echo "[Running] Run case study (cloud-safe MPI defaults)..."

# # ----------------------------
# # Helper: identify head node
# # ----------------------------
# NODE_RANK="${ADVISER_NODE_RANK:-0}"

# ----------------------------
# Environment setup (head + workers)
# ----------------------------
# If you truly need these deps on *every* node (e.g., if mpirun launches python on workers),
# it's okay to install on every node. Otherwise you can restrict installs to rank 0.
source /home/firedrake/firedrake/bin/activate

python -m pip install --upgrade pip
python -m pip install gmsh

pip install -e ./icepack
pip install -e ./modelfunc

# Make sure mpirun exists (needed on head at minimum; on some systems also needed on workers)
# if ! command -v mpirun >/dev/null 2>&1; then
#   echo "[setup] mpirun not found; installing Open MPI..."
#   sudo apt-get update
#   sudo apt-get install -y openmpi-bin
# fi
# sudo apt-get update
# sudo apt-get install -y openmpi-bin



cd RWArchive/inversionInputs
bash runInv > inversion_$(date +%Y%m%d_%H%M%S).log 2>&1
