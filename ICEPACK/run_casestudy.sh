#!/usr/bin/env bash
set -euo pipefail

echo "[Running] Run case study..."

source /home/firedrake/firedrake/bin/activate

python -m pip install --upgrade pip
python -m pip install gmsh

pip install -e ./icepack
pip install -e ./modelfunc

sudo apt-get update
sudo apt-get install -y openmpi-bin

cd casestudy
for a in 0.3 0.5 0.7 1.0 1.3 1.6 2.0; do
  python -m experiments.run_forward \
    --out ../adviser_output/sweep_A_${a} \
    --fluidity-scale "$a" \
    --dx 5000
done


for dx in 8000 6000 5000 4000 3000; do
  python -m experiments.run_forward \
    --out adviser_output/scale_dx_${dx} \
    --fluidity-scale 1.0 \
    --dx $dx
done



# for np in 1 2 4 8 16; do
#   mpirun -np $np \
#     python -m experiments.run_forward \
#       --out adviser_output/mpi_np_${np} \
#       --fluidity-scale 1.0 \
#       --dx 4000
# done

for np in 1 2 4 8 16; do
  mpiexec -n $np \
    python -m experiments.run_forward \
      --out adviser_output/mpi_np_${np} \
      --fluidity-scale 1.0 \
      --dx 4000
done