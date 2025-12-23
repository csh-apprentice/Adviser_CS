DEBUG mode: we first create the cluster, then run the simulation from the existed cluster:

./adviser cluster create \
  --short-name pism \
  --cloud aws \
  --region us-west-2 \
  --num-nodes 4 \
  --instance-type t2.medium


./adviser cluster create \
  --short-name icepack \
  --cloud aws \
  --region us-west-2 \
  --instance-type c6i.4xlarge



./adviser run \
  --cluster 1155 \
   --container-image-uri docker.io/firedrakeproject/firedrake-vanilla:2025-01 \
  "
    set -euo pipefail
    rm -rf Adviser_CS
    git clone --recurse-submodules https://github.com/csh-apprentice/Adviser_CS.git
    cd Adviser_CS/ICEPACK
    chmod +x run_inverse.sh
    ./run_inverse.sh
  "


./adviser run \
  --cluster 1154 \
  "
    git clone --recurse-submodules https://github.com/csh-apprentice/Adviser_CS.git
    cd Adviser_CS/PISM
    chmod +x setup_pism.sh
    ./setup_pism.sh
  "

./adviser run \
  --cluster 1146 \
  --num-nodes 4 \
  "
    cd Adviser_CS/PISM
    chmod +x first_run_new.sh
    ./first_run_new.sh
  "


./adviser run \
  --cluster 1146 \
  --num-nodes 4 \
  "
    chmod +x debug.sh
    ./debug.sh
  "



./adviser run \
  --cluster 1153 \
  "
    git clone --recurse-submodules https://github.com/csh-apprentice/Adviser_CS.git
    cd Adviser_CS/PISM
    chmod +x setup_pism.sh
    ./setup_pism.sh
  "


./adviser run \
  --cluster 1155 \
  "
    chmod +x debug_new.sh
    ./debug_new.sh
  "

./adviser run \
  --cluster 1153 \
  --num-nodes 4 \
  "
    chmod +x debug_new.sh
    ./debug_new.sh
  "

./adviser run \
  --cluster 1154 \
  "
    ls
    cp first_run_debug.sh Adviser_CS/PISM
    cd Adviser_CS/PISM
    chmod +x first_run_debug.sh
    ./first_run_debug.sh
  "

./adviser run \
  --cluster 1154 \
  --num-nodes 4 \
  "
    pwd
    cd Adviser_CS/PISM
    chmod +x first_run_new.sh
    ./first_run_new.sh
  "

