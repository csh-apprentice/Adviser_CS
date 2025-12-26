DEBUG mode: we first create the cluster, then run the simulation from the existed cluster:

./adviser cluster create \
  --short-name pism \
  --cloud aws \
  --region us-west-2 \
  --num-nodes 4 \
  --instance-type t2.medium

./adviser cluster create \
  --short-name pism \
  --cloud aws \
  --region us-west-2 \
  --instance-type t2.medium


./adviser cluster create \
  --short-name intel \
  --cloud aws \
  --region us-west-2 \
  --instance-type c6i.2xlarge

./adviser cluster create \
  --short-name intel \
  --cloud aws \
  --region us-west-2 \
  --instance-type c6i.2xlarge


./adviser cluster create \
  --short-name pism \
  --cloud aws \
  --region us-west-2 \
  --num-nodes 4 \
  --instance-type c7a.2xlarge

./adviser run \
  --cluster 1197 \
  --cpu 8 \
   --container-image-uri docker.io/firedrakeproject/firedrake-vanilla-default:2025.10.2 \
  "
    set -euo pipefail
    rm -rf Adviser_CS
    git clone --recurse-submodules https://github.com/csh-apprentice/Adviser_CS.git
    cd Adviser_CS/ICEPACK
    chmod +x test_parallel_scaling.sh
    bash test_parallel_scaling.sh
  "

./adviser cluster create \
  --short-name pism \
  --cloud aws \
  --region us-west-2 \
  --num-nodes 4 \
  --instance-type c7g.2xlarge


./adviser cluster create \
  --short-name icepack \
  --cloud aws \
  --region us-west-2 \
  --instance-type c6i.4xlarge



./adviser run \
  --cluster 1157 \
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
  --cluster 1157 \
   --container-image-uri docker.io/firedrakeproject/firedrake-vanilla:2025-01 \
  "
    set -euo pipefail
    cd Adviser_CS/ICEPACK
    chmod +x run_forward.sh
    ./run_forward.sh
  "

./adviser run \
  --cluster 1157 \
   --container-image-uri docker.io/firedrakeproject/firedrake-vanilla:2025-01 \
  "
    cp icepack_debug.sh Adviser_CS/ICEPACK
    cd Adviser_CS/ICEPACK
    chmod +x icepack_debug.sh
    ./icepack_debug.sh
  "


./adviser run \
  --cluster 1156 \
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
  --cluster 1156 \
  "
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


./adviser run \
  --cluster 1156 \
  --num-nodes 4 \
  "
    cd Adviser_CS/PISM
    chmod +x first_run_new.sh
    ./first_run_new.sh
  "