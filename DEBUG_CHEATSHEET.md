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
  --cluster 1148 \
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
  --cluster 1146 \
  "
    git clone --recurse-submodules https://github.com/csh-apprentice/Adviser_CS.git
    cd Adviser_CS/PISM
    chmod +x setup_pism.sh
    ./setup_pism.sh
  "