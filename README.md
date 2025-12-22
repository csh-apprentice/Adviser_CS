# ICEPACK
## Icepack Case Study (Synthetic Ice Shelf)
ICEPACK/casestuy contains a small, headless (non-notebook) Icepack/Firedrake workload that is designed for PEARC-style platform evaluation. The main entry point is the script:
```
run_casestudy.sh
```
This script sets up the environment, installs required dependencies, and runs a sequence of forward-model experiments to study parameter sensitivity, resolution scaling, and single-node, multi-core scaling.

*Running the Demonstration*


```
./adviser run \
  --region us-west-2 \
  --instance-type c6i.4xlarge \
  --container-image-uri docker.io/firedrakeproject/firedrake-vanilla:2025-01 \
  -- \
  bash -lc '
    set -euxo pipefail
    git clone --recurse-submodules https://github.com/csh-apprentice/Adviser_CS.git
    cd Adviser_CS/ICEPACK
    chmod +x run_casestudy.sh
    ./run_casestudy.sh
  '

```

## Pine Island Basin Scale Model

We follow the instructions in [icesheetModels](https://github.com/fastice/icesheetModels), but fixing some code issues caused by the dependcies version mismatch. 

*Running the Demonstration*
```
./adviser run \
  --region us-west-2 \
  --instance-type c6i.4xlarge \
  --container-image-uri docker.io/firedrakeproject/firedrake-vanilla:2025-01 \
  -- \
  bash -lc '
    set -euxo pipefail
    git clone --recurse-submodules https://github.com/csh-apprentice/Adviser_CS.git
    cd Adviser_CS/ICEPACK
    chmod +x run_casestudy.sh
    ./run_casestudy.sh
  '

```