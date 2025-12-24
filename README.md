# ICEPACK & PISM CASESTUDY
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

*Running the Demonstration [INVERSE]*
```
./adviser run \
  --region us-west-2 \
  --instance-type c6i.4xlarge \
  --num-nodes 4 \
  --container-image-uri docker.io/firedrakeproject/firedrake-vanilla:2025-01 \
  -- \
  bash -lc '
    set -euxo pipefail
    git clone --recurse-submodules https://github.com/csh-apprentice/Adviser_CS.git
    cd Adviser_CS/ICEPACK
    chmod +x run_inverse.sh
    ./run_inverse.sh
    chmod +x run_forward.sh
    ./run_forward.sh
  '

```

*Running the Demonstration [FORWARD]*
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
    chmod +x run_inverse.sh
    ./run_inverse.sh
    chmod +x run_forward.sh
    ./run_forward.sh
  '

```

*Running the Demonstration [INVERSE+FORWARD]*
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
    chmod +x run_full.sh
    ./run_full.sh
  '

```




# PISM
We follow the instructions in [PISM first run](https://www.pism.io/docs/manual/std-greenland/run-1.html).
```
./adviser run \
  --region us-west-2 \
  --num-nodes 4 \
  -- \
  bash -lc '
    git clone --recurse-submodules https://github.com/csh-apprentice/Adviser_CS.git
    cd Adviser_CS/PISM
    chmod +x setup_pism.sh
    ./setup_pism.sh
    chmod +x first_run_new.sh
    ./first_run_new.sh
  '

```


# BENCHMARK STUDY
We also studies how different instance type affects the solver speed, and use the icepack forward problem as a benchmark problem.

## Intel (baseline)
```
c6i.2xlarge   (8 vCPU)
```
## AMD (comparison)
```
c7a.2xlarge   (8 vCPU) 0.41
```
## Graviton (ARM)
```
c7g.2xlarge   (8 vCPU) 0.29$ 
```
## Intel “new gen”
```
c7i.2xlarge (8 vCPU)
```
## Useful “cost/perf” Intel variant: C7i-flex
```
c7i-flex.2xlarge
```

*Running the Benchmark [1 warm up and 100 experiements]*

```
./adviser run \
  --container-image-uri docker.io/firedrakeproject/firedrake-vanilla-default:2025.10.2 \
  --instance-type c7a.2xlarge \
  "
    set -euo pipefail
    rm -rf Adviser_CS
    git clone --recurse-submodules https://github.com/csh-apprentice/Adviser_CS.git
    cd Adviser_CS/ICEPACK
    chmod +x run_benchmark.sh
    bash run_benchmark.sh 100 2000
  "


./adviser run \
   --container-image-uri docker.io/firedrakeproject/firedrake:2025.10.2 \
  --instance-type c7g.2xlarge \
  "
    set -euo pipefail
    rm -rf Adviser_CS
    git clone --recurse-submodules https://github.com/csh-apprentice/Adviser_CS.git
    cd Adviser_CS/ICEPACK
    chmod +x run_benchmark.sh
    bash run_benchmark.sh 100 2000
  "

```
