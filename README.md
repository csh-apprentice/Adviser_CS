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
  -- \
  bash -lc '
    set -euxo pipefail
    git clone https://github.com/csh-apprentice/Adviser_CS.git
    cd Adviser_CS/ICEPACK
    chmod +x run_casestudy.sh
    ./run_casestudy.sh
  '

```