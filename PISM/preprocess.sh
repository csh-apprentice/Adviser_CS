#!/usr/bin/env bash
set -euo pipefail


## First Run


cd $HOME/pism-stable/examples/std-greenland
export PATH="$HOME/pism/bin:$PATH"
./preprocess.sh