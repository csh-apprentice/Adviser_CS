#!/usr/bin/env bash
set -euo pipefail




cd $HOME/pism-stable/examples/std-greenland
export PATH="$HOME/pism/bin:$PATH"
./preprocess.sh