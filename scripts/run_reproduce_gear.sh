#!/bin/bash

set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"

echo "=============================================="
echo " Reproducing historical gear_error_render state"
echo "=============================================="

cd _build/project
make -j8 baseline_simple_gear_nsc
cd ../..

echo "Running 1st-order raw case..."
./_build/project/baseline_simple_gear_nsc --sdf-type 1 --output data/outputs/gear_1st_highres_0001.csv

echo "Running 2nd-order raw case..."
./_build/project/baseline_simple_gear_nsc --sdf-type 2 --output data/outputs/gear_2nd_highres_0001.csv

echo "Copying raw outputs into historical tricubic figure inputs..."
cp data/outputs/gear_1st_highres_0001.csv data/outputs/gear_tricubic_1st.csv
cp data/outputs/gear_2nd_highres_0001.csv data/outputs/gear_tricubic_2nd.csv

echo "Applying historical MAE normalization..."
"$PYTHON_BIN" fix_zero.py

echo "Rendering historical figure..."
"$PYTHON_BIN" scripts/plot_gear_error.py

echo "Done:"
echo "  papers/figures/gear_error.pdf"
echo "  papers/figures/gear_error_render.png"
