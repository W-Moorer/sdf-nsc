#!/bin/bash
set -e

echo "=== Phase 1: Bootstrapping Chrono (Download) ==="
bash scripts/bootstrap_chrono.sh

echo ""
echo "=== Phase 1.5: Bootstrapping OpenVDB/NanoVDB ==="
bash scripts/bootstrap_vdb.sh

echo ""
echo "=== Phase 2: Building Chrono (Minimal Setup) ==="
bash scripts/build_chrono_minimal.sh

echo ""
echo "=== Phase 3: Building Main Project ==="
bash scripts/build_project.sh

echo ""
echo "=== All baseline bootstrap phases completed. ==="
