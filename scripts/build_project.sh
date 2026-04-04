#!/bin/bash
set -e

ROOT_DIR="$(pwd)"
PROJECT_BUILD_DIR="$ROOT_DIR/_build/project"
CHRONO_INSTALL_DIR="$ROOT_DIR/_deps/chrono-install"

echo "Configuring the main project..."
mkdir -p "$PROJECT_BUILD_DIR"
cd "$PROJECT_BUILD_DIR"

# Pass the Chrono installation path so find_package(Chrono) can locate it.
cmake ../.. \
    -DCMAKE_BUILD_TYPE=Release \
    -DChrono_DIR="$CHRONO_INSTALL_DIR/lib/cmake/Chrono" \
    -DSPCC_ENABLE_VDB=ON \
    -DNANOVDB_ROOT="$ROOT_DIR/third_party/openvdb/nanovdb" \
    -DOPENVDB_INCLUDE_DIR="/usr/include" \
    -DOPENVDB_LIBRARY="/usr/lib/x86_64-linux-gnu/libopenvdb.so"

echo "Building the main project..."
CORES=4
if command -v nproc >/dev/null 2>&1; then CORES=$(nproc); fi
cmake --build . --config Release -j $CORES

echo "Project build completed successfully."
