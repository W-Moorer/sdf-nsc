#!/bin/bash
set -e

# Define directories
ROOT_DIR="$(pwd)"
CHRONO_SRC_DIR="$ROOT_DIR/third_party/chrono"
CHRONO_BUILD_DIR="$ROOT_DIR/_build/chrono"
CHRONO_INSTALL_DIR="$ROOT_DIR/_deps/chrono-install"

if [ ! -d "$CHRONO_SRC_DIR" ]; then
    echo "Error: third_party/chrono does not exist. Please run bootstrap_chrono.sh first."
    exit 1
fi

echo "Configuring minimal Chrono build (Out-of-source)..."
mkdir -p "$CHRONO_BUILD_DIR"
cd "$CHRONO_BUILD_DIR"

# TODO: If any of these CMake flags break due to upstream Chrono updates, 
# please adjust the ENABLE_MODULE_* flags accordingly.
cmake "$CHRONO_SRC_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$CHRONO_INSTALL_DIR" \
    -DENABLE_MODULE_VEHICLE=OFF \
    -DENABLE_MODULE_SENSOR=OFF \
    -DENABLE_MODULE_FSI=OFF \
    -DENABLE_MODULE_ROS=OFF \
    -DENABLE_MODULE_PYTHON=OFF \
    -DENABLE_MODULE_MATLAB=OFF \
    -DENABLE_MODULE_OPENGL=OFF \
    -DENABLE_MODULE_IRRLICHT=OFF \
    -DENABLE_MODULE_VSG=OFF \
    -DENABLE_MODULE_GPU=OFF \
    -DENABLE_MODULE_POSTPROCESS=OFF

echo "Building Chrono (Minimal)..."
# Using cross-platform nproc, defaulting to 4 if unavailable
CORES=4
if command -v nproc >/dev/null 2>&1; then CORES=$(nproc); fi
cmake --build . --config Release -j $CORES

echo "Installing Chrono locally to $CHRONO_INSTALL_DIR..."
cmake --install . --config Release

echo "Chrono minimal build and install completed successfully."
