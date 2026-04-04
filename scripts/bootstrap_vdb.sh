#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== VDB Bootstrap: installing OpenVDB system package ==="
if command -v apt-get >/dev/null 2>&1; then
    sudo apt-get update
    sudo apt-get install -y libopenvdb-dev libtbb-dev
else
    echo "apt-get not found. Please install OpenVDB manually for your platform."
    exit 1
fi

echo "=== VDB Bootstrap: fetching NanoVDB headers ==="
mkdir -p "$ROOT_DIR/third_party"

if [ -f "$ROOT_DIR/third_party/openvdb/nanovdb/nanovdb/NanoVDB.h" ]; then
    echo "NanoVDB headers already present at third_party/openvdb/nanovdb"
else
    if [ ! -d "$ROOT_DIR/third_party/openvdb/.git" ]; then
        git clone --depth 1 --filter=blob:none --sparse \
            https://github.com/AcademySoftwareFoundation/openvdb.git \
            "$ROOT_DIR/third_party/openvdb"
        cd "$ROOT_DIR/third_party/openvdb"
        git sparse-checkout set nanovdb LICENSE
    else
        cd "$ROOT_DIR/third_party/openvdb"
        git sparse-checkout set nanovdb LICENSE || true
        git pull --ff-only || true
    fi

    if [ ! -f "$ROOT_DIR/third_party/openvdb/nanovdb/nanovdb/NanoVDB.h" ]; then
        echo "Failed to fetch NanoVDB headers."
        exit 1
    fi
fi

echo "=== VDB Bootstrap completed ==="
echo "OpenVDB header : /usr/include/openvdb/openvdb.h"
echo "NanoVDB header : $ROOT_DIR/third_party/openvdb/nanovdb/nanovdb/NanoVDB.h"
