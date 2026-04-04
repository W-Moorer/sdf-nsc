#!/bin/bash
set -e

echo "Checking prerequisites..."
command -v git >/dev/null 2>&1 || { echo >&2 "Git is required but it's not installed. Aborting."; exit 1; }
command -v cmake >/dev/null 2>&1 || { echo >&2 "CMake is required but it's not installed. Aborting."; exit 1; }

# Basic compiler check
if ! command -v c++ >/dev/null 2>&1 && ! command -v g++ >/dev/null 2>&1 && ! command -v clang++ >/dev/null 2>&1; then
    echo "Warning: No unified C++ compiler command found. Ensure MSVC or relevant compiler is in PATH."
fi

mkdir -p third_party

if [ -d "third_party/chrono/.git" ]; then
    echo "Project Chrono is already cloned in third_party/chrono. Skipping download."
else
    echo "Cloning official Project Chrono repository..."
    git clone -b main https://github.com/projectchrono/chrono.git third_party/chrono
    echo "Chrono downloaded successfully."
fi
