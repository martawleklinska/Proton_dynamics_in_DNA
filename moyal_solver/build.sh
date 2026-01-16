#!/bin/bash

echo "=== Building Moyal Solver ==="

if [ ! -d "build" ]; then
    mkdir build
    echo "Created build directory"
fi

cd build

echo "Configuring with CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
echo "Building..."
make -j$(nproc)

if [ $? -eq 0 ]; then
    echo "Build successful!"
    echo "Executable: ./main"
else
    echo "Build failed!"
    exit 1
fi
