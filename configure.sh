#!/bin/bash

set -e

cd "$(dirname "$0")"
mkdir -p build
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
