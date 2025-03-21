#!/bin/bash

mkdir build
cd build

cmake -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=${PREFIX} \
    -DCMAKE_C_COMPILER=clang \
    -DCMAKE_CXX_COMPILER=clang++ \
    -DCMAKE_CXX_FLAGS=-Wunused-command-line-argument \
    -DRUN_CHECKS=ON \
    -DPHYSICA_MKL=ON \
    -DPHYSICA_HDF5=ON \
    -DPHYSICA_CUDA=ON ..
cmake --build .

cd test
ctest --parallel
cd ..

cmake --install .
