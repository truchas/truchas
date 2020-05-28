#!/bin/bash

set -ex

export INTELDIR=/opt/intel/
export PATH=$INTELDIR/bin/:$PATH
export LD_LIBRARY_PATH=$INTELDIR/lib/intel64
export FC=ifort
export CC=icc
export CXX=icpc

# Add openmpi executables into path
export PATH=$HOME/ext/bin/:$PATH
export LD_LIBRARY_PATH=$HOME/ext/lib:$LD_LIBRARY_PATH

mkdir build
cd build
cmake \
    -C ../config/intel-opt.cmake \
    -DTRUCHAS_TPL_DIR=$HOME/ext \
    -DCMAKE_INSTALL_PREFIX=inst \
    ..
make -j8
make install

cd inst/bin/
../../../ci/make_dist.py
cd ../..

ctest --output-on-failure
