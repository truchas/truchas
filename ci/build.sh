#!/bin/bash

set -ex

if [[ $1 == "Debug" ]]; then
    build_type="Debug"
else
    build_type="Release"
fi

# Set Intel compilers env variables
export INTELDIR=/opt/intel/oneapi/2024.0
export PATH=$INTELDIR/bin/:$PATH
export LD_LIBRARY_PATH=$INTELDIR/lib
export FC=ifx
export CC=icx
export CXX=icpx

# Add openmpi executables into path
export PATH=$HOME/ext/bin/:$PATH
export LD_LIBRARY_PATH=$HOME/ext/lib:$LD_LIBRARY_PATH

# Prepend our `python` executable into path
export PATH="${HOME}/ext/python-install/bin:$PATH"


# Install Truchas
mkdir build
cd build
cmake --version
cmake \
    -C ../config/linux-intel-llvm.cmake \
    -DCMAKE_BUILD_TYPE=${build_type} \
    -DTRUCHAS_TPL_DIR=$HOME/ext \
    -DCMAKE_INSTALL_PREFIX=inst \
    -DUSE_PORTAGE=Yes \
    ..
make -j8
make install
