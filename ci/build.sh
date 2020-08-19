#!/bin/bash

set -ex

if [[ $1 == "Debug" ]]; then
    build_type="Debug"
else
    build_type="Release"
fi

# Set Intel compilers env variables
export INTELDIR=/opt/intel/
export PATH="$INTELDIR/bin/:$PATH"
export LD_LIBRARY_PATH=$INTELDIR/lib/intel64
export FC=ifort
export CC=icc
export CXX=icpc

# Add openmpi executables into path
export PATH=$HOME/ext/bin/:$PATH
export LD_LIBRARY_PATH=$HOME/ext/lib:$LD_LIBRARY_PATH

sudo rm -f /usr/bin/python
sudo ln -s ${HOME}/ext/python-install/bin/python /usr/bin/python
sudo ln -s ${HOME}/ext/python-install/bin/python3 /usr/bin/python3


# Install Truchas
mkdir build
cd build
cmake \
    -C ../config/linux-intel.cmake \
    -DCMAKE_BUILD_TYPE=${build_type} \
    -DTRUCHAS_TPL_DIR=$HOME/ext \
    -DCMAKE_INSTALL_PREFIX=inst \
    ..
make -j8
make install
