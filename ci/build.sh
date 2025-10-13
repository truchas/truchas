#!/bin/bash

set -ex

build_type="$1"

if [[ $CC == "icx" ]]; then
    compiler="intel"
else
    compiler="gcc"
fi

cmake --version
cmake -S. -Bbuild \
    -C config/linux-${compiler}.cmake \
    -DCMAKE_BUILD_TYPE=${build_type} \
    -DTRUCHAS_TPL_DIR=$HOME/ext \
    -DCMAKE_INSTALL_PREFIX=inst \
    -DUSE_PORTAGE=Yes
cd build
make -j8
make install
