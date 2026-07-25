#!/usr/bin/env bash
set -euo pipefail

jobs=${JOBS:-$(nproc)}
prefix=${HOME}/ext
work=${HOME}/build

set +u
source /opt/intel/oneapi/setvars.sh --force >/dev/null
set -u
# This image deliberately uses OpenBLAS.  MUMPS treats a defined MKLROOT as a
# request to configure against MKL, even when the MKL package is not installed.
unset MKLROOT
export PATH="${prefix}/cmake/bin:${prefix}/bin:${prefix}/python-install/bin:${PATH}"
export LD_LIBRARY_PATH="${prefix}/lib:${prefix}/lib64:${LD_LIBRARY_PATH:-}"
export CC=icx CXX=icpx FC=ifx

mkdir -p "${work}" "${prefix}"
cd "${work}"

curl -fsSLO https://github.com/Kitware/CMake/releases/download/v3.30.5/cmake-3.30.5-linux-x86_64.tar.gz
mkdir -p "${prefix}/cmake"
tar -xzf cmake-3.30.5-linux-x86_64.tar.gz --strip-components=1 -C "${prefix}/cmake"

git clone --depth=1 https://github.com/python-cmake-buildsystem/python-cmake-buildsystem.git python-cmake-buildsystem
cmake -S python-cmake-buildsystem -B python-build \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_INSTALL_PREFIX="${prefix}/python-install"
cmake --build python-build --parallel "${jobs}"
cmake --install python-build
"${prefix}/python-install/bin/python" -m ensurepip
"${prefix}/python-install/bin/python" -m pip install --no-cache-dir numpy h5py scipy matplotlib fypp

curl -fsSLO https://www.mpich.org/static/downloads/5.0.1/mpich-5.0.1.tar.gz
tar -xzf mpich-5.0.1.tar.gz
cd mpich-5.0.1
./configure --prefix="${prefix}" --with-pm=hydra
make -j"${jobs}"
make install
cd "${work}"

curl -fsSLO https://github.com/OpenMathLib/OpenBLAS/archive/refs/tags/v0.3.26.tar.gz
tar -xzf v0.3.26.tar.gz
cd OpenBLAS-0.3.26
make -j"${jobs}" MAKE_NB_JOBS="${jobs}" CC=gcc FC=gfortran DYNAMIC_ARCH=1 USE_OPENMP=0
make PREFIX="${prefix}" install
cd "${work}"

git clone https://gitlab.com/truchas/truchas-tpl.git
cd truchas-tpl
git checkout fb5ef5e38
git apply /usr/local/share/truchas-package/mumps-openblas.patch
git apply /usr/local/share/truchas-package/fvtkhdf-prefix.patch
cmake -S . -B build -C config/linux-intel.cmake \
  -DCMAKE_INSTALL_PREFIX="${prefix}" \
  -DSEARCH_FOR_HDF5=no \
  -DSEARCH_FOR_NETCDF=no \
  -DBUILD_PORTAGE=OFF
if ! cmake --build build --parallel "${jobs}"; then
  find build/mumps/src/mumps-stamp -type f -name 'mumps-configure-*.log' -print -exec cat {} \; 2>/dev/null || true
  exit 1
fi

rm -rf "${work}"
