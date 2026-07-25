#!/usr/bin/env bash
set -euo pipefail

prefix=${HOME}/ext
jobs=${JOBS:-32}
state=${prefix}/.rhel8-24.04-state
work=${prefix}/sources
cmake_dir=${prefix}/cmake
marker=${prefix}/.truchas-rhel8-24.04-dependencies

set +u
source /opt/intel/bin/compilervars.sh intel64
set -u
export FC=ifort CC=icc CXX=icpc

mkdir -p "${state}" "${work}"

if [[ ! -f ${state}/mpich ]]; then
  if [[ -x ${prefix}/bin/mpifort && -x ${prefix}/bin/mpiexec.hydra && -f ${prefix}/lib/libmpi.so ]]; then
    touch "${state}/mpich"
  else
    cd "${work}"
    if [[ ! -d mpich-3.3.2 ]]; then
      curl -fsSLO https://www.mpich.org/static/downloads/3.3.2/mpich-3.3.2.tar.gz
      echo '4bfaf8837a54771d3e4922c84071ef80ffebddbb6971a006038d91ee7ef959b9  mpich-3.3.2.tar.gz' | sha256sum --check
      tar -xzf mpich-3.3.2.tar.gz
    fi
    cd mpich-3.3.2
    [[ -f Makefile ]] || ./configure --prefix="${prefix}"
    make -j"${jobs}"
    make install
    touch "${state}/mpich"
  fi
fi

if [[ ! -f ${state}/cmake ]]; then
  mkdir -p "${cmake_dir}"
  if [[ ! -x ${cmake_dir}/bin/cmake ]]; then
    cd "${work}"
    if [[ ! -f cmake-3.22.2-linux-x86_64.tar.gz ]]; then
      curl -fsSLO https://github.com/Kitware/CMake/releases/download/v3.22.2/cmake-3.22.2-linux-x86_64.tar.gz
    fi
    echo '38b3befdee8fd2bac06954e2a77cb3072e6833c69d8cc013c0a3b26f1cfdfe37  cmake-3.22.2-linux-x86_64.tar.gz' | sha256sum --check
    rm -rf "${cmake_dir}"
    mkdir -p "${cmake_dir}"
    tar -xzf cmake-3.22.2-linux-x86_64.tar.gz --strip-components=1 -C "${cmake_dir}"
  fi
  touch "${state}/cmake"
fi
export PATH="${cmake_dir}/bin:${prefix}/bin:${PATH}"
export LD_LIBRARY_PATH="${prefix}/lib:${LD_LIBRARY_PATH:-}"

if [[ ! -f ${state}/python ]]; then
  cd "${work}"
  if [[ ! -d python-cmake-buildsystem/.git ]]; then
    rm -rf python-cmake-buildsystem
    git clone https://github.com/python-cmake-buildsystem/python-cmake-buildsystem.git
  fi
  cd python-cmake-buildsystem
  git checkout 11a369e29b1be82e3b3108259ab89b7095854b47
  # CPython 3.9's CMake wrapper crashes during sysconfig generation with icc.
  # Python is a packaged runtime dependency; the Fortran/C++ toolchain remains Intel.
  if [[ -f build/CMakeCache.txt ]] && ! grep -q '^CMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc$' build/CMakeCache.txt; then
    rm -rf build
  fi
  cmake -S . -B build -DCMAKE_C_COMPILER=gcc -DCMAKE_INSTALL_PREFIX:PATH="${prefix}/python-install"
  cmake --build build -j"${jobs}"
  cmake --install build
  "${prefix}/python-install/bin/python" -m ensurepip
  "${prefix}/python-install/bin/python" -m pip install --no-cache-dir numpy==1.26.4 h5py==3.10.0
  touch "${state}/python"
fi

if [[ ! -f ${state}/tpl ]]; then
  cd "${work}"
  if [[ ! -d truchas-tpl/.git ]]; then
    rm -rf truchas-tpl
    git clone https://gitlab.com/truchas/truchas-tpl.git
  fi
  cd truchas-tpl
  git checkout v22
  cmake -S . -B build -C config/linux-intel.cmake \
    -DCMAKE_INSTALL_PREFIX="${prefix}" \
    -DSEARCH_FOR_HDF5=no -DSEARCH_FOR_NETCDF=no -DBUILD_PORTAGE=OFF
  cmake --build build -j"${jobs}"
  touch "${state}/tpl"
fi

cat > "${marker}" <<EOF
MPICH=3.3.2
BLAS_LAPACK=Intel-MKL-2020.4.304
CMAKE=3.22.2
PYTHON_CMAKE_BUILDSYSTEM=11a369e29b1be82e3b3108259ab89b7095854b47
TRUCHAS_TPL=v22
EOF
