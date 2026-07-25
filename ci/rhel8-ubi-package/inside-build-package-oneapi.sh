#!/usr/bin/env bash
set -euo pipefail

jobs=${JOBS:-32}
export PATH=/opt/intel/oneapi/compiler/2026.0/bin:${HOME}/ext/cmake/bin:${HOME}/ext/bin:${HOME}/ext/python-install/bin:${PATH}
export LD_LIBRARY_PATH=/opt/intel/oneapi/compiler/2026.0/lib:/opt/intel/oneapi/compiler/2026.0/opt/compiler/lib:${HOME}/ext/lib:${LD_LIBRARY_PATH:-}
export FC=ifx CC=icx CXX=icpx

work=${HOME}/truchas
rm -rf "${work}"
mkdir -p "${work}"
cp -a /source/. "${work}/"
cd "${work}"

version=${TRUCHAS_VERSION:?TRUCHAS_VERSION must be set}
rm -f .git
git init --quiet
git config user.name "RHEL package build"
git config user.email "package-build@localhost"
git add --all
git commit --quiet -m "Package source snapshot"
git tag "${version}"
printf '%s\n' "${version}" > version

mkdir build
cd build
cmake -C ../config/linux-intel.cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=icx \
  -DCMAKE_CXX_COMPILER=icpx \
  -DCMAKE_Fortran_COMPILER=ifx \
  -DTRUCHAS_TPL_DIR="${HOME}/ext" \
  -DCMAKE_INSTALL_PREFIX=inst \
  -DBLA_VENDOR=OpenBLAS \
  -DUSE_PORTAGE=No \
  ..
make -j"${jobs}"

# Keep the complete suite result with the package artifacts.  Some known
# oneAPI/UBI8 numerical differences are accepted, so do not block packaging.
set +e
# Run the complete suite serially.  MPICH unit tests can hang when scheduled
# concurrently with the broader MPI workload.
ctest --test-dir . --output-on-failure 2>&1 | tee /output/ctest.log
ctest_status=${PIPESTATUS[0]}
set -e

if (( ctest_status == 0 )); then
  printf 'All CTest tests passed.\n' > /output/ctest-summary.txt
else
  {
    printf 'CTest exited with status %d. Failed tests:\n' "${ctest_status}"
    sed -n -E 's/^[[:space:]]*[0-9]+ - (.*) \\(Failed\\)$/  - \\1/p' /output/ctest.log
  } | tee /output/ctest-summary.txt
fi

make install
cd inst/bin
ln -sf t-linux.x86_64.intelllvm t-linux.x86_64.intel
"${HOME}/ext/python-install/bin/python" ../../../ci/make_dist.py
cp truchas-*.tar.bz2 /output/
cp ../../../version /output/
