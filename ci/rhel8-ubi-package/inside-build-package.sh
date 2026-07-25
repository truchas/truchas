#!/usr/bin/env bash
set -euo pipefail

jobs=${JOBS:-32}

set +u
source /opt/intel/bin/compilervars.sh intel64
set -u
export FC=ifort CC=icc CXX=icpc
export PATH="${HOME}/ext/cmake/bin:${HOME}/ext/bin:${HOME}/ext/python-install/bin:${PATH}"
export LD_LIBRARY_PATH="${HOME}/ext/lib:${LD_LIBRARY_PATH:-}"

work=${HOME}/truchas
rm -rf "${work}"
mkdir -p "${work}"
cp -a /source/. "${work}/"
cd "${work}"

# The mounted source can be a linked Git worktree whose Git metadata is not
# available in the container.  Recreate only the version metadata locally.
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
cmake \
  -C ../config/linux-intel.cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DTRUCHAS_TPL_DIR="${HOME}/ext" \
  -DCMAKE_INSTALL_PREFIX=inst \
  -DUSE_PORTAGE=No \
  ..
make -j"${jobs}"
make install
cd inst/bin
"${HOME}/ext/python-install/bin/python" ../../../ci/make_dist.py
cp truchas-*.tar.bz2 /output/
cp ../../../version /output/
