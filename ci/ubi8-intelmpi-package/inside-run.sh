#!/usr/bin/env bash
set -euo pipefail

usage() {
  printf 'Usage: %s {package|test}\n' "${0##*/}" >&2
}

command=${1:-}
case "${command}" in
  package|test) ;;
  *) usage; exit 2 ;;
esac

jobs=${JOBS:-$(nproc)}
work=${HOME}/truchas-${command}
rm -rf "${work}"
mkdir -p "${work}"
cp -a /source/. "${work}/"
cd "${work}"

if [[ "${command}" == package ]]; then
  version=${TRUCHAS_VERSION:?TRUCHAS_VERSION must be set}
  rm -f .git
  git init --quiet
  git config user.name 'UBI package build'
  git config user.email 'package-build@localhost'
  git add --all
  git commit --quiet -m 'Package source snapshot'
  git tag "${version}"
  printf '%s\n' "${version}" > version
fi

cmake -S . -B build -C config/linux-intel.cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DTRUCHAS_TPL_DIR="${TRUCHAS_TPL_DIR}" \
  -DCMAKE_INSTALL_PREFIX=inst \
  -DBLA_VENDOR=OpenBLAS \
  -DUSE_PORTAGE=No
cmake --build build --parallel "${jobs}"

if [[ "${command}" == test ]]; then
  ctest --test-dir build --repeat until-pass:2 --output-on-failure 2>&1 | tee /output/ctest.log
  exit 0
fi

cmake --install build
cd inst/bin
# The current master distribution script and launcher use the GNU binary name.
ln -sf t-linux.x86_64.intelllvm t-linux.x86_64.gnu
TRUCHAS_SOURCE_DIR="${work}" "${TRUCHAS_TPL_DIR}/python-install/bin/python" /tools/make_dist.py
cp truchas-*.tar.bz2 /output/
cp ../../version /output/
cp /opt/truchas-package/versions /output/toolchain-versions.txt
cp /opt/truchas-package/toolchain-info.txt /output/
