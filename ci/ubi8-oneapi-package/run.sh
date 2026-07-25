#!/usr/bin/env bash
set -euo pipefail

usage() {
  printf 'Usage: %s {image|package|test|smoke}\n' "${0##*/}" >&2
}

command=${1:-}
case "${command}" in
  image|package|test|smoke) ;;
  *) usage; exit 2 ;;
esac

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_dir=${SOURCE_DIR:-$(cd "${script_dir}/../.." && pwd)}
image=${IMAGE_TAG:-truchas-ubi8-oneapi-build:2026.0}
jobs=${JOBS:-$(nproc)}

if [[ "${command}" == image ]]; then
  podman build --build-arg JOBS="${jobs}" --tag "${image}" \
    --file "${script_dir}/Containerfile" "${script_dir}"
  exit 0
fi

source_snapshot=$(mktemp -d)
trap 'rm -rf "${source_snapshot}"' EXIT
git -C "${source_dir}" archive --format=tar HEAD | tar -xf - -C "${source_snapshot}"

case "${command}" in
  package)
    output_dir=${OUTPUT_DIR:-${PWD}/ubi8-oneapi-package-output}
    version=$(git -C "${source_dir}" describe --tags --dirty)
    mkdir -p "${output_dir}"
    podman run --rm --user 0 --env HOME=/home/swuser --env JOBS="${jobs}" \
      --env TRUCHAS_VERSION="${version}" \
      --volume "${script_dir}:/tools:ro,Z" \
      --volume "${source_snapshot}:/source:ro,Z" \
      --volume "${output_dir}:/output:Z" \
      "${image}" bash -o pipefail -c \
      'bash /tools/inside-run.sh package 2>&1 | tee /output/package-build.log'
    ;;
  test)
    output_dir=${OUTPUT_DIR:-${PWD}/ubi8-oneapi-ctest-output}
    mkdir -p "${output_dir}"
    podman run --rm --user 0 --env HOME=/home/swuser --env JOBS="${jobs}" \
      --volume "${script_dir}:/tools:ro,Z" \
      --volume "${source_snapshot}:/source:ro,Z" \
      --volume "${output_dir}:/output:Z" \
      "${image}" bash -o pipefail -c \
      'bash /tools/inside-run.sh test 2>&1 | tee /output/ctest-build.log'
    ;;
  smoke)
    package_dir=${PACKAGE_DIR:-${PWD}/ubi8-oneapi-package-output}
    test -f "${package_dir}"/truchas-*.tar.bz2
    podman run --rm \
      --volume "${source_snapshot}:/source:ro,Z" \
      --volume "${package_dir}:/package:ro,Z" \
      registry.access.redhat.com/ubi8/ubi:8.10 bash -lc \
      'dnf -y install bzip2 && mkdir -p /work/dist && cp /package/truchas-*.tar.bz2 /work/dist/ && cd /work && /source/ci/test_tarball.sh'
    ;;
esac
