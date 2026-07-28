#!/usr/bin/env bash
set -euo pipefail

truchas_dir=$(dirname "$(readlink -f "$0")")

# Keep the documented Intel MPI regular-mode default while allowing a site to
# select a different fabric combination before invoking the package.
: "${I_MPI_FABRICS:=shm:ofi}"
export I_MPI_FABRICS

# libfabric loads providers dynamically. Prefer bundled providers while
# retaining a site-supplied search path for its fabric.
export FI_PROVIDER_PATH="${truchas_dir}/../lib/prov${FI_PROVIDER_PATH:+:${FI_PROVIDER_PATH}}"

ulimit -s unlimited
exec "${truchas_dir}/t-linux.x86_64.gnu" "$@"
