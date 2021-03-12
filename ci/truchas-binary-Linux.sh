#!/bin/bash

# Exit on any errors
set -e

# Determine the directory that this script is installed in (not the directory it
# is executed from), which also contains the Truchas binary
TRUCHAS_DIR="$(dirname "$(readlink -f "$0")")"

# Set unlimited stack size, which is necessary for Intel Fortran compiled
# binaries, otherwise they will exceed the stack size and segfault for moderate
# to big problems
ulimit -s unlimited

# Launch the Truchas binary itself, forwarding all arguments
${TRUCHAS_DIR}/t-linux.x86_64.intel $@
