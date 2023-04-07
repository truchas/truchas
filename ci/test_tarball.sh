#!/bin/bash

# This tests that the binary tarball contains all the files that we require and
# that all executables / Python scripts work as expected. The overall idea is
# that unless a feature from the tarball is tested in this script, we can assume
# it is broken otherwise. When you add a new feature to the tarball, please add
# a test here for it.

set -ex

tar xaf dist/truchas-*.tar.bz2
cd truchas-*/
cat README.md
./bin/t-linux.x86_64.gnu -h
./bin/truchas -h
./bin/python -c "import numpy, h5py; print(h5py.__file__)"
./bin/python bin/write-restart.py -h
cd examples/broken-dam
../../bin/truchas broken-dam.inp
../../bin/python ../../bin/write-restart.py broken-dam_output/broken-dam.h5
sed -i s/0.075/0.005/ broken-dam.inp
../../bin/mpiexec -n 2 ../../bin/truchas -f broken-dam.inp
