# Truchas Binary Tarball (Linux)

This tarball contains pre-compiled Truchas and all libraries it depends on that
should work on any Linux distribution.

Unpack the tarball using:

    tar xjf truchas-*-Linux.tar.bz2
    cd truchas-*-Linux

## Example: Broken Dam

You can run the example as follows in serial:

    cd examples/broken-dam
    ../../bin/truchas broken-dam.inp

or parallel (4 cores):

    cd examples/broken-dam
    ../../bin/mpiexec -n 4 ../../bin/truchas broken-dam.inp
