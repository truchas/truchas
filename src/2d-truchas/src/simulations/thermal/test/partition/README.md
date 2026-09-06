# VTKHDF partition test

This test uses a deliberately simple problem: a uniform material, a linear
initial temperature, and matching Dirichlet conditions on every boundary. It
is intended to isolate output and ordering behavior rather than test thermal
physics. The nonconstant temperature field is important: a constant field
would not reveal a cell-data permutation error.

The test runs the same input with one and four MPI processes, reads both
VTKHDF files through the Python reader, and compares:

- external cell IDs and external node IDs;
- cell-center coordinates;
- final temperature and enthalpy values.

It therefore exercises the VTKHDF global-ID and ghost metadata and the reader's
reconstruction of external order. It also verifies that the reported solution
is independent of the MPI partition. It is not an analytic heat-transport
test; those are kept under `../linear/`.
