# Freezing flow

This is a two-dimensional adaptation of the truchas-pbf
`test/freezing-flow` problem. An initially superheated phase-change material
moves through an inviscid channel. Cooling at the upper wall advances a
horizontal solidification front downward, reducing the available liquid
channel width while the flow remains approximately plug-like.

The test compares the initial state and outputs at t=0.5 and t=3 with the
four-process NAG reference in `reference/out.vtkhdf`. It checks phase
fractions, temperature, enthalpy, and accepted-step counts against the
reference. On active cells, pressure and velocity are checked against the
analytic solution p=0, u=(1,0). Pressure and velocity are expected to be
NaN in solid cells.
