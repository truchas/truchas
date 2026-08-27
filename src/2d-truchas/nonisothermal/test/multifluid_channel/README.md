# Multifluid channel tests

These tests exercise coupled material transport, enthalpy transport, and
incompressible flow for two liquid materials in a translating channel.  They
use BDF1 and compare each output against an analytic solution rather than
against stored field data.

- `input_horizontal.json`: stationary horizontal layers with different
  material densities.  It checks that a discontinuous material composition can
  remain at rest without generating spurious flow.
- `input_vertical.json`: a horizontal interface translating vertically with a
  uniform velocity.  It checks interface transport, temperature, enthalpy, and
  velocity for unequal material densities.
- `input_density.json`: a vertical interface translating horizontally with a
  uniform velocity.  It is a focused regression for the conservative
  variable-density momentum update: old momentum must use the old cell density,
  while the predictor uses the new density.
