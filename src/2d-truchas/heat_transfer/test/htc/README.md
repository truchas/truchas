# Transient HTC cooling test

This test exercises an external heat-transfer-coefficient boundary condition
in a transient single-material calculation. It uses an 8-by-8 quad mesh on a
unit square with density, specific heat, and conductivity all equal to one.
The initial temperature is uniformly one, and all four sides use `h = 1` with
ambient temperature zero.

The Python checker reads the final VTKHDF output and verifies that every cell
has cooled below its initial temperature while remaining positive. It also
checks that enthalpy equals temperature for the unit density and specific-heat
values. The case is run in both serial and four-process configurations.

This is a transient response smoke test. The rotated linear HTC cases in
`../linear` provide the analytic accuracy test for HTC fluxes and spatially
varying ambient temperatures.
