# Transient radiation tests

The `transient-kelvin` and `transient-celsius` cases use the annular-sector
mesh and output at times `1e-4` and `1e-3`. Both tests run with four processes
against a serial Kelvin reference. The Kelvin case compares both non-initial
states with reference data. The Celsius case checks that its final temperature
and enthalpy fields, shifted by 273.15, agree with the Kelvin reference. The
comparison tolerance is `1e-4`, allowing the expected small
partition-dependent differences while remaining tighter than the input
absolute integration tolerance of `1e-3`. The absolute integration tolerances
make the Kelvin/Celsius comparison invariant under the temperature-scale
shift.
