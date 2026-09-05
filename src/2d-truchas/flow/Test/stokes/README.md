# Stokes tests

These integrated tests run the `truchas-2d` executable with `--simulation flow`
and inertia disabled:

- pressure-driven channel flow, aligned and rotated, compared with the
  analytic Poiseuille profile;
- two-fluid Couette flow with a mesh-aligned interface and contrasting
  viscosities, compared with the analytic piecewise-linear profile;
- the same channel test with mesh perturbations;
- body-force balance and a closed-box hydrostatic state;
- `lid_driven_cavity.json`, a developer example that is not registered with
  CTest.

The Python tests run the executable in a temporary directory and inspect its
VTKHDF output.
