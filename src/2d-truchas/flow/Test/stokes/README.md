# Stokes tests

These integrated tests run the `flow_2d` executable:

- pressure-driven channel flow, aligned and rotated, compared with the
  analytic Poiseuille profile;
- the same channel test with mesh perturbations;
- body-force balance and a closed-box hydrostatic state;
- `lid_driven_cavity.json`, a developer example that is not registered with
  CTest.

The Python tests run the executable in a temporary directory and inspect its
VTKHDF output.
