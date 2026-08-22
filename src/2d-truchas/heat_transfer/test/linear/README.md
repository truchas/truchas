# Linear conduction tests

This directory contains the same linear conduction problem on two element
types and in two orientations:

- `input_quad_x.json` and `test_quad_x.py` impose `T = 1 - x` on quads;
- `input_quad_y.json` and `test_quad_y.py` impose `T = 1 - y` on quads;
- `input_tri_x.json` and `test_tri_x.py` impose `T = 1 - x` on triangles;
- `input_tri_y.json` and `test_tri_y.py` impose `T = 1 - y` on triangles.

The `*_xy` cases use both coordinates:

- `input_quad_xy.json` and `test_quad_xy.py` use quads;
- `input_tri_xy.json` and `test_tri_xy.py` use triangles;
- `input_mixed_xy.json` and `test_mixed_xy.py` use mixed quads and triangles.

They impose and initialize `T = 1 + x + 2*y` and check both temperature and
enthalpy (`H = 2*T`).

Two additional quad cases exercise geometry independently of topology:

- `input_quad_xy_rotated.json` rotates an orthogonal mesh by 17 degrees;
- `input_quad_xy_noisy.json` perturbs the mesh with `noise-factor: 0.1`.

These use analytic Dirichlet conditions on all four sides, since constant
normal flux values would not remain consistent with the linear solution after
rotating or perturbing the boundary geometry.

The analytic checkers compare the final cell temperatures with the appropriate
cell-center profile. Each case is run in both serial and four-process
configurations. All analytic checks use a `1e-10` tolerance. The triangle and
mixed inputs use tighter solver and initial-condition tolerances to achieve
this accuracy.
