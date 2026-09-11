# He et al. single-mode Rayleigh–Taylor benchmark

This example is the two-dimensional single-mode Rayleigh–Taylor instability
benchmark introduced by He et al. (1999). It uses the nondimensional
parameters

- Atwood number `At = 0.5`;
- Reynolds number `Re = 256`;
- a unit-width, height-four domain;
- a 128-by-512 quadrilateral mesh.

The light fluid initially lies below the heavy fluid. The interface is

```text
y = 2 + 0.1 cos(2 pi x)
```

and gravity acts downward. The side and bottom walls are free-slip, while the
top boundary has prescribed pressure zero. The simulation integrates to
`t=4.5` and writes 18 equally spaced output intervals. The resulting rising
bubbles and falling spikes provide a visual and qualitative benchmark for
multifluid free-surface flow.

## Interface function

The active input uses a built-in polynomial that closely approximates the
cosine interface. The exact interface function is retained in
`rt-interface.f90` and can be dynamically loaded instead.

From this directory, compile the shared library with GNU Fortran:

```sh
gfortran -O2 -fPIC -shared -o libregions.so rt-interface.f90
```

The exported function is `rt_interface`, and the commented input block already
contains the corresponding `library-path`, `library-symbol`, and parameter
values. Uncomment that block and comment out the polynomial block in
`input.json` to use the exact cosine interface.

Because the loader interprets a library name without a directory component
relative to the current working directory, run the case from this example
directory when using `"library-path": "libregions.so"`:

```sh
cd examples/flow/rayleigh-taylor
mpiexec -n 16 truchas-2d --simulation flow \
  --output-dir output --force input.json
```
