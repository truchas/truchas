These are some auxiliary results from the drag-1d test problem. The *.agr files
are save files from the 2D plotting program xmgrace.

* u-vel-probe is a plot of x-velocity at the mid-channel probe location using
  a time step of 1/16 for the 3 standard choices of drag_implicitness (0, 1/2,
  and 1) compared with the exact solution.

* convergence is a plot of the errors at t=1 for the 3 standard choices of
  drag_implicitness and a sequence of time step sizes. This shows the 1st
  order convergence of forward and backward Euler and 2nd order convergence
  of trapezoid.
