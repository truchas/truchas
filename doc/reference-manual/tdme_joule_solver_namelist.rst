.. _TDME_JOULE_SOLVER_Namelist:

.. toctree::
   :maxdepth: 1

TDME_JOULE_SOLVER Namelist
===========================
This namelist specifies the parameters that pertain to the time-domain
Joule heat solver, which integrates the periodically-forced Maxwell equations
in time until a periodic steady state solution is attained, and then computes
the time-averaged Joule heat over one cycle of the forcing. This solver is used
by the induction heating physics model; see the
:ref:`INDUCTION_HEATING<INDUCTION_HEATING_Namelist>` namelist.

Namelist Variables
------------------

.. contents::
   :local:

steps_per_cycle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The number of time steps per cycle of the periodic forcing used to integrate
the EM field equations.

:Type: integer
:Default: 20
:Valid Values: > 0
:Notes: Increasing the number of time steps per cycle increases the accuracy
   and stability of the Joule heat computation, while generally increasing
   the execution time. A reasonable range of values is [10,40]; anything
   less than 10 is severely discouraged.


steady_state_tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Convergence to the periodic steady state is measured by comparing the
computed Joule heat field averaged over the last source field cycle,
:math:`q_\text{last}`, with the result from the previous cycle,
:math:`q_\text{prev}`. The Joule heat computation is considered converged
when :math:`\lVert{q_\text{last}-q_\text{prev}}\rVert_\infty < \delta\,
\lVert{q_\text{last}}\rVert_\infty` where :math:`\delta` is the value of
this parameter.

:Type: real
:Default: :math:`10^{-2}`
:Valid Values: > 0
:Notes: Depending on the accuracy of the other physics, :math:`10^{-2}` or
        :math:`10^{-3}` are adequate values. The measured error is more
        properly the error in :math:`q_\text{prev}`; the actual error in
        :math:`q_\text{last}` will be significantly less.


max_source_cycles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Specifies the time limit, measured in cycles of the periodic forcing, allowed
for convergence to the periodic steady state. If convergence is not attained
within this limit, the last result is returned and a warning issued, but
execution continues.

:Type: integer
:Default: 5
:Valid Values: > 0
:Notes: To avoid ringing, the amplitude of the external source field is phased
        in and is not at full strength until after approximately two cycles
        have passed. Consequently this input variable should normally be at
        least 3. Convergence to a periodic steady state is commonly attained
        within 5 cycles; see `steady_state_tol`_.


rel_tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The convergence tolerance :math:`\eta` for the conjugate gradient (CG) solution
of the linear time step system. The CG iteration is considered converged when
the residual :math:`r` satisfies :math:`\lVert{r}\rVert<\eta\lVert{r_0}\rVert`.

:Type: real
:Default: :math:`10^{-8}`
:Valid Values: :math:`(0, 0.1)`
:Notes: The numerical characteristics of the EM system require that the linear
        systems be solved to significantly greater accuracy than might
        otherwise be required. Too loose a tolerance will manifest itself in
        a significant build-up of noise in the solution of the electric field
        over the course of the simulation.


max_iter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The maximum number of conjugate gradient (CG) iterations allowed. It is a
fatal error if convergence is not attained within this number of iterations.

:Type: integer
:Default: 500
:Valid Values: > 0


c_ratio
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In the quasi-magnetostatic regime, Maxwell's equations can become become very
ill-conditioned and difficult to solve in regions with zero electrical
conductivity, such as free space. Artificially reducing the speed of light
in those regions helps ameliorate that difficulty, and this numerical parameter
is that reduction factor :math:`r`. A value of 1 means no reduction and is the
default.

:Type: real
:Default: 1
:Valid Values: :math:`\in(0,1]`
:Notes: This is implemented by altering the permittivity,
        :math:`\epsilon\to\epsilon/r^2`,
        which only modifies the displacement current term
        :math:`\partial\epsilon\vec{E}/\partial t`.
        In the magnetostatic regime, this term functions as a very small
        perturbation in Maxwell's equations.

        For common induction heating applications, this parameter can be taken
        fairly small, :math:`10^{-3}` or :math:`10^{-4}`, without any signficant
        effect on the computed Joule heat.

        More generally, :math:`\omega\epsilon\ll\sigma` in the magnetostatic
        regime, where :math:`\omega` is the frequency of the external driving
        field, and one wants to ensure that this separation of scales is
        maintained when choosing the reduction factor,
        :math:`\omega\epsilon/r^2\ll\sigma`.


graphics_output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A flag that enables output of the EM solution for graphics visualization.

:Type: logical
:Default: ``.false.``
:Notes: When enabled, the solver generates its own graphics data file in
   VTKHDF format, which is compatible with the ParaView visualization software.
   The solver produces a separate sequentially-numbered file for each use,
   identifiable by the ``.vtkhdf`` suffix. This setting does not affect the
   standard graphics output generated by Truchas, which is configured
   elsewhere. Since these files are simply HDF5 files structured to the
   VTKHDF specification, they can also be used with generic HDF5 tools.


print_level
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Controls the verbosity of the time-domain Joule heat solver.

:Type: integer
:Default: 1
:Valid Values: 1, 2, or 3
:Notes: At the default level, 1, a status message is output at the end of
  each source field cycle showing the progress toward steady state. Level 2
  adds a summary of the CG iteration for each time step. Level 3 adds
  convergence info for each CG iterate. Levels 1 and 2 are typical.


relax_type (experimental)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Specifies the type of relaxation used within the linear system preconditioner.

:Type: integer
:Default: 0
:Valid values: 0, Gauss-Seidel; 1, Hypre Boomer AMG
:Notes: Boomer is significantly more expensive without any obvious benefit.
