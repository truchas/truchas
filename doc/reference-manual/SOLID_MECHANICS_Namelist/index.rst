SOLID_MECHANICS Namelist
========================

The ``SOLID_MECHANICS`` namelist sets parameters that are specific to the solid
mechanics model and algorithm. This namelist is read whenever the :ref:`PHYSICS
Namelist<PHYSICS_Namelist/index:PHYSICS Namelist>` option
:ref:`solid_mechanics<PHYSICS_Namelist/index:solid_mechanics>` is enabled. Solid
mechanics boundary conditions are defined using :ref:`SOLID_MECHANICS_BC
Namelists<SOLID_MECHANICS_BC_Namelist/index:SOLID_MECHANICS_BC Namelist>`.

:Required/Optional: Required when solid mechanics physics is enabled.
:Single/Multiple Instances: Single

.. contents:: Components
   :local:

Physics Options
---------------

contact_distance
^^^^^^^^^^^^^^^^

A length scale parameter :math:`\beta` for the contact function

:math:`\lambda = \lambda_s * \lambda_{\tau}`

.. math::
   :label: poly_lambda_s
   :nowrap:

   \[
   \lambda_s =
   \left\{
   \begin{array}{
     @{}% no padding
     l@{\quad}% some padding
     r@{}% no padding
     >{{}}r@{}% no padding
     >{{}}l@{}% no padding
   }
     1 & \mathrm{if} &s \le 0 \\
     0 & \mathrm{if} &s \ge \beta \\
     2*(\frac{s}{\beta} - 1)^3 + 3(\frac{s}{\beta} - 1)^2 ,& \mathrm{if} & 0 \le s \le \beta
   \end{array}
   \right.
   \]

and

.. math::
   :label: poly_lambda_tau
   :nowrap:

   \[
   \lambda_{\tau} =
   \left\{
   \begin{array}{
     @{}% no padding
     l@{\quad}% some padding
     r@{}% no padding
     >{{}}r@{}% no padding
     >{{}}l@{}% no padding
   }
     1 & \mathrm{if} &\tau_n \le 0 \\
     0 & \mathrm{if} &\tau_n \ge \tau^{*} \\
     2*(\frac{\tau_n}{\tau^{*}} - 1)^3 + 3(\frac{\tau_n}{\tau^{*}} - 1)^2 & \mathrm{if} & 0 \le \tau_n \le \tau^{*}
   \end{array}
   \right.
   \]

:Physical dimension: :math:`L`
:Type: real
:Default: 1.0e-7
:Valid Values: (0, :math:`\infty`]

.. note::
   The default value is usually a good value for mesh cell sizes in the 1 - 10
   mm size range.


contact_norm_trac
^^^^^^^^^^^^^^^^^

A parameter :math:`\tau^{*}` for the contact function

:math:`\lambda = \lambda_s * \lambda_{\tau}`

where

.. math::
   :label: poly_lambda_s1
   :nowrap:

   \[
   \lambda_s =
   \left\{
   \begin{array}{
     @{}% no padding
     l@{\quad}% some padding
     r@{}% no padding
     >{{}}r@{}% no padding
     >{{}}l@{}% no padding
   }
     1 & \mathrm{if} &s \le 0 \\
     0 & \mathrm{if} &s \ge \beta \\
     2*(\frac{s}{\beta} - 1)^3 + 3(\frac{s}{\beta} - 1)^2 ,& \mathrm{if} & 0 \le s \le \beta
   \end{array}
   \right.
   \]

and

.. math::
   :label: poly_lambda_tau1
   :nowrap:

   \[
   \lambda_{\tau} =
   \left\{
   \begin{array}{
     @{}% no padding
     l@{\quad}% some padding
     r@{}% no padding
     >{{}}r@{}% no padding
     >{{}}l@{}% no padding
   }
     1 & \mathrm{if} &\tau_n \le 0 \\
     0 & \mathrm{if} &\tau_n \ge \tau^{*} \\
     2*(\frac{\tau_n}{\tau^{*}} - 1)^3 + 3(\frac{\tau_n}{\tau^{*}} - 1)^2 & \mathrm{if} & 0 \le \tau_n \le \tau^{*}
   \end{array}
   \right.
   \]

:math:`\tau_n` is the normal traction at the interface where a positive value
corresponds to a tensile force normal to the surface.

:Physical dimension: :math:`F/L^2`
:Type: real
:Default: 1.0e4
:Valid Values: [0, :math:`\infty`]

.. note::
   The default value is probably appropriate for materials with elastic
   constants in the range :math:`10^9 - 10^{11}`. This parameter should probably
   be scaled proportionately for elastic constants that differ from this range.

Numerical Parameters
--------------------

maximum_iterations
^^^^^^^^^^^^^^^^^^

Maximum allowed number of iterations of the nonlinear solver.

:Type: integer
:Default: 500
:Valid Values: :math:`[0,\infty)`


maximum_outer_iterations
^^^^^^^^^^^^^^^^^^^^^^^^

Maximum allowed number of iterations around the nonlinear solver if contact is
present. Contact is a nonlinear effect, and convergence depends strongly on a
good preconditioner. It may be necessary to restart the solver a number of times
with a recomputed preconditioner for convergence. When contact is not present,
this is always set to 1.

:Type: integer
:Default: 5
:Valid Values: :math:`[0,\infty)`


nlk_vector_tolerance
^^^^^^^^^^^^^^^^^^^^

The vector drop tolerance for the NLK method. When assembling the acceleration
subspace vector by vector, a vector is dropped when the sine of the angle
between the vector and the subspace less than this value.

:Type: real
:Default: 0.01
:Valid Values: :math:`(0,1)`


nlk_max_vectors
^^^^^^^^^^^^^^^

For the NLK method, the maximum number of acceleration vectors to be used.

:Type: integer
:Default: 20
:Valid Values: :math:`[0,\infty)`


rel_stress_tol
^^^^^^^^^^^^^^

The tolerance for the relative error of the residual used by the solid mechanics
solver.

:Physical Dimension: same as the ‘stress’ variable
:Type: real
:Default: 1e-10
:Valid Values: :math:`\gt 0`


pc_abs_lame_tol
^^^^^^^^^^^^^^^^

The tolerance for the absolute change in the Lame parameters, above which a

The tolerance :math:`\epsilon` for the absolute error component of the
displacement error norm used by the nonlinear solver. If :math:`\delta u` is a
displacement field increment with reference displacement field :math:`u`,then
this error norm is

.. math::
   |||\delta u||| \equiv \mathop{{max}_j} |\delta u_j|/(\epsilon + \eta |u_j|)

The relative error tolerance :math:`\eta` is given by `rel_displ_tol`_.

:Physical Dimension: :math:`\Theta`
:Type: real
:Default: 1e-10
:Valid Values: :math:`\geq 0`

.. note::
   The error norm is dimensionless and normalized.

.. note::
   For :math:`u_j` sufficiently small the norm approximates an absolute norm
   with tolerance :math:`\epsilon`, and for :math:`u_j` sufficiently large the
   norm approximates a relative norm with tolerance :math:`\eta`. If
   :math:`\epsilon = 0` then the norm is a pure relative norm and the
   displacement must be bounded away from 0.


pc_rel_lame_tol
^^^^^^^^^^^^^^^^

The tolerance :math:`\eta` for the relative error component of the displacement
error norm used by the nonlinear solver. If :math:`\delta u` is a displacement
field increment with reference displacement field :math:`u`, then this error
norm is

.. math::
   |||\delta u||| \equiv \mathop{{max}_j} |\delta u_j|/(\epsilon + \eta |u_j|)

The absolute error tolerance :math:`\epsilon` is given by `abs_displ_tol`_.

:Physcial Dimension: dimensionless
:Type: real
:Default: 1e-10
:Valid Values: (0, 1)

.. note::
   See the notes for `abs_displ_tol`_.


abs_displ_tol (expert)
^^^^^^^^^^^^^^^^^^^^^^

The tolerance :math:`\epsilon` for the absolute error component of the
displacement error norm used by the nonlinear solver. If :math:`\delta u` is a
displacement field increment with reference displacement field :math:`u`,then
this error norm is

.. math::
   |||\delta u||| \equiv \mathop{{max}_j} |\delta u_j|/(\epsilon + \eta |u_j|)

The relative error tolerance :math:`\eta` is given by `rel_displ_tol`_.

:Physical Dimension: :math:`\Theta`
:Type: real
:Default: 1e100
:Valid Values: :math:`\geq 0`

.. note::
   The error norm is dimensionless and normalized.

.. note::
   For :math:`u_j` sufficiently small the norm approximates an absolute norm
   with tolerance :math:`\epsilon`, and for :math:`u_j` sufficiently large the
   norm approximates a relative norm with tolerance :math:`\eta`. If
   :math:`\epsilon = 0` then the norm is a pure relative norm and the
   displacement must be bounded away from 0.


rel_displ_tol (expert)
^^^^^^^^^^^^^^^^^^^^^^

The tolerance :math:`\eta` for the relative error component of the displacement
error norm used by the nonlinear solver. If :math:`\delta u` is a displacement
field increment with reference displacement field :math:`u`, then this error
norm is

.. math::
   |||\delta u||| \equiv \mathop{{max}_j} |\delta u_j|/(\epsilon + \eta |u_j|)

The absolute error tolerance :math:`\epsilon` is given by `abs_displ_tol`_.

:Physcial Dimension: dimensionless
:Type: real
:Default: 1e100
:Valid Values: :math:`\geq 0`

.. note::
   See the notes for `abs_displ_tol`_.


nlk_tol (expert)
^^^^^^^^^^^^^^^^

The convergence tolerance for the NLK nonlinear solver. The nonlinear system is
considered solved by the current iterate if the norm of the last solution
correction is less than this value.

:Type: real
:Default: 1.0
:Valid Values: (0, 1]

.. note::
   This tolerance is relative to the dimensionless and normalized BDF2
   integrator norm; see `abs_displ_tol`_, for example. The nonlinear system only
   needs to be solved to an accuracy equal to the acceptable local truncation
   error for the step, which is roughly 1. Solving to a greater accuracy is
   wasted effort.


preconditioner_method (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The preconditioning method to use. The current default is `"boomeramg"`, which
performs well and generally shouldn't need to be changed. Other options include
`"ssor"`, and `"ds"` for diagonal scaling (Jacobi).

:Type: string
:Default: `"boomeramg"`
:Valid Values: `"boomeramg"`, `"ssor"`, or `"ds"`

relaxation_parameter (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The relaxation parameter for the diagonal scaling preconditioner.

:Type: real
:Default: 1.0
:Valid Values: (0, 1]


stress_relaxation_parameter (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The relaxation parameter for the stress part of the diagonal scaling
preconditioner. For reference, the original solid mechanics solver used a value of 16 / 9.

:Type: real
:Default: 1.0
:Valid Values: :math:`\gt 0`


preconditioning_steps (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Number of iterations for the preconditioner. In the current version, a diagonal
scaling preconditioner is used, such that there is no sense in changing this
value unless the `relaxation_parameter (expert)`_ variable is also changed from the
default.

:Type: integer
:Default: 2 when `preconditioner_method = "boomeramg"` (default) or `preconditioner_method = "ssor"`, and 1 when `preconditioner_method = "ds"`.
:Valid Values: :math:`\geq 1`


contact_penalty (expert)
^^^^^^^^^^^^^^^^^^^^^^^^

A penalty factor for the penetration constraint in the contact algorithm.
Changing this is probably not a good idea in the current version.

:Physical Dimension: dimensionless
:Type: real
:Default: 1e3
:Valid Values: [0, :math:`\infty`)
