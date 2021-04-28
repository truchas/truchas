.. _SOLID_MECHANICS_Namelist:

.. toctree::
   :maxdepth: 1

SOLID_MECHANICS Namelist
============================

Overview
----------
The :ref:`SOLID_MECHANICS<SOLID_MECHANICS_Namelist>` namelist sets parameters that are specific to the solid mechanics model and algorithm. This namelist is read whenever the :ref:`PHYSICS<PHYSICS_Namelist>` namelist option :ref:`Solid_Mechanics<PHYSICS_SM>` is enabled. Solid mechanics boundary conditions are defined using :ref:`SOLID_MECHANICS_BC<SM_BC_Namelist>` namelists.

SOLID_MECHANICS Namelist Features
-----------------------------------
| **Required/Optional        :** Required when solid mechanics physics is enabled.
| **Single/Multiple Instances:** Single

Components
------------

Physics Options
^^^^^^^^^^^^^^^^^^
* :ref:`contact_distance <SM_CD>`
* :ref:`contact_norm_trac <SM_CNT>`

Numerical Parameters
^^^^^^^^^^^^^^^^^^^^^^^
* :ref:`maximum_iterations <SM_MI>`
* :ref:`nlk_vector_tolerance <SM_NVT>`
* :ref:`nlk_max_vectors <SM_NMV>`
* :ref:`abs_stress_tol <SM_AST>`
* :ref:`abs_displ_tol <SM_ADT>`
* :ref:`rel_displ_tol <SM_RDT>`
* :ref:`nlk_tol <SM_NT>` (expert)
* :ref:`relaxation_parameter <SM_RP>` (expert)
* :ref:`preconditioning_steps <SM_PS>` (expert)
* :ref:`contact_penalty <SM_CP>` (expert)

.. _SM_CD:

contact_distance
^^^^^^^^^^^^^^^^^^^

| **Description** : A length scale parameter :math:`\beta` for the contact function

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
     1 & if &s \le 0 \\
     0 & if &s \ge \beta \\
     2*(\frac{s}{\beta} - 1)^3 + 3(\frac{s}{\beta} - 1)^2 ,& if & 0 \le s \le \beta
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
     1 & if &\tau_n \le 0 \\
     0 & if &\tau_n \ge \tau^{*} \\
     2*(\frac{\tau_n}{\tau^{*}} - 1)^3 + 3(\frac{\tau_n}{\tau^{*}} - 1)^2 & if & 0 \le \tau_n \le \tau^{*}
   \end{array}
   \right.
   \]

| **Physical dimension**: :math:`L`
| **Type**        : real
| **Default**     : 1.0e-7
| **Valid Values**: (0, :math:`\infty`]
| **Notes**       : The default value is usually a good value for mesh cell sizes in the 1 - 10 mm size range.

.. _SM_CNT:

contact_norm_trac
^^^^^^^^^^^^^^^^^^^

| **Description** : A parameter :math:`\tau^{*}` for the contact function

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
     1 & if &s \le 0 \\
     0 & if &s \ge \beta \\
     2*(\frac{s}{\beta} - 1)^3 + 3(\frac{s}{\beta} - 1)^2 ,& if & 0 \le s \le \beta
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
     1 & if &\tau_n \le 0 \\
     0 & if &\tau_n \ge \tau^{*} \\
     2*(\frac{\tau_n}{\tau^{*}} - 1)^3 + 3(\frac{\tau_n}{\tau^{*}} - 1)^2 & if & 0 \le \tau_n \le \tau^{*}
   \end{array}
   \right.
   \]

:math:`\tau_n` is the normal traction at the interface where a positive value corresponds to a tensile force normal to the surface.

| **Physical dimension**: :math:`F/L^2`
| **Type**        : real
| **Default**     : 1.0e4
| **Valid Values**: [0, :math:`\infty`]
| **Notes**       : The default value is probably appropriate for materials with elastic constants in the range :math:`10^9 - 10^{11}`. This parameter should probably be scaled proportionately for elastic constants that differ from this range.

.. _SM_MI:

maximum_iterations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : Maximum allowed number of iterations of the nonlinear solver.
| **Type**         : integer
| **Default**      : :math:`100`
| **Valid Values** : :math:`[0,\infty)`

.. _SM_NVT:

nlk_vector_tolerance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : The vector drop tolerance for the NLK method. When assembling the acceleration subspace vector by vector, a vector is dropped when the sine of the angle between the vector and the subspace less than this value.
| **Type**         : real
| **Default**      : :math:`0.01`
| **Valid Values** : :math:`(0,1)`

.. _SM_NMV:

nlk_max_vectors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : For the NLK method, the maximum number of acceleration vectors to be used.
| **Type**         : integer
| **Default**      : :math:`20`
| **Valid Values** : :math:`[0,\infty)`

.. _SM_AST:

abs_stress_tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The tolerance for the absolute error of the residual used by the solid mechanics solver.

| **Physical Dimension**: same as the ‘stress’ variable
| **Type**        : real
| **Default**     : :math:`1e-10`
| **Valid Values**: :math:`\gt 0`

.. _SM_ADT:

abs_displ_tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The tolerance :math:`\epsilon`  for the absolute error component of the temperature error norm used by the nonlinear solver. If :math:`\delta u` is a displacement field increment with reference displacement field :math:`u`,then this error norm is

.. math::
   |||\delta u||| \equiv \mathop{{max}_j} |\delta u_j|/(\epsilon + \eta |u_j|)

The relative error tolerance :math:`\eta` is given by :ref:`rel_displ_tol<SM_RDT>`.

| **Physical Dimension**: :math:`\Theta`
| **Type**        : real
| **Default**     : :math:`1e-10`
| **Valid Values**: :math:`\geq 0`
| **Notes** : The error norm is dimensionless and normalized.

For :math:`u_j` sufficiently small the norm approximates an absolute norm with tolerance  :math:`\epsilon`  , and for :math:`u_j` sufficiently large the norm approximates a relative norm with tolerance :math:`\eta`. If :math:`\epsilon = 0`  then the norm is a pure relative norm and the displacement must be bounded away from 0.

.. _SM_RDT:

rel_displ_tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The tolerance  :math:`\eta` for the relative error component of the displacement error norm used by the nonlinear solver. If :math:`\delta u` is a displacement field increment with reference displacement field :math:`u`, then this error norm is

.. math::
   |||\delta T||| \equiv \mathop{{max}_j} |\delta T_j|/(\epsilon + \eta |T_j|)

The absolute error tolerance :math:`\epsilon` is given by :ref:`abs_temp_tol<SM_ADT>`.

| **Physcial Dimension**: dimensionless
| **Type**        : real
| **Default**     : :math:`1e-10`
| **Valid Values**: (0, 1)
| **Notes** : See the notes for :ref:`abs_displ_tol<SM_ADT>`.

.. _SM_NT:

nlk_tol (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The convergence tolerance for the NLK nonlinear solver. The nonlinear system is considered solved by the current iterate if the norm of the last solution correction is less than this value.
| **Type**        : real
| **Default**     : 1.0
| **Valid Values**: (0, 1]
| **Notes** : This tolerance is relative to the dimensionless and normalized BDF2 integrator norm; see :ref:`abs_displ_tol<SM_ADT>`, for example. The nonlinear system only needs to be solved to an accuracy equal to the acceptable local truncation error for the step, which is roughly 1. Solving to a greater accuracy is wasted effort.

.. _SM_RP:

relaxation_parameter (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The relaxation parameter for the diagonal scaling preconditioner.
| **Type**        : real
| **Default**     : 1.0
| **Valid Values**: (0, 1]

.. _SM_PS:

preconditioning_steps (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Number of iterations for the preconditioner. In the current version, a diagonal scaling preconditioner is used, such that there is no sense in changing this value unless the :ref:`relaxation_parameter<SM_RP>` variable is also changed from the default.
| **Type**        : integer
| **Default**     : 1
| **Valid Values**: :math:`\geq 1`

.. _SM_CP:

contact_penalty (expert)
^^^^^^^^^^^^^^^^^^^

| **Description** : A penalty factor for the penetration constraint in the contact algorithm. Changing this is probably not a good idea in the current version.
| **Physical Dimension**: dimensionless
| **Type**         : real
| **Default**      : :math:`1e3`
| **Valid Values** : [0, :math:`\infty`)
