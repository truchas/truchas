.. _SOLID_MECHANICS_Namelist:

.. toctree::
   :maxdepth: 1

SOLID_MECHANICS Namelist
============================

Overview
----------
The :ref:`SOLID_MECHANICS<SOLID_MECHANICS_Namelist>` namelist sets parameters that are specific to the solid mechanics model and algorithm. This namelist is read whenever the :ref:`PHYSICS<PHYSICS_Namelist>` namelist option :ref:`Solid_Mechanics<PHYSICS_SM>` is enabled. Parameters for the nonlinear solver used by the algorithm and its preconditioner are specified in :ref:`NONLINEAR_SOLVER<NONLINEAR_SOLVER_Namelist>` and :ref:`LINEAR_SOLVER<LINEAR_SOLVER_Namelist>` namelists. Optional material viscoplasticity models are defined in :ref:`VISCOPLASTIC_MODEL<VISCOPLASTIC_MODEL_Namelist>` namelists.

SOLID_MECHANICS Namelist Features
-----------------------------------
| **Required/Optional        :** Required when solid mechanics physics is enabled.
| **Single/Multiple Instances:** Single

Components
------------
* :ref:`Contact_Distance <SM_CD>`
* :ref:`Contact_Norm_Trac <SM_CNT>`
* :ref:`Contact_Penalty <SM_CP>`
* :ref:`Displacement_Nonlinear_Solution <SM_DNS>`
* :ref:`Solid_Mechanics_Body_Force <SM_SMBF>`
* :ref:`Stress_Reduced_Integration <SM_SRI>`
* :ref:`Strain_Limit <SM_SL>`
* :ref:`Convergence_Criterion <SM_CC>`
* :ref:`Maximum_Iteration <SM_MI>`
* :ref:`NLK_Vector_Tolerance <SM_NVT>`
* :ref:`NLK_Max_Vectors <SM_NMV>`

.. _SM_CD:

Contact_Distance
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

Contact_Norm_Trac
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

.. _SM_CP:

Contact_Penalty
^^^^^^^^^^^^^^^^^^^

| **Description** : A penalty factor for the penetration constraint in the contact algorithm. Changing this is probably not a good idea in the current version.
| **Physical Dimension**: dimensionless
| **Type**         : real
| **Default**      : :math:`1.0e3`
| **Valid Values** : [0, :math:`\infty`]

.. _SM_DNS:

Displacement_Nonlinear_Solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : A character string pointer to the nonlinear solution algorithm parameters to be used in a Newton-Krylov solution of the nonlinear thermo-elastic viscoplastic equations. This string “points” to a particular :ref:`NONLINEAR_SOLVER<NONLINEAR_SOLVER_Namelist>` namelist if it matches the :ref:`Name<NL_Name>` input variable string in the :ref:`NONLINEAR_SOLVER<NONLINEAR_SOLVER_Namelist>` namelist.
| **Type**         : string
| **Default**      : "default"
| **Valid Values** : arbitrary string
| **Notes**        : If this string does not match a :ref:`Name<NL_Name>` input variable string specified in a :ref:`NONLINEAR_SOLVER<<NONLINEAR_SOLVER_Namelist>` namelist, then the default set of nonlinear solution algorithm parameters is used for the thermo-elastic viscoplastic equations.

.. _SM_SMBF:

Solid_Mechanics_Body_Force
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : Body forces will be included in the solid mechanics calculation.
| **Type**         : logical
| **Default**      : .false.

.. _SM_SL:

Strain_Limit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : This parameter controls the use of the ODE integrator in the plastic strain calculation. It should be set to the minimum significant value of the plastic strain increment for a time step. If convergence seems poor when a viscoplastic material model is used, it may help to reduce the value.
| **Type**         : real
| **Default**      : :math:`1.0e-10`
| **Valid Values** : :math:`\geq 0`
| **Notes**        : This parameter can not be currently used to control the time step. It may be used for such purposes in future releases.

.. _SM_CC:

Convergence_Criterion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : Tolerance used to determine when nonlinear convergence has been reached. 
| **Type**         : real
| **Default**      : :math:`1.0e-12`
| **Valid Values** : (0,0.1)
| **Note**         : We refer to the input value of Convergence_Criterion as:math:`\epsilon`. :math:`F(x) = 0` is the nonlinear system being solved.

The nonlinear iteration is stopped when either of the following two conditions are met:

* reduction in the 2-norm of the nonlinear residual meets the criterion, i.e.:

:math:`\frac{||F(x_{k+1})||_2}{||F(x_0)||_2} \lt \gamma` 

* The relative change in the max-norm of the solution meets the criterion, i.e.:

:math:`\frac{||\gamma_x||_{\infty}}{||x_{k+1}||_{\infty}} \lt \gamma` 

where :math:`\gamma` is the input desired tolerance, modified using an estimate of the convergence rate, i.e.:

:math:`\gamma = (1-\rho)\epsilon`

and:

:math:`\frac{||x_{k+1} - x_k||_{\infty} / ||x_{k+1}||_{\infty}}{||x_{k} - x_{k-1}||_{\infty} / ||x_{k}||_{\infty}}`

This is an attempt to prevent false convergence if the solution stagnates, but allow iteration to stop if the solution is acceptable.

.. _SM_MI:

Maximum_Iteration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : Maximum allowed number of iterations of the nonlinear solver.
| **Type**         : integer
| **Default**      : :math:`100`
| **Valid Values** : :math:`[0,\infty)`

.. _SM_NVT:

NLK_Vector_Tolerance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : The vector drop tolerance for the NLK method. When assembling the acceleration subspace vector by vector, a vector is dropped when the sine of the angle between the vector and the subspace less than this value.
| **Type**         : real
| **Default**      : :math:`0.01`
| **Valid Values** : :math:`(0,1)`

.. _SM_NMV:

NLK_Vector_Tolerance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : For the NLK method, the maximum number of acceleration vectors to be used.
| **Type**         : integer
| **Default**      : :math:`20`
| **Valid Values** : :math:`[0,\infty)`