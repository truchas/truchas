.. _DIFFUSION_Solver_Namelist:

.. toctree::
   :maxdepth: 1

DIFFUSION Solver Namelist
==========================

Overview
------------
The DIFFUSION_SOLVER namelist sets the parameters that are specific to the
heat and species transport solver. The namelist is read when either of the
:ref:`PHYSICS<PHYSICS_Namelist>` namelist options
:ref:`heat_transport<physics-ht>` or :ref:`species_transport<physics-st>`
are enabled.

The solver has two time integration methods which are selected by the variable
:ref:`Stepping_Method<DIFF_SOL_SM>`. The default is a variable step-size,
implicit second-order BDF2 method that controls the local truncation error of
each step to a user-defined tolerance by adaptively adjusting the step size.
The step size is chosen so that an a priori estimate of the error will be
within tolerance, and steps are rejected when the actual erroris too large.
A failed step may be retried with successively smaller step sizes.

The other integration method is a non-adaptive, implicit first-order BDF1
method specifically designed to handle the exceptional difficulties that
arise when heat transfer is coupled to a fluid flow system that includes
void. In this context the heat transfer domain changes from one step to the
next because of the moving void region, and mesh cells may only be partially
filled with material. For this method the timestep is controlled by flow or
other physics models.

Both methods share a common nonlinear solver and preconditioning options.

The initial step size and upper and lower bounds for the step size are set
in the :ref:`NUMERICS<NUMERICS_Namelist>` namelist. In addition, the step
size selected by the adaptive solver may be further limited by other physics
models or by the :ref:`NUMERICS<NUMERICS_Namelist>` variables
:ref:`Dt_Grow<NUMERICS_DTG>` and :ref:`Dt_Constant<NUMERICS_DTC>`. When only
diffusion solver physics are enabled, it isimportant that these variables be
set appropriately so as not to unnecessarily impede the normal functioning of
the diffusion solver.

.. admonition:: Namelist Usage

   :Required/Optional: Required when :ref:`heat_transport<physics-ht>` and/or
      :ref:`species_transport<physics-st>` physics is enabled.
   :Single/Multiple Instances: Single

Components
------------
* :ref:`Abs_Conc_Tol<DIFF_SOL_ACT>`
* :ref:`Abs_Enthalpy_Tol<DIFF_SOL_AET>`
* :ref:`Abs_Temp_Tol<DIFF_SOL_ATT>`
* :ref:`Cond_Vfrac_Threshold<DIFF_SOL_CVT>`
* :ref:`Cutvof<DIFF_SOL_CV>`
* :ref:`Hypre_AMG_Debug<DIFF_SOL_HAD>`
* :ref:`Hypre_AMG_Logging_Level<DIFF_SOL_HALL>`
* :ref:`Hypre_AMG_Print_Level<DIFF_SOL_HAPL>`
* :ref:`Max_NLK_Itr<DIFF_SOL_MNI>`
* :ref:`Max_NLK_Vec<DIFF_SOL_MNV>`
* :ref:`Max_Step_Tries<DIFF_SOL_MST>`
* :ref:`NLK_Preconditioner<DIFF_SOL_NP>`
* :ref:`NLK_Tol<DIFF_SOL_NT>`
* :ref:`NLK_Vec_Tol<DIFF_SOL_NVT>`
* :ref:`PC_AMG_Cycles<DIFF_SOL_PAC>`
* :ref:`PC_Freq<DIFF_SOL_PF>`
* :ref:`PC_SSOR_Relax<DIFF_SOL_PSR>`
* :ref:`PC_SSOR_Sweeps<DIFF_SOL_PSS>`
* :ref:`Rel_Conc_Tol<DIFF_SOL_RCT>`
* :ref:`Rel_Enthalpy_Tol<DIFF_SOL_RET>`
* :ref:`Rel_Temp_Tol<DIFF_SOL_RTT>`
* :ref:`Residual_Atol<DIFF_SOL_RA>`
* :ref:`Residual_Rtol<DIFF_SOL_RR>`
* :ref:`Stepping_Method<DIFF_SOL_SM>`
* :ref:`Verbose_Stepping<DIFF_SOL_VS>`
* :ref:`Void_Temperature<DIFF_SOL_VT>`

.. _DIFF_SOL_ACT:

Abs_Conc_Tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The tolerance :math:`\epsilon` for the absolute error component of the concentration error norm used by the BDF2 integrator. If :math:`\delta c` is a concentration field increment with reference concentration field :math:`c`,then this error norm is

.. math::
   |||\delta c||| \equiv \mathop{{max}_j} |\delta c_j|/(\epsilon + \eta |c_j|)

The relative error tolerance :math:`\eta` is given by :ref:`Rel_Conc_Tol<DIFF_SOL_RCT>`. This variable is only relevant to the adaptive integrator and to diffusion systems that include concentration as a dependent variable.

| **Physical Dimension**: same as the ‘concentration’ variable
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\geq 0` 
| **Notes** : The error norm is dimensionless and normalized. The BDF2 integrator will accept time steps where the estimated truncation error is less than 2, and chooses the next suggested time step so that its prediction of the next truncation error is :math:`\frac{1}{2}`.

For :math:`c_j` sufficiently small the norm approximates an absolute norm with tolerance :math:`\epsilon`, and for :math:`c_j` sufficiently large the norm approximates a relative norm with tolerance :math:`\eta`. If :math:`\epsilon` = 0 then the norm is a pure relative norm and the concentration must be bounded away from 0. The same tolerance is used for all concentration components.

.. _DIFF_SOL_AET:

Abs_Enthalpy_Tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The tolerance :math:`\epsilon` for the absolute error component of the enthalpy error norm used by the BDF2 integrator. If :math:`\delta H` is a enthalpy field increment with reference enthalpy field :math:`H`, then this error norm is

.. math::
   |||\delta H||| \equiv \mathop{{max}_j} |\delta H_j|/(\epsilon + \eta |H_j|)

The relative error tolerance :math:`\eta` is given by :ref:`Rel_Enthalpy_Tol<DIFF_SOL_RET>`. This variable is only relevant to the adaptive integrator and to diffusion systems that include enthalpy as a dependent variable.

| **Physical Dimension**: :math:`E/(\Theta * L^3)`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\geq 0` 
| **Notes** : The error norm is dimensionless and normalized. The BDF2 integrator will accept time steps where the estimated truncation error is less than 2, and chooses the next suggested time step so that its prediction of the next truncation error is :math:`\frac{1}{2}`.

For :math:`H_j` sufficiently small the norm approximates an absolute norm with tolerance  :math:`\epsilon` , and for :math:`H_j` sufficiently large the norm approximates a relative norm with tolerance :math:`\eta`. If  :math:`\epsilon = 0` then the norm is a pure relative norm and the enthalpy must be bounded away from 0.

The same tolerance is used for all concentration components.

.. _DIFF_SOL_ATT:

Abs_Temp_Tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The tolerance :math:`\epsilon`  for the absolute error component of the temperature error norm used by the BDF2 integrator. If :math:`\delta T` is a temperature field increment with reference temperature field :math:`T`,then this error norm is

.. math::
   |||\delta T||| \equiv \mathop{{max}_j} |\delta T_j|/(\epsilon + \eta |T_j|)

The relative error tolerance :math:`\eta` is given by :ref:`Rel_Temp_Tol<DIFF_SOL_RTT>`. This variable is only relevant to the adaptive integrator and to diffusion systems that include temperature as a dependent variable.

| **Physical Dimension**: :math:`\Theta`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\geq 0` 
| **Notes** : The error norm is dimensionless and normalized. The BDF2 integrator will accept time steps where the estimated truncation error is less than 2, and chooses the next suggested time step so that its prediction of the next truncation error is :math:`\frac{1}{2}`.

For :math:`T_j` sufficiently small the norm approximates an absolute norm with tolerance  :math:`\epsilon`  , and for :math:`T_j` sufficiently large the norm approximates a relative norm with tolerance :math:`\eta`. If :math:`\epsilon = 0`  then the norm is a pure relative norm and the temperature must be bounded away from 0.

The same tolerance is used for all concentration components.

.. _DIFF_SOL_CVT:

Cond_Vfrac_Threshold
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Material volume fraction threshold for inclusion in heat conduction when using the non-adaptive integrator.
| **Type**        : real
| **Default**     : 0.001
| **Valid Values**: (0,1)
| **Notes** : Fluid flow systems that include void will result in partially filled cells, often times with only a tiny fragment of material. Including such cells in the heat conduction problem can cause severe numerical difficulties. By excluding cells with a material volume fraction less than this threshold from participation in heat conduction we can obtain a much better conditioned system. Note that we continue to track enthalpy for such cells, including enthalpy that may be advected into or out of the cell; we just do not consider diffusive transport of enthalpy.

.. _DIFF_SOL_CV:

Cutvof (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description**       : When flow with void is enabled, the heat transfer solver will remove
                          all material from any cell whose non-void volume fraction falls below
                          the value of this threshold.
| **Physical dimension**: dimensionless
| **Type**              : real
| **Default**           : 1e-8 
| **Valid Values**      : [0.0, 1.0)
| **Notes**             : See also :ref:`vol_frac_cutoff<FLOW_VFC>`.

.. _DIFF_SOL_HAD:

Hypre_AMG_Debug
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Enable debugging output from Hypre’s **BoomerAMG** solver. Only relevant when :ref:`NLK_Preconditioner<DIFF_SOL_NP>` is set to **’Hypre_AMG’**.
| **Type**        : logical
| **Default**     : .false. (off)
| **Notes** : See **HYPRE_BoomerAMGSetDebugFlag** in the Hypre Reference Manual.

.. _DIFF_SOL_HALL:

Hypre_AMG_Logging_Level
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Enable additional diagnostic computation by Hypre’s BoomerAMG solver. Only relevant when :ref:`NLK_Preconditioner<DIFF_SOL_NP>` is set to **’Hypre_AMG’**.
| **Type**        : integer
| **Default**     : 0 (none)
| **Valid Values**: 0, none ; >0, varying amounts. 
| **Notes**       :Refer to the Hypre Reference Manual description of **HYPRE_BoomerAMGSetLogging** for details.

.. _DIFF_SOL_HAPL:

Hypre_AMG_Print_Level
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The diagnostic output verbosity level of Hypre’s BoomerAMG solver. Only relevant when :ref:`NLK_Preconditioner<DIFF_SOL_NP>` is set to **’Hypre_AMG’**.

| **Type**        : integer
| **Default**     : 0 
| **Valid Values**: 

.. list-table::
   :widths: 10 50
   :header-rows: 1 
   :class: tight-table   
 
   * - Option
     - Description
   * - 0
     - no output
   * - 1
     - write set up information
   * - 2
     - write solve information
   * - 3
     - write both set up and solve information

| **Notes**       :See **HYPRE_BoomerAMGSetPrintLevel** in the Hypre Reference Manual.


.. _DIFF_SOL_MNI:

Max_NLK_Itr
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The maximum number of NLK nonlinear solver iterations allowed.
| **Type**        : integer
| **Default**     : 5
| **Valid Values**: :math:`\geq 2`
| **Notes** : This variable is used by both the adaptive and non-adaptive integrators, though the appropriate values differ significantly.

For the adaptive integrator, the failure of a nonlinear iteration to convergeis not necessarily fatal;the BDF2 integration procedure expects that this will occur, using it as an indication that the preconditioner for the nonlinear system needs to be updated. If still unsuccessful, the step may be retried with a halved time step size, perhaps repeatedly. Therefore it is important that the maximum number of iterations not be set too high, as this merely delays the recognition that some recovery strategy needs to be taken, and can result in much wasted effort.

By contrast, a nonlinear solver convergence failure is fatal for the non-adaptive solver. Thus the maximum number of iterations should be set to some suitably large value; if the number of iterations ever exceeds this value the simulation is terminated.

The default value is appropriate for the adaptive integrator.

.. _DIFF_SOL_MNV:

Max_NLK_Vec
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The maximum number of acceleration vectors to use in the NLK nonlinear solver.
| **Type**        : integer
| **Default**     : Max_NLK_Itr − 1
| **Valid Values**: :math:`\gt 0`
| **Notes** : The acceleration vectors are derived from the difference of successive nonlinear function iterates accumulated over the course of a nonlinear solve. Thus the maximum possible number of acceleration vectors available is one less than the maximum number of NLK iterations, and so specifying a larger number merely wastes memory. If a large number of NLK iterations is allowed (as when using the non-adaptive integrator) then it may be appropriate to use a smaller value for this parameter, otherwise the default value is fine.

.. _DIFF_SOL_MST:

Max_Step_Tries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The maximum number of attempts to successfully take a time step before giving up. The step size is reduced between each try. This is only relevant to the adaptive solver.
| **Type**        : integer
| **Default**     : 10
| **Valid Values**: :math:`\geq 1`
| **Notes** : If other physics is enabled then this variable is effectively assigned the value 1, overriding the input value. This is required for compatibility with the other physics solvers which currently have no way of recovering from a failed step.

.. _DIFF_SOL_NP:

NLK_Preconditioner
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The choice of preconditioner for the NLK iteration. There are currently two preconditioners to choose from: symmetric successive over-relaxation (SSOR) and the BoomerAMG implementation of algebraic multigrid (AMG) from the HYPRE library. AMG is the default preconditioning method and should normally be used; SSOR is intended mainly for developer use.
| **Type**        : string
| **Default**     : `"hypre_amg"`
| **Valid Values**: `"ssor"` or `"hypre_amg"`
| **Notes** : When `"hypre_amg"` is selected, use the variable :ref:`PC_AMG_Cycles<DIFF_SOL_PAC>` to specify the number of AMG V-cycles per preconditioning step. The BoomerAMG implementation supports a great many configuration options, some which may be set using this namelist. The default configuration is generally suitable, but may be modified using the following variables. Refer to the `HYPRE documentation <https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html>`_ on BoomerAMG for details.

* **hypre_amg_strong_threshold** (real, default 0.5)
* **hypre_amg_coarsen_type** (integer, default 10)
* **hypre_amg_interp_type** (integer, default 6)
* **hypre_amg_relax_down_type** (integer, default 13)
* **hypre_amg_relax_up_type** (integer, default 14)

Note that the variables correspond to similarly named `HYPRE library functions <https://hypre.readthedocs.io/en/latest/api-sol-parcsr.html#>`_ and not actual HYPRE variables.

When `"ssor"` is selected, use the variable :ref:`PC_SSOR_Relax<DIFF_SOL_PSR>` to specify the overrelaxation parameter and :ref:`PC_SSOR_Sweeps<DIFF_SOL_PSS>` to specify the number of SSOR sweeps per preconditioning step.


.. _DIFF_SOL_NT:

NLK_Tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The convergence tolerance for the NLK nonlinear solver. The nonlinear system is considered solved by the current iterate if the BDF2 integrator norm of the last solution correction is less than this value. This variable is only relevant to the adaptive integrator.
| **Type**        : real
| **Default**     : 0.1
| **Valid Values**: (0, 1)
| **Notes** : This tolerance is relative to the dimensionless and normalized BDF2 integrator norm; see :ref:`Abs_Conc_Tol<DIFF_SOL_ACT>`, for example. The nonlinear system only needs to be solved to an accuracy equal to the acceptable local truncation error for the step, which is roughly 1. Solving to a greater accuracy is wasted effort. Using a tolerance in the range (0.01,0.1) is generally adequate to ensure a sufficently converged nonlinear iterate.

.. _DIFF_SOL_NVT:

NLK_Vec_Tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The NLK vector drop tolerance. When assembling the acceleration subspace vector by vector, a vector is dropped when the sine of the angle between the vector and the subspace less than this value.
| **Type**        : real
| **Default**     : 0.001
| **Valid Values**: :math:`\gt 0`

.. _DIFF_SOL_PAC:

PC_AMG_Cycles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The number of V-cycles to take per preconditioning step of the nonlinear iteration.
| **Physcial Dimension**: dimensionless
| **Type**        : integer
| **Default**     : 2
| **Valid Values**: :math:`\gt 1`
| **Notes** : We use standard **V(1,1)** cycles. Parameters other than the number of V cycles cannot be controlled by the user.


.. _DIFF_SOL_PF:

PC_Freq
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : This controls how frequently the preconditioner is updated in the adaptive BDF2 integrator. A value of **N** will allow a preconditioner to be used for as many as N consecutive time steps before being updated, although it may be updated more frequently based on other criteria. A value of 1 causes the preconditioner to be updated every time step. The default behavior is to not require any minimum update frequency.
| **Type**        : integer
| **Default**     : :math:`\infty`
| **Valid Values**: :math:`\geq 1`
| **Notes** : A basic strategy of the adaptive BDF2 integrator is to use a preconditioner for as many time steps as possible, and only update it when a nonlinear time step iteration fails to converge. This generally works quite well. But if you find that the integrator is thrashing — evidenced by the number of times a step failed with an old preconditioner and was retried (this is the NNR diagnostic value in the terminal output) being a significant fraction of the number of time steps — it may be more cost effective to set this value to 1, for example.

.. _DIFF_SOL_PSR:

PC_SSOR_Relax
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The relaxation parameter used in the SSOR preconditioning of the nonlinear system.
| **Physcial Dimension**: dimensionless
| **Type**        : real
| **Default**     : 1.4
| **Valid Values**: (0, 2)
| **Notes** : A value less than 1 gives under-relaxation and a value greater than 1 over-relaxation.

.. _DIFF_SOL_PSS:

PC_SSOR_Sweeps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The number of sweeps used in the SSOR preconditioning of the nonlinear system.
| **Type**        : integer
| **Default**     : 4
| **Valid Values**: :math:`\geq 1`
| **Notes** : The effectiveness of the SSOR preconditioner (measured by the convergence rate of the nonlinear iteration) improves as the number of sweeps increases, though at increasing cost. For especially large systems where the effectiveness of SSOR deteriorates, a somewhat larger value than the default 4 sweeps may be required. Using fewer than 4 sweeps is generally not recommended.

.. _DIFF_SOL_RCT:

Rel_Conc_Tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The tolerance  :math:`\eta`  for the relative error component of the concentration error norm used by the BDF2 integrator. If :math:`\delta c`  is a concentration field increment with reference concentration field :math:`c`, then this error norm is

.. math::
   |||\delta c||| \equiv \mathop{{max}_j} |\delta c_j|/(\epsilon + \eta |c_j|)

The absolute error tolerance :math:`\epsilon`  is given by :ref:`Abs_Conc_Tol<DIFF_SOL_ACT>`. This variable is only relevant to the adaptive solver and to diffusion systems that include concentration as a dependent variable.

| **Physcial Dimension**: dimensionless
| **Type**        : real
| **Default**     : 0.0
| **Valid Values**: (0, 1)
| **Notes** : See the notes for :ref:`Abs_Conc_Tol<DIFF_SOL_ACT>`.

.. _DIFF_SOL_RET:

Rel_Enthalpy_Tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The tolerance  :math:`\eta`  for the relative error component of the enthalpy error norm used by the BDF2 integrator. If :math:`\delta H`  is a concentration field increment with reference concentration field :math:`H`, then this error norm is

.. math::
   |||\delta H||| \equiv \mathop{{max}_j} |\delta H_j|/(\epsilon + \eta |H_j|)

The absolute error tolerance :math:`\epsilon`  is given by :ref:`Abs_Enthalpy_Tol<DIFF_SOL_AET>`. This variable is only relevant to the adaptive solver and to diffusion systems that include enthalpy as a dependent variable.

| **Physcial Dimension**: dimensionless
| **Type**        : real
| **Default**     : 0.0
| **Valid Values**: (0, 1)
| **Notes** : See the notes for :ref:`Abs_Enthalpy_Tol<DIFF_SOL_AET>`.

.. _DIFF_SOL_RTT:

Rel_Temp_Tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The tolerance  :math:`\eta`  for the relative error component of the temperature error norm used by the BDF2 integrator. If :math:`\delta T`  is a concentration field increment with reference concentration field :math:`T`, then this error norm is

.. math::
   |||\delta T||| \equiv \mathop{{max}_j} |\delta T_j|/(\epsilon + \eta |T_j|)

The absolute error tolerance :math:`\epsilon`  is given by :ref:`Abs_Temp_Tol<DIFF_SOL_ATT>`. This variable is only relevant to the adaptive solver and to diffusion systems that include temperature as a dependent variable.

| **Physcial Dimension**: dimensionless
| **Type**        : real
| **Default**     : 0.0
| **Valid Values**: (0, 1)
| **Notes** : See the notes for :ref:`Abs_Temp_Tol<DIFF_SOL_ATT>`.

.. _DIFF_SOL_RA:

Residual_Atol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The absolute residual tolerance :math:`\epsilon_1`  used by the iterative nonlinear solver of the non-adaptive integrator. If :math:`r_o` denotes the initial nonlinear residual, iteration stops when the current residual **r** satisfies

.. math::

   ||r||_2 \leq max[\epsilon_1, \epsilon_2||r_o||_2]

| **Type**        : real
| **Default**     : 0
| **Valid Values**: :math:`\geq 0`
| **Notes** : Ideally this tolerance should be set to 0, but in some circumstances, especially at the start of a simulation, the initial residual may be so small that it is impossible to reduce it by the factor2 due to finite precision arithmetic. In such cases it is necessary to provide this absolute tolerance. It is impossible, however, to say what a suitable value would be, as this depends on the nature of the particular nonlinear system. Some guidance can be obtained through trial-and-error by enabling :ref:`Verbose_Stepping<DIFF_SOL_VS>` and observing the magnitude of the residual norms in the resulting diagnostic output file.

.. _DIFF_SOL_RR:

Residual_Rtol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The relative residual tolerance :math:`\epsilon_2`  used by the iterative nonlinear solver of the non-adaptive integrator. If :math:`r_o` denotes the initial nonlinear residual, iteration stops when the current residual **r** satisfies

.. math::

   ||r||_2 \leq max[\epsilon_1, \epsilon_2||r_o||_2]

| **Type**        : real
| **Default**     : none
| **Valid Values**: (0, 1)

.. _DIFF_SOL_SM:

Stepping_Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The choice of time integration method.
| **Type**        : string
| **Default**     : `Adaptive BDF2`
| **Valid Values**: `Adaptive BDF2` or `Non-adaptive BDF1`
| **Notes** : The non-adaptive integrator must be selected when fluid flow is enabled and void material is present. Otherwise use the default adaptive integrator.

.. _DIFF_SOL_VS:

Verbose_Stepping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A flag that enables the output of detailed BDF2 time stepping information. The human-readable information is written to a file with the suffix ``.bdf2.out`` in the output directory.
| **Type**        : logical
| **Default**     : false

.. _DIFF_SOL_VT:

Void_Temperature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : An arbitrary temperature assigned to cells that contain only the void material. The value has no effect on the simulation, and is only significant to visualization.
| **Type**        : real
| **Default**     : 0
