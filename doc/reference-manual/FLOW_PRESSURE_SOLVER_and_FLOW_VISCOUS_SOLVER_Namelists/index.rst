.. _FLOW_PRESSURE_SOLVER_and_FLOW_VISCOUS_SOLVER_Namelists:

.. toctree::
   :maxdepth: 1

FLOW_PRESSURE_SOLVER and FLOW_VISCOUS_SOLVER Namelist 
========================================================

Overview
-------------------

The flow algorithm requires the solution of two linear systems at each time step: the implicit viscous velocity update system and the pressure Poisson system. Truchas uses the hybrid solver from the `HYPRE software library <https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods>`_ to solve these systems. 

The hybrid solver first uses a diagonally-scaled iterative Krylov solver. If it determines that convergence is too slow, the solver switches to a more expensive but more effective preconditioned Krylov solver that uses an algebraic multigrid (AMG) preconditioner (BoomerAMG).

The **FLOW_VISCOUS_SOLVER** namelist sets the HYPRE hybrid solver parameters for the solution of the implicit viscous velocity update system, and the **FLOW_PRESSURE_SOLVER** namelist sets the solver parameters for the solution of the pressure Poisson system. The same variables are used in both namelists.

FLOW_VISCOUS_SOLVER Namelist Features
---------------------------------------
| **Required/Optional        :** Required only for viscous flow with :ref:`viscous_implicitness<FLOW_VI>` > 0.
| **Single/Multiple Instances:** Single

FLOW_VISCOUS_SOLVER Namelist Features
---------------------------------------
| **Required/Optional        :** Required 
| **Single/Multiple Instances:** Single


Components
--------------

* :ref:`krylov_method<FPS_FVS_KM>`
* :ref:`krylov_dim<FPS_FVS_KD>`
* :ref:`conv_rate_tol<FPS_FVS_CRT>`
* :ref:`abs_tol<FPS_FVS_AT>`
* :ref:`rel_tol<FPS_FVS_AT>`
* :ref:`max_ds_iter<FPS_FVS_MDI>`
* :ref:`max_amg_iter<FPS_FVS_MAI>`
* :ref:`print_level<FPS_FVS_PL>`

.. _FPS_FVS_KM:

krylov_method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Selects the Krylov method used by the HYPRE hybrid solver.
| **Type**        : string
| **Default**     : cg
| **Valid Values**: cg, gmres, bicgstab

.. _FPS_FVS_KD:

krylov_dim
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The Krylov subspace dimension for the restarted GMRES method.
| **Type**        : integer
| **Default**     : :math:`5`
| **Valid Values**: :math:`\gt 0`

.. _FPS_FVS_CRT:


conv_rate_tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The convergence rate tolerance :math:`\theta` where the hybrid solver switches to the more expensive **AMG** preconditioned **Krylov** solver. The average convergence rate after :math:`n` iterations of the diagonally-scaled Krylov solver is :math:`\rho_n = (‖r_n‖/‖r_0‖)^{1/n}`, where :math:`r_n = Ax_n − b` is the residual of the linear system, and its convergence is considered too slow when

:math:`[1-\frac{|\rho_n-\rho_{n-1}|}{max(\rho_n,\rho_{n-1})}]\rho_n \gt \theta`

| **Type**        : real
| **Default**     : :math:`0.9`
| **Valid Values**: :math:`(0,1)`

.. _FPS_FVS_AT:


abs_tol, rel_tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The absolute and relative error tolerances :math:`\epsilon_1` and :math:`\epsilon_2` for the solution of the linear system. The test for convergence is :math:`‖r‖ \leq max\{\epsilon_1,\epsilon_2‖b‖\}`, where :math:`r = Ax − b` is the residual of the linear system.

| **Type**        : real
| **Default**     : none

.. _FPS_FVS_MDI:

max_ds_iter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The maximum number of diagonally scaled Krylov iterations allowed. If convergence is not achieved within this number of iterations the hybrid solver will switch to the preconditioned Krylov solver.
| **Type**        : integer

.. _FPS_FVS_MAI:

max_amg_iter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The maximum number of preconditioned Krylov iterations allowed. 
| **Type**        : integer

.. _FPS_FVS_PL:

print_level
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Set this parameter to 2 to have HYPRE write diagnostic data to the terminal for each solver iteration. This is only useful in debugging situations.
| **Default**        : 0 (no output)

Additional HYPRE parameters (Expert)
---------------------------------------
Some additional HYPRE solver parameters and options can be set using these namelists. Nearly all of these are associated with the BoomerAMG preconditioner, and all have reasonable defaults set by HYPRE. See the ParCSR Hybrid Solver section in the `HYPRE reference manual <https://computing.llnl.gov/sites/default/files/public/hypre-2.11.1_ref_manual.pdf>`_  for details. The `HYPRE documentation <https://hypre.readthedocs.io/en/latest/ch-intro.html>`_ has some additional information. The variables that can be set are listed below. Note that the variables correspond to similarly-named HYPRE library functions and not actual HYPRE variables. Also note that there are many parameters and options that cannot currently be set by the namelists.

* cg_use_two_norm(logical)
* amg_strong_threshold(real)
* amg_max_levels(integer)
* amg_coarsen_method(integer)
* amg_smoothing_sweeps(integer)
* amg_smoothing_method(integer)
* amg_interp_method(integer)




  