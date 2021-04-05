.. _FLOW_Namelist:

.. toctree::
   :maxdepth: 1

FLOW Namelist 
==============================

Overview
------------
The :ref:`FLOW<FLOW_Namelist>` namelist specifies the parameters for the fluid flow model and algorithm. This namelist is read whenever the :ref:`PHYSICS<PHYSICS_Namelist>` namelist option :ref:`Flow<PHYSICS_F>` is enabled. Parameters for the linear solvers employed by the algorithm are specified using :ref:`FLOW_VISCOUS_SOLVER<FLOW_PRESSURE_SOLVER_and_FLOW_VISCOUS_SOLVER_Namelists>` and :ref:`FLOW_PRESSURE_SOLVER<FLOW_PRESSURE_SOLVER_and_FLOW_VISCOUS_SOLVER_Namelists>` namelists. Flow boundary conditions are defined using :ref:`FLOW_BC<FLOW_BC_Namelist>` namelists.

FLOW Namelist Features
----------------------------
| **Required/Optional        :** Required when flow physics is enabled.
| **Single/Multiple Instances:** Single

Components
------------
Physics Options
^^^^^^^^^^^^^^^^^^
* :ref:`inviscid<FLOW_inviscid>`

Numerical Parameters
^^^^^^^^^^^^^^^^^^^^^^^
* :ref:`courant_number<FLOW_CN>`
* :ref:`viscous_number<FLOW_VN>`
* :ref:`viscous_implicitness<FLOW_VI>`
* :ref:`track_interfaces<FLOW_TI>`
* :ref:`material_priority<FLOW_MP>`
* :ref:`vol_track_subcycles<FLOW_VTS>`
* :ref:`nested_dissection<FLOW_ND>` (expert)
* :ref:`vol_frac_cutoff<FLOW_VFC>` (expert)
* :ref:`fischer_dim<FLOW_FD>` (expert)
* :ref:`fluid_frac_threshold<FLOW_FFT>` (expert)
* :ref:`min_face_fraction<FLOW_MFF>` (expert)
* :ref:`void_collapse<FLOW_VC>` (experimental)
* :ref:`void_collapse_relaxation<FLOW_VCR>` (experimental)
* :ref:`wisp_redistribution<FLOW_WR>` (experimental)
* :ref:`wisp_cutoff<FLOW_WC>` (experimental)
* :ref:`wisp_neighbor_cutoff<FLOW_WNC>` (experimental)

.. _FLOW_inviscid:

inviscid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : This option omits the viscous forces from the flow equations. In this case there is no viscous system to solve and the :ref:`FLOW_VISCOUS_SOLVER<FLOW_PRESSURE_SOLVER_and_FLOW_VISCOUS_SOLVER_Namelists>` namelist is not required.
| **Type**        : logical
| **Default**     : false

.. _FLOW_CN:

courant_number
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : This parameter sets an upper bound on the time step that is associated with stability of the explicit fluid advection algorithm. A value of 1 corresponds (roughly) to the stability limit, with smaller values resulting in proportionally smaller allowed time steps. Truchas uses the largest timestep possible subject to this and other limits.
| **Type**        : real
| **Default**     : 0.5
| **Valid values**: (0.0,1.0]
| **Notes**       : The Courant number for a cell is the dimensionless value :math:`C_i = u_i \Delta t/\Delta x_i` where :math:`\Delta t` is the timestep, :math:`u_i` the fluid velocity magnitude on the cell, and :math:`\Delta x_i` a measure of the cell size. The time steplimit is the largest :math:`\Delta t` such that :math:`max\{C_i\}` equals the value of :ref:`courant_number<FLOW_CN>`. The interpretation of :math:`u_i` and :math:`\Delta x_i` for a general cell is somewhat sticky. Currently a ratio :math:`u_f/h_f` is computed for each face of a cell and the maximum taken for the value of :math:`u_i/\Delta x_i`. Here :math:`u_f` is the normal fluxing velocity on the face, and :math:`h_f` is the inscribed cell height at the face; that is, the minimum normal distance between the face and cell nodes not belonging to the face.

.. _FLOW_VN:

viscous_number
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : This parameter sets an upper bound on the time step that is associated with stability of an explicit treatment of viscous flow stress tensor. A value of 1 corresponds roughly to the stability limit, with smaller values resulting in proportionally smaller allowed time steps. Truchas uses the largest time step possible subject to this and other limits. For an implicit treatment of the viscous flow stress tensor with :ref:`viscous_implicitness<FLOW_VI>` at least 12, which is always strongly recommended, the viscous discretization is unconditionally stable and no :math:`time \enspace step \enspace limit` is needed. In this case, the parameter can still be used to limit the time step for accuracy. A value of :math:`0` will disable this limit entirely.
| **Type**        : real
| **Default**     : :math:`0`
| **Valid values**: :math:`\geq 0`
| **Notes**       : The viscous number for a cell is the dimensionless value :math:`V_i = \nu_i \Delta t / \Delta x_i^2`, where :math:`\Delta_t` is the timestep, :math:`\nu_i` the kinematic viscosity :math:`(\mu/\rho)` on the cell, and :math:`\Delta x_i` a measure of the cell size. The time steplimit is the largest :math:`\Delta t` such that :math:`max\{V_i\}` equals the value of :ref:`viscous_number<FLOW_VN>`. Currently the measure of cell size mirrors that used in the definition of the :ref:`courant_number<FLOW_CN>`, namely that :math:`\Delta x_i` is taken as the minimum of the inscribed heights :math:`h_f` for the faces of the cell.

.. _FLOW_VI:

viscous_implicitness
^^^^^^^^^^^^^^^^^^^^^^

| **Description** : The degree of time implicitness :math:`\theta` used for the velocity field in the discretization of the viscous flow stress tensor in the fluid momentum conservation equation. The velocity is given by the :math:`\theta`-method, :math:`u_\theta = (1−\theta)u_n − \theta u_{n+1}:\theta = 0` gives an explicit discretization and :math:`\theta = 1` a fully implicit discretization. In practice only the values :math:`1, \frac{1}{2}` (trapezoid method), and 0 are useful, and use of the latter explicit discretization is generally not recommended. Note that an implicit discretization, :math:`\theta > 0`, will require the solution of a linear system; see :ref:`FLOW_VISCOUS_SOLVER<FLOW_PRESSURE_SOLVER_and_FLOW_VISCOUS_SOLVER_Namelists>`. This parameter is not relevant to :ref:`inviscid<FLOW_inviscid>` flow problems
| **Type**    : real
| **Default** : 1
| **Valid Values**: [0,1]
| **Notes**   : The discretization is first order except for the trapezoid method (:math:`\theta` = 12) which is second order. However note that the flow algorithm overall is only first order irrespective of the choice of :math:`\theta`. The advanced velocity :math:`u_{n+1}` is actually the predicted velocity :math:`u_{n+1}^*` from the predictor stage of the flow algorithm. 

.. _FLOW_TI:

track_interfaces
^^^^^^^^^^^^^^^^^^^^^^

| **Description** : This option enables the tracking of material interfaces. The default is to track interfaces whenever the problem involves more than one material. If the problem involves a single fluid and it is known a priori that there will never be any mixed material cells containing fluid, then this option can be set to false to short-circuit some unnecessary work, but otherwise the default should be used.
| **Type**    : logical
| **Default** : true

.. _FLOW_MP:

material_priority
^^^^^^^^^^^^^^^^^^^^^^

| **Description** : A list of material names that defines the priority order in which material interfaces are reconstructed within a cell for volume tracking. All fluid material names must be listed, and if the problem includes any non-fluid materials, this list must include the case-sensitive keyword "SOLID", which stands for all non-fluid materials lumped together. The default is the list of fluid materials in input file order, followed by "SOLID" for the lumped non-fluids.
| **Type**    : string list
| **Notes**   : Different priorities will result in somewhat different results. Unfortunately there are no hard and fast rules for selecting the priorities.
  
.. _FLOW_VTS:

vol_track_subcycles
^^^^^^^^^^^^^^^^^^^^^^

| **Description** : The number of sub-time steps :math:`n` taken by the volume tracker for every time step of the flow algorithm. If the flow time step size is :math:`\Delta t` then the volume tracker will taken time steps of size :math:`\Delta t/n` to compute the net flux volumes and advance the volume fractions for the flow step. 
| **Type**   : integer
| **Default**: 2
| **Notes**  : With the current unsplit advection algorithm :footcite:`rider1998reconstructing` it is necessary to sub-cycle the the volume tracking time integration method in order to obtain good “corner coupling” of the volume flux terms.

.. _FLOW_ND:

nested_dissection (Expert Parameter)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : This option enables use of the nested dissection algorithm to reconstruct material interfaces in cells containing 3 or more materials. If set false the less accurate and less expensive onion skin algorithm will be used.
| **Type**   : logical
| **Default**: True

.. _FLOW_VFC:

vol_frac_cutoff (Expert Parameter)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : The smallest material volume fraction allowed. If a material volume fraction drops below this cutoff, the material is removed entirely from the cell, and its volume fraction replaced by proportional increases to the volume fractions of the remaining materials, or if the cell contains void, by increasing the void volume fraction alone.
| **Type**   : real
| **Default**: :math:`10^{-8}`
| **Valid Values**  : (0,1)

.. _FLOW_FD:

fischer_dim (Expert Parameter)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : The dimension :math:`d` of the subspace used in Fischer’s projection method :footcite:`fischer1998projection` for computing aninitial guess for pressure projection system based on previous solutions. Memory requirements for the method are :math:`2(d+ 1)` cell-based vectors. Set this variable to 0 to disable use of this method.
| **Type**   : integer
| **Default**: :math:`6`

.. _FLOW_FFT:

fluid_frac_threshold (Expert Parameter)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : Cells with a total fluid volume fraction less than this threshold are ignored by the flow solver, being regarded as ‘solid’ cells.
| **Type**   : real
| **Default**: :math:`10^{-2}`
| **Valid values**: (0,1)


.. _FLOW_MFF:

min_face_fraction (Expert Parameter)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : The variable sets the minimum value of the fluid density associated with a face for the pressure projection system. It is specified as a fraction of the minimum fluid density (excluding void) of any fluid in the problem.
| **Type**   : real
| **Default**: :math:`10^{-3}`

.. _FLOW_VC:

void_collapse (Experimental)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : The volume-of-fluid algorithm effectively treats small fragments of void entrained in fluid as incompressible, resulting in unphysical void “bubbles” that persist in the flow. A model that drives the collapse of these void fragments will be enabled when this variable is set to true. See :ref:`void_collapse_relaxation<FLOW_VCR>` for a model parameter.
| **Type**   : logical
| **Default**: False

.. _FLOW_VCR:

void_collapse_relaxation (Experimental)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : The relaxation parameter in the void collapse model. See :ref:`void_collapse<FLOW_VC>`.
| **Type**   : real
| **Default**: 0.1
| **Valid Values**: [0,1]
| **Notes**   : The relaxation factor is roughly inversely proportional to the number of timesteps required for all the void in a cell to collapse as dictated by inertial forces. Thus a relaxation factor of 1 would allowfor all the void in a cell to collapse over a single timestep. A factor of 0.1 would allow for the void in a cell to collapse over the course of 10 timesteps. Larger values tend to cause more mass loss (on theorder of 0.5%), although the results do improve with increased subcycling.

.. _FLOW_WR:

wisp_redistribution (Experimental)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : A cell containing a small amount of fluid (wisps) can sometimes trigger pathological behavior, manifesting as a perpetual acceleration that drives the timestep to 0. A model that redistributes these wisps to other other fluid cells will be enabled when this variable is set to true. See :ref:`wisp_cutoff<FLOW_WC>` and :ref:`wisp_neighbor_cutoff<FLOW_WNC>` model parameters.
| **Type**   : logical
| **Default**: false

.. _FLOW_WC:

wisp_cutoff (Experimental)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : Fluid cells with a fluid volume fraction below this value may be treated as wisps and havetheir fluid material moved around to other fluid cells. Generally, moving fluid around the domain can have an undesirable effect on the flow physics. It is therefore advisable to keep this value as small as possible. Based on numerical experimentation, the recommended value is 0.05. Smaller values, such as, 0.01 did not result in robust simulations. See :ref:`wisp_redistribution<FLOW_WR>`.
| **Type**   : real
| **Default**: 0.05

.. _FLOW_WNC:

wisp_neighbor_cutoff (Experimental)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : Fluid cells with a fluid volume fraction below wisp_cutoff can only be considered a wisp if the amount of fluid in the neighboring cells is also "small". This definition of "small" is controlled by wisp_neighbor_cutoff. In addition, for a a cell to receive wisp material, the fluid volume fraction of it and its neighbors must be larger than wisp_neighbor_cutoff. See :ref:`wisp_redistribution<FLOW_WR>`.
| **Type**   : real
| **Default**: 0.25




.. footbibliography::
   
