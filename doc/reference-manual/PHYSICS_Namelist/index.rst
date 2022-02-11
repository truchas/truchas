.. _PHYSICS_Namelist:

.. toctree::
   :maxdepth: 1

PHYSICS Namelist
==================

Overview
-----------

The PHYSICS namelist specifies which physics models are active in the simulation. The models are implemented by the four primary physics kernels — fluid flow, heat/species transport, induction heating, and solid mechanics — which are weakly coupled using time splitting. A brief overview of the physics kernels follows; see Truchas Physics and Algorithms for more details.

**Fluid Flow**. The fluid flow physics model simulates multi-material, incompressible flow with interface tracking. A gravitational body force is defined using the :ref:`Body_Force_Density<PHYSICS_BFD>` variable. See the :ref:`MATERIAL<MATERIAL_and_PHASE_Namelists>` namelist for a description of the material properties required by the fluid flow model.

**Heat and Species Transport**.  The heat and species transport physics kernel models both heat conduction with thermal (view factor) radiation, and solutal species diffusion and thermodiffusion. These (primarily) diffusive transport processes are fully coupled; advection of enthalpy and solutal species are handled by the fluid flow physics kernel and incorporated as loosely-coupled source terms. Heat transport is enabled using the :ref:`Heat_Transport<PHYSICS_HT>` flag, and solves the heat equation

.. math::
  \frac{\partial H}{\partial T} = \Delta . K \Delta T + Q + Q_{joule} + Q_{adv} 
 :label: eq_ht

with dependent variables temperature :math:`T` and enthalpy density :math:`H`. The enthalpy density is algebraically related to temperature as :math:`H=f(T)` where :math:`f′(T) =\rho c_p`is the volumetric heat capacity. See the :ref:`MATERIAL<MATERIAL_and_PHASE_Namelists>` namelist for a description of the material properties required by the heat equation. The optional volumetric heat source :math:`Q` is defined through the :ref:`DS_SOURCE<DS_SOURCE_Namelist>` namelist using "temperature" as the equation name. The Joule heating source :math:`Q_{joule}` is computed by the induction heating kernel, and the advected heat :math:`Q_{adv}` by the flow kernel. The boundary conditions on :math:`T` are defined through the :ref:`THERMAL_BC<THERMAL_BC_Namelist>` namelists. The initial value of :math:`T` are defined through the Temperature variable of the :ref:`BODY<BODY_Namelist>` namelists. View factor radiation systems which couple to the heat equation are defined using :ref:`ENCLOSURE_RADIATION<ENCLOSURE_RADIATION_Namelist>` namelists. Solutal species transport is enabled using the :ref:`Species_Transport<PHYSICS_ST>` flag, which solves the :math:`n` coupled equations

.. math::
  \frac{\partial \phi_i}{\partial t} = \Delta . D_i (\Delta \phi_i[+S_i\Delta T]) + Q_i + Q_{{i},{adv}} 
 :label: eq_st

for species concentrations :math:`\phi_i`. The number of components :math:`n` is defined by :ref:`Number_of_Species <PHYSICS_NOS>`. The thermodiffusion term in [·] is only included when coupled with heat transport. See the :ref:`MATERIAL<MATERIAL_and_PHASE_Namelists>` namelist for defining the diffusivities :math:`D_i` and Soret coefficients :math:`S_i`. The optional volumetric source :math:`Q_i` is defined through the :ref:`DS_SOURCE<DS_SOURCE_Namelist>` namelist using "concentration i" as the equation name. The advected species source :math:`Q_{i,adv}` is computed by the flow kernel. Boundary conditions on :math:`\phi_i` are defined through the :ref:`SPECIES_BC<SPECIES_BC_Namelist>` namelists. The initial value of the :math:`\phi_i` are defined through the :ref:`Phi<B_phi>` variable of the :ref:`BODY<BODY_Namelist>` namelists.

**Induction Heating**.  The induction heating physics kernel solves for the Joule heat that is used as a source in heat transport. It is enabled using the :ref:`Electromagnetics<PHYSICS_EM>` flag. See the :ref:`MATERIAL<MATERIAL_and_PHASE_Namelists>` namelist for a description of the material properties required by the electromagnetics solver. The :ref:`Electromagnetics<PHYSICS_EM>` namelist is used to describe the induction heating problem.

**Solid Mechanics**.  The solid mechanics physics kernel models small strain elastic deformation of solid material phases, including deformations induced by temperature changes. It is enabled using the :ref:`Solid_Mechanics<PHYSICS_SM>` flag. See the :ref:`MATERIAL<MATERIAL_and_PHASE_Namelists>` namelist for a description of the material properties required by the solid mechanics kernel. Displacement and traction boundary conditions are defined using :ref:`SM_BC<SM_BC_Namelist>` namelists. A gravitational body force may be defined using the :ref:`Body_Force_Density<PHYSICS_BFD>` variable.

**Legacy Solid Mechanics**.  The solid mechanics physics kernel models small strain elastic and plastic deformation of solid material phases, including deformations induced by temperature changes and solid state phase changes. It is enabled using the :ref:`Legacy_Solid_Mechanics<PHYSICS_LSM>` flag. See the :ref:`MATERIAL<MATERIAL_and_PHASE_Namelists>` namelist for a description of the material properties required by the solid mechanics kernel. Parameters which define the plasticity model are defined using the :ref:`VISCOPLASTIC_MODEL<VISCOPLASTIC_MODEL_Namelist>` namelist. Displacement and traction boundary conditions are defined using the :ref:`BC<BC_Namelist>` namelist. The effect of the gravitational body force defined by :ref:`Body_Force_Density<PHYSICS_BFD>` can be included by enabling the :ref:`Solid_Mechanics_Body_Force<LSM_SMBF>` flag in the :ref:`LEGACY_SOLID_MECHANICS<LEGACY_SOLID_MECHANICS_Namelist>` namelist.

PHYSICS Namelist Features
---------------------------
| **Required/Optional        :** Required
| **Single/Multiple Instances:** Single

Components
------------
* :ref:`Body_Force_Density<PHYSICS_BFD>`
* :ref:`Electromagentics<PHYSICS_EM>`
* :ref:`Flow<PHYSICS_F>`
* :ref:`Heat_Transport<PHYSICS_HT>`
* :ref:`Legacy_Solid_Mechanics<PHYSICS_LSM>`
* :ref:`Materials<PHYSICS_M>`
* :ref:`Number_of_Species<PHYSICS_NOS>`
* :ref:`Solid_Mechanics<PHYSICS_SM>`
* :ref:`Species_Transport<PHYSICS_ST>`

.. _PHYSICS_BFD:

Body_Force_Density
^^^^^^^^^^^^^^^^^^^

| **Description** : A constant force per unit mass, **g**, that acts throughout material volumes. The net force ona volume is the integral of its density times **g** over the volume. Typically **g** is the gravitational acceleration.
| **Physical dimension**: :math:`L/T^2`
| **Type**        : real 3-vector
| **Default**     : (0.0, 0.0, 0.0)
| **Note**: The fluid flow and solid mechanics models always include this body force. The legacy solid mechanics model has the option of including this body force or not; see :ref:`Solid_Mechanics_Body_Force<LSM_SMBF>`.

.. _PHYSICS_EM:

Electromagentics
^^^^^^^^^^^^^^^^^^

| **Description** : Enables the calculation of Joule heating.
| **Type**        : logical
| **Default**     : false

.. _PHYSICS_HT:

Heat_Transport
^^^^^^^^^^^^^^^^^^

| **Description** : Enables the calculation of heat conduction, advection, and radiation using the heat/species transport physics kernel.
| **Type**        : logical
| **Default**     : false

.. _PHYSICS_F:

Flow
^^^^^^^^^^^^^^^^^^

| **Description** : Enables the simulation of fluid flow.
| **Type**        : logical
| **Default**     : false

.. _PHYSICS_LSM:

Legacy_Solid_Mechanics
^^^^^^^^^^^^^^^^^^^^^^^

| **Description** : Enables the legacy package for calculation of solid material stresses and strains.
| **Type**        : logical
| **Default**     : false

.. _PHYSICS_M:

Materials
^^^^^^^^^^^^^^^^^^

| **Description** : A list of materials to include in the simulation. These are material names defined in :ref:`MATERIAL<MATERIAL_and_PHASE_Namelists>` namelists. The list must include all materials assigned to a region in a :ref:`BODY<BODY_Namelist>` namelist, or specified as an :ref:`inflow_material<FLOW_BC_IM>` in a fluid flow boundary condition, but it need not include all materials defined in the input file. Use the reserved name "VOID" to refer to the built-in void pseudo-material.
| **Type**        : string list

.. _PHYSICS_NOS:

Number_of_Species
^^^^^^^^^^^^^^^^^^

| **Description** : The number of species components. Required when :ref:`Species_Transport<PHYSICS_ST>` is enabled.
| **Type**        : integer
| **Default**     : 0
| **Valid Values**: > 0

.. _PHYSICS_SM:

Solid_Mechanics
^^^^^^^^^^^^^^^^^^

| **Description** : Enables the calculation of solid material stresses and strains.
| **Type**        : logical
| **Default**     : false

.. _PHYSICS_ST:

Species_Transport
^^^^^^^^^^^^^^^^^^

| **Description** : Enables the calculation of species diffusion and advection using the heat/species transport physics kernel. The number of species components must be specified using :ref:`Number_of_Species<PHYSICS_NOS>`.
| **Type**        : logical
| **Default**     : false
