.. _PHYSICS_Namelist:

.. toctree::
   :maxdepth: 1

PHYSICS Namelist
==================
The PHYSICS namelist specifies which physics models are active in the
simulation. The models are implemented by the four primary physics solvers --
fluid flow, heat/species transport, electromagnetics, and solid mechanics --
which are weakly coupled using time splitting. A brief overview of the physics
solvers follows; see Truchas Physics and Algorithms for more details.

**Fluid Flow**. The fluid flow physics model simulates multi-material,
incompressible flow with interface tracking. A gravitational body force is
defined using the `body_force_density`_ variable. See the
:ref:`MATERIAL<MATERIAL_and_PHASE_Namelists>` namelist for a description
of the material properties required by the fluid flow model.

**Heat and Species Transport**.  The heat and species transport physics model
simulates both heat conduction with thermal (view factor) radiation, and
solutal species diffusion and thermodiffusion. These (primarily) diffusive
transport processes are fully coupled; advection of enthalpy and solutal
species are handled by the fluid flow physics solver and incorporated as
explicit source terms. Heat transport is enabled using the `heat_transport`_
flag, and solves the heat equation

.. math::
  :label: eq_ht
  
  \frac{\partial H}{\partial T} = \nabla\cdot K \nabla T + Q + Q_\text{adv} 

with dependent variables temperature :math:`T` and volumetric enthalpy density
:math:`H`. The enthalpy density is algebraically related to temperature as
:math:`H=f(T)` where :math:`f'(T)=\rho c_p` is the volumetric heat capacity.
See the :ref:`MATERIAL<MATERIAL_and_PHASE_Namelists>` namelist for a description
of the material properties required by the heat equation. The optional
volumetric heat source :math:`Q` is defined using the
:ref:`THERMAL_SOURCE<THERMAL_SOURCE_Namelist>` namelist,
and the advected heat :math:`Q_\text{adv}` by the flow solver. The boundary
conditions on :math:`T` are defined using the
:ref:`THERMAL_BC<THERMAL_BC_Namelist>` namelists. The initial value of
:math:`T` is defined by the `temperature` variable of the
:ref:`BODY<BODY_Namelist>` namelists. View factor radiation systems which
couple to the heat equation are defined using
:ref:`ENCLOSURE_RADIATION<ENCLOSURE_RADIATION_Namelist>` namelists. Solutal
species transport is enabled using the `species_transport`_ flag, which solves
the :math:`n` coupled equations

.. math::
  :label: eq_st
  
  \frac{\partial\phi_i}{\partial t} = \nabla\cdot D_i
      (\nabla\phi_i~[{}+S_i\nabla T]~) + Q_i + Q_{i,\text{adv}} 

for species concentrations :math:`\phi_i`. The number of components :math:`n`
is defined by `number_of_species`_. The thermodiffusion
term in :eq:`eq_st` is only included when coupled with heat transport.
See the :ref:`MATERIAL<MATERIAL_and_PHASE_Namelists>` namelist for defining
the diffusivities :math:`D_i` and Soret coefficients :math:`S_i`. The optional
volumetric source :math:`Q_i` is defined using the
:ref:`SPECIES_SOURCE<SPECIES_SOURCE_Namelist>` namelist. The advected species
source :math:`Q_{i,\text{adv}}` is computed
by the flow solver. Boundary conditions on :math:`\phi_i` are defined through
the :ref:`SPECIES_BC<SPECIES_BC_Namelist>` namelists. The initial value of the
:math:`\phi_i` are defined through the `conc` variable of the
:ref:`BODY<BODY_Namelist>` namelists.

**Induction Heating**. The induction heating model solves the heat transport
model above with an additional Joule heat source computed by an auxiliary
electromagnetics problem defined by the
:ref:`INDUCTION_HEATING<INDUCTION_HEATING_Namelist>` namelist that models the
eddy currents induced by an external low-frequency magnetic field. It is
enabled by the `induction_heating`_ flag. See the
:ref:`MATERIAL<MATERIAL_and_PHASE_Namelists>` namelist for a
description of the material properties required by the electromagnetics solver.

**Microwave Heating**. The microwave heating model solves the heat transport
model with an additional dielectric heat source computed by an auxiliary
microwave-frequency electromagnetics problem defined by the
:ref:`MICROWAVE_HEATING<MICROWAVE_HEATING_Namelist>` namelist. It is enabled
by the `microwave_heating`_ flag. See the
:ref:`MATERIAL<MATERIAL_and_PHASE_Namelists>` namelist for a description of
the material properties required by the electromagnetics solver. 

**Solid Mechanics**.  The solid mechanics physics kernel models small strain
elastic deformation of solid material phases, including deformations induced
by temperature changes. It is enabled using the `solid_mechanics`_ flag. See
the :ref:`MATERIAL<MATERIAL_and_PHASE_Namelists>` namelist for a description
of the material properties required by the solid mechanics kernel.
Displacement and traction boundary conditions are defined using
:ref:`SM_BC<SM_BC_Namelist>` namelists. A gravitational body force may be
defined using the `body_force_density`_ variable. Parameters which define
the plasticity model are defined using the
:ref:`VISCOPLASTIC_MODEL<VISCOPLASTIC_MODEL_Namelist>` namelist.


.. admonition:: Namelist usage

   :Required/Optional: Required
   :Single/Multiple Instances: Single


Namelist Variables
------------------
.. contents::
   :local:

body_force_density
^^^^^^^^^^^^^^^^^^^
A constant force per unit mass, **g**, that acts throughout material volumes.
The net force on a volume is the integral of its density times **g** over the
volume. Typically **g** is the gravitational acceleration.

:Type: real 3-vector
:Default: (0.0, 0.0, 0.0)
:Note: The fluid flow and solid mechanics models always include this body force.

.. _physics-ih:

induction_heating
^^^^^^^^^^^^^^^^^^
Enables the simulation of induction heating.

:Type: logical
:Default: false


heat_transport
^^^^^^^^^^^^^^^^^^
Enables the simulation of heat conduction, advection, and radiation using the
heat/species transport physics solver.

:Type: logical
:Default: false


flow
^^^^^^^^^^^^^^^^^^
Enables the simulation of fluid flow.

:Type: logical
:Default: false


materials
^^^^^^^^^^^^^^^^^^
A list of materials to include in the simulation. These are material names
defined in :ref:`MATERIAL<MATERIAL_and_PHASE_Namelists>` namelists. The list
must include all materials assigned to a region in a :ref:`BODY<BODY_Namelist>`
namelist, or specified as an :ref:`inflow_material<FLOW_BC_IM>` in a fluid flow
boundary condition, but it need not include all materials defined in the input
file. Use the reserved name **"VOID"** to refer to the built-in void
pseudo-material.

:Type: string list
:Default: none

.. _physics-mwh:

microwave_heating
^^^^^^^^^^^^^^^^^^
Enables the simulation of microwave heating.

:Type: logical
:Default: false


number_of_species
^^^^^^^^^^^^^^^^^^
The number of species components. Required when `species_transport`_ is enabled.

:Type: integer
:Default: 0
:Valid Values: > 0


solid_mechanics
^^^^^^^^^^^^^^^^^^
Enables the calculation of solid material stresses and strains.

:Type: logical
:Default: false


species_transport
^^^^^^^^^^^^^^^^^^
Enables the calculation of species diffusion and advection using the
heat/species transport physics kernel. The number of species components must
be specified using `number_of_species`_.

:Type: logical
:Default: false
