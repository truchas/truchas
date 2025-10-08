.. _THERMAL_BC_Namelist:

.. toctree::
   :maxdepth: 2

THERMAL_BC Namelist
====================
The THERMAL_BC namelist is used to define boundary conditions for the heat
transfer model at external boundaries and internal interfaces. Each instance
of the namelist defines a particular condition to impose over a subset of the
domain boundary. The boundary subset :math:`\Gamma` is specified using mesh
face sets. The namelist variable `face_set_ids`_ takes a list of face set IDs,
and the boundary condition is imposed on all faces belonging to those face
sets. Note that ExodusII mesh side sets are imported into Truchas as face sets
with the same IDs.

.. note::

   :Required/Optional: Required
   :Single/Multiple Instances: Multiple


Conditions on External Boundaries
---------------------------------
The following types of external boundary conditions can be defined. The outward
unit normal to the boundary :math:`\Gamma` is denoted :math:`\hat{n}`.

Temperature
+++++++++++
A temperature Dirichlet condition,

.. math::
   T = T_b\text{ on $\Gamma$},

is defined by setting `type`_ to **"temperature"**. The boundary value
:math:`T_b` is specified using either `temp`_ for a constant value, or
`temp_func`_ for a function.

Heat Flux
+++++++++
A heat flux condition,

.. math::
   -\kappa\nabla T\cdot\hat{n} = q_b \text{ on $\Gamma$},

is defined by setting `type`_ to **"flux"**. The heat flux :math:`q_b` is
specified using either `flux`_ for a constant value, or `flux_func`_ for a
function.

Oriented Heat Flux
++++++++++++++++++
An oriented heat flux condition,

.. math::
   -\kappa\nabla T\cdot\hat{n}=\alpha\,\vec{q_b}\cdot\hat{n}\text{ on $\Gamma$},

is defined by setting `type`_ to **"oriented-flux"**. The 3-vector heat flux
:math:`\vec{q_b}` is specified using either `vflux`_ for a constant value, or
`vflux_func`_ for a function, and the absorptivity coefficient :math:`\alpha`
is specified using either `absorptivity`_ for a constant value, or
`absorptivity_func`_ for a function.

Heat Transfer
+++++++++++++
An external heat transfer flux condition,

.. math::
   -\kappa\nabla T\cdot\hat{n} = \alpha(T - T_{\infty}) \text{ on $\Gamma$},

is defined by setting `type`_ to **"htc"**. The heat transfer coefficient
:math:`\alpha` is specified using either `htc`_ for a constant value, or
`htc_func`_ for a function, and the ambient temperature :math:`T_{\infty}`
is specified using either `ambient_temp`_ for a constant value, or
`ambient_temp_func`_ for a function.

Ambient Radiation
+++++++++++++++++
A simple ambient thermal radiation condition,

.. math::
   -\kappa\nabla T\cdot\hat{n} = \epsilon\sigma\bigl((T-T_0)^4 -
     (T_{\infty}-T_0)^4\bigr) \text{ on $\Gamma$},

is defined by setting `type`_ to **"radiation"**. The emissivity
:math:`\epsilon` is specified using either `emissivity`_ for a constant value
or `emissivity_func`_ for a function, and the temperature of the ambient
environment :math:`T_{\infty}` is specified using either `ambient_temp`_ for
a constant value, or `ambient_temp_func`_ for a function. Here :math:`\sigma`
is the Stefan-Boltzmann constant and :math:`T_0` is the absolute zero
temperature. They must be redefined if the problem units differ from
the default SI units using the :ref:`Stefan_Boltzmann<PhyCo_SB>` and
:ref:`Absolute_Zero<PhyCo_AZ>` variables of the
:ref:`PHYSICAL_CONSTANTS<PHYSICAL_CONSTANTS_Namelist>` namelist.


Conditions on Internal Interfaces
---------------------------------
Internal interfaces are merely coincident pairs of conforming external mesh
boundaries. These are modifications to the mesh created by Truchas and are
defined using the :ref:`Interface_Side_Sets<M_ISS>` variable from the
:ref:`MESH<MESH_Namelist>` namelist. Only the face set IDs referenced there
can be used in the definition of the following interface conditions.

Interface Heat Transfer
+++++++++++++++++++++++
An interface heat transfer condition models heat transfer across an imperfect
contact between two bodies or across a thin subscale material layer lying along
an interface :math:`\Gamma`. It imposes continuity of the heat flux
:math:`-\kappa\nabla T\cdot\hat{n}` across the interface :math:`\Gamma` and
gives this flux as

.. math::

   -\kappa\nabla T\cdot\hat{n} = -\alpha [T] \text{ on $\Gamma$},

where :math:`[T]` is the jump in :math:`T` across :math:`\Gamma` in the
direction :math:`\hat{n}`. It is defined by setting `type`_
to **"interface-htc"**. The heat transfer coefficient :math:`\alpha` is
specified using either `htc`_ for a constant value, or
`htc_func`_ for a function.

Gap Radiation
+++++++++++++
A gap radiation condition models radiative heat transfer
across a thin open gap lying along an interface :math:`\Gamma`. It imposes
continuity of the heat flux :math:`-\kappa\nabla T\cdot\hat{n}` across
:math:`\Gamma` and gives the flux as

.. math::
   -\kappa\nabla T\cdot\hat{n} = \epsilon_{\Gamma}\,\sigma\bigl(
  (T_{-}-T_0)^4-(T_{+}-T_0)^4\bigr) \text{ on $\Gamma$},

where :math:`T_{-}` and :math:`T_{+}` denote the values of :math:`T` on the
inside and outside gap surfaces with respect to the normal :math:`\hat{n}`
to :math:`\Gamma`. It is defined by setting `type`_ to **"gap-radiation"**.
The gap emissivity :math:`\epsilon_{\Gamma}` is specified using either
`emissivity`_ for a constant value, or `emissivity_func`_ for a function.
The effective gap emissivity :math:`\epsilon_{\Gamma}` depends on the
emissivities :math:`\epsilon_{-}` and :math:`\epsilon_{+}` of the surfaces
on either side of the gap and is given by

.. math::
   \epsilon_{\Gamma} = \frac{\epsilon_{-}\epsilon_{+}}{\epsilon_{-} +
       \epsilon_{+} - \epsilon_{-}\epsilon_{+}}.

The value of the Stefan-Boltzmann constant :math:`\sigma` and the absolute
zero temperature :math:`T_0` must be redefined if the problem units differ
from the default SI units using the :ref:`Stefan_Boltzmann<PhyCo_SB>` and
:ref:`Absolute_Zero<PhyCo_AZ>` variables of the
:ref:`PHYSICAL_CONSTANTS<PHYSICAL_CONSTANTS_Namelist>` namelist.


Usage Limitations and Requirements
----------------------------------
Use of the THERMAL_BC namelist is subject to a number of limitations:

* Instances of the same `type`_ may not *overlap* (i.e., apply to intersecting
  portions of the boundary), with the exception of **"oriented-flux"** where
  the net heat flux is the sum of the fluxes from the individual instances.

* A **"temperature"** type instance may not overlap with instances of any
  other `type`_.

* Instances with different external flux types (**"flux"**, **"htc"**, etc.)
  are allowed to overlap. The net flux will be the sum of the individual fluxes.

* Instances with different interface flux types (**"interface-htc"**,
  **"gap-radiation"**) are allowed to overlap. The net flux will be the sum of
  the individual fluxes.

* Interface flux instances may not overlap with external flux instances.

It is generally required that the specified thermal boundary conditions
completely cover the domain boundary; *there are no default thermal
boundary conditions.* However when enclosure radiation systems are present,
no boundary condition need be imposed on any part of the boundary that
belongs to an enclosure because the enclosure radiation system defines the
heat flux. Nevertheless, either **"temperature"** or flux-type conditions may
still be imposed there, and in the latter case the net heat flux is the sum
of boundary condition fluxes and radiative (from enclosure radiation) fluxes.

Namelist Variables
------------------

.. contents::
   :local:

name
+++++++++++++++++++++++++++++++++
A unique name used to identify a particular instance of this namelist.

:Type: string
:Default: none

.. _tbc-fsi:

face_set_ids
+++++++++++++++++++++++++++++++++
A list of face set IDs that define the portion of the boundary where the
boundary condition will be imposed.

:Type: integer list (32 max)
:Default: none


type
+++++++++++++++++++++++++++++++++
The type of boundary condition.

:Type: string
:Default: none
:Valid values:

.. csv-table::
   :header: "Value", "Description", "Associated variables"
   :class: tight-table
   :widths: 1 2 2

   **"temperature"**, "Temperature Dirichlet condition", "`temp`_, `temp_func`_"
   **"flux"**, "Heat flux condition", "`flux`_, `flux_func`_"
   **"oriented-flux"**, "Oriented heat flux condition", "`vflux`_, `vflux_func`_"
   **"htc"**, "External heat transfer condition", "`htc`_, `htc_func`_, `ambient_temp`_, `ambient_temp_func`_"
   **"radiation"**,"Simple thermal radiation to ambient environment", "`emissivity`_, `emissivity_func`_, `ambient_temp`_, `ambient_temp_func`_"
   **"interface-htc"**, "Internal interface heat transfer condition", "`htc`_, `htc_func`_"
   **"gap-radiation"**, "A gap thermal radiation condition", "`emissivity`_, `emissivity_func`_"


temp
+++++++++++++++++++++++++++++++++
The constant value of boundary temperature for a
:ref:`"temperature"<thermal_bc_namelist/index:temperature>` boundary
condition. To specify a function, use `temp_func`_ instead.

:Type: real
:Default: none


temp_func
+++++++++++++++++++++++++++++++++
The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function
that gives the boundary temperature for a
:ref:`"temperature"<thermal_bc_namelist/index:temperature>` boundary condition.
The function is expected to be a function of :math:`(t,x,y,z)`.

:Type: string
:Default: none


flux
+++++++++++++++++++++++++++++++++
The constant value of the outward heat flux for a
:ref:`"flux"<thermal_bc_namelist/index:heat flux>` boundary
condition. To specify a function, use `flux_func`_ instead.

:Type: real
:Default: none


flux_func
+++++++++++++++++++++++++++++++++
The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function
that gives the outward boundary heat flux for a
:ref:`"flux"<thermal_bc_namelist/index:heat flux>` boundary condition.
The function is expected to be a function of :math:`(t,x,y,z)`.

:Type: string
:Default: none


vflux
+++++++++++++++++++++++++++++++++
The constant value of the heat flux for an
:ref:`"oriented-flux"<thermal_bc_namelist/index:oriented heat flux>` boundary
condition. To specify a function, use `vflux_func`_ instead.

:Type: real 3-vector
:Default: none


vflux_func
+++++++++++++++++++++++++++++++++
The name of a :ref:`VFUNCTION<VFUNCTION_Namelist>` namelist defining a function
that gives the heat flux for an
:ref:`"oriented-flux"<thermal_bc_namelist/index:oriented heat flux>` boundary
condition. The function is expected to be a 3-vector function of
:math:`(t,x,y,z)`.

:Type: string
:Default: none


htc
+++++++++++++++++++++++++++++++++
The constant value of the heat transfer coefficient for
:ref:`"htc"<thermal_bc_namelist/index:heat transfer>` or
:ref:`"interface-htc"<thermal_bc_namelist/index:interface heat transfer>` type
boundary conditions. To specify a function, use `htc_func`_ instead.

:Type: real
:Default: none


htc_func
+++++++++++++++++++++++++++++++++
The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function
that gives the heat transfer coefficient for
:ref:`"htc"<thermal_bc_namelist/index:heat transfer>` or
:ref:`"interface-htc"<thermal_bc_namelist/index:interface heat transfer>` type
boundary conditions. The function is expected to be a function of
:math:`(t,x,y,z)` for an "htc" condition, and a function of
:math:`(T,t,x,y,z)` for an "interface-htc" condition. In the latter case
:math:`T` is taken to be the maximum of the two temperatures on either side of
the interface.

:Type: string
:Default: none


ambient_temp
+++++++++++++++++++++++++++++++++
The constant value of the ambient temperature for
:ref:`"htc"<thermal_bc_namelist/index:heat transfer>` or
:ref:`"radiation"<thermal_bc_namelist/index:ambient radiation>` type boundary
conditions. To specify a function, use `ambient_temp_func`_ instead.

:Type: real
:Default: none


ambient_temp_func
+++++++++++++++++++++++++++++++++
The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function
that gives the ambient temperature for
:ref:`"htc"<thermal_bc_namelist/index:heat transfer>` or
:ref:`"radiation"<thermal_bc_namelist/index:ambient radiation>` type boundary
boundary conditions. The function is expected to be a function of
:math:`(t,x,y,z)`.

:Type: string
:Default: none


absorptivity
+++++++++++++++++++++++++++++++++
The constant value of absorptivity for the
:ref:`"oriented-flux"<thermal_bc_namelist/index:oriented heat flux>` type
boundary condition. To specify a function, use `absorptivity_func`_ instead.

:Type: real
:Default: none
:Valid values: :math:`[0.0,1.0]`


absorptivity_func
+++++++++++++++++++++++++++++++++
The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function
that gives the absorptivity for the
:ref:`"oriented-flux"<thermal_bc_namelist/index:oriented heat flux>` type
boundary condition. This is expected to be a function of temperature :math:`T`
alone.

:Type: string
:Default: none


emissivity
+++++++++++++++++++++++++++++++++
The constant value of emissivity for a
:ref:`"radiation"<thermal_bc_namelist/index:ambient radiation>` or
:ref:`"gap-radiation"<thermal_bc_namelist/index:gap radiation>` type
boundary condition. To specify a function, use `emissivity_func`_ instead.

:Type: real
:Default: none
:Valid values: :math:`[0.0,1.0]`


emissivity_func
+++++++++++++++++++++++++++++++++
The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function
that gives the emissivity for a
:ref:`"radiation"<thermal_bc_namelist/index:ambient radiation>` or
:ref:`"gap-radiation"<thermal_bc_namelist/index:gap radiation>` type
boundary condition. The function is expected to be a function of
:math:`(T,t,x,y,z).` For a "gap-radiation" condition :math:`T` is taken to be
the maximum of the two temperatures on either side of the interface.

:Type: string
:Default: none
