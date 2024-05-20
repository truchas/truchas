.. _FLOW_BC_Namelist:

.. toctree::
   :maxdepth: 1

FLOW_BC Namelist
==============================

The FLOW_BC namelist is used to define boundary conditions for the fluid flow
model at external boundaries. At inflow boundaries it also specifies the value
of certain intensive material quantities, like temperature, that may be
associated with other physics models.

Each instance of the namelist defines a particular condition to impose over
a subset :math:`\Gamma` of the domain boundary. The boundary subset
:math:`\Gamma` is specified using mesh face sets. The namelist variable
`face_set_ids`_ takes a list of face set IDs, and the boundary condition is
imposed on all faces belonging to those face sets. Note that ExodusII mesh
side sets are imported into Truchas as face sets with the same IDs.


.. admonition:: Namelist Usage

   :Required/Optional: Required
   :Single/Multiple Instances: Multiple


Boundary Conditions
-------------------
The following common types of boundary conditions can be defined:

* **Pressure**. A pressure Dirichlet condition :math:`p=p_b` on :math:`\Gamma`
  is defined by setting `type`_ to `"pressure"`. The boundary value :math:`p_b`
  is specified using either `pressure`_ for a constant value or
  `pressure_func`_ for a function.

* **Velocity**. A velocity Dirichlet condition :math:`u=u_b` on :math:`\Gamma`
  is defined by setting `type`_ to `"velocity"`. The boundary value :math:`u_b`
  is specified using either `velocity`_ for a constant value, or
  `velocity_func`_ for a function.

* **No Slip**. The special velocity Dirichlet condition :math:`u = 0` on
  :math:`\Gamma` is defined by setting `type`_ to `"no-slip"`.

* **Free slip**. A free-slip condition where fluid is not permitted to
  penetrate the boundary, :math:`\hat{n}\cdot u = 0` on :math:`\Gamma`, where
  :math:`\hat{n}` is the unit normal to :math:`\Gamma`, but is otherwise free
  to slide along the boundary (no tangential traction) is defined by setting
  `type`_ to `"free-slip"`. This is the default boundary condition.

* **Tangential surface tension**

These boundary condition types are mutually exclusive: namely, no two types
may be defined on overlapping subsets of the boundary. Any subset of the
boundary not explicitly assigned a boundary condition will be implicitly
assigned a free-slip condition.

Currently it is only possible to assign boundary conditions on the external
mesh boundary. However in many multiphysics applications the boundary of the
fluid flow domain will not coincide with the boundary of the larger problem
mesh. In some cases the boundary will coincide with an internal mesh-conforming
interface that separates fluid cells and solid (non-fluid) cells, where a
boundary condition could conceivably be assigned. In other cases, typically
those involving phase change, the boundary is only implicit, passing through
mixed fluid/solid cells, and will not conform to the mesh. In either case,
the flow algorithm aims to impose an effective no-slip condition for viscous
flows, or a free-slip condition for inviscid flows. A possible modeling
approach in the former mesh-conforming case is to define an internal mesh
interface using the :ref:`MESH<MESH_Namelist>` namelist variable
:ref:`interface_side_sets<M_ISS>`. This effectively creates new external
mesh boundary where flow boundary conditions can be assigned.


Namelist Variables
------------------

.. contents::
   :local:


name
^^^^^^^^^^^^^^^^^^
A unique name used to identify a particular instance of this namelist.

:Type: string
:Default: none

face_set_ids
^^^^^^^^^^^^^^^^^^
A list of face set IDs that define the portion of the boundary :math:`\Gamma`
where the boundary condition will be imposed.

type
^^^^^^^^^^^^^^^^^^
The type of boundary condition. The available values are:

.. list-table::
   :class: tight-table
   :header-rows: 1
   :widths: 1 5 2

   * - Value
     - Description
     - Associated Variables
   * - **"pressure"**
     - Pressure is prescribed on the boundary.
     - `pressure`_ or `pressure_func`_
   * - **"velocity"**
     - Velocity is prescribed on the boundary.
     - `velocity`_ or `velocity_func`_
   * - **"no-slip"**
     - Zero velocity is imposed on the boundary.
       This is incompatible with inviscid flow.
     - (none)
   * - **"free-slip"**
     - Zero velocity normal to the boundary, but the tangential velocity
       is otherwise free (no traction forces).
     - (none)
   * - **"marangoni"**
     - Like "free-slip" except a tangential traction is applied that is due to
       temperature dependence of surface tension. This is incompatible with
       inviscid flow.
     - `dsigma`_


pressure
^^^^^^^^^
The constant value of boundary pressure for a pressure-type boundary condition.
To specify a function, use `pressure_func`_ instead.

:Default: none
:Type: real


pressure_func
^^^^^^^^^^^^^^^
The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function
that gives the boundary pressure for a pressure-type boundary condition. The
function is expected to be a function of :math:`(t,x,y,z)`.

:Default: none
:Type: string


velocity
^^^^^^^^^
The constant value of boundary velocity for a velocity-type boundary condition.
To specify a function,  use `velocity_func`_ instead.

:Default: none
:Type: real 3-vector


velocity_func
^^^^^^^^^^^^^^^
The name of a :ref:`VFUNCTION<VFUNCTION_Namelist>` namelist defining a function
that gives the boundary velocity for a velocity-type boundary condition. The
function is expected to be a function of :math:`(t,x,y,z)`.

:Default: none
:Type: string


dsigma
^^^^^^^
The constant value of :math:`d\sigma/dT` for the marangoni-type condition.
Here :math:`\sigma(T)` is the temperature-dependent surface tension
coefficient.

:Default: none
:Type: real


inflow_material
^^^^^^^^^^^^^^^^
Velocity and pressure boundary conditions may result in fluid flow into
the domain across the boundary. This optional variable specifies the name
of the fluid material to flux in. If not specified, materials are fluxed
into a cell through a boundary face in the same proportion as the material
volume fractions present in the cell. To specify a function, use
`inflow_material_func`_ instead.

:Type: string


inflow_material_func
^^^^^^^^^^^^^^^^^^^^
The name of a :ref:`VFUNCTION<VFUNCTION_Namelist>` namelist defining a
function that gives the inflow material volume fractions for an inflow
boundary condition. The function is expected to be a function of
:math:`(t,x,y,z)`, and return a vector of length equal to the number of
liquid phases, plus void (if present). The order of the output components
is the same as the input file ordering provided by :ref:`materials<PHYSICS_M>`.

:Type: string


inflow_temperature
^^^^^^^^^^^^^^^^^^^^^^
Velocity and pressure boundary conditions may result in fluid flow into the
domain across the boundary. This optional variable specifies the temperature
of the material fluxed in. If not specified, materials are fluxed into a call
through a boundary face at the same temperature as the cell.

:Type: real


inflow_conc
^^^^^^^^^^^^^^^^^^^^^^
Velocity and pressure boundary conditions may result in fluid flow into the
domain across the boundary. This optional variable specifies the concentration
of the passively advected scalar fluxed into a cell through a boundary face.
If not specified, the fluxed concentration is the same as that of the cell.
To specify a function, use `inflow_conc_func`_ instead. This is a vector-valued
variable with elements that correspond to the different species components of
the species advection-diffusion model.

:Type: real


inflow_conc_func
^^^^^^^^^^^^^^^^^^^^^^
The names of :ref:`FUNCTION<FUNCTION_Namelist>` namelists defining functions
that give the inflow concentrations for an inflow boundary condition. The
functions are expected to be functions of :math:`(t,x,y,z)`. This is a
vector-valued variable with elements that correspond to the different species
components of the species advection-diffusion model.

:Type: string
