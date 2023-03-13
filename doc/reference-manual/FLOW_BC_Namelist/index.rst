.. _FLOW_BC_Namelist:

.. toctree::
   :maxdepth: 1

FLOW_BC Namelist
==============================

Overview
------------
The :ref:`FLOW_BC<FLOW_BC_Namelist>` namelist is used to define boundary conditions for the fluid flow model at external boundaries. At inflow boundaries it also specifies the value of certain intensive material quantities, like temperature, that may be associated with other physics models.

Each instance of the namelist defines a particular condition to impose over a subset :math:`\Gamma` of the domain boundary. The boundary subset :math:`\Gamma` is specified using mesh face sets. The namelist variable :ref:`face_set_ids<FLOW_BC_FSI>` takes a list of face set IDs, and the boundary condition is imposed on all faces belonging to those face sets.Note that ExodusII mesh sides sets are imported into Truchas as face sets with the same IDs. The following common types of boundary conditions can be defined:

* **Pressure** A pressure Dirichlet condition :math:`p=p_b` on :math:`\Gamma` is defined by setting :ref:`type<FLOW_BC_Type>` to **pressure**. The boundary value :math:`p_b` is specified using either :ref:`pressure<FLOW_BC_P>` for a constant value, or :ref:`pressure_func<FLOW_BC_PF>` for a function.

* **Velocity** A velocity Dirichlet condition :math:`u=u_b` on :math:`\Gamma` is defined by setting :ref:`type<FLOW_BC_Type>` to **velocity**. The boundary value :math:`u_b` is specified using either :ref:`velocity<FLOW_BC_V>` for a constant value, or :ref:`velocity_func<FLOW_BC_VF>` for a function.

* **No Slip** The special velocity Dirichlet condition :math:`u = 0` on :math:`\Gamma` is defined by setting :ref:`type<FLOW_BC_Type>` to **no-slip**

* **Free slip** A free-slip condition where fluid is not permitted to penetrate the boundary,:math:`\hat{n}.u = 0` on :math:`\Gamma`, where :math:`\hat{n}` is the unit normal to :math:`\Gamma`, but is otherwise free to slide along the boundary (no tangential traction) is defined by setting :ref:`type<FLOW_BC_Type>` to **free-slip**.

* **Tangential surface tension**  These boundary condition types are mutually exclusive: namely, no two types may be defined onoverlapping subsets of the boundary. Any subset of the boundary not explicitly assigned a boundary condition will be implicitly assigned a free-slip condition.

Currently it is only possible to assign boundary conditions on the external mesh boundary. However in many multiphysics applications the boundary of the fluid flow domain will not coincide with the boundary of the larger problem mesh. In some cases the boundary will coincide with an internal mesh-conforming interface that separates fluid cells and solid (non-fluid) cells, where a boundary condition could conceivably be assigned. In other cases, typically those involving phase change, the boundary is only implicit, passing through mixed fluid/solid cells, and will not conform to the mesh. In either case, the flow algorithm aims to impose an effective no-slip condition for viscous flows, or a free-slip condition for inviscid flows. A possible modeling approach in the former mesh-conforming case is to define an internal mesh interface using the :ref:`MESH<MESH_Namelist>` namelist variable :ref:`interface_side_sets<M_ISS>`. This effectively creates new external mesh boundary where flow boundary conditions can be assigned.

FLOW_BC Namelist Features
---------------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Multiple

Components
-------------
* :ref:`name<FLOW_BC_N>`
* :ref:`face_set_ids<FLOW_BC_FSI>`
* :ref:`type<FLOW_BC_Type>`
* :ref:`pressure<FLOW_BC_P>`
* :ref:`pressure_func<FLOW_BC_PF>`
* :ref:`velocity<FLOW_BC_V>`
* :ref:`velocity_func<FLOW_BC_VF>`
* :ref:`dsigma<FLOW_BC_D>`
* :ref:`inflow_material<FLOW_BC_IM>`
* :ref:`inflow_material_func<FLOW_BC_IMF>`
* :ref:`inflow_temperature<FLOW_BC_IT>`

.. _FLOW_BC_N:

name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A unique name used to identify a particular instance of this namelist
| **Type**        : string (31 characters max)
| **Default**     : none

.. _FLOW_BC_FSI:

face_set_ids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A list of face set IDs that define the portion of the boundary where the boundary condition will be imposed.
| **Type**        : integer list (32 max)
| **Default**     : none

.. _FLOW_BC_Type:

type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The type of boundary condition. The available options are:

.. _flow_bc_options:

.. csv-table::
   :header: "Option", "Description"
   :class: tight-table
   :widths: 1 5

   "**pressure**", "Pressure is prescribed on the boundary. Use :ref:`pressure<FLOW_BC_P>` or :ref:`pressure_func<FLOW_BC_PF>` to specify its value."
   "**velocity**", "Velocity is prescribed on the boundary. Use :ref:`velocity<FLOW_BC_V>` or :ref:`velocity_func<FLOW_BC_VF>` to specify its value."
   "**no-slip**", "0-velocity is imposed on the boundary. This is incompatible with inviscid flow."
   "**free-slip**", "No velocity normal to the boundary, but the tangential velocity is otherwise free (notraction forces)."
   "**marangoni**", "Like **free-slip** except a tangential traction is applied that is due to temperature dependence of surface tension. Use :ref:`dsigma<FLOW_BC_D>` to specify the value of :math:`d\sigma/dT`. This is incompatible with inviscid flow."

| **Type**        : string
| **Default**     : none
| **Notes**       :The different boundary condition types are mutually exclusive; no two can be specified on a common portion of the boundary.

.. _FLOW_BC_P:

pressure
^^^^^^^^^
| **Description** : The constant value of boundary pressure for a pressure-type boundary condition. To specify a function, use :ref:`pressure_func<FLOW_BC_PF>` instead.
| **Default**     : none
| **Type**        : real

.. _FLOW_BC_PF:

pressure_func
^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function that gives the boundary pressure for a pressure-type boundary condition. The function is expected to be a function of :math:`(t,x,y,z)`.
| **Default**     : none
| **Type**        : string

.. _FLOW_BC_V:

velocity
^^^^^^^^^
| **Description** : The constant value of boundary velocity for a velocity-type boundary condition. To specify a function,  use :ref:`velocity_func<FLOW_BC_VF>` instead.
| **Default**     : none
| **Type**        : real 3-vector

.. _FLOW_BC_VF:

velocity_func
^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`VFUNCTION<VFUNCTION_Namelist>` namelist defining a function that gives the boundary velocity for a velocity-type boundary condition. The function is expected to be a function of :math:`(t,x,y,z)`.
| **Default**     : none
| **Type**        : string

.. _FLOW_BC_D:

dsigma
^^^^^^^
| **Description** : The constant value of :math:`d\sigma/dT` for the marangoni-type condition. Here :math:`\sigma(T)` is the temperature dependent surface tension coefficient.
| **Default**     : none
| **Type**        : real

.. _FLOW_BC_IM:

inflow_material
^^^^^^^^^^^^^^^^
| **Description** : Velocity and pressure boundary conditions may result in fluid flow into the domain across the boundary. This parameter specifies the name of the fluid material to flux in. If not specified, materials are fluxed into a cell through a boundary face in the same proportion as the material volume fractions present in the cell. To specify a function, use :ref:`inflow_material_func<FLOW_BC_IMF>` instead.
| **Default**     : none

.. _FLOW_BC_IMF:

inflow_material_func
^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`VFUNCTION<VFUNCTION_Namelist>` namelist defining a function that gives the inflow material volume fractions for an inflow boundary condition. The function is expected to be a function of :math:`(t,x,y,z)`, and return a vector of length equal to the number of liquid phases, plus void (if present). The order of the output components is the same as the input file ordering provided by :ref:`materials<PHYSICS_M>`.
| **Default**     : none
| **Type**        : string

.. _FLOW_BC_IT:

inflow_temperature
^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Velocity and pressure boundary conditions may result in fluid flow into the domain across the boundary. This parameter specifies the temperature of the material fluxed in. If not specified, materials are fluxed into a call through a boundary face at the same temperature as the cell.
| **Default**     : none
