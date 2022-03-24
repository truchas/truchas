.. _SM_BC_Namelist:

.. toctree::
   :maxdepth: 1

SOLID_MECHANICS_BC Namelist
==============================

Overview
------------
The :ref:`SOLID_MECHANICS_BC<SM_BC_Namelist>` namelist is used to define boundary conditions for the solid mechanics model at external and internal boundaries. Each instance of the namelist defines a particular condition to impose over a subset :math:`\Gamma` of the domain boundary. The boundary subset :math:`\Gamma` is specified using mesh face sets. The namelist variable :ref:`face_set_ids<SM_BC_FSI>` takes a list of face set IDs, and the boundary condition is imposed on all faces belonging to those face sets. Note that ExodusII mesh sides sets are imported into Truchas as face sets with the same IDs. Displacement x/y/z boundary conditions may alternatively be applied to node sets with the namelist variable :ref:`node_set_ids<SM_BC_NSI>`.

External boundaries
^^^^^^^^^^^^^^^^^^^^^^
The following types of external boundary conditions can be defined. The outward unit normal to the boundary :math:`\Gamma` is denoted :math:`\hat{n}`. Note the solid mechanics package defaults to a traction-free surface with no displacement constraints at external boundaries.

* **Traction** An applied traction condition :math:`(\sigma \cdot \hat{n}) \cdot \hat{m}=\tau_b` for a given vector :math:`\hat{m}` on :math:`\Gamma` is defined by setting :ref:`type<SM_BC_Type>` to **traction-n**, **traction-x**, **traction-y**, or **traction-z**. The boundary value :math:`\tau_b` is specified using either :ref:`traction<SM_BC_T>` for a constant value, or :ref:`traction_func<SM_BC_TF>` for a function.

* **Displacement** A dirichlet displacement condition :math:`u=u_b` on :math:`\Gamma` is defined by setting :ref:`type<SM_BC_Type>` to **displacement-n**, **displacement-x**, **displacement-y**, or **displacement-z**. The boundary value :math:`\tau_b` is specified using either :ref:`displacement<SM_BC_D>` for a constant value, or :ref:`displacement_func<SM_BC_DF>` for a function.

Internal interfaces
^^^^^^^^^^^^^^^^^^^^^
Internal interfaces are merely coincident pairs of conforming external mesh boundaries. These are modifications to the mesh created by Truchas and are defined using the :ref:`Interface_Side_Sets<M_ISS>` parameter from the :ref:`MESH<MESH_Namelist>` namelist. Only the face set IDs referenced there can be used in the definition of the following interface conditions. The following types of internal interface conditions can be defined.

* **Gap Contact**. A gap contact condition models no-penetration sliding contact along internal surfaces :math:`\Gamma`. It is defined by setting :ref:`type<SM_BC_Type>` to **gap-contact**.

SOLID_MECHANICS_BC Namelist Features
------------------------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Multiple

Components
-------------
* :ref:`name<SM_BC_N>`
* :ref:`face_set_ids<SM_BC_FSI>`
* :ref:`node_set_ids<SM_BC_NSI>`
* :ref:`type<SM_BC_Type>`
* :ref:`displacement<SM_BC_D>`
* :ref:`displacement_func<SM_BC_DF>`
* :ref:`traction<SM_BC_T>`
* :ref:`traction_func<SM_BC_TF>`

.. _SM_BC_N:

name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A unique name used to identify a particular instance of this namelist
| **Type**        : string (31 characters max)
| **Default**     : none

.. _SM_BC_FSI:

face_set_ids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A list of face set IDs that define the portion of the boundary where the boundary condition will be imposed.
| **Type**        : integer list (32 max)
| **Default**     : none

.. _SM_BC_NSI:

node_set_ids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A list of node set IDs that define the portion of the boundary where the boundary condition will be imposed. It is used to remove the null-space from the solution which is present when a solid part is otherwise free to move or rotate uniformly through space. This option is only valid if :ref:`type<SM_BC_Type>` is set to "**displacement-x**", "**displacement-y**", "**displacement-z**".
| **Type**        : integer list (32 max)
| **Default**     : none

.. _SM_BC_Type:

type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The type of boundary condition. The available options are:

.. _flow_bc_options:

.. csv-table::
   :header: "Option", "Description"
   :class: tight-table
   :widths: 1 5

   "**displacement-n**", "Normal displacement is prescribed on the boundary. Use :ref:`displacement<SM_BC_D>` or :ref:`displacement_func<SM_BC_DF>` to specify its value."
   "**displacement-x**", "X-direction displacement is prescribed on the boundary. Use :ref:`displacement<SM_BC_D>` or :ref:`displacement_func<SM_BC_DF>` to specify its value."
   "**displacement-y**", "Y-direction displacement is prescribed on the boundary. Use :ref:`displacement<SM_BC_D>` or :ref:`displacement_func<SM_BC_DF>` to specify its value."
   "**displacement-z**", "Z-direction displacement is prescribed on the boundary. Use :ref:`displacement<SM_BC_D>` or :ref:`displacement_func<SM_BC_DF>` to specify its value."
   "**traction-n**", "Normal traction is applied on the boundary. Use :ref:`traction<SM_BC_T>` or :ref:`traction_func<SM_BC_DF>` to specify its value."
   "**traction-x**", "X-direction traction is applied on the boundary. Use :ref:`traction<SM_BC_T>` or :ref:`traction_func<SM_BC_DF>` to specify its value."
   "**traction-y**", "Y-direction traction is applied on the boundary. Use :ref:`traction<SM_BC_T>` or :ref:`traction_func<SM_BC_DF>` to specify its value."
   "**traction-z**", "Z-direction traction is applied on the boundary. Use :ref:`traction<SM_BC_T>` or :ref:`traction_func<SM_BC_DF>` to specify its value."
   "**gap-contact**", "An internal interface sliding contact condition."

| **Type**        : string
| **Default**     : none
| **Notes**       : X/Y/Z variants can be applied together, e.g. **displacement-x**, **displacement-y**, and **displacement-z** to zero for a no-displacement boundary. Orthogonal conditions of different types may also be applied on the same surfaces, e.g. **traction-z** and **displacement-y**. External boundaries with unspecified BCs are given a zero-traction constraint with free displacement.

.. _SM_BC_D:

displacement
^^^^^^^^^^^^^^^^
| **Description** : The constant value of boundary displacement for a displacement-type boundary condition. To specify a function, use :ref:`displacement_func<SM_BC_DF>` instead.
| **Default**     : none
| **Type**        : real

.. _SM_BC_DF:

displacement_func
^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function that gives the boundary displacement for a displacement-type boundary condition. The function is expected to be a function of :math:`(t,x,y,z)`.
| **Default**     : none
| **Type**        : string

.. _SM_BC_T:

traction
^^^^^^^^^
| **Description** : The constant value of boundary traction for a traction-type boundary condition. To specify a function, use :ref:`traction_func<SM_BC_TF>` instead.
| **Default**     : none
| **Type**        : real

.. _SM_BC_TF:

traction_func
^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist defining a function that gives the boundary traction for a traction-type boundary condition. The function is expected to be a function of :math:`(t,x,y,z)`.
| **Default**     : none
| **Type**        : string
