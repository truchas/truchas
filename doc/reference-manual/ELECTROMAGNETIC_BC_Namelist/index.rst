.. _ELECTROMAGNETIC_BC_Namelist:

.. toctree::
   :maxdepth: 1

ELECTROMAGNETIC_BC Namelist
===========================

The ELECTROMAGNETIC_BC namelist is used to define boundary conditions for the
electromagnetic (EM) models at external boundaries. Each instance of the
namelist defines a particular condition to impose over a subset of the domain
boundary. The boundary subset :math:`\Gamma` is specified using mesh face sets.
The namelist variable `face_set_ids`_ takes a list of face set IDs, and the
boundary condition is applied at all faces belonging to those face sets. Note
that the mesh in view here is the one specified by the :ref:`EM_MESH
<EM_Mesh_Namelist>` namelist. Also note that ExodusII mesh side sets are
imported into Truchas as face sets with the same IDs.

.. admonishment:: Namelist Usage

   :Required/Optional: Required
   :Single/Multiple Instances: Multiple


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

:Type: integer list (32 max)
:Default: none

type
^^^^^^^^^^^^^^^^^^
The type of boundary condition. The available values are:

.. list-table::
   :class: tight-table
   :header-rows: 1
   :widths: 1 5

   * - Value
     - Description
   * - **"pec"**
     - Perfect electric conductor/no tangential E-field:
       :math:`\hat{n}\times\vec{E}=0` on :math:`\Gamma`.
   * - **"ih-hfield"**
     - A tangential H-field condition:
       :math:`\hat{n}\times\vec{H} = \hat{n}\times\vec{H}_\text{ext}` on
       :math:`\Gamma`, where :math:`\vec{H}_\text{ext}` is the external
       magnetic field defined by the
       :ref:`INDUCTION_SOURCE_FIELD<INDUCTION_SOURCE_FIELD_Namelist>` namelist.

:Notes: The default boundary condition is :math:`\hat{n}\times\vec{H}=0`.
