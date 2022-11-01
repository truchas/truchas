.. _ALTMESH_Namelist:

.. toctree::
   :maxdepth: 1

ALTMESH Namelist
==================

Overview
----------
The **ALTMESH** namelist specifies the alternate mesh used by the induction heating solver. This is a 3D tetrahedral mesh imported from an ExodusII format disk file.

ALTMESH Namelist Features
---------------------------
| **Required/Optional        :** Required when :ref:`Electromagnetics<ELECTROMAGNETICS_Namelist>` is true.
| **Single/Multiple Instances:** Single

Components
------------
* :ref:`Altmesh_Coordinate_Scale_Factor <AM_ACSF>`
* :ref:`Altmesh_File <AM_AF>`
* :ref:`First_Partition <AM_FiPa>`
* :ref:`Grid_Transfer_File <AM_GTF>`
* :ref:`Partitioner <AM_P>`
* :ref:`Partition_File <AM_PF>`
* :ref:`rotation_angles <AM_RA>`
* :ref:`data_mapper_kind <AM_DMK>`

.. _AM_ACSF:

Altmesh_Coordinate_Scale_Factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : An optional factor by which to scale all mesh node coordinates.
| **Type**        : real
| **Default**     : 1.0
| **Valid Values**: > 0

.. _AM_AF:

Altmesh_File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Specifies the path to the ExodusII mesh file. If not an absolute path, it will be interpreted relative the the Truchas input file directory.
| **Type**        : case-sensitive string
| **Default**     : none

.. _AM_FiPa:

First_Partition 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Specifies the number given the first partition in the numbering convention used in the partition file. Either 0-based or 1-based numbering is allowed.
| **Type**        : integer
| **Default**     : 0
| **Valid Values**: 0 or 1

.. _AM_GTF:

Grid_Transfer_File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Certain fields must be mapped between the main and alternative meshes during the course of a simulation. These mappings are accomplished using some fixed grid-mapping data that depend only on the two meshes. This optional variable specifies the path of a file containing this grid mapping data. If specified, and if the file exists, it will be read and its mapping data checked to ensure that it corresponds to the two meshes being used. If it corresponds, the data will be used for the calculation. Otherwise, the grid mapping data is computed and written to the file **altmesh_mapping_data.bin** in the output directory for use in future calculations, avoiding a needless and potentially costly recomputation of the same data.
| **Type**        : case-sensitive string
| **Default**     : none
| **Note**        : The mapping data depends on the internal ordering of the nodes and cells of each mesh, in addition to the meshes themselves. Thus it is recommended that this file be named in such a way that reflects the identity of the two meshes and the number of processors used to compute the mapping data; mapping data computed with one number of processors will not be usable in a calculation with a different number of processors, even when the same pair of meshes is used.

.. _AM_P:

Partitioner
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The partitioning method used to generate the parallel decomposition of the EM mesh.
| **Type**        : case-insensitive string
| **Default**     : "metis"
| **Valid Values**: "chaco", "metis", "file", "block"
| **Notes**       : See the :ref:`MESH<MESH_Namelist>` namelist variable :ref:`Partitioner <M_P>` for a description of the options and their associated input variables.

.. _AM_PF:

Partition_File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Specifies the path to the EM mesh cell partition file, and is required when :ref:`Partitioner<AM_P>` is "file". If not an absolute path, it will be interpreted as a path relative to the Truchas input file directory.
| **Type**        : case-sensitive string
| **Default**     : none
| **Notes**       : See the Notes for the :ref:`MESH<MESH_Namelist>` namelist variable :ref:`Partition_File <M_PF>`

.. _AM_FP:

First_Partitioner
^^^^^^^^^^^^^^^^^^^
| **Description** : Specifies the number given the first partition in the numbering convention used in the partition file. Either 0-based or 1-based numbering is allowed.
| **Type**        : integer
| **Default**     : 0
| **Valid Values**: 0 or 1

.. _AM_RA:

rotation_angles
^^^^^^^^^^^^^^^^^^^
| **Description** :  list of 3 angles, given in degrees, specifying the amount of counter-clockwise rotation about the x, y, and z-axes to apply to the mesh. The rotations are done sequentially in that order. A negative angle is a clockwise rotation, and a zero angle naturally implies no rotation.
| **Type**        : real 3-vector
| **Default**     : (0.0, 0.0, 0.0)

.. _AM_DMK:

data_mapper_kind (Experimental)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : This specifies the tool that will be used to map fields between the main heat transfer mesh and this alternative mesh. If "portage" is selected, an experimental data mapper based on the Portage toolkit, https://laristra.github.io/portage will be used. This data mapper is capable of handling main meshes containing prism and pyramid cells, but it has not yet been thoroughly vetted. Otherwise the normal data mapping tool will be used by default.
| **Type**        : string
| **Default**     : "default"
| **Valid Values**: "default", "portage"

