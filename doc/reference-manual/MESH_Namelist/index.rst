.. _MESH_Namelist:

.. toctree::
   :maxdepth: 1

MESH Namelist
===============

Overview
++++++++
The :ref:`MESH<MESH_Namelist>` namelist specifies the common mesh used by all physics models other than the induction heating model, which uses a separate tetrahedral mesh specified by the :ref:`EM_MESH<EM_MESH_Namelist>` namelist. For simple demonstration problems, a rectilinear hexahedral mesh of a brick domain can be defined, but for most applications the mesh will need to be generated beforehand by some third party tool or tools and saved asa file that Truchas will read. At this time Exodus II :footcite:`sjaardema2006exodus` is the only supported mesh format (also sometimes known as Genesis). This well-known format is used by some mesh generation tools (`Cubit <https://cubit.sandia.gov/>`_, for example) and utilities exist for translating from other formats to Exodus II. The unstructured 3D mesh may be a general mixed-element mesh consisting of non-degenerate hexehedral, tetrahedral, pyramid, and wedge/prism elements. The Exodus II format supports a partitioning of the elements into element blocks. It also supports the definition of side sets, which are collections of oriented element faces that describe mesh surfaces, either internal or boundary, and node sets which are collections of mesh nodes. Extensive use is made of this additional mesh metadata in assigning materials, initial conditions, boundary conditions, etc., to the mesh.

MESH Namelist Features
----------------------------
| **Required/Optional        :** Required
| **Single/Multiple Instances:** Single


Components
++++++++++

.. contents::
   :local:


External Mesh file
----------------------
In typical usage, the mesh will be read from a specified Exodus II mesh file. Other input variables that follow specify optional modifications that can be made to the mesh after it is read.

.. _M_MF:

mesh_file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Specifies the path to the Exodus II mesh file. If not an absolute path, it will be interpreted as a path relative to the Truchas input file directory.
| **Type**        : case-sensitive string
| **Default**     : none 

.. _M_ISS:

Interface_Side_Sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A list of side set IDs from the ExodusII mesh identifying internal mesh surfaces that will be treated specially by the heat/species transport solver.
| **Type**        : integer list
| **Default**     : An empty list of side set IDs.
| **Valid Values**: Any side set ID whose faces are internal to the mesh.
| **Notes**       : The heat/species transport solver requires that boundary conditions are imposed along the specified surface. Typically these will be interface conditions defined by :ref:`THERMAL_BC<THERMAL_BC_Namelist>` namelists, but in unusual use cases they could also be external boundary conditions defined by the same namelists. In the latter case it is necessary to understand that the solver views the mesh as having been sliced open along the specified internal surfaces creating matching pairs of additional external boundary and, where interface conditions are not imposed, boundary conditions must be imposed on both sides of the interface.

.. _M_GEB:

gap_element_blocks (deprecated)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A list of element block IDs from an Exodus II mesh that are to be treated as gap elements.
| **Type**        : integer list
| **Default**     : An empty list of element block IDs.
| **Valid Values**: Any element block ID.
| **Notes**       : Any element block ID in the mesh file can be specified, but elements that are not connected such that they can function as gap elements or are not consistent with side set definitions will almost certainly result in incorrect behavior. The code does not check for these inconsistencies. The heat/species transport solver drops these elements from its view of the mesh and treats them instead as an internal interface; see the notes to :ref:`interface_side_sets<M_ISS>`. The block IDs specified here can be used as values for :ref:`face_set_ids<TB_FSI>` from the :ref:`THERMAL_BC<THERMAL_BC_Namelist>` namelist.

.. _M_EBM:

exodus_block_modulus
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : When importing an Exodus II mesh, the element block IDs are replaced by their value modulo this parameter. Set the parameter to 0 to disable this procedure.
| **Type**        : integer 
| **Default**     : 10000
| **Valid Values**: :math:`\geq 0`
| **Notes**       : This parameter helps solve a problem posed by mixed-element meshes created by Cubit and Trelis. In those tools a user may define an element block comprising multiple element types. But when exported in the Exodus II format, which doesn’t support blocks with mixed element types, the element block will be written as multiple Exodus II blocks, one for each type of element. One of the blocks will retain the user-specified ID of the original block. The IDs of the others will be that ID plus an offset specific to the element type. For example, if the original block ID was 1, hexahedra in the block will be written to a block with ID 1, tetrahedra to a block with ID :math:`10001`, pyramids to a block with ID 100001, and wedges to a block with ID :math:`200001`. These are the default offset values, and they can be set in Cubit/Trelis; see their documentation for details on how the IDs are generated. It is important to note that this reorganization of element blocks occurs silently and so the user may be unaware that it has happened. In order to reduce the potential for input errors, Truchas will by default convert the block IDs to congruent values modulo :math:`N` in the interval :math:`[1,N−1]` where :math:`N` is the value of this parameter. The default value :math:`10000` is appropriate for the default configuration of Cubit/Trellis, and restores the original user-specified block IDs. Note that this effectively limits the range of element block IDs to :math:`[1,N−1]`.

The element block IDs are modified immediately after reading the file. Any input parameters that refer to block IDs must refer to the modified IDs.

Internally Generated Mesh
----------------------------
A rectilinear hexahedral mesh for a brick domain :math:`[x_{min},x_{max}] × [y_{min},y_{max}] × [z_{min},z_{max}]` can be generated internally as part of a Truchas simulation using the following input variables. The mesh is the tensor product of 1D grids in each of the coordinate directions. Each coordinate grid is defined by a coarse grid whose intervals are subdivided into subintervals, optionally with biased sizes. The generated Exodus II mesh consists of a single element block with ID 1, and a side set is defined for each of the six sides of the domain with IDs 1 through 6 for the :math:`x=x_{min}, x=x_{max}, y=y_{min}, y=y_{max}, z=z_{min},` and :math:`z=z_{max}` sides, respectively. In addition, a different node set is defined for each of the 8 nodes at the corners of the domain, with IDs 1 through 8. The first node set is the :math:`(x_{min},y_{min},z_{min})` corner. It is followed by the remaining corners on the :math:`z=z_{min}` side in a counter-clockwise order with respect to the :math:`z` axis, and then the corners on the :math:`z=z_{max}` side in the analogous manner. Note that while the mesh is formally structured, it is represented internally as a general unstructured mesh.

.. _M_XYZA:

x_axis, y_axis, z_axis
^^^^^^^^^^^^^^^^^^^^^^^
Data that describes the grid in each of the coordinate directions. The tensor product of these grids define the nodes of the 3D mesh. The data for each coordinate grid consists of these three component arrays:

.. _grid_mesh_options:
.. csv-table:: 
   :header: "Option", "Description"
   :class: tight-table
   :widths: 1 5

   "**%coarse_grid**","A strictly increasing list of two or more real values that define the points of the coarse grid for the coordinate direction. The first and last values define the extent of the domain in this direction."
   "**%intervals**","A list of postive integers defining the number of subintervals into which each corresponding coarse grid interval should be subdivided. The number of values must be one less than the number of coarse grid points."
   "**%ratio**","An optional list of positive real values that define the ratio of the lengths of successive subintervals for each coarse grid interval. The default is to subdivide into equal length subintervals. If specified, the number of values must be one less than the number of coarse grid points."

See :numref:`Figure %s<fig_mesh_internal_schematic>` for an example. That is a mesh generated by the following input.

.. _internal_mesh_example_script:

::

   x_axis%coarse_grid = 0.0, 0.67, 1.33, 2.0
   x_axis%intervals   = 3, 1, 4
   x_axis%ratio       = 1.3, 1.0, 0.7
   y_axis%coarse_grid = 0.0, 0.5, 1.0
   y_axis%intervals   = 1, 3
   y_axis%ratio       = 0.7
   z_axis%coarse_grid = 0.0, 1.0
   z_axis%intervals   = 1

.. _fig_mesh_internal_schematic:
.. figure:: images/grid-plot.png
   :width: 650px
   :align: center
   
   Top xy surface of the rectilinear mesh generated by the :ref:`example input<internal_mesh_example_script>` shown.


.. _M_NF:

noise_factor (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : If specified with a positive value, the coordinates of each mesh node will be perturbed by uniformly distributed random amount whose magnitude will not exceed this value times the local cell size at the node. Nodes on the boundary are not perturbed in directions normal to the boundary. This is only useful for testing.
| **Default**     : 0
| **Valid Values**: :math:`\in [0,0.3]`

Common variables
--------------------

The following variables apply to both types of meshes.

.. _M_CSF:

coordinate_scale_factor 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : An optional factor by which to scale all mesh node coordinates.
| **Type**        : real
| **Default**     : 1.0
| **Valid Values**: :math:`\gt 0`

.. _M_RA:

rotation_angles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A list of 3 angles, given in degrees, specifying the amount of counter-clockwise rotation about the x, y, and z-axes to apply to the mesh. The rotations are done sequentially in that order. A negative angle is a clockwise rotation, and a zero angle naturally implies no rotation.
| **Type**        : real 3-vector
| **Default**     : :math:`(0.0,\:0.0,\:0.0)`


.. _M_P:

partitioner
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The partitioning method used to generate the parallel decomposition of the mesh.
| **Type**        : case-insensitive string
| **Default**     : "metis"
| **Valid Values**: "chaco", "metis", "file", "block"
| **Notes**       :

.. _partitioner_options:
.. csv-table:: 
   :header: "Option", "Description"
   :class: tight-table
   :widths: 1 5

   "**chaco**","uses a graph partitioning method from the Chaco library :footcite:`leland1995chaco` to compute the mesh decomposition at run time. Support for this legacy library may be removed in the future."
   "**metis**","uses the well-known METIS library :footcite:`karypis1998fast` to partition the dual graph of the mesh at runtime. This method has a number of options which are described below."
   "**file**","reads the partitioning of the mesh cells from a disk file; see :ref:`partition_file<M_PF>`."
   "**block**","partitions the mesh cells into nearly equal-sized blocks of consecutively numbered cells according their numbering in the mesh file. The quality of this naive decomposition entirely depends on the given ordering of mesh cells, and thus this option is not generally recommended."

.. _M_PF:

partition_file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Specifies the path to the mesh cell partition file, and is required when :ref:`partitioner<M_P>` is **file**. If not an absolute path, it will be interpreted as a path relative to the Truchas input file directory.
| **Type**        : case-sensitive string
| **Default**     : none
| **Notes**       : The format of this text file consists of a sequence of integer values, one value or multiple values per line. The first value is the partition number of the first cell, the second value the partition number of the second cell, and so forth. The number of values must equal the number of mesh cells. The file may use either a 0-based or 1-based numbering convention for the partitions. Popular mesh partitioning tools typically use 0-based partition numbering, and so the default is to assume 0-based numbering; use :ref:`first_partition<M_FP>` to specify 1-based numbering.

.. _M_FP:

first_partition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Specifies the number given the first partition in the numbering convention used in the partition file. Either 0-based or 1-based numbering is allowed.
| **Type**        : integer
| **Default**     : 0
| **Valid Values**: 0 or 1

METIS Mesh Partitioning
-------------------------
When **metis** is specified for **partitioner**, a graph partitioning procedure from the METIS library isused to partition the dual graph of the mesh. This is the graph whose nodes are the mesh cells and edges are the faces shared by cells. The partitioning procedures have the following integer-valued options that may be specified, though all have reasonable defaults so that none must be specified. See the METIS documentation :footcite:`karypis1998fast` for more details on these options.

metis_ptype
^^^^^^^^^^^^
| **Description** : Specifies the partitioning method. 
| **Default**     : 0
| **Valid Values**: 
| 0 - Multilevel recursive bisection
| 1 - Multilevel :math:`k`-way partitioning

metis_iptype
^^^^^^^^^^^^
| **Description** : specifies the algorithm used during initial partitioning (recursive bisection only). Possible values are:
| **Default**     : 0
| **Valid Values**:
| 0 - Grows a bisection using a greedy strategy
| 1 - Computes a bisection at random followed by a refinement

metis_ctype
^^^^^^^^^^^^
| **Description** : specifies the matching scheme to be used during coarsening. 
| **Default**     : 1
| **Valid Values**:
| 0 - Random matching
| 1 - Sorted heavy-edge matching

metis_ncuts
^^^^^^^^^^^^
| **Description** : specifies the number of different partitionings that will be computed. The final partitioning will be the one that achieves the best edgecut.
| **Default**     : 1

metis_niter
^^^^^^^^^^^^
| **Description** : specifies the number of iterations of the refinement algorithm at each stage of the uncoarsening process.
| **Default**     : 10

metis_ufactor
^^^^^^^^^^^^^^
| **Description** : specifies the maximum allowed load imbalance among the partitions. A value of :math:`n` indicates that the allowed load imbalance is :math:`(1 +n)/1000`. The default is :math:`1` for recursive bisection (i.e., an imbalance of :math:`1.001`) and the default value is :math:`30` for :math:`k`-way partitioning (i.e., an imbalance of :math:`1.03`).

metis_minconn
^^^^^^^^^^^^^^
| **Description** : specifies whether the partitioning procedure should seek to minimize the maximum degree of the subdomain graph (1); or not (0, default). The subdomain graph is the graph in which each partition is a node, and edges connect subdomains with a shared interface.
| **Default** : 0

metis_contig
^^^^^^^^^^^^^^
| **Description**  : specifies whether the partitioning procedure should produce partitions that are contiguous (1); or not (0, default). If the dual graph of the mesh is not connected this option is ignored.
| **Default**: :math:`0`

metis_seed
^^^^^^^^^^^
| **Description** : specifies the seed for the random number generator.

metis_dbglvl
^^^^^^^^^^^^^
| **Description** : specifies the amount and type of diagnostic information that will be written to **stderr** by the partitioning procedure. The default is :math:`0`, no output. Use :math:`1` to write some basic information. Refer to the METIS documentation for the many other possible values and the output they generate.
| **Default** : :math:`0`

                  
.. footbibliography::
