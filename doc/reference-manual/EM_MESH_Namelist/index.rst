.. _EM_MESH_Namelist:

.. toctree::
   :maxdepth: 1

EM_MESH Namelist
==================
The EM_MESH namelist specifies the mesh used by the induction heating solver.
This is a 3D tetrahedral mesh imported from an ExodusII format disk file.

.. admonition:: Namelist Usage

   :Required/Optional: Required when :ref:`induction_heating<physics-ih>` or
      :ref:`microwave_heating<physics-mwh>` PHYSICS options are enabled.
   :Single/Multiple Instances: Single

Namelist Variables
------------------

.. contents::
   :local:

.. _AM_AF:

mesh_file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Specifies the path to the ExodusII mesh file. If not an absolute path, it
will be interpreted as a path relative the the Truchas input file directory.

:Type: case-sensitive string
:Default: none

.. _AM_ACSF:

coord_scale_factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
An optional factor by which to scale all mesh node coordinates.

:Type: real
:Default: 1.0
:Valid Values: > 0

.. _AM_RA:

rotation_angles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
An optional list of 3 angles, given in degrees, specifying a counter-clockwise
rotation about the x, y, and z-axes to apply to the mesh. The rotations are
done sequentially in that order. A negative angle is a clockwise rotation, and
a zero angle naturally implies no rotation.

:Type: real 3-vector
:Default: (0.0, 0.0, 0.0)

.. _AM_P:

partitioner
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The partitioning method used to generate the parallel decomposition of the mesh.

:Type: case-insensitive string
:Default: "metis"
:Valid Values: "metis", "file", "block"
:Notes: See the :ref:`MESH<MESH_Namelist>` namelist variable
        :ref:`partitioner <M_P>` for a description of the options and their associated
        input variables. In particular, the METIS options listed there may also be
        specified in this namelist.

.. _AM_PF:

partition_file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Specifies the path to the mesh cell partition file, and is required when
:ref:`partitioner<AM_P>` is "file". If not an absolute path, it will be
interpreted as a path relative to the Truchas input file directory.

:Type: case-sensitive string
:Default: none
:Note: See the Notes for the :ref:`MESH<MESH_Namelist>` namelist variable
       :ref:`partition_file <M_PF>`

.. _AM_FiPa:

first_partition 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Specifies the number given the first partition in the numbering convention
used in the partition file. Either 0-based or 1-based numbering is allowed.

:Type: integer
:Default: 0
:Valid Values: 0 or 1
