ENCLOSURE Namelist
==================

The **ENCLOSURE** namelist specifies the mesh used by genre, along with the side sets where view factors are to be computed.

:Required/Optional: Required
:Single/Multiple Instances: Single

.. contents:: Components
   :local:


name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A user-supplied name for the enclosure.

:Type: string
:Default: ``"RadE"``


mesh_file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Path to the Exodus/Genesis mesh file.

:Type: case-sensitive string
:Default: none


ignore_block_IDs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

List of element blocks to mask-off from the mesh.

:Type: integer list
:Default: none


side_set_ids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

List of side sets specifying the surface.

:Type: integer list
:Default: none


symmetries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

List of up to 3 symmetry operations

:Type: string list (3 max)
:Default: none
:Valid Values:
   - ``"Mirror<a>"``, a = X, Y, Z; e.g. ``"MirrorZ"``
   - ``"Rot<a><n>"``, a = X, Y, Z; n integer; e.g. ``"RotZ3"``, ``"RotX16"``


displace_side_set_ids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

List of side sets whose surfaces will be displaced.

:Type: integer list
:Default: none


displacement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(x,y,z) displacement vector.

:Type: real 3-vector
:Default: none
