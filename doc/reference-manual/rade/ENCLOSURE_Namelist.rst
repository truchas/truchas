ENCLOSURE Namelist
==================

The **ENCLOSURE** namelist specifies the mesh used by genre, along with the side
sets where view factors are to be computed.

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

List of side sets specifying the surface. All surface faces must belong to the
boundary of the mesh less any element blocks specified by `ignore_block_IDs`.
By default it is a fatal error if not satisfied; see `warn_non_boundary`.

:Type: integer list
:Default: none


warn_non_boundary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If any surface faces specified by `side_set_ids` are non-boundary faces, this
normally fatal error is reduced to a warning when this parameter is true. Such
faces are *not included in the surface*, and the number of such faces for each
side set are written to the terminal.

:Type: logical
:Default: F


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

List of side sets whose surfaces will be displaced by the vector specified by
`displacement`_.

:Type: integer list
:Default: none

.. note::
   Displacing surfaces is optional. Any displaced surfaces must be totally
   disconnected from the remaining enclosure surfaces (this is checked). The
   constant displacement must be compatibile with the specified symmetries (this
   is not checked). This is useful when the mesh consists of disconnected parts
   that will be shifted relative to one another. Instead of generating multiple
   3D meshes, a single 3D mesh plus use of this displacement option is all that
   is needed.


displacement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(x,y,z) displacement vector.

:Type: real 3-vector
:Default: none
