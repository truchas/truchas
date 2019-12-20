.. sectionauthor:: David Neill-Asanza <dhna@lanl.gov>

PATCHES Namelist
================

.. code-block:: console

  &PATCHES
    patch_algorithm = 'PAVE'
    verbosity_level = 3
    max_angle = 30.0
    pave_split_patch_size = 4
  /

:superscript:`Example PATCHES namelist`

The `PATCHES` namelist defines the parameters used by the patching algorithms. The namelist supports
ten parameters, but not all parameters are used by all algorithms. Parameters only used by a
particular algorithm are prefixed with the algorithm's name.

These parameters are described in detail below.

.. contents:: PATCHES namelist parameters
   :local:
   :backlinks: none


General Parameters
------------------

PATCH_ALGORITHM
+++++++++++++++
Selects one of the three available algorithms, or disables patching.

   **Type:** ``STRING``

   **Domain:** Must be one of ``'NONE'``, ``'PAVE'``, ``'VAC'``, or ``'VSA'``

   **Default:** `patch_algorithm =` ``'PAVE'``

Each option selects a different patch algorithm:

#. **NONE:** No patches will be generated. All other parameters are ignored. This is equivalent to
   an absent `PATCHES` namelist.
#. **PAVE:** Generate patches with the :doc:`PAVE algorithm <pave>`.
#. **VAC:** Generate patches with the :doc:`VAC algorithm <vac>`.
#. **VSA:** Generate patches with the :doc:`VSA algorithm <vsa>`.


VERBOSITY_LEVEL
+++++++++++++++
Defines the verbosity level for all console output of the patch algorithm.

   **Type:** ``INTEGER``

   **Domain:** `verbosity_level >= 0`

   **Default:** `verbosity_level = 1`

The verbosity levels are defined as follows:

.. list-table::
   :widths: 15 30
   :header-rows: 1

   * - Value
     - Description
   * - verbosity_level = 0
     - Suppress all output.
   * - verbosity_level = 1
     - Print a summary of the run when algorithm finishes.
   * - verbosity_level > 1
     - Print detailed run information, used for debugging.


MAX_ANGLE
+++++++++
Defines the maximum allowable angle (in degrees) between adjacent faces.

   **Type:** ``REAL``

   **Domain:** `0.0 <= max_angle <= 180.0`

   **Default:** `max_angle = 20.0`

All the patch algorithms construct the `adjacency matrix
<http://mathworld.wolfram.com/AdjacencyMatrix.html>`_ of the enclosure faces to efficiently
determine which faces are adjacent to others. If the normals of two 'topologically adjacent' faces
exceed *max_angle*, then the faces will not be neighbors in the internal adjacency matrix.

.. figure:: images/connected_components.png
   :figwidth: 45%
   :align: center

   The connected components of the outer surface of a furnace funnel. MAX_ANGLE is set to 20
   degrees. Each component is a different color. The face edges are omitted for clarity.

The patch algorithms guarantee that patches will be *connected sets* of faces. Therefore,
``max_angle`` divides the enclosure into connected components of faces wherever there are 'sharp'
edges whose angle exceeds the parameter. Patches will never span more than one component.

.. note::
  ``max_angle`` only applies to *pairs of adjacent faces*, so two faces within a patch may be at an
  angle greater than ``max_angle`` if the faces between them are at sufficiently large angles. This
  is unlikely in practice, given a reasonably smooth enclosure and small ``max_angle``.

.. seealso::
   The effects of ``max_angle`` vary by algorithm. Refer to the documentation of the :doc:`PAVE
   <pave>`, :doc:`VAC <vac>`, and :doc:`VSA <vsa>` algorithms for more details.



PAVE Parameters
---------------
The following namelist parameters apply only to the PAVE algorithm. For more
information, refer to the :doc:`PAVE algorithm documentation <pave>`.


PAVE_MERGE_LEVEL
++++++++++++++++
Controls the aggressiveness of patch merging for the :doc:`PAVE algorithm <pave>`.

   **Type:** ``INTEGER``

   **Domain:** `pave_merge_level >= 0`

   **Default:** `pave_merge_level = 3`

After paving is complete, there will be a valid patching of the enclosure. The algorithm then
attempts to merge patches in order to reduce the patch count.

The merge levels are defined as follows:

.. list-table::
   :widths: 15 30
   :header-rows: 1

   * - Value
     - Description
   * - pave_merge_level = 0
     - No merging.
   * - pave_merge_level = 1
     - Merge patches that are within the faces of a vertex.
   * - pave_merge_level = 2
     - Same as 1. Additionally, merge patches that are within the faces of pairs
       of adjacent vertices. The old patches are requeued with their original
       weight so that a merge is only performed if the merge candidate has a
       lower weight than any of its consituent patches.
   * - pave_merge_level >= 3
     - Same as 2. Additionally, merge patches within the faces of pairs of
       adjacent vertices, but add a large weight to the requeued old patches.
       This ensures that the merge is always performed.


PAVE_SPLIT_PATCH_SIZE
+++++++++++++++++++++
Defines the maximum size of patches to be split during patch merging for the :doc:`PAVE algorithm <pave>`.

   **Type:** ``INTEGER``

   **Domain:** `pave_split_patch_size > 1`

   **Default:** `pave_split_patch_size = 3`

Before merging patches, all :ref:`merge methods
<tools/RadE/patches/patches_namelist:PAVE_MERGE_LEVEL>` find patches with less than
``pave_split_patch_size`` faces and 'split' them into 1-face patches. The original patches aren't
actually modified, rather they are re-queued along with their constituent faces. This allows the
algorithm to find more merge candidates and then 'fill in the gaps' with the 1-face patches.

The 1-face patches have a large weight, so they will only be used after all other patches are set.
Therefore, the enclosure will tend retain the same patches as before the split, unless this is not
possible due to a merge.

.. note::
   For best results, set ``pave_split_patch_size`` to 3 for quadrilateral meshes
   and to 5 for triangular meshes. This avoids splitting too many patches.



VAC Parameters
--------------
The following namelist parameters apply only to the VAC algorithm. For more
information, refer to the :doc:`VAC algorithm documentation <vac>`.


VAC_MERGE_LEVEL
+++++++++++++++
Controls the aggressiveness of patch merging for the :doc:`VAC algorithm <vac>`.

   **Type:** ``INTEGER``

   **Domain:** `vac_merge_level >= 0`

   **Default:** `vac_merge_level = 3`

After the main stage of the VAC algorithm, there will be a valid patching of the enclosure. The
algorithm then attempts to merge patches in order to reduce the patch count.

The merge levels are defined as follows:

.. list-table::
   :widths: 15 30
   :header-rows: 1

   * - Value
     - Description
   * - vac_merge_level = 0
     - No merging.
   * - vac_merge_level = 1
     - Merge patches that are within the faces of a vertex.
   * - vac_merge_level = 2
     - Same as 1. Additionally, merge patches that are within the faces of pairs
       of adjacent vertices. The old patches are requeued with their original
       weight so that a merge is only performed if the merge candidate has a
       lower weight than any of its consituent patches.
   * - vac_merge_level >= 3
     - Same as 2. Additionally, merge patches within the faces of pairs of
       adjacent vertices, but add a large weight to the requeued old patches.
       This ensures that the merge is always performed.


VAC_SPLIT_PATCH_SIZE
++++++++++++++++++++
Defines the maximum size of patches to be split during patch merging for the :doc:`VAC algorithm <vac>`.

   **Type:** ``INTEGER``

   **Domain:** `vac_split_patch_size > 1`

   **Default:** `vac_split_patch_size = 3`

Before merging patches, all :ref:`merge methods
<tools/RadE/patches/patches_namelist:VAC_MERGE_LEVEL>` find patches with less than
``vac_split_patch_size`` faces and 'split' them into 1-face patches. The original patches aren't
actually modified, rather they are re-queued along with their constituent faces. This allows the
algorithm to find more merge candidates and then 'fill in the gaps' with the 1-face patches.

The 1-face patches have a large weight, so they will only be used after all other patches are set.
Therefore, the enclosure will tend retain the same patches as before the split, unless this is not
possible due to a merge.

.. note::
   For best results, set ``vac_split_patch_size`` to 3 for quadrilateral meshes
   and to 5 for triangular meshes. This avoids splitting too many patches.



VSA Parameters
--------------
The following namelist parameters apply only to the VSA algorithm. For more
information, refer to the :doc:`VSA algorithm documentation <vsa>`.


VSA_MAX_ITER
++++++++++++
Defines the maximum number of iterations for the :doc:`VSA algorithm <vsa>`.

   **Type:** ``INTEGER``

   **Domain:** `vsa_max_iter >= 1`

   **Default:** `vsa_max_iter = 1000`

The algorithm stops when ``vsa_max_iter`` is reached, regardless of other
terminating conditions.


VSA_MIN_DELTA
+++++++++++++
Defines the minimum allowable change in patch proxies between successive
iterations of the :doc:`VSA algorithm <vsa>`.

   **Type:** ``REAL``

   **Domain:** `vsa_min_delta >= 0.0`

   **Default:** `vsa_min_delta = 1.0E-6`

At the end of each iteration, the new patch proxies for the next iteration are computed and compared
against the old proxies. The algorithm keeps track of the *minimum* change between the old and new
proxies. This change is computed as the sum of the squares of the difference between the old and new
proxy vectors. If the minimum change in patch proxies is less than ``vsa_min_delta``, the algorithm
stops at that iteration.


VSA_AVG_FACES_PER_PATCH
+++++++++++++++++++++++
Defines the average faces per patch, and by extension the total number of patches.

   **Type:** ``REAL``

   **Domain:** vsa_avg_faces_per_patch >= 1.0

   **Default:** vsa_avg_faces_per_patch = 4.0

The average faces per patch is given by

.. math::
   \text{(Total Faces)}/\text{(Total Patches)}

Since the number of faces is fixed, this parameter determines the total number of patches in the
final configuration:

.. math::
   \text{(Total Patches)} = \text{(Total Faces)} *
   \text{vsa_avg_faces_per_patch}

Rather than set the number of patches explicitly, which is mesh dependent, expressing this parameter
as an average allows the same value to apply to a variety of meshes.
