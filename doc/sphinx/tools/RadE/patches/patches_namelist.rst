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
many parameters, but not all parameters are used by all algorithms. Parameters only used by a
particular algorithm are prefixed with the algorithm's name.

These parameters are described in detail below.

.. contents:: PATCHES namelist parameters
   :local:
   :backlinks: none


General Parameters
------------------

PATCH_ALGORITHM
+++++++++++++++
Selects one of the available algorithms, or disables patching.

.. namelist_parameter::
   :type: STRING
   :domain: Must be one of ``'NONE'``, ``'PAVE'``, ``'VAC'``, ``'VSA'``, ``'METIS'``, or ``'FILE'``
   :default: patch_algorithm = ``'PAVE'``

Each option selects a different patch algorithm:

#. **NONE:** No patches will be generated. All other parameters are ignored. This is equivalent to
   an absent `PATCHES` namelist.
#. **PAVE:** Generate patches with the :doc:`PAVE algorithm <pave>`.
#. **VAC:** Generate patches with the :doc:`VAC algorithm <vac>`.
#. **VSA:** Generate patches with the :doc:`VSA algorithm <vsa>`.
#. **METIS:** Generate patches with the :doc:`METIS algorithm <metis>`.
#. **FILE:** Patches will be read from a file. Because the cost of computing
   patches can be quite substantial for very large enclosure meshes, this
   pseudo-algorithm is provided to enable the use of previously computed
   patches.


VERBOSITY_LEVEL
+++++++++++++++
Defines the verbosity level for all console output of the patch algorithm.

.. namelist_parameter::
   :type: INTEGER
   :domain: verbosity_level >= 0
   :default: verbosity_level = 1

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

.. namelist_parameter::
   :type: REAL
   :domain: 0.0 <= max_angle <= 180.0
   :default: max_angle = 20.0

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


FILE Parameters
---------------
The following namelist parameter applies only to the FILE algorithm.

PATCH_FILE
++++++++++
The path to an existing radiation enclosure file containing patch information.

.. namelist_parameter::
   :type: STRING
   :domain: patch_file must be a valid path
   :default: patch_file = ``''``

The enclosure defined by the file must be identical to the current enclosure. This may be an
absolute path or a relative path.


PAVE Parameters
---------------
The following namelist parameters apply only to the PAVE algorithm. For more
information, refer to the :doc:`PAVE algorithm documentation <pave>`.


PAVE_MERGE_LEVEL
++++++++++++++++
Controls the aggressiveness of patch merging.

.. namelist_parameter::
   :type: INTEGER
   :domain: pave_merge_level >= 0
   :default: pave_merge_level = 3

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
Defines the maximum size of patches to be split during patch merging.

.. namelist_parameter::
   :type: INTEGER
   :domain: pave_split_patch_size > 1
   :default: pave_split_patch_size = 3

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


PAVE_RANDOM_SEED
++++++++++++++++
Defines the seed for the random number generator used to pick the initial seed patches.

.. namelist_parameter::
   :type: INTEGER
   :domain: pave_random_seed > 0
   :default: ``NONE``, the seed is taken from the system clock.

The PAVE algorithm begins by creating a 'seed patch' in each connected component of the enclosure.
Each component is then 'paved' or 'tiled' with patches, starting from the seed patch. The seed
patches are chosen randomly from a set of patches determined to produce optimal results. Refer to
the :ref:`seed patches section <tools/RadE/patches/pave:Choosing Seed Patches>` of the PAVE
documentation for more information on how the seed patches are selected.

This parameter sets the seed for the random number generator used to pick the seed patches.
Therefore, runs with the same value for this parameter will produce identical results. If this
parameter is not specified, then the seed is taken from the system clock and results will likely
vary from run to run.


VAC Parameters
--------------
The following namelist parameters apply only to the VAC algorithm. For more
information, refer to the :doc:`VAC algorithm documentation <vac>`.


VAC_MERGE_LEVEL
+++++++++++++++
Controls the aggressiveness of patch merging.

.. namelist_parameter::
   :type: INTEGER
   :domain: vac_merge_level >= 0
   :default: vac_merge_level = 3

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
Defines the maximum size of patches to be split during patch merging.

.. namelist_parameter::
   :type: INTEGER
   :domain: vac_split_patch_size > 1
   :default: vac_split_patch_size = 3

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
Defines the maximum number of iterations.

.. namelist_parameter::
   :type: Integer
   :domain: vsa_max_iter >= 1
   :default: vsa_max_iter = 1000

The algorithm stops when ``vsa_max_iter`` is reached, regardless of other
terminating conditions.


VSA_MIN_DELTA
+++++++++++++
Defines the minimum allowable change in patch proxies between successive iterations.

.. namelist_parameter::
   :type: REAL
   :domain: vsa_min_delta >= 0.0
   :default: vsa_min_delta = 1.0E-6

At the end of each iteration, the new patch proxies for the next iteration are computed and compared
against the old proxies. The algorithm keeps track of the *minimum* change between the old and new
proxies. This change is computed as the sum of the squares of the difference between the old and new
proxy vectors. If the minimum change in patch proxies is less than ``vsa_min_delta``, the algorithm
stops at that iteration.


VSA_FACE_PATCH_RATIO
++++++++++++++++++++
Defines the ratio of total faces to total patches, and by extension the total number of patches.

.. namelist_parameter::
   :type: REAL
   :domain: vsa_face_patch_ratio >= 1.0
   :default: vsa_face_patch_ratio = 4.0

Since the number of faces is fixed, this parameter determines the total number of patches in the
final configuration:

.. math::
   \text{(Total Patches)} = \text{(Total Faces)}\ /\ \text{vsa_face_patch_ratio}

Rather than set the number of patches explicitly, which is mesh dependent, expressing this
parameter as a ratio allows the same value to apply to a variety of meshes.


VSA_MAX_PATCH_RADIUS
++++++++++++++++++++
Defines the desired maximum radius for a patch.

.. namelist_parameter::
   :type: REAL
   :domain: vsa_max_patch_radius > 0.0
   :default: vsa_max_patch_radius = sqrt(huge(0.0_r8))

This parameter is used to compute the *size bias* term of the weight of a face relative to
a patch proxy. Refer to the :ref:`size bias section <tools/RadE/patches/vsa:Size Bias>` of the
VSA documentation for more information on how the parameter affects the face weight computation.

Note that the default value of this parameter is :fortran:`sqrt(huge(0.0_r8))` because it is squared
in the face weight computation. By taking the root of :fortran:`huge(0.0_r8)` we prevent floating
point overflow errors. Numerically, the default value on the order of `1.34*10^{154}`.


VSA_NORMALIZE_DIST
++++++++++++++++++
Determines whether to normalize the distance bias.

.. namelist_parameter::
   :type: LOGICAL
   :domain: Must be ``.true.`` or ``.false.``
   :default: vsa_normalize_dist = ``.true.``

This parameter affects the computation of the *distance bias* term of the weight of a face relative
to a patch proxy. Broadly speaking, enabling normalization tends to produce patches with a similar
number of faces, regardless of the physical size of each patch. Conversely, disabling normalization
tends to make all patches about the same physical size, regardless of the number of faces in each
patch.

Refer to the :ref:`distance bias section <tools/RadE/patches/vsa:Distance Bias>`
of the VSA documentation for more information on how the parameter affects the face weight
computation.


VSA_RANDOM_SEED
+++++++++++++++
Defines the seed for the random number generator used to pick the initial seed patches.

.. namelist_parameter::
   :type: INTEGER
   :domain: pave_random_seed > 0
   :default: ``NONE``, the seed is taken from the system clock.

The VSA algorithm uses a 'farthest-point' initialization method to choose the seed patches for the
first iteration. To start, a random face in each connected component of the enclosure is chosen as a
seed patch. Then, seed patches are added one at a time by performing a :ref:`partitioning
<tools/RadE/patches/vsa:Geometry Partitioning>` and then choosing the face with highest total
distortion as the new seed patch.

This parameter sets the seed for the random number generator used to pick the first seed patch in
each connected component. Therefore, runs with the same value for this parameter will produce
identical results. If this parameter is not specified, then the seed is taken from the system clock
and results will likely vary from run to run.



METIS Parameters
----------------
The following namelist parameters apply only to the METIS algorithm. For more
information, refer to the :doc:`METIS algorithm documentation <metis>`.

The METIS algorithm constructs the weighted dual graph of the enclosure and passes it to the METIS
library :cite:`patches-nml-Karypis:1998:METIS` to partition the dual graph. The METIS namelist
parameters are thus divided into two: those that are used to construct the dual graph, and those
that are passed directly to the METIS graph partitioner.

We first discuss the three parameters used during initialization, and then briefly present the 12
METIS library parameters passed to the graph partitioner.


Initialization Parameters
+++++++++++++++++++++++++

METIS_FACE_PATCH_RATIO
^^^^^^^^^^^^^^^^^^^^^^
Defines the ratio of total faces to total desired patches, and by extension the final number of
patches generated.

.. namelist_parameter::
   :type: REAL
   :domain: metis_face_patch_ratio >= 1.0
   :default: meti_face_patch_ratio = 4.0

This parameter determines the number of partitions NPART passed to the METIS graph partitioner:

.. math::
   \text{NPART} = \frac{\text{NFACE}}{\text{METIS_FACE_PATCH_RATIO}}

where NFACE is the total number of faces. Since the METIS library is free to produce less partitions
than requested, NPART is not necessarily the final number of patches.

The METIS library must ensure that the constraints on the objective function are satisfied (see
:ref:`partitioning objective <tools/RadE/patches/metis:Partitioning objective>`), and can thus
produce a drastically different number of partitions than requested. In particular, when
:ref:`METIS_FACE_WEIGHT <tools/RadE/patches/metis:METIS_FACE_WEIGHT>` is enabled for an enclosure
with faces of vastly different sizes, the requirement to evenly divide the total enclosure surface
area among the patches might produce significantly fewer partitions than requested.

Moreover, after the METIS library partitions the dual graph the patch splitting step breaks up
disconnected patches which may increase the final patch count. In short, NPART is only a suggestion
for the final patch count. Consider tweaking other parameters if an exact patch count is desired.


METIS_EDGE_WEIGHT
^^^^^^^^^^^^^^^^^
Determines whether to weight the edges of the dual graph by the corresponding enclosure edge lengths.

.. namelist_parameter::
   :type: LOGICAL
   :domain: Must be ``.true.`` or ``.false.``
   :default: metis_edge_weight = ``.true.``

This parameter determines whether the Euclidean length of the enclosure edges are assigned as edge
weights in the dual graph passed to the METIS library. If the parameter is false, then the dual
graph edges are assigned a weight of 1.

Refer to the :ref:`edge weight section <tools/RadE/patches/metis:Edge Weight>` of the METIS
algorithm documentation for more information on how the parameter affects the final patch
configuration.


METIS_FACE_WEIGHT
^^^^^^^^^^^^^^^^^
Determines whether to weight the vertices of the dual graph by the corresponding enclosure face
areas.

.. namelist_parameter::
   :type: LOGICAL
   :domain: Must be ``.true.`` or ``.false.``
   :default: metis_face_weight = ``.true.``

This parameter determines whether the area of the enclosure faces are assigned as vertex weights in
the dual graph passed to the METIS library. If the parameter is false, then the dual graph vertices
are assigned a weight of 1.

Refer to the :ref:`face weight section <tools/RadE/patches/metis:Face Weight>` of the METIS
algorithm documentation for more information on how the parameter affects the final patch
configuration.



METIS library parameters
++++++++++++++++++++++++
The METIS graph partitioning routine admits the following integer-valued options that may be
specified, though all have reasonable defaults so that none must be specified. See the METIS
documentation :cite:`patches-nml-Karypis:1998:METIS` for more details on these options.

METIS_PTYPE
^^^^^^^^^^^
Specifies the partitioning method.

.. namelist_parameter::
   :type: INTEGER
   :domain: metis_ptype `\in` {0,1}
   :default: metis_ptype = 0

The partitioning methods are encoded as follows:

.. list-table::
   :widths: 15 30
   :header-rows: 1

   * - Value
     - Description
   * - metis_ptype = 0
     - Multilevel recursive bisection
   * - metis_ptype = 1
     - Multilevel `k`-way partitioning


METIS_OBJTYPE
^^^^^^^^^^^^^
Specifies the type of objective.

.. namelist_parameter::
   :type: INTEGER
   :domain: metis_objtype `\in` {0,1}
   :default: metis_objtype = 0

The objective types are encoded as follows:

.. list-table::
   :widths: 15 30
   :header-rows: 1

   * - Value
     - Description
   * - metis_objtype = 0
     - Edge-cut minimization.
   * - metis_objtype = 1
     - Total communication volume minimization.


METIS_CTYPE
^^^^^^^^^^^
Specifies the matching scheme to be used during coarsening.

.. namelist_parameter::
   :type: INTEGER
   :domain: metis_ctype `\in` {0,1}
   :default: metis_ctype = 1

The matching schemes are encoded as follows:

.. list-table::
   :widths: 15 30
   :header-rows: 1

   * - Value
     - Description
   * - metis_ctype = 0
     - Random matching
   * - metis_ctype = 1
     - Sorted heavy-edge matching


METIS_IPTYPE
^^^^^^^^^^^^
Specifies the algorithm used during initial partitioning (recursive bisection only).

.. namelist_parameter::
   :type: INTEGER
   :domain: metis_iptype `\in` {0,1,2,3}
   :default: metis_iptype = 0

The partitioning algorithms are encoded as follows:

.. list-table::
   :widths: 15 30
   :header-rows: 1

   * - Value
     - Description
   * - metis_iptype = 0
     - Grows a bisection using a greedy strategy
   * - metis_iptype = 1
     - Computes a bisection at random followed by a refinement
   * - metis_iptype = 2
     - Derives a separator from an edge cut.
   * - metis_iptype = 3
     - Grow a bisection using a greedy node-based strategy


METIS_NCUTS
^^^^^^^^^^^
Specifies the number of different partitionings that will be computed. The final partitioning will
be the one that achieves the best edge-cut or communication volume.

.. namelist_parameter::
   :type: INTEGER
   :domain: metis_ncuts >= 1
   :default: metis_ncuts = 1


METIS_NITER
^^^^^^^^^^^
Specifies the number of iterations of the refinement algorithm at each stage of the uncoarsening
process.

.. namelist_parameter::
   :type: INTEGER
   :domain: metis_niter >= 1
   :default: metis_niter = 10


METIS_SEED
^^^^^^^^^^
Specifies the seed for the random number generator.

.. namelist_parameter::
   :type: INTEGER
   :domain: metis_seed `\in \mathbb{Z}`
   :default: metis_seed = -1


METIS_MINCONN
^^^^^^^^^^^^^
Specifies whether the partitioning procedure should seek to minimize the maximum degree of the
subdomain graph.

.. namelist_parameter::
   :type: INTEGER
   :domain: metis_minconn `\in` {0,1}
   :default: metis_minconn = 0

The subdomain graph is the graph in which each partition is a node, and edges connect subdomains
with a shared interface. This parameter is encoded as follows:

.. list-table::
   :widths: 15 30
   :header-rows: 1

   * - Value
     - Description
   * - metis_minconn = 0
     - Does not explicitly minimize the maximum connectivity.
   * - metis_minconn = 1
     - Explicitly minimize the maximum connectivity.


METIS_NO2HOP
^^^^^^^^^^^^
Specifies that the coarsening will not perform any 2–hop matchings when the standard matching
approach fails to sufficiently coarsen the graph.

.. namelist_parameter::
   :type: INTEGER
   :domain: metis_no2hop `\in` {0,1}
   :default: metis_no2hop = 1

The 2–hop matching is very effective for graphs with power-law degree distributions. This parameter
is encoded as follows:

.. list-table::
   :widths: 15 30
   :header-rows: 1

   * - Value
     - Description
   * - metis_no2hop = 0
     - Performs a 2–hop matching.
   * - metis_no2hop = 1
     - Does not perform a 2–hop matching.


METIS_CONTIG
^^^^^^^^^^^^
Specifies whether the partitioning procedure should produce partitions that are contiguous.

.. namelist_parameter::
   :type: INTEGER
   :domain: metis_contig `\in` {0,1}
   :default: metis_contig = 0

If the dual graph of the mesh is not connected this option is ignored. This parameter is encoded as
follows:

.. list-table::
   :widths: 15 30
   :header-rows: 1

   * - Value
     - Description
   * - metis_contig = 0
     - Does not force contiguous partitions.
   * - metis_contig = 1
     - Forces contiguous partitions.


METIS_UFACTOR
^^^^^^^^^^^^^
Specifies the maximum allowed load imbalance among the partitions.

.. namelist_parameter::
   :type: INTEGER
   :domain: metis_ufactor >= 1
   :default: metis_ufactor = 1

A value of `n` indicates that the allowed load imbalance is `(1+n)/1000`. The default is `1` for
recursive bisection (i.e., an imbalance of `1.001`) and the default value is `30` for `k`-way
partitioning (i.e., an imbalance of `1.03`).


METIS_DBGLVL
^^^^^^^^^^^^
Specifies the amount and type of diagnostic information that will be written to **stderr** by the
partitioning procedure.

.. namelist_parameter::
   :type: INTEGER
   :domain: metis_dbglvl >= 1
   :default: metis_dbglvl = 0

The default `0` means no output. Use `1` to write some basic information. Refer to the METIS
documentation :cite:`patches-nml-Karypis:1998:METIS` for the many other possible values and the
output they generate.



References
----------
.. bibliography:: references.bib
   :style: unsrt
   :keyprefix: patches-nml-
