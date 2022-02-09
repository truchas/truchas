.. sectionauthor:: David Neill-Asanza <dhna@lanl.gov>

PATCHES Namelist
==================

The `PATCHES` namelist defines the parameters used by the patching algorithms.
The namelist supports many parameters, but not all parameters are used by all
algorithms. Parameters only used by a particular algorithm are prefixed with the
algorithm's name.

These parameters are described in detail below.

:Required/Optional: Optional
:Single/Multiple Instances: Single

.. code-block::

  &PATCHES
    patch_algorithm = 'METIS'
    metis_face_patch_ratio = 2
    metis_face_weight = F
    metis_edge_weight = T
  /

:superscript:`Example PATCHES namelist`

.. contents:: Components
   :local:

General Parameters
------------------

patch_algorithm
^^^^^^^^^^^^^^^
Selects one of the available algorithms, or disables patching.

:Type: string
:Default: ``'PAVE'``
:Valid Values: - ``'NONE'``: No patches will be generated. All other parameters are ignored. This is equivalent to an absent `PATCHES` namelist.
               - ``'PAVE'``: Generate patches with the `PAVE algorithm <https://www.truchas.org/docs/sphinx/tools/RadE/patches/pave.html>`_.
               - ``'VAC'``: Generate patches with the `VAC algorithm <https://www.truchas.org/docs/sphinx/tools/RadE/patches/vac.html>`_.
               - ``'VSA'``: Generate patches with the `VSA algorithm <https://www.truchas.org/docs/sphinx/tools/RadE/patches/vsa.html>`_.
               - ``'METIS'``: Generate patches with the `METIS algorithm <https://www.truchas.org/docs/sphinx/tools/RadE/patches/metis.html>`_.
               - ``'FILE'``: Patches will be read from a file. Because the cost of computing patches can be quite substantial for very large enclosure meshes, this pseudo-algorithm is provided to enable the use of previously computed patches.


verbosity_level
^^^^^^^^^^^^^^^
Defines the verbosity level for all console output of the patch algorithm.

:Type: integer
:Default: 1
:Valid Values: -  0: Suppress all output.
               -  1: Print a summary of the run when algorithm finishes.
               - >1: Print detailed run information, used for debugging.


max_angle
^^^^^^^^^
Defines the maximum allowable angle (in degrees) between adjacent faces.

:Type: real
:Default: 20.0
:Valid Values: [0, 180]

All the patch algorithms construct the `adjacency matrix
<http://mathworld.wolfram.com/AdjacencyMatrix.html>`_ of the enclosure faces to
efficiently determine which faces are adjacent to others. If the normals of two
'topologically adjacent' faces exceed *max_angle*, then the faces will not be
neighbors in the internal adjacency matrix.

.. figure:: images/connected_components.png
   :figwidth: 45%
   :align: center

   The connected components of the outer surface of a furnace funnel. MAX_ANGLE
   is set to 20 degrees. Each component is a different color. The face edges are
   omitted for clarity.

The patch algorithms guarantee that patches will be *connected sets* of faces.
Therefore, ``max_angle`` divides the enclosure into connected components of
faces wherever there are 'sharp' edges whose angle exceeds the parameter.
Patches will never span more than one component.

.. note::
  ``max_angle`` only applies to *pairs of adjacent faces*, so two faces within a
  patch may be at an angle greater than ``max_angle`` if the faces between them
  are at sufficiently large angles. This is unlikely in practice, given a
  reasonably smooth enclosure and small ``max_angle``.

.. seealso::
   The effects of ``max_angle`` vary by algorithm. Refer to the documentation of
   the `PAVE <https://www.truchas.org/docs/sphinx/tools/RadE/patches/pave.html>`_, `VAC <https://www.truchas.org/docs/sphinx/tools/RadE/patches/vac.html>`_, and `VSA <https://www.truchas.org/docs/sphinx/tools/RadE/patches/vsa.html>`_
   algorithms for more details.


FILE Parameters
---------------
The following namelist parameter applies only to the FILE algorithm.

patch_file
^^^^^^^^^^
The path to an existing radiation enclosure file containing patch information.
The enclosure defined by the file must be identical to the current enclosure.
This may be an absolute path or a relative path.

:Type: case-sensitive string
:Default: ``''``
:Valid Values: must be a valid path


PAVE Parameters
---------------
The following namelist parameters apply only to the PAVE algorithm. For more
information, refer to the `PAVE algorithm documentation <https://www.truchas.org/docs/sphinx/tools/RadE/patches/pave.html>`_.


pave_merge_level
^^^^^^^^^^^^^^^^
Controls the aggressiveness of patch merging. After paving is complete, there
will be a valid patching of the enclosure. The algorithm then attempts to merge
patches in order to reduce the patch count.

:Type: integer
:Default: 3
:Valid Values: - 0: No merging.
               - 1: Merge patches that are within the faces of a vertex.
               - 2: Same as 1. Additionally, merge patches that are within the faces of pairs of adjacent vertices. The old patches are requeued with their original weight so that a merge is only performed if the merge candidate has a lower weight than any of its consituent patches.
               - :math:`\geq 3`: Same as 2. Additionally, merge patches within the faces of pairs of adjacent vertices, but add a large weight to the requeued old patches. This ensures that the merge is always performed.


pave_split_patch_size
^^^^^^^^^^^^^^^^^^^^^
Defines the maximum size of patches to be split during patch merging.

:Type: integer
:Default: 3
:Valid Values: :math:`\gt 1`

Before merging patches, all :ref:`merge methods
<rade/PATCHES_Namelist:pave_merge_level>` find patches with less than
``pave_split_patch_size`` faces and 'split' them into 1-face patches. The
original patches aren't actually modified, rather they are re-queued along with
their constituent faces. This allows the algorithm to find more merge candidates
and then 'fill in the gaps' with the 1-face patches.

The 1-face patches have a large weight, so they will only be used after all
other patches are set. Therefore, the enclosure will tend retain the same
patches as before the split, unless this is not possible due to a merge.

.. note::
   For best results, set ``pave_split_patch_size`` to 3 for quadrilateral meshes
   and to 5 for triangular meshes. This avoids splitting too many patches.


pave_random_seed
^^^^^^^^^^^^^^^^
Defines the seed for the random number generator used to pick the initial seed
patches.

:Type: integer
:Default: ``NONE``, the seed is taken from the system clock.
:Valid Values: :math:`\gt 0`

The PAVE algorithm begins by creating a 'seed patch' in each connected component
of the enclosure. Each component is then 'paved' or 'tiled' with patches,
starting from the seed patch. The seed patches are chosen randomly from a set of
patches determined to produce optimal results. Refer to the `seed patches
section
<https://www.truchas.org/docs/sphinx/tools/RadE/patches/pave.html#choosing-seed-patches>`_
of the PAVE documentation for more information on how the seed patches are
selected.

This parameter sets the seed for the random number generator used to pick the
seed patches. Therefore, runs with the same value for this parameter will
produce identical results. If this parameter is not specified, then the seed is
taken from the system clock and results will likely vary from run to run.


VAC Parameters
--------------
The following namelist parameters apply only to the VAC algorithm. For more
information, refer to the `VAC algorithm documentation <https://www.truchas.org/docs/sphinx/tools/RadE/patches/vac.html>`_.


vac_merge_level
^^^^^^^^^^^^^^^
Controls the aggressiveness of patch merging. After the main stage of the VAC
algorithm, there will be a valid patching of the enclosure. The algorithm then
attempts to merge patches in order to reduce the patch count.

:Type: integer
:Default: 3
:Valid Values: - 0: No merging.
               - 1: Merge patches that are within the faces of a vertex.
               - 2: Same as 1. Additionally, merge patches that are within the faces of pairs of adjacent vertices. The old patches are requeued with their original weight so that a merge is only performed if the merge candidate has a lower weight than any of its consituent patches.
               - :math:`\geq 3`: Same as 2. Additionally, merge patches within
                 the faces of pairs of adjacent vertices, but add a large weight
                 to the requeued old patches. This ensures that the merge is
                 always performed.


vac_split_patch_size
^^^^^^^^^^^^^^^^^^^^
Defines the maximum size of patches to be split during patch merging.

:Type: integer
:Default: 3
:Valid Values: :math:`\gt 1`

Before merging patches, all :ref:`merge methods
<rade/PATCHES_Namelist:vac_merge_level>` find patches with less than
``vac_split_patch_size`` faces and 'split' them into 1-face patches. The
original patches aren't actually modified, rather they are re-queued along with
their constituent faces. This allows the algorithm to find more merge candidates
and then 'fill in the gaps' with the 1-face patches.

The 1-face patches have a large weight, so they will only be used after all
other patches are set. Therefore, the enclosure will tend retain the same
patches as before the split, unless this is not possible due to a merge.

.. note::
   For best results, set ``vac_split_patch_size`` to 3 for quadrilateral meshes
   and to 5 for triangular meshes. This avoids splitting too many patches.



VSA Parameters
--------------
The following namelist parameters apply only to the VSA algorithm. For more
information, refer to the `VSA algorithm documentation <https://www.truchas.org/docs/sphinx/tools/RadE/patches/vsa.html>`_.


vsa_max_iter
^^^^^^^^^^^^
Defines the maximum number of iterations.

:Type: integer
:Default: 1000
:Valid Values: :math:`\geq 1`

The algorithm stops when ``vsa_max_iter`` is reached, regardless of other
terminating conditions.


vsa_min_delta
^^^^^^^^^^^^^
Defines the minimum allowable change in patch proxies between successive
iterations.

:Type: real
:Default: 1.0E-6
:Valid Values: :math:`\geq 0.0`

At the end of each iteration, the new patch proxies for the next iteration are
computed and compared against the old proxies. The algorithm keeps track of the
*minimum* change between the old and new proxies. This change is computed as the
sum of the squares of the difference between the old and new proxy vectors. If
the minimum change in patch proxies is less than ``vsa_min_delta``, the
algorithm stops at that iteration.


vsa_face_patch_ratio
^^^^^^^^^^^^^^^^^^^^
Defines the ratio of total faces to total patches, and by extension the total
number of patches.

:Type: real
:Default: 4.0
:Valid Values: :math:`\geq 1.0`

Since the number of faces is fixed, this parameter determines the total number
of patches in the final configuration:

.. math::
   \text{(Total Patches)} = \text{(Total Faces)}\ /\ \text{vsa_face_patch_ratio}

Rather than set the number of patches explicitly, which is mesh dependent,
expressing this parameter as a ratio allows the same value to apply to a variety
of meshes.


vsa_max_patch_radius
^^^^^^^^^^^^^^^^^^^^
Defines the desired maximum radius for a patch.

:Type: real
:Default: ``sqrt(huge(0.0_r8))``
:Valid Values: :math:`\gt 0.0`

This parameter is used to compute the *size bias* term of the weight of a face
relative to a patch proxy. Refer to the `size bias section
<https://www.truchas.org/docs/sphinx/tools/RadE/patches/vsa.html#size-bias>`_ of
the VSA documentation for more information on how the parameter affects the face
weight computation.

Note that the default value of this parameter is ``sqrt(huge(0.0_r8))``
because it is squared in the face weight computation. By taking the root of
``huge(0.0_r8)`` we prevent floating point overflow errors. Numerically,
the default value on the order of :math:`1.34\times 10^{154}`.


vsa_normalize_dist
^^^^^^^^^^^^^^^^^^
Determines whether to normalize the distance bias.

:Type: logical
:Default: True

This parameter affects the computation of the *distance bias* term of the weight
of a face relative to a patch proxy. Broadly speaking, enabling normalization
tends to produce patches with a similar number of faces, regardless of the
physical size of each patch. Conversely, disabling normalization tends to make
all patches about the same physical size, regardless of the number of faces in
each patch.

Refer to the `distance bias section
<https://www.truchas.org/docs/sphinx/tools/RadE/patches/vsa.html#distance-bias>`_
of the VSA documentation for more information on how the parameter affects the
face weight computation.


vsa_random_seed
^^^^^^^^^^^^^^^
Defines the seed for the random number generator used to pick the initial seed
patches.

:Type: integer
:Default: ``NONE``, the seed is taken from the system clock.
:Valid Values: :math:`\gt 0`

The VSA algorithm uses a 'farthest-point' initialization method to choose the
seed patches for the first iteration. To start, a random face in each connected
component of the enclosure is chosen as a seed patch. Then, seed patches are
added one at a time by performing a `partitioning
<https://www.truchas.org/docs/sphinx/tools/RadE/patches/vsa.html#geometry-partitioning>`_
and then choosing the face with highest total distortion as the new seed patch.

This parameter sets the seed for the random number generator used to pick the
first seed patch in each connected component. Therefore, runs with the same
value for this parameter will produce identical results. If this parameter is
not specified, then the seed is taken from the system clock and results will
likely vary from run to run.



METIS Parameters
----------------
The following namelist parameters apply only to the METIS algorithm. For more
information, refer to the `METIS algorithm documentation <https://www.truchas.org/docs/sphinx/tools/RadE/patches/metis.html>`_.

The METIS algorithm constructs the weighted dual graph of the enclosure and
passes it to the METIS library :footcite:`Karypis:1998:METIS` to
partition the dual graph. The METIS namelist parameters are thus divided into
two: those that are used to construct the dual graph, and those that are passed
directly to the METIS graph partitioner.

We first discuss the three parameters used during initialization, and then
briefly present the 12 METIS library parameters passed to the graph partitioner.


Initialization Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

metis_face_patch_ratio
''''''''''''''''''''''
Defines the ratio of total faces to total desired patches, and by extension the
final number of patches generated.

:Type: real
:Default: 4.0
:Valid Values: :math:`\geq 1.0`

This parameter determines the number of partitions :math:`N_p` passed to the METIS
graph partitioner:

.. math::
   N_p = \frac{N_f}{R}

where :math:`N_f` is the total number of faces. Since the METIS library is free to
produce less partitions than requested, :math:`N_p` is not necessarily the final
number of patches.

The METIS library must ensure that the constraints on the objective function are
satisfied (see `partitioning objective
<https://www.truchas.org/docs/sphinx/tools/RadE/patches/metis.html#partitioning-objective>`_),
and can thus produce a drastically different number of partitions than
requested. In particular, when `metis_face_weight`_ is enabled for an enclosure
with faces of vastly different sizes, the requirement to evenly divide the total
enclosure surface area among the patches might produce significantly fewer
partitions than requested.

Moreover, after the METIS library partitions the dual graph the patch splitting
step breaks up disconnected patches which may increase the final patch count. In
short, :math:`N_p` is only a suggestion for the final patch count. Consider tweaking
other parameters if an exact patch count is desired.


metis_edge_weight
'''''''''''''''''
Determines whether to weight the edges of the dual graph by the corresponding
enclosure edge lengths.

:Type: logical
:Default: True

This parameter determines whether the Euclidean length of the enclosure edges
are assigned as edge weights in the dual graph passed to the METIS library. If
the parameter is false, then the dual graph edges are assigned a weight of 1.

Refer to the `edge weight section
<https://www.truchas.org/docs/sphinx/tools/RadE/patches/metis.html#edge-weight>`_
of the METIS algorithm documentation for more information on how the parameter
affects the final patch configuration.


metis_face_weight
'''''''''''''''''
Determines whether to weight the vertices of the dual graph by the corresponding
enclosure face areas.

:Type: logical
:Default: True

This parameter determines whether the area of the enclosure faces are assigned
as vertex weights in the dual graph passed to the METIS library. If the
parameter is false, then the dual graph vertices are assigned a weight of 1.

Refer to the `face weight section
<https://www.truchas.org/docs/sphinx/tools/RadE/patches/metis.html#face-weight>`_
of the METIS algorithm documentation for more information on how the parameter
affects the final patch configuration.



METIS library parameters
^^^^^^^^^^^^^^^^^^^^^^^^^
The METIS graph partitioning routine admits the following integer-valued options
that may be specified, though all have reasonable defaults so that none must be
specified. See the METIS documentation :footcite:`Karypis:1998:METIS`
for more details on these options.

metis_ptype
'''''''''''
Specifies the partitioning method.

:Type: integer
:Default: 0
:Valid Values: - 0: Multilevel recursive bisection
               - 1: Multilevel :math:`k`-way partitioning


metis_objtype
'''''''''''''
Specifies the type of objective.

:Type: integer
:Default: metis_objtype = 0
:Valid Values: - 0: Edge-cut minimization.
               - 1: Total communication volume minimization.


metis_ctype
'''''''''''
Specifies the matching scheme to be used during coarsening.

:Type: integer
:Default: 1
:Valid Values: - 0: Random matching
               - 1: Sorted heavy-edge matching


metis_iptype
''''''''''''
Specifies the algorithm used during initial partitioning (recursive bisection
only).

:Type: integer
:Default: metis_iptype = 0
:Valid Values: - 0: Grows a bisection using a greedy strategy
               - 1: Computes a bisection at random followed by a refinement
               - 2: Derives a separator from an edge cut.
               - 3: Grow a bisection using a greedy node-based strategy
                   

metis_ncuts
'''''''''''
Specifies the number of different partitionings that will be computed. The final
partitioning will be the one that achieves the best edge-cut or communication
volume.

:Type: integer
:Default: 1
:Valid Values: :math:`\geq 1`


metis_niter
'''''''''''
Specifies the number of iterations of the refinement algorithm at each stage of
the uncoarsening process.

:Type: integer
:Default: 10
:Valid Values: :math:`\geq 1`


metis_seed
''''''''''
Specifies the seed for the random number generator.

:Type: integer
:Default: -1


metis_minconn
'''''''''''''
Specifies whether the partitioning procedure should seek to minimize the maximum
degree of the subdomain graph. The subdomain graph is the graph in which each
partition is a node, and edges connect subdomains with a shared interface.

:Type: integer
:Default: 0
:Valid Values: - 0: Does not explicitly minimize the maximum connectivity.
               - 1: Explicitly minimize the maximum connectivity.


metis_no2hop
''''''''''''
Specifies that the coarsening will not perform any 2–hop matchings when the
standard matching approach fails to sufficiently coarsen the graph.

:Type: integer
:Default: 1
:Valid Values: - 0: Performs a 2–hop matching.
               - 1: Does not perform a 2–hop matching.

.. note::
   The 2–hop matching is very effective for graphs with power-law degree distributions.               


metis_contig
''''''''''''
Specifies whether the partitioning procedure should produce partitions that are
contiguous. If the dual graph of the mesh is not connected this option is
ignored.

:Type: integer
:Default: 0
:Valid Values: - 0: Does not force contiguous partitions.
               - 1: Forces contiguous partitions.


metis_ufactor
'''''''''''''
Specifies the maximum allowed load imbalance among the partitions. A value of
:math:`n` indicates that the allowed load imbalance is :math:`(1+n)/1000`.

:Type: integer
:Default: 1 for recursive bisection (i.e., an imbalance of 1.001); 30 for
          :math:`k`-way partitioning (i.e., an imbalance of 1.03).
:Valid Values: :math:`\geq 1`


metis_dbglvl
''''''''''''
Specifies the amount and type of diagnostic information that will be written to
**stderr** by the partitioning procedure.

:Type: integer
:Default: 0
:Valid Values: :math:`\geq 1`

The default `0` means no output. Use `1` to write some basic information. Refer
to the METIS documentation :footcite:`Karypis:1998:METIS` for the many
other possible values and the output they generate.



.. footbibliography::
