.. sectionauthor:: David Neill-Asanza <dhna@lanl.gov>

PATCHES Namelist
==================

The **PATCHES** namelist defines the parameters used by the patching algorithms.
The namelist supports many parameters, but not all parameters are used by all
algorithms. Parameters only used by a particular algorithm are prefixed with the
algorithm's name. These parameters are described in detail below.

:Required/Optional: Optional
:Single/Multiple Instances: Single


General Guidance
----------------

Genre's patching allows the Chaparral library to group faces together into
patches. This has two main purposes, both especially useful in large problems:

1. To reduce the output file size. In many scenarios, the radiation enclosure
   output produced by Genre is prohibitively large (on the order of 100s of GB).
   It can take a very long time to write these files to disk, or it may be
   impossible for Truchas to fit into memory.

2. To reduce the need for very high :ref:`hemicube_resolution
   <rade/CHAPARRAL_Namelist:hemicube_resolution>` on fine meshes. As face mesh
   resolution increases, increasingly fine ``hemicube_resolution`` is needed to
   resolve view factors. This quadratically increases the genre runtime. By
   grouping faces into patches, the need for high ``hemicube_resolution`` is
   less.

Genre offers several patching algorithms, but for most situations, the
``'METIS'`` algorithm (default) is the best choice for speed, accuracy, and
simplicity. The `metis_face_patch_ratio`_ can be tuned as needed. A sample input
is given below, and is considered a good starting point for most problems. One
should always view the patches using :ref:`Vizre <rade/index:Invoking Vizre>` to
assess their quality.

.. code-block::

  &PATCHES
    metis_face_patch_ratio = 4
  /

:superscript:`Example PATCHES namelist`

.. image:: images/basic_hemi_metis_1.png
   :width: 40%
   :align: center

:superscript:`Result of patching on the 'basic hemi' enclosure.`


Components
----------

.. contents::
   :local:

.. toctree::
   :maxdepth: 1

   Expert Options <PATCHES_Namelist_Expert_Options>


METIS Parameters
^^^^^^^^^^^^^^^^

The following namelist parameters apply only to the METIS algorithm. For more
information, refer to the `METIS algorithm documentation
<https://www.truchas.org/docs/sphinx/tools/RadE/patches/metis.html>`_.

The METIS algorithm constructs the weighted dual graph of the enclosure and
passes it to the METIS library :footcite:`Karypis:1998:METIS` to
partition the dual graph. The METIS namelist parameters are thus divided into
two: those that are used to construct the dual graph, and those that are passed
directly to the METIS graph partitioner.

We first discuss the three parameters used during initialization, and then
briefly present the 12 METIS library parameters passed to the graph partitioner.


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

where :math:`N_f` is the total number of faces and :math:`R` is the
``metis_face_patch_ratio``. Since the METIS library is free to produce less
partitions than requested, :math:`N_p` is not necessarily the final number of
patches.

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


metis_ncuts
'''''''''''
Specifies the number of different partitionings that will be computed. The final
partitioning will be the one that achieves the best edge-cut or communication
volume.

:Type: integer
:Default: 1
:Valid Values: :math:`\geq 1`


metis_seed
''''''''''
Specifies the seed for the random number generator.

:Type: integer
:Default: -1

.. note::
   With the default option, patches are not deterministic (depending on the
   platform). If it is necessary for two runs with identical input to produce
   identical output, set this parameter.


metis_ufactor
'''''''''''''
Specifies the maximum allowed load imbalance among the partitions. A value of
:math:`n` indicates that the allowed load imbalance is :math:`(1+n)/1000`.

:Type: integer
:Default: 1 for recursive bisection (i.e., an imbalance of 1.001); 30 for
          :math:`k`-way partitioning (i.e., an imbalance of 1.03).
:Valid Values: :math:`\geq 1`

.. note::
   For large enclosures, this parameter can be useful to give METIS more leeway
   in optimizing the patches.

.. footbibliography::
