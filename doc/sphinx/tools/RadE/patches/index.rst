.. sectionauthor:: David Neill-Asanza <dhna@lanl.gov>

Enclosure Patches
=================

.. TODO:: finish RadE patches page

Patch Algorithms
----------------
The :program:`genre` program supports the following patching algorithms.

.. toctree::
   :titlesonly:

   pave
   vsa
   metis
   vac

.. |pave_patches| image:: images/basic_hemi_pave_1.png
   :width: 100%
   :align: middle

.. |vsa_patches| image:: images/basic_hemi_vsa_1.png
   :width: 100%
   :align: middle

.. |metis_patches| image:: images/basic_hemi_metis_1.png
   :width: 100%
   :align: middle

.. table::
   :align: center
   :width: 100%

   +-----------------+-----------------+-----------------+
   | |pave_patches|  | |vsa_patches|   | |metis_patches| |
   +-----------------+-----------------+-----------------+
   | Result of running **PAVE** (left), **VSA**          |
   | (middle), and **METIS** (right) on                  |
   | the 'basic hemi' enclosure.                         |
   +-----------------------------------------------------+

General guidance
++++++++++++++++
We provide general guidance for choosing an algorithm given performance constraints and desired
properties for the resulting patches.

- PAVE

  - **Pros:** produces neat tilings that account for patch planarity and irregularity. The algorithm
    is fast (linear in the number of faces).

  - **Cons:** uniform coarsening of the mesh. Only supports faces-per-patch ratios of roughly 4 or
    6, depending on whether the mesh is hexahedral or tetrahedral.

- VSA

  - **Pros:** works well for arbitrary faces-per-patch ratios. Provides flexible control of the
    geometric patch diameters. Maximizes patch planarity. Provides choice between patches with
    roughly the same number of faces, or roughly the same surface area.

  - **Cons:** can be extremely slow for large meshes (scales badly), due to a quadratic dependence on
    the number of faces. However, once computed, the result can be reused with the :ref:`FILE method
    <tools/RadE/patches/patches_namelist:PATCH_ALGORITHM>` allowing experimentation with Chaparral
    parameters without rerunning the algorithm.


- METIS

  - **Pros:** works well for arbitrary faces-per-patch ratios. The algorithm is fast by leveraging
    the METIS library :cite:`patches-index-Karypis:1998:METIS`. Provides choice between patches with
    roughly the same number of faces, or roughly the same surface area. In the second case, the
    constraint is tight, so the faces-per-patch ratio is largely ignored and need only be
    approximate.

  - **Cons:** New and not extensively tested. Patch planarity is not explicitly maximized, but at
    least patches will not straddle "sharp edges" (see :ref:`here <tools/RadE/patches/metis:Dual
    Graph>` for more details). Note that it is not clear whether patch planarity, or lack thereof,
    affects the quality of the thermal radiation model.

- VAC

  - Retained for developer use and not recommended for regular users. It has similar limitations to
    PAVE, but produces lower quality patches.

Refer to each algorithm's documentation for extensive details on its properties, limitations, and
supported inputs.


PATCHES Namelist
----------------
.. toctree::
   :maxdepth: 2
   :hidden:

   patches_namelist

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

Refer to the :doc:`PATCES namelist documentation <patches_namelist>` for detailed information on
these parameters.


References
----------
.. bibliography:: references.bib
   :style: unsrt
   :keyprefix: patches-index-
