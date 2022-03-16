CHAPARRAL Namelist
==================

The **CHAPARRAL** namelist specifies solver parameters for the Chaparral view
factor library :footcite:`glass_chaparral_1995`.

:Required/Optional: Optional
:Single/Multiple Instances: Single

.. note::
   The CHAPARRAL namelist is optional, but Genre will perform different tasks
   depending on whether one is supplied. If there is a CHAPARRAL namelist, the
   specified enclosure surface is generated and written to the enclosure file
   along with the calculated view factors. If there is *no* CHAPARRAL namelist,
   just the enclosure surface is written; this is useful for examining the
   surface for correctness prior to performing the expensive view factor
   calculation.

General Guidance
----------------

- Use `blocking_enclosure`_, `partial_enclosure`_, and `partial_area`_ to define
  the problem.
- Use `hemicube_resolution`_ and `min_separation`_ to set appropriate solver
  parameters :footcite:`neill-asanza_genre_2019`. The defaults are a good
  starting point, and may be sufficient.
- The defaults for the other parameters are sufficient for most problems.
- It is a good idea to consider :ref:`patching <rade/PATCHES_Namelist:PATCHES
  Namelist>`, especially for high-resolution problems.
- Read the genre output, noting:

  - The "Maximum desired surface subdivision" is less than or equal to the
    "Actual maximum surface subdivision". If not, you may need to change either
    `min_separation`_ or :ref:`max_subdivisions
    <rade/CHAPARRAL_Namelist:max_subdivisions (expert)>`, or address a problem
    in the mesh.
  - The number of "Nonzero lower/upper triangular entries" changes by a few % at
    most during the smoothing step. If the change is excessive, e.g. 30%,
    consider the tip in the `hemicube_resolution`_ section.
  - The message "WARNING!! Solution is not converged!" should not appear in the
    smoothing stage. If it does, follow the suggestions in the
    `hemicube_resolution`_ tip. Increasing the maximum smoothing iterations is
    not advised, this failure generally indicates an underlying issue.
       
The best way to tell if your view factor calculation is trustworthy is by
looking at the output for the smoothing stage. A well-resolved simulation does
not require significantly changing the number of non-zero entries, and should
require no more than tens of iterations. For instance, consider this smoothing
output:

.. code-block::

  ******************************************************************
   V I E W F A C T O R    M A T R I X    S M O O T H I N G
  ******************************************************************

   Smoothing of viewfactor matrix for enclosure <OUTER>
     PCG Solver, wt = 2.0
     Max iterations = 200
     Tolerance      = 1e-08

     Enforcing reciprocity by averaging...
       Nonzero lower triangular entries = 14783488 changed to 14889500
       Nonzero upper triangular entries = 14715955 changed to 14889500
       Elapsed time                     = 1.67
     Number of passes     = 1
     Number of iterations = 9
     Elapsed time         = 4.51

The number of nonzero entries changed from ~14.7m to ~14.9m, roughly a 1%
increase. Furthermore, only 9 iterations were performed. If a smoothing step
increases the number of nonzero entries significantly (e.g. by 30%), or takes
closer to the 50-iteration limit, then something is likely wrong and the results
should not be trusted. Steps to improve the simulation should be:

1. Ensure the mesh does not have fundamental issues, such as poorly shaped
   elements or vanishingly tiny elements.

2. Ensure the problem is set up correctly, with appropriate values for
   `blocking_enclosure`_, `partial_enclosure`_, and `partial_area`_.

3. Improve solver parameters. Either the `hemicube_resolution`_ needs to be
   increased, `min_separation`_ needs to be increased, :ref:`patching
   <rade/PATCHES_Namelist:PATCHES Namelist>` needs to be introduced, or
   :ref:`metis_face_patch_ratio <rade/PATCHES_Namelist:metis_face_patch_ratio>`
   needs to be increased (assuming METIS patching).

.. contents:: Components
   :local:


blocking_enclosure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Defines whether this is a blocking problem; i.e., whether the path between two
surfaces might be blocked by a third surface. If the shape produced by the
enabled radiating surfaces is *entirely convex*, this parameter may be set to
false.

:Type: logical
:Default: true

.. warning::
   If ``blocking_enclosure = F`` is provided to a blocking problem, the view
   factors will be incorrect. If ``blocking_enclosure = T`` is provided to a
   non-blocking problem, the view factor calculation will run slower than
   necessary.


partial_enclosure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Defines whether this is a partial enclosure problem (i.e., one with gaps to an
ambient temperature).

:Type: logical
:Default: false

.. note::
   When true, ``partial_area`` must be given as well.

.. note::
   This parameter informs Genre whether a given geometry contains gaps *after*
   all symmetries are applied. For example, if computing view factors on the
   interior of a complete sphere, with quarter-symmetry applied, one ought to
   supply ``partial_enclosure = F`` (the default).

.. tip::
   If an enclosure is water-tight, it is not a partial enclosure. If an
   enclosure might leak when filled with water, it is a partial enclosure.


partial_area
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Area of the gaps in a partial enclosure problem, *before* symmetries are
applied. Required when ``partial_enclosure`` is true.

:Type: real
:Default: none
:Valid Values: :math:`\gt 0`

.. note::
   There are many possible surfaces which will fill gaps in a surface. For best
   results, this area ought to be near the minimum partial enclosure area. We
   suggest the provided ``partial_area`` be greater than or equal to the minimum
   partial enclosure area, and within a factor of two of the minimum partial
   enclosure area. Genre computes this and reports the value after computing
   view factors (for example, see below). It may be necessary to run genre
   twice; once to compute the minimum partial enclosure area, and a second time
   with that value provided as the ``partial_area`` value.
   

   .. code-block::

        ******************************************************************
         V I E W F A C T O R    C A L C U L A T I O N
        ******************************************************************

         Calculating viewfactors for enclosure <OUTER>
           enclosure geometry:    3D
           enclosure type:        partial (area=0.001257), blocking

         <snip>

         Minimum Partial Enclosure Area = 0.00122526


hemicube_resolution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The number of 1D subdivisions for the hemicube over each face. Given a
``hemicube_resolution`` of :math:`n`, there will be a total of :math:`n^2` total
subdivisions per hemicube.

:Type: integer
:Default: 500
:Valid Values: :math:`\geq 4`

.. note::
   This is one of the most significant solver factors. Increasing
   ``hemicube_resolution`` will both increase runtime and improve result
   quality. The given default is a good starting point. This tends to be more
   important than ``min_separation`` on *fine meshes*. This is because on fine
   meshes, faces are less likely to need subdivision, while higher resolution
   hemicubes are needed to accurately hit faces.

.. tip::
   Read the `General Guidance`_ section for tips on how to decide whether
   ``hemicube_resolution`` ought to be increased.


min_separation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The minimum ratio of distance to diameter between any two faces.

A face :math:`j` is subdivided until the following condition is satisfied for
all subfaces :math:`k`:

.. math::
   \delta \le \frac{\Delta x^\mathrm{min}_j}{d_{kj}}

where :math:`\delta` is the ``min_separation``, :math:`d_{kj}` is the diameter
of subface :math:`k` of face :math:`j`, and :math:`\Delta x^\mathrm{min}_j` is
the minimum distance between face :math:`j` and all other faces.

:Type: real
:Default: 20
:Valid Values: :math:`\geq 0`

.. note::
   This is one of the most significant solver factors. Increasing
   ``min_separation`` will both increase runtime and improve result quality. The
   given default is a good starting point. This tends to be more important than
   ``hemicube_resolution`` on *coarse meshes*. This is because on coarse meshes,
   low-resolution hemicubes tend to hit most faces, while coarse faces are more
   likely to need subdivision.

.. tip::
   A range of 10 - 40 is most useful for most problems.

.. warning:: 
   The number of subdivisions is limited by :ref:`max_subdivisions
   <rade/CHAPARRAL_Namelist:max_subdivisions (expert)>`.


verbosity_level
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Determines the detail and frequency of terminal output.

:Type: integer
:Default: 2
:Valid Values: :math:`\geq 0`


max_subdivisions (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The maximum face subdivisions in 1D. Given a ``max_subdivisions`` of :math:`n`,
a total of :math:`n^2` subdivisions per face will be allowed.

:Type: integer
:Default: 100
:Valid Values: :math:`\geq 0`

.. note::
   The default is set such that, in most scenarios, the minimum separation will
   be reached.


.. warning::
   This limit will not be exceeded, regardless of ``min_separation``. Genre
   prints the maximum number of subdivisions needed to satisfy the given
   ``min_separation``, shown below. If the maximum desired surface subdivision
   exceeds the actual maximum allowed by ``max_subdivisions``, the computed view
   factors may be of poor quality. Below is an example of a good output, where
   the maximum desired surface subdivision did not reach the maximum allowed.

   .. code-block::

      ******************************************************************
       V I E W F A C T O R    C A L C U L A T I O N
      ******************************************************************

         <snip>
         Maximum desired surface subdivision = 38, 37
         Actual maximum surface subdivision  = 100


BSP_max_tree_depth (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a Chaparral-internal parameter. :footcite:`glass_chaparral_1995`

:Type: integer
:Default: 50
:Valid Values: :math:`\geq 1`


BSP_min_leaf_length (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a Chaparral-internal parameter. :footcite:`glass_chaparral_1995`

:Type: integer
:Default: 25
:Valid Values: :math:`\geq 1`


spatial_tolerance (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a Chaparral-internal parameter. :footcite:`glass_chaparral_1995`

:Type: real
:Default: :math:`10^{-8}`
:Valid Values: :math:`\gt 0`


smoothing_max_iter (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is an upper limit to the number of smoothing iterations permitted

:Type: integer
:Default: 50
:Valid Values: :math:`\geq 0`

.. warning::
   This upper limit should not be increased. If a simulation is failing to
   converge at the smoothing step, the results ought not be trusted. Some other
   issue is likely present, for instance insufficient ``hemicube_resolution`` or
   insufficient :ref:`patching <rade/PATCHES_Namelist:PATCHES Namelist>`. The
   number of smoothing iterations used is printed by genre, see below:

   .. code-block::

     ******************************************************************
      V I E W F A C T O R    M A T R I X    S M O O T H I N G
     ******************************************************************
        <snip>
        Number of passes     = 1
        Number of iterations = 17

   If smoothing failed to converge in the maximum permitted steps, the following
   message will appear:

   ``WARNING!!  Solution is not converged!``


smoothing_tolerance (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a Chaparral-internal parameter. :footcite:`glass_chaparral_1995`

:Type: real
:Default: :math:`10^{-8}`
:Valid Values: :math:`\gt 0`


smoothing_weight (expert)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a Chaparral-internal parameter. :footcite:`glass_chaparral_1995`

:Type: real
:Default: 2.0
:Valid Values: :math:`\gt 0`


.. footbibliography::
