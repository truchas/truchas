CHAPARRAL Namelist
==================

The **CHAPARRAL** namelist specifies solver parameters for the Chaparral view factor library.

:Required/Optional: Optional
:Single/Multiple Instances: Single

.. note::

   The CHAPARRAL namelist is optional, but Genre will perform different tasks depending on
   whether one is supplied. If there is a CHAPARRAL namelist, the specified enclosure surface is
   generated and written to the enclosure file along with the calculated view factors. If there is
   *no* CHAPARRAL namelist, just the enclosure surface is written; this is useful for examining the
   surface for correctness prior to performing the expensive view factor calculation.

.. contents:: Components
   :local:


blocking_enclosure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

:Type: logical
:Default: true


partial_enclosure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

:Type: logical
:Default: false


partial_area
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

:Type: real
:Default: none
:Valid Values: :math:`\gt 0`


BSP_max_tree_depth
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

:Type: integer
:Default: 15
:Valid Values: :math:`\geq 1`


BSP_min_leaf_length
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

:Type: integer
:Default: 25
:Valid Values: :math:`\geq 1`


spatial_tolerance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

:Type: real
:Default: none
:Valid Values: :math:`\gt 0`


hemicube_resolution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

:Type: integer
:Default: none
:Valid Values: :math:`\geq 4`


min_separation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

:Type: real
:Default: none
:Valid Values: :math:`\geq 0`


max_subdivisions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

:Type: integer
:Default: none
:Valid Values: :math:`\geq 0`


smoothing_tolerance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

:Type: real
:Default: none
:Valid Values: :math:`\gt 0`


smoothing_weight
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

:Type: real
:Default: 2.0
:Valid Values: :math:`\gt 0`


smoothing_max_iter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

:Type: integer
:Default: none
:Valid Values: :math:`\geq 0`


verbosity_level
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

:Type: integer
:Default: 2
:Valid Values: :math:`\geq 0`
