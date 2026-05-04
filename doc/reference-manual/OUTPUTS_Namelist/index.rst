.. _OUTPUTS_Namelist:

.. toctree::
   :maxdepth: 1

OUTPUTS Namelist
===================

Overview
------------

The OUTPUTS namelist defines the problem end time and various output options.

OUTPUTS Namelist Features
----------------------------
| **Required/Optional        :** Required
| **Single/Multiple Instances:** Single

Components
------------
* :ref:`Times<O_Times>`
* :ref:`Subintervals<O_Subintervals>`
* :ref:`Fields<O_Fields>`
* :ref:`Move_Block_IDs<O_MBI>`
* :ref:`Move_Toolpath_Name<O_MTN>`

.. _O_Times:

Times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A sequence of time values for defining time spans that have distinct output subdivision counts. The last time is the problem end time.
| **Physical Dimension**: T
| **Type**        : real array
| **Default**     : none
| **Valid Values**: strictly increasing sequence of two or more values.

.. _O_Subintervals:

Subintervals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Number of equal subintervals for each time span. Output is written at the endpoint of each subinterval.
| **Type**        : integer array
| **Default**     : 1
| **Valid Values**: :math:`\gt 0`

.. _O_Fields:

Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A list of visualization fields to write to the VTKHDF output file.
| **Type**        : string array
| **Default**     : all available visualization fields
| **Valid Values**: ``'all'`` or one or more field names from the list below. Field names are case sensitive.
| **Notes**       : Omitting this variable currently writes all available fields for compatibility with earlier output behavior. To explicitly request that behavior, specify ``fields = 'all'``. The value ``'all'`` may not be combined with other field names. A requested field that is unavailable for the active physics is omitted with a warning. In parallel runs, the default and ``'all'`` selections include ``process_id``.

The following field names are available, depending on the active physics and
model configuration:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Field name
     - Description
   * - ``temperature``
     - Cell temperature.
   * - ``enthalpy``
     - Cell enthalpy.
   * - ``phi_N``
     - Species field ``N``, for example ``phi_1`` or ``phi_2``.
   * - ``vfrac``
     - Shorthand for all material volume-fraction fields.
   * - ``vfrac_<phase-name>``
     - Volume fraction for one material phase, for example ``vfrac_water``.
   * - ``process_id``
     - MPI process rank that owns each cell. Written only for parallel runs.
   * - ``em_heat``
     - Electromagnetic heat generation.
   * - ``velocity``
     - Cell-centered flow velocity.
   * - ``pressure``
     - Cell-centered flow pressure.
   * - ``displacement``
     - Solid-mechanics displacement.
   * - ``strain``
     - Solid-mechanics strain tensor.
   * - ``thermal_strain``
     - Solid-mechanics thermal strain tensor.
   * - ``stress``
     - Solid-mechanics stress tensor.
   * - ``rotation``
     - Solid-mechanics rotation.
   * - ``contact_normal_displ``
     - Normal displacement on solid-mechanics contact interface blocks.
   * - ``contact_pressure``
     - Contact pressure on solid-mechanics contact interface blocks.
   * - ``plastic_strain``
     - Solid-mechanics plastic strain tensor; available when viscoplasticity is enabled.
   * - ``plastic_strain_rate``
     - Solid-mechanics plastic strain rate; available when viscoplasticity is enabled.
   * - ``ustruc``
     - Shorthand for all configured microstructure fields.

The ``vfrac`` and ``ustruc`` entries are selectors rather than dataset names.
The corresponding VTKHDF datasets are named after the expanded fields, such as
``vfrac_water`` or ``ustruc1_gl_t_sol``.

.. _O_MBI:

Move_Block_IDs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A list of element block IDs that are associated with a translation written to the output file. Use :ref:`Move_Toolpath_Name<O_MTN>` to specify the translation.
| **Type**        : a list of up to 32 integers
| **Default**     : none
| **Notes**       : Use of this feature does not alter the mesh data that is written to the HDF5 output file. It merely adds some additional data that associates a time-dependent translation with element blocks. Use of the data, if any, is left to users of the file. At this time the Paraview Truchas output reader (postversion 5.2) uses this information to translate the mesh blocks for visualization.

.. _O_MTN:

Move_Toolpath_Name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`TOOLPATH<TOOLPATH_Namelist>` namelist that defines the translation to apply to the element blocks given by :ref:`Move_Block_IDs<O_MBI>`.
| **Type**        : string
| **Default**     : none
