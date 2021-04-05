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
* :ref:`Int_Output_Dt_Multiplier<O_IODM>`
* :ref:`Output_Dt<O_OD>`
* :ref:`Output_Dt_Multiplier<O_ODM>`
* :ref:`Output_T<O_OT>`
* :ref:`Probe_Output_Cycle_Multiplier<O_POCM>`
* :ref:`Short_Output_Dt_Multiplier<O_SODM>`
* :ref:`Move_Block_IDs<O_MBI>`
* :ref:`Move_Toolpath_Name<O_MTN>`

.. _O_IODM:

Int_Output_Dt_Multiplier
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Factor multiplying :ref:`Output_Dt<O_OD>` for time interval to write interface output data.
| **Type**        : integer array
| **Default**     : none
| **Valid Values**: :math:`\geq 0`

.. _O_OD:

Output_Dt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Output time interval for each output time span.
| **Physical Dimension**: T
| **Type**        : real array
| **Default**     : none
| **Valid Values**: :math:`\gt 0`

.. _O_OT:

Output_T
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A sequence of time values for defining time spans that have distinct output time intervals. The last time is the problem end time.
| **Physical Dimension**: T
| **Type**        : real array
| **Default**     : none
| **Valid Values**: strictly increasing sequence of two or more values.

.. _O_POCM:

Probe_Output_Cycle_Multiplier
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Factor multiplying truchas cycle to determine frequency of writing probe output.
| **Type**        : integer
| **Default**     : 1
| **Valid Values**: :math:`\gt 0`

.. _O_SODM:

Short_Output_Dt_Multiplier
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Factor multiplying :ref:`Output_Dt<O_OD>` for time interval to write short edits.
| **Type**        : integer array
| **Default**     : 0
| **Valid Values**: :math:`\geq 0`

.. _O_ODM:

Output_Dt_Multiplier
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Factor multiplying :ref:`Output_Dt<O_OD>` for time interval to write output.
| **Type**        : integer array
| **Default**     : 0
| **Valid Values**: :math:`\geq 0`

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