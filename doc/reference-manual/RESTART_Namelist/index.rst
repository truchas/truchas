.. _RESTART_Namelist:

.. toctree::
   :maxdepth: 1

RESTART Namelist
============================

Overview
------------

Truchas is able to use data from a previous calculation to initialize a new calculation. Such a restart calculation is invoked by using the ``-r`` command line argument (:numref:`Table %s <truchas_command_line_option>`) to the executable with the path name of the restart data file. By default, all appropriate data from the restart file is used. This optional namelist provides variables to limit the restart data that will be used.

RESTART Namelist Features
----------------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Single

Components
------------
* :ref:`Ignore_T<R_IT>`
* :ref:`Ignore_Dt<R_ID>`
* :ref:`Ignore_EM_Heat<R_IJH>`

.. _R_IT:

Ignore_T
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : When restarting, the initial time and starting cycle count are normally extracted from the restart file. If this flag is true, then those values are ignored and the first value of the :ref:`Output_T<O_OT>` array is used as the initial time and the cycle count starts at 0, as happens with a non-restart run.
| **Type**        : logical
| **Default**     : false

.. _R_ID:

Ignore_Dt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : When restarting, the initial time step size is normally extracted from the restart file. If this flag is true, then that value is ignored and the value specified by :ref:`Dt_Init<NUMERICS_DTI>` from the :ref:`NUMERICS<NUMERICS_Namelist>` namelist is used instead. Note that if :ref:`Dt_Constant<NUMERICS_DTC>` in the :ref:`NUMERICS<NUMERICS_Namelist>` namelist is specified then its value is used regardless.
| **Type**        : logical
| **Default**     : false

.. _R_IJH:

Ignore_EM_Heat
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : If this flag is true, the electromagnetic heat data in the
   restart file (if any) will be ignored when initializing the code. This variable
   is only relevant for restart calculations with
   :ref:`induction_heating<physics-ih>` or :ref:`microwave_heating<physics-mwh>`
   enabled in the PHYSICS namelist.
| **Type**        : logical
| **Default**     : false
