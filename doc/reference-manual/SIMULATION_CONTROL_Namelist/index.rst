.. _SIMULATION_CONTROL_Namelist:

.. toctree::
   :maxdepth: 1

SIMULATION_CONTROL Namelist (Experimental)
============================================

Overview
------------

There may be points in time during a simulation when something changes abruptly; a boundary condition or source turns on/off, or a physics model is enabled/disabled, for example. In such circumstances it is best to hit these times precisely with a time step and then continue from that point with a reduced step size appropriate to resolving the time transients that result from the impulsive forcing of the model—in essence, to split the simulation seamlessly into a sequence of phases where each phase is a new simulation whose initial state is the final state of the preceding phase. This experimental namelist provides a means for achieving this. The start time of each additional phase subsequent to the initial phase is specified using the :ref:`Phase_Start_Times<SC_PST>` array, and the initial step size using either the :ref:`Phase_Init_Dt<SC_PID>` or :ref:`Phase_Init_Dt_Factor<SC_PIDF>` variables. Truchas will hit those times precisely with a time step, smoothly adjusting the step size in advance to avoid abrupt step size changes, and then effectively **restart** the timestepping. Currently this only effects the second-order diffusion solver which maintains a (smooth) history of states at recent time steps. When restarting that history is deleted and time stepping begins fresh using only the current state.

SIMULATION_CONTROL Namelist Features
---------------------------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Single

Components
------------
* :ref:`Event_Lookahead<SC_EL>`
* :ref:`Phase_Init_Dt<SC_PID>`
* :ref:`Phase_Init_Dt_Factor<SC_PIDF>`
* :ref:`Phase_Start_Times<SC_PST>`

.. _SC_EL:

Equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : When approaching an event time, such as an output time, the time step sizes are gradually adjusted to hit the time precisely. This variable specifies the number of steps over which this occurs.
| **Type**        : integer
| **Default**     : 5
| **Valid Values**: :math:`\geq 2`
| **Note**        : A value of 1 would mean no step size adjustment until a step would go beyond the event time. This can result in needing to take an arbitrarily small time step and one smaller than the minimum allowed, and thus this is not allowed. The greater the lookahead, the more gradually the time step will be reduced.

.. _SC_PID:

Phase_Init_Dt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The initial time step size to use for each of the simulation phases. This is the analog of :ref:`Dt_Init<NUMERICS_DTI>`, which is used for the initial (and default) phase of the simulation. Either :ref:`Phase_Init_Dt<SC_PID>` or :ref:`Phase_Init_Dt_Factor<SC_PIDF>` must be specified, but not both.
| **Type**        : real
| **Default**     : none

.. _SC_PIDF:

Phase_Init_Dt_Factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The initial time step size used for each of the simulation phases is this factor times the last step size of the preceding phase. Either :ref:`Phase_Init_Dt<SC_PID>` or :ref:`Phase_Init_Dt_Factor<SC_PIDF>` must be specified, but not both.
| **Type**        : real
| **Default**     : none

.. _SC_PST:

Phase_Start_Times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The list of starting times of each of the phases.
| **Type**        : real array
| **Default**     : none
| **Note**        : The initial simulation phase, which is otherwise the only phase, need not be included in this list, though it may be. The provided list of times is sorted, and the first time greater than the initial time is taken as the start of the first phase following the default initial phase; earlier times are ignored.
