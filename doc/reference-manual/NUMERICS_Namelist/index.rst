.. _NUMERICS_Namelist:

.. toctree::
   :maxdepth: 1

NUMERICS Namelist
===================

The NUMERICS namelist specifies general numerical parameters not specific to
any particular physics, especially those controlling the overall time stepping
of the Truchas model.


.. admonition:: Namelist Usage

   :Required/Optional: Required
   :Single/Multiple Instances: Single


Namelist Variables
------------------

.. contents::
   :local:


.. _NUMERICS_CM:

Cycle_Max
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description**       : The maximum cycle number allowed in the given simulation, where one “cycle” corresponds to one integration time step over all physical model equations. The simulation will terminate gracefully when the cycle number reaches Cycle_Max.
| **Type**              : integer
| **Default**           : 1000000 
| **Valid Values**      : :math:`(0,\infty)`
| **Notes**             : A simulation will also terminate gracefully if the last entry in the :ref:`Output_T<O_OT>` input variable (real) array (in the :ref:`OUTPUTS<OUTPUTS_Namelist>` namelist) is exceeded by the current simulation time :ref:`t<NUMERICS_T>`.

.. _NUMERICS_CN:

Cycle_Number
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description**       : The cycle number to be used just prior to the first actual cycle (time step) taken in the given simulation. The first simulation cycle number is then taken to be cycle_number+ 1.
| **Type**              : integer
| **Default**           : 0
| **Valid Values**      : [1, :math:`\infty`)
| **Notes**             : The default value of 0 results in the first computational cycle taken to be 1. This input variable is most useful when starting simulations from a restart file that already contains a restart cycle number > 0.

.. _NUMERICS_DTC:

Dt_Constant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description**       : A constant integration time step value to be used for all time steps in the simulation.
| **Physical dimension**: T
| **Type**              : real
| **Default**           : none 
| **Valid Values**      : :math:`(0,\infty)`
| **Notes**             : This :ref:`Dt_Constant<NUMERICS_DTC>` input variable should be used with extreme caution, as its specification overrules all other time step choices that are controlled by linear stability and accuracy considerations. In particular, :ref:`NUMERICS<NUMERICS_Namelist>` namelist input variables :ref:`Dt_Grow<NUMERICS_DTG>`, :ref:`Dt_Init<NUMERICS_DTI>`, :ref:`Dt_Max<NUMERICS_DTMax>`, and :ref:`Dt_Min<NUMERICS_DTMin>` are all ignored if :ref:`Dt_Constant<NUMERICS_DTC>` is specified. Use of :ref:`Dt_Constant<NUMERICS_DTC>` is therefore advised only for controlled numerical algorithm experiments, and not for simulations intended for applications analysis and validation.

.. _NUMERICS_DTG:

Dt_Grow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description**       : A factor to multiply the current integration time step :math:`(\delta t)` for the purpose of estimating the time step used during the next cycle :math:`(\delta t_g =` Dt_Grow :math:`∗\delta t)`. This candidate time step :math:`(\delta t_g)` is chosen only if all other currently active time step restrictions have not already limited its value below :math:`(\delta t_g)`.
| **Physical dimension**: dimensionless
| **Type**              : real
| **Default**           : 1.05 
| **Valid Values**      : [1, :math:`\infty`)
| **Notes**             : :ref:`Dt_Grow<NUMERICS_DTG>` is ignored if :ref:`Dt_Constant<NUMERICS_DTC>` is specified.

.. _NUMERICS_DTI:

Dt_Init
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description**       : A Integration time step value used for the first computational cycle.
| **Physical dimension**: T
| **Type**              : real
| **Default**           : :math:`10^{-6}`
| **Valid Values**      : (1, :math:`\infty`)
| **Notes**             : The default value of :ref:`Dt_Init<NUMERICS_DTI>` is completely arbitrary, as it is always problem dependent. Unless aconstant time step is desired (via specification of :ref:`Dt_Constant<NUMERICS_DTC>`), then, a value for :ref:`Dt_Init<NUMERICS_DTI>` should be specified (instead of relying on the default) that is consistent with the time scales of interest for the simulation at hand. A general rule of thumb is to set :ref:`Dt_Init<NUMERICS_DTI>` smaller than the time step ultimately desired, thereby allowing the time step to grow gradually and steadily (say, over 10-20 cycles) to the desired time step value.

For restart calculations, the initial time step size is extracted from the restart file and used instead of :ref:`Dt_Init<NUMERICS_DTI>`, unless :ref:`Ignore_Dt<R_ID>` in the :ref:`RESTART<RESTART_Namelist>` namelist has been set to true. :ref:`Dt_Init<NUMERICS_DTI>` is ignored if :ref:`Dt_Constant<NUMERICS_DTC>` is specified.

.. _NUMERICS_DTMax:

Dt_Max
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description**       : Maximum allowable value for the time step
| **Physical dimension**: T
| **Type**              : real
| **Default**           : 10
| **Valid Values**      : (0, :math:`\infty`)
| **Notes**             : The time step is not allowed to exceed this value, even if other accuracy- or stability-based criteria would permit it. The default value of :ref:`Dt_Max<NUMERICS_DTMax>` is completely arbitrary, as it is always problem dependent. Unless a constant time step is desired (via specification of :ref:`Dt_Constant<NUMERICS_DTC>`), then, a value for :ref:`Dt_Max<NUMERICS_DTMax>` should be specified (instead of relying on the default) that is consistent with the timescales of interest for the simulation at hand.

:ref:`Dt_Max<NUMERICS_DTMax>` is ignored if :ref:`Dt_Constant<NUMERICS_DTC>` is specified.

.. _NUMERICS_DTMin:

Dt_Min
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description**       : Minimum allowable value for the time step
| **Physical dimension**: T
| **Type**              : real
| **Default**           : :math:`10^{-6}`
| **Valid Values**      : (0, :math:`\infty`)
| **Notes**             : If the time step falls below this value, the user is informed as such and the simulation terminates gracefully. Termination occurs at this condition because is it very probable that either numerical algorithm or physical model problems have occurred and the simulation is unable to recover. The default value of :ref:`Dt_Min<NUMERICS_DTMin>` is completely arbitrary, as it is always problem dependent. Unless a constant time step is desired (via specification of :ref:`Dt_Constant<NUMERICS_DTC>`), then, a value for :ref:`Dt_Min<NUMERICS_DTMin>` should be specified (instead of relying on the default) that is consistent with the time scales of interest for the simulation at hand. A good rule of thumb to use for :ref:`Dt_Min<NUMERICS_DTMin>` is to use a value several orders of magnitude below the time scale of interest.

:ref:`Dt_Min<NUMERICS_DTMin>` is ignored if :ref:`Dt_Constant<NUMERICS_DTC>` is specified.

.. _NUMERICS_T:

t
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description**       : The simulation time to be used at the beginning of the first computational cycle.
| **Physical dimension**: T
| **Type**              : real
| **Default**           : 0
