.. _DED_HEAD_Namelist:

.. toctree::
   :maxdepth: 1

DED_HEAD Namelist (Experimental)
================================

Overview
------------
The variables in an :ref:`DED_HEAD<DED_HEAD_Namelist>` namelist specify the physical characteristics of heat source for additive manufacturing and welding simulations.


DED_HEAD Namelist Features
-----------------------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Single

Components
------------
* :ref:`Toolpath<TP_DH>`
* :ref:`Laser_type<LT_DH>`
* :ref:`Laser_time_const<LRC_DH>`
* :ref:`Laser_absorptivity<LA_DH>`
* :ref:`Laser_power<LP_DH>`
* :ref:`Laser_power_func<LPF_DH>`
* :ref:`Laser_sigma<LS_DH>`

.. _TP_DH:

Toolpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Unique name of the vector function for the laser tool path referenced to :ref:`name<VF_Name>` in the :ref:`VFUNCTION<VFUNCTION_Namelist>`.
| **Type**        : A case-sensitive string of up to 31 characters
| **Default**     : none

.. _LT_DH:

Laser_type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Type of the laser heat source can be specified using this variable. A simple, constant gaussian profile is specified using `guassian` and `beam` is used to specify a complex laser profile.

The mathematical definition of `gaussian` is as follows
.. math::

        P(x,y) = P_o/{2*\pi*{\sigma}^2}.exp{-{dx^2+dy^2}/{2*{\sigma}^2}}

where :math:`P_o` is the power specified by :ref:`Laser_power<LP_DH>` or :ref:`Laser_power_func<LPF_DH>`, :math:`\sigma` is the value specified by :ref:`laser_sigma<LS_DH>`

The mathematical definition of `beam` is as follows
.. math::
        P(x,y,z) = P_o/{\pi*c_o}.exp{-{dx^2 + dy^2}/c_o}
        
        c_o = w_r^2/2*[1 + ({dz/{\pi*w_r^2}}*{\lambda * msq}^2)]

where :math:`P_o` is the power specified by :ref:`Laser_power<LP_DH>` or :ref:`Laser_power_func<LPF_DH>`, :math:`w_r` is the waist radius, :math:`\lambda` is the wavelength and math:`msq` is ??????? 


| **Type**        : string
| **Default**     : none
| **Valid Values**: `gaussian`, `beam` 

.. _LRC_DH:

Laser_time_const
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : This variable is used to specify the ramp up time for laser from :math:`0` to the value specified by :ref:`Laser_power<LP_DH>` or :ref:`Laser_power_func<LPF_DH>`. This helps in introducing transition time for turning on and off of the laser power thereby preventing abrupt transition which is undesirable for the ODE solver. 
| **Physical dimension**: :math:`T`
| **Type**        : real
| **Default**     : none

.. _LA_DH:

Laser_absorptivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Absorptivity (percentage of laser power absorbed by the interacting material) is spceified using this variable.
| **Type**        : real
| **Default**     : none
| **Valid Values**: [0,1]

.. _LP_DH:

Laser_power
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description**         : A constant value of laser power is specified using this variable. 
| **Physical dimension**  : :math:`M.L^2.T^{-3}`
| **Type**                : real
| **Default**             : none
| **Valid Values**        : :math:`(0, \infty)`

.. _LPF_DH:

Laser_power_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A unique function name referenced to the variable :ref:`name<FUNC_name>` in :ref:`FUNCTION<FUNCTION_Namelist>` namelist. The laser power can be defined as a function of time using this variable.
| **Type**        : A case-sensitive string of up to 31 characters. 
| **Default**     : none

.. _LS_DH:

Laser_sigma
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Laser sigma is the following function of laser diameter :math:`laser sigma = {diameter/2}*{1/sqrt{-2*log{1-fraction power target}}`. `Fraction power target` is the fraction of total power that falls within the laser diameter.
| **Physical dimension**: :math:`L`
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`(0, \infty)`
