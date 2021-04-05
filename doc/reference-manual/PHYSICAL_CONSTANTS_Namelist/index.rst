.. _PHYSICAL_CONSTANTS_Namelist:

.. toctree::
   :maxdepth: 1

PHYSICAL_CONSTANTS Namelist
=============================

Overview
------------
The values of physical constants used in Truchas’ physics models are set through this namelist. The default for all these constants is their value in SI units. If a different system of units is used, these may need to be assigned the appropriate values.

PHYSICAL_CONSTANTS Namelist Features
--------------------------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Single

Components
------------
* :ref:`Absolute_Zero<PhyCo_AZ>`
* :ref:`Stefan_Boltzmann<PhyCo_SB>`
* :ref:`Vacuum_Permeability<PhyCo_VPm>`
* :ref:`Vacuum_Permittivity<PhyCo_VPt>`

.. _PhyCo_AZ:

Absolute_Zero
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The value of absolute-zero in the temperature scale used. The default value is 0 for Kelvin. Centigrade would use the value −273.15, for example. This constant is used by thermal radiation boundary conditions and enclosure (view factor) radiation.
| **Physical dimension** :math:`\Theta`
| **Type**        : real
| **Default**     : 0 K
| **Valid Values**: any value

.. _PhyCo_SB:

Stefan_Boltzmann
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Stefan-Boltzmann constant for thermal radiation. The default is its value in SI units. This constant is used by thermal radiation boundary condititions and enclosure (view factor) radiation.
| **Physical dimension** :math:`E/{T.L^2.\Theta^4}`
| **Type**        : real
| **Default**     ::math:`5.67 . 10^{-8} W/{m^2.K}`
| **Valid Values**: any positive value

.. _PhyCo_VPm:

Vacuum_Permeability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The magnetic permeability of free space. The default is its value in SI units. This parameter is used by the electromagnetics solver.
| **Physical dimension** :math:`M.L.T^{-2}.I^{-2}`
| **Type**        : real
| **Default**     ::math:`4\pi . 10^{-7} H/m`
| **Valid Values**: any positive value

.. _PhyCo_VPt:

Vacuum_Permittivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The electric permittivity of free space. The default is its value in SI units. This parameter is used by the electromagnetics solver.
| **Physical dimension** :math:`M^{-1}.L^{-3}.T^{4}.I^{2}`
| **Type**        : real
| **Default**     ::math:`8.854188 x 10^{-12} F/m`
| **Valid Values**: any positive value