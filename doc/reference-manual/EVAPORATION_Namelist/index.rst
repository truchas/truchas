.. _EVAPORATION_Namelist:

.. toctree::
   :maxdepth: 1

EVAPORATION Namelist (Experimental)
=====================================

Overview
------------
This namelist defines a special heat flux boundary condition that models heat loss due to the evaporation of material. Its intended use is in the simulation of additive manufacturing or welding processes where the laser heat source can produce localized surface temperatures approaching and exceeding the boiling temperature. The form of the heat flux is an Arrhenius-type function.

.. math::

   f(T) = A T^\beta e^{E_a/{RT}}

where **T** is temperature, **R** the gas constant, and **A**, :math:`\beta`, and :math:`E_a` are model parameters defined by input. When using this model, the simulation should be using Kelvin for temperature.   

EVAPORATION Namelist Features
--------------------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Single

Components
------------
* :ref:`Face_Set_IDs<Evap_FSD>`
* :ref:`Prefactor<Evap_P>`
* :ref:`Temp_Exponent<Evap_TE>`
* :ref:`Activation_Energy<Evap_AE>`

.. _Evap_FSD:

Face_Set_IDs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A list of face set IDs that define the subset of the boundary where the evaporation boundary condition model will be imposed.
| **Type**        : a list of up to 32 integers
| **Default**     : none
| **Valid Values**: any valid mesh face set ID
| **Notes**       : Exodus II mesh side sets are interpreted by Truchas as face sets having the same IDs.

.. _Evap_P:

Prefactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The prefactor A.
| **Type**        : real
| **Default**     : none

.. _Evap_TE:

Temp_Exponent
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The temperature exponent :math:`\beta`
| **Type**        : real
| **Default**     : none

.. _Evap_AE:

Activation_Energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The activation energy :math:`E_a`. This value must be specified in Joules per mole units, regardless of the units used elsewhere in the simulation.
| **Type**        : real
| **Default**     : none

