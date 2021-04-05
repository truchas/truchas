.. _DS_SOURCE_Namelist:

.. toctree::
   :maxdepth: 1

DS_SOURCE Namelist
====================

The DS_SOURCE namelist is used to define external volumetric sources for species transport model.

Overview
------------

DS_SOURCE Namelist Features
----------------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Multiple

Components
------------
* :ref:`Equation<DS_E>`
* :ref:`Cell_Set_IDs<DS_CSI>`
* :ref:`Source_Constant<DS_SC>`
* :ref:`Source_Function<DS_SF>`

.. _DS_E:

Equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of the equation this source applies to.
| **Type**        : string
| **Default**     : none
| **Valid Values**: "concentration1","concentration2", . . .
| **Note**        : Any name may be specified, but Truchas will only look for and use DS_SOURCE namelists with the indicated equation names; any others are silently ignored.

.. _DS_CSI:

Cell_Set_IDs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A list of cell set IDs that define the subdomain where the source is applied.
| **Type**        : a list of up to 32 integers
| **Default**     : none
| **Valid Values**: any valid mesh cell set ID
| **Note**        : Different instances of this namelist with a given :ref:`Equation<DS_E>` value must apply to disjoint subdomains; overlapping of source functions is not allowed. Exodus II mesh element blocks are interpreted by Truchas as cell sets having the same IDs.

.. _DS_SC:

Source_Constant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The constant value of the source function.
| **Type**        : real
| **Default**     : none
| **Note**        : Either Source_Constant or :ref:`Source_Function<DS_SF>` must be specified, but not both.

.. _DS_SF:

Source_Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist that defines the source function. That function is expected to be a function of (t,x,y,z).
| **Type**        : string
| **Default**     : none
| **Note**        : Either Source_Function or :ref:`Source_Constant<DS_SC>` must be specified, but not both.

