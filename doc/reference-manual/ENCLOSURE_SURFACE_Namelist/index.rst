.. _ENCLOSURE_SURFACE_Namelist:

.. toctree::
   :maxdepth: 1

ENCLOSURE_SURFACE Namelist 
==============================

Overview
------------

ENCLOSURE_SURFACE Namelist Features
---------------------------------------
| **Required/Optional        :** Required when :ref:`ENCLOSURE_RADIATION<ENCLOSURE_RADIATION_Namelist>` namelists are active.
| **Single/Multiple Instances:** Multiple

Components
------------
* :ref:`Name<ES_N>`
* :ref:`Enclosure_Name<ES_EN>`
* :ref:`Face_Block_IDs<ES_FBI>`
* :ref:`Emmisivity_Constant<ES_EC>`
* :ref:`Emmisivity_Function<ES_EF>`

.. _ES_N:

Equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A unique name for this enclosure surface.
| **Type**        : string (31 characters max)
| **Default**     : none

.. _ES_EN:

Enclosure_Name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of the :ref:`ENCLOSURE_RADIATION<ENCLOSURE_RADIATION_Namelist>` namelist that defines the enclosure radiation system to which this surface belongs.
| **Type**        : string
| **Default**     : none

.. _ES_FBI:

Face_Block_IDs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A list of face block IDs that define this enclosure surface
| **Type**        : a list of up to 32 integers
| **Default**     : none
| **Valid Values**: any valid face block ID from the enclosure file
| **Notes**       : The surface faces in an enclosure file are divided into blocks. When the **genre** program from the **RadE** tool suite is used to create this file, these blocks are automatically generated; each block corresponds to the Exodus II mesh side set used to defined it, and the side set ID is assigned as the face block ID. In this case, then, the IDs that can be specified here will be certain side set IDs from the Truchas Exodus II mesh file.


.. _ES_EC:

Emmisivity_Constant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The constant emissivity value for this enclosure surface.
| **Physical dimensions**: dimensionless
| **Type**        : real
| **Default**     : none
| **Notes**       : Either Emissivity_Constant or :ref:`Emissivity_Function<ES_EF>` must be specified, but not both.

.. _ES_EF:

Emmisivity_Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist that defines the emissivity function. That function is expected to be a function of (t,x,y,z).
| **Type**        : string
| **Default**     : none
| **Notes**       : Either Emissivity_Function or :ref:`Emissivity_Constant<ES_EC>` must be specified, but not both.