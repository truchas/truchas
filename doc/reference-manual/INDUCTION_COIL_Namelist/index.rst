.. _INDUCTION_COIL_Namelist:

.. toctree::
   :maxdepth: 1

INDUCTION_COIL Namelist
============================

Overview
------------
The variables in an :ref:`INDUCTION_COIL<INDUCTION_COIL_Namelist>` namelist specify the physical characteristics of an induction coil that is to produce an external magnetic field to drive the electromagnetic Joule heat calculation. :numref:`Figure %s<fig_induc_coil_schematic>` shows the idealized model of a coil that is used to analytically evaluate the driving field. The coil axis is assumed to be oriented with the problem :ref:`symmetry axis<EM_SA>` as defined in the :ref:`ELECTROMAGNETICS<ELECTROMAGNETICS_Namelist>` namelist. Multiple coils may be specified; the net driving field is the superposition of the fields due to the individual coils, and a spatially :ref:`uniform field<EM_US>` that can be specified in the :ref:`ELECTROMAGNETICS<ELECTROMAGNETICS_Namelist>` namelist. The coils carry a sinusoidally-varying current with a common frequency and phase. The :ref:`frequency<EM_SF>` is specified in the :ref:`ELECTROMAGNETICS<ELECTROMAGNETICS_Namelist>` namelist, while the phase value is irrelevant due to the time averaging of the calculated Joule heat. Each coil, however, has an independent current amplitude which is specified here. In addition, the current and frequency may be piecewise constant functions of time.

.. _fig_induc_coil_schematic:
.. figure:: images/ic-coil.jpg
   :width: 650px
   :align: center
   
   Physical 4-turn helical coil with extended wire cross section (left), and its idealized model as a stacked array of circular current loops (right).

INDUCTION_COIL Namelist Features
-----------------------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Multiple

Components
------------
* :ref:`Center<IC_Cen>`
* :ref:`Current<IC_Cur>`
* :ref:`Length<IC_Len>`
* :ref:`NTurns<IC_NT>`
* :ref:`Radius<IC_Rad>`

.. _IC_Cen:

Center
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A 3-vector :math:`x_0` giving the position of the center of the coil; cf. :numref:`Figure %s<fig_induc_coil_schematic>`.
| **Physical dimension**: L
| **Type**        : real
| **Default**     : :math:`(0, 0, 0)`
| **Valid Values**: any 3-vector

.. _IC_Cur:

Current
^^^^^^^^^
| **Description** : Amplitude of the sinusoidally-varying current in the coil.
| **Physical dimension**: I
| **Type**        : real
| **Default**     : none
| **Valid Values**: Any single value, or any sequence of values.
| **Notes**       : A sequence of up to 32 values may be assigned to this variable in order to specify a time-dependent current amplitude; see :ref:`Source_Times<EM_ST>` in the :ref:`ELECTROMAGNETICS<ELECTROMAGNETICS_Namelist>` namelist and :numref:`Figure %s<fig_rm_em_pc>`.

.. _IC_Len:

Length
^^^^^^^^^
| **Description** : Length :math:`l` of the coil; cf. :numref:`Figure %s<fig_induc_coil_schematic>`
| **Physical dimension**: L
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`(0, \infty)`
| **Notes**       : Length is not required, nor meaningful, if :ref:`NTurns<IC_NT>` is :math:`1`.

.. _IC_NT:

NTurns
^^^^^^^^^
| **Description** : Number of turns of the coil; cf. :numref:`Figure %s<fig_induc_coil_schematic>`
| **Type**        : integer
| **Default**     : none
| **Valid Values**: Any positive integer.

.. _IC_Rad:

Radius
^^^^^^^^^
| **Description** : Radius :math:`r` of the coil; cf. :numref:`Figure %s<fig_induc_coil_schematic>`
| **Physical dimension**: L
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`(0, \infty)`