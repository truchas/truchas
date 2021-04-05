.. _TURBULENCE_Namelist:

.. toctree::
   :maxdepth: 1

TURBULENCE Namelist
============================

Overview
------------
The presence of the TURBULENCE namelist enables a simple algebraic turbulence model for viscous flow problems. The turbulent kinetic viscosity :math:`\nu_t = (\mu_t/\rho)` is taken to be

.. math::
   \nu_t = c_{\mu}k^{\frac{1}{2}}l 
   :label: turb_eq_1
   

where :math:`c_{\mu}` is the proportionality constant, :math:`l` is a length scale corresponding to the eddy size, and

.. math::
   k = f.\frac{1}{2}.u^2 
   :label: turb_eq_2
   

is the local turbulent kinetic energy per unit mass, modeled as a fraction :math:`f` of the mean kinetic energy. The namelist variables give values for the model parameters. See Truchas Physics and Algorithms for more details.


TURBULENCE Namelist Features
----------------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Single

Components
------------
* :ref:`CMU<Turb_CMU>`
* :ref:`KE_fraction<Turb_KEf>`
* :ref:`Length<Turb_Length>`

.. _Turb_CMU:

CMU
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Value of the parameter :math:`c_{\mu}` in :eq:`turb_eq_1`.
| **Type**        : real
| **Default**     : 0.05
| **Valid Values**: :math:`\gt 0`
| **Notes**       : The default value is appropriate in most situations.

.. _Turb_KEF:

KE_fraction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Value of the parameter :math:`f` in :eq:`turb_eq_2`.
| **Type**        : real
| **Default**     : 0.1
| **Valid Values**: (0,1)
| **Notes**       : The default value is appropriate in most situations.

.. _Turb_Length:

Length
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : Value of the parameter :math:`l` in :eq:`turb_eq_1`.
| **Physical dimension**: L
| **Type**        : real
| **Default**     : none
| **Valid Values**: :math:`\gt 0`
| **Notes**       : The value should correspond to the size of the turbulent eddies. In turbulent pipe flow, for example, this would be one third the pipe radius.

