.. _PHYSICAL_CONSTANTS_Namelist:

.. toctree::
   :maxdepth: 1

PHYSICAL_CONSTANTS Namelist
=============================
The value of physical constants used in Truchas physics models are set using
this namelist. The default for all these constants is their value in SI units.
When the input for a problem uses a different system of units, these may need
to be assigned the corresponding value in that system of units.

.. admonition:: Namelist Usage

   :Required/Optional: Optional
   :Single/Multiple Instances: Single

Namelist Variables
------------------

.. contents::
   :local:
   
.. _PhyCo_AZ:

absolute_zero
^^^^^^^^^^^^^
The value of absolute zero temperature. This constant is used by thermal
radiation boundary conditions and the enclosure radiation model.

:Default: :math:`0~\text{K}`
:Note: Use :math:`-273.15` for temperature in centigrade.

.. _PhyCo_UGC:

gas_constant
^^^^^^^^^^^^
The universal gas constant.
This constant is used by the evaporation heat flux model.

:Default: :math:`8.31446261815324~\text{J}/\text{mol}\cdot\text{K}`

.. _PhyCo_SB:

stefan_boltzmann
^^^^^^^^^^^^^^^^
The Stefan-Boltzmann constant :math:`\sigma`.
This constant is used by thermal radiation boundary condititions and the
enclosure radiation model.

:Default: :math:`5.67\times10^{-8}~\text{W}/\text{m}^2\cdot\text{K}^4`

.. _PhyCo_VPm:

vacuum_permeability
^^^^^^^^^^^^^^^^^^^
The magnetic permeability of free space.
This constant is used by the electromagnetics solver.

:Default: :math:`4\pi\times10^{-7}~\text{N}/\text{A}^2`

.. _PhyCo_VPt:

vacuum_permittivity
^^^^^^^^^^^^^^^^^^^
The electric permittivity of free space.
This constant is used by the electromagnetics solver.

:Default: :math:`8.854188\times10^{-12}~\text{F}/\text{m}`
