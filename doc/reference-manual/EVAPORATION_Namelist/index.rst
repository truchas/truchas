.. _EVAPORATION_Namelist:

.. toctree::
   :maxdepth: 1

EVAPORATION Namelist (Experimental)
=====================================
This namelist applies a special boundary heat flux that models heat loss due
to the evaporation of a material. Its intended use is in the simulation of
additive manufacturing or welding processes where the heat source can produce
localized surface temperatures approaching and exceeding the vaporization
temperature.

The heat flux due to evaporation,

.. math::

   q_{\text{evap}} = L \dot{m},

is expressed in terms of the evaporated mass flux :math:`\dot{m}` and
latent heat of vaporization :math:`L`. The mass flux (mass per unit
area per unit time) is determined from the Hertz-Knudsen equation
:footcite:`he2003alloying`

.. math::
   \dot{m} = \lambda p_v(T)\sqrt{\frac{M}{2\pi R T}}

where :math:`M` is the molar mass of the evaporating material, :math:`R` is
the universal gas constant, :math:`p_v` is the vapor pressure of the material
at temperature :math:`T`, and :math:`\lambda` is an empirical constant that
corrects for the condensation of a portion of the vaporized atoms. The vapor
pressure is given by the Clausius-Claperyon equation :footcite:`velasco2009clausius`

.. math::
   p_v(T) = p_0 \exp\Bigl[\frac{ML}{R}\Bigl(\frac1{T}-\frac1{T_0}\Bigr)\Bigr]

where :math:`p_0` is the pressure of the ambient atmosphere, and
:math:`T_0` is the vaporization temperature at that pressure.

.. caution::

   This model is only compatible with problems using Kelvin for temperature.
   
   If any other problem units differ from SI units, the value of the universal
   gas constant :math:`R` must be redefined accordingly using the
   :ref:`gas_constant<PhyCo_UGC>` variable of the
   :ref:`PHYSICAL_CONSTANTS<PHYSICAL_CONSTANTS_Namelist>` namelist.

.. admonition:: Namelist Usage

   :Required/Optional: Optional
   :Single/Multiple Instances: Single


Namelist Variables
------------------

.. contents::
   :local:
   

face_set_ids
^^^^^^^^^^^^
A list of face set IDs that define the portion of the boundary where the
evaporation model will be applied.

:Default: none

vaporization_heat
^^^^^^^^^^^^^^^^^
The latent heat of vaporization :math:`L` (energy per unit mass).

:Default: none

vaporization_temp
^^^^^^^^^^^^^^^^^
The vaporization temperature :math:`T_0` (K) at the ambient pressure.

:Default: none

molar_mass
^^^^^^^^^^
The molar mass :math:`M` of the evaporating material (mass per mole).

:Default: none

ambient_pressure
^^^^^^^^^^^^^^^^
The pressure of the ambient atmosphere :math:`p_0` (force per unit area).

:Default: :math:`1.01325\times10^5~\text{N}/\text{m}^2` (1 atm)

condensation_factor
^^^^^^^^^^^^^^^^^^^
The empirical constant :math:`\lambda` that accounts for the condensation
of vaporized atoms.

:Default: :math:`0.1`
:Note: For evaporation in a vacuum :math:`\lambda=1`, but at 1 atm of
       pressure the Hertz-Knudsen equation overestimates the mass flux by an
       order of magnitude. :footcite:`he2003alloying`

----

.. footbibliography::
