.. _PHASE_CHANGE_Namelist:

.. toctree::
   :maxdepth: 1

PHASE_CHANGE Namelist 
==============================

Overview
-----------
The :ref:`PHASE_CHANGE<PHASE_CHANGE_Namelist>` namelist defines the phase change model used for the transformation between two phases of a multi-phase material. An instance of this namelist is required for each consecutive pair of phases specified by the :ref:`phases<MP_phase_pm>` variable of a :ref:`MATERIAL<MATERIAL_and_PHASE_Namelists>` namelist.

Truchas models the transformation from low temperature phase ("solid") to high temperature phase ("liquid") using a "solid fraction" function :math:`f_{\text{sol}}(T)` that gives the fraction of low temperature phase as a function of temperature. There are two temperatures :math:`T_{\text{sol}} \lt T_{\text{liq}}` such that :math:`f_{\text{sol}} = 1` for :math:`T \leq T_{\text{sol}},f_{\text{sol}} = 0` for :math:`T \geq T_{\text{liq}}`, and :math:`0 \lt f_{\text{sol}}  \lt 1` for :math:`T_{\text{sol}} \lt T \lt T_{\text{liq}}`. There are two alternative methods for specifying this function.

The first simple, but non-physical, method uses a smooth Hermite cubic polynomial to interpolate :math:`f_{\text{sol}}` between :math:`T_{\text{sol}}` and :math:`T_{\text{liq}}`. This is pictured below. The second method uses a table of :math:`(T,f_{\text{sol}})` values to define :math:`f_{\text{sol}}(T)`. This allows physics-based phase change models like lever and Scheil to be defined, and data from CALPHAD-type tools to be used directly.

Note that while described in terms of solid and liquid, these phase change models also apply to solid-solid transformations with the obvious reinterpretation of notation.

Namelist Variables
----------------------

**low_temp_phase**, **high_temp_phase**
      These are the names of the two material phases.

**solidus_temp**, **liquidus_temp**
      The solidus and liquidus temperatures, :math:`T_{\text{sol}}` and :math:`T_{\text{liq}}`. If these variables are defined, the simple solid fraction model is used, which interpolates the solid fraction over the interval :math:`[T_{\text{sol}}, T_{\text{liq}}]` using a smooth Hermite cubic polynomial. These variables are incompatible with **solid_frac_table**.

      .. _fig_phase_change:
      .. figure:: images/fs-smooth.jpg
         :width: 300px
         :align: center

         Phase Change Schematic

**solid_frac_table**
      A table of temperature-solid fraction values. The format is the same as described for the :ref:`FUNCTION<FUNCTION_Namelist>` namelist :ref:`Tabular_Data<FUNC_TD>` variable. However the table may be given in either increasing or decreasing temperature order, and the solid fraction values must be strictly monotone. The table endpoints must specify the 1 and 0 solid fractions, and the corresponding temperatures are taken to be :math:`T_{\text{sol}}` and :math:`T_{\text{liq}}`, respectively. The solid fraction function is interpolated from this table using Akima smoothing; see the :ref:`FUNCTION<FUNCTION_Namelist>` namelist :ref:`Tabular_Interp<FUNC_TI>`. Additional data points outside the interval :math:`[T_{\text{sol}}, T_{\text{liq}}]` are automatically added to ensure that the solid fraction is constant outside the transformation interval. A plotfile of the resulting smoothed solid fraction function is written to the output directory. This is a standard 2-column text file of :math:`(T,f_{\text{sol}}(T))` values over the transition interval, which can be imported into your favorite plotting tool. This variable is incompatible with **solidus_temp** and **liquidus_temp**.

      ::

         solid_frac_table = 930.0 1.00
                            931.0 0.95
                            ...
                            940.0 0.00

**latent_heat**
      The energy per unit mass :math:`L` absorbed during the transformation of the material from the low temperature phase to the high temperature phase at the liquidus temperature :math:`T_{\text{liq}}`. This property is required when the specific enthalpy of the high temperature phase is defined from its specific heat. Otherwise it is ignored.

      .. _fig_latent_heat:
      .. figure:: images/latent-heat.jpg
         :width: 300px
         :align: center

         Latent Heat Schematic

