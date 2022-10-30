.. _MICROSTRUCTURE_Namelist:

MICROSTRUCTURE Namelist
=======================

Truchas is able to capture and analyze data about the solidification process
as it evolves within each mesh cell, and output the results for later analysis.
This data is generally used by analytic 0-dimensional models to predict the
local microstructural characteristics and defects of the solidified material.
Use the MICROSTRUCTURE namelist to enable this feature and set the parameters
that control it.

.. note::

   :Required/Optional: Optional
   :Single/Multiple Instances: Single

Common Namelist Variables
-------------------------
The following variables are common to all analysis modules and are required.

**material**

      The name of the material for which solidification data will be collected.
      The material must have a single high-temperature liquid phase and one
      or more lower-temperature solid phases. The analysis is applied to its
      liquid-solid phase transformation.

**cell_set_ids**

      A list of cell set IDs that specify the part of the domain where
      solidification data will be collected. Cell-based data is output for the
      entire mesh, using dummy values for cells that are not included in the
      analysis and included cells without valid data.

Basic GL Analysis Module
------------------------
By default, the thermal gradient :math:`G=\nabla T` and cooling rate
rate :math:`L=-\partial T/\partial t` are recorded at the onset of
of solidification, and the local solidification time (time spent in the mushy
zone) is recorded at the completion of solidification. These are written to
the output file as the cell-based fields "uStruc-G", "uStruc-L", and
"uStruc-t-sol", respectively. A dummy value of 0 is written for cells not
included in the analysis and cells without valid GL or solidification time
data.

.. tip::
   Many, if not most, analytic models for microstructure and defect prediction
   can be defined in terms of the fundamental G, L, and solidification time
   data produced by this basic analysis module. These models can be easily
   computed and visualized within Paraview using its calculator tool. For
   example, the (vector) solidification front velocity is
   :math:`V = (L/\Vert{G}\Vert^2)G`, and the solidification
   front speed is :math:`\Vert{V}\Vert = L / \Vert{G}\Vert`.

There are two methods of specifying the parameters that control how the
data is collected, which are pictured below: one based on solid fraction
and the other based on temperature. They are mutually exclusive. The
collection procedure is designed to only report data for those cells that
have ultimately passed monotically from liquid to solid (or are in the
process of doing so). Thus some cells included in the analysis may report
invalid data even after all solidification is complete.

   .. figure:: images/fig.png
      :width: 600px
      :align: left

Namelist Variables for the Temperature Based Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**begin_temp**

   Solidification is regarded as having started when the temperature first
   falls below this threshold; it need not be the actual liquidus temperature.
   This is a required variable.

**end_temp**

   Solidification is regarded as complete when the temperature first drops
   below this threshold; it need not be the actual solidus temperature.
   This is a required variable.

**gl_temp**

   When the temperature first falls below this threshold, G and L are set
   to the current value of the thermal gradient and cooling rate. The value
   must belong to the interval [begin_temp, end_temp). This is an optional
   variable; its value defaults to the value of begin_temp.

**begin_temp_reset**

   If the temperature rises back above this threshold while solidifying, it
   is regarded as having returned to a liquid state and any G and L data that
   may have been computed is erased. This is an optional variable; its value
   defaults to the value of begin_temp, and if specified it must be greater
   than or equal to that value.

**end_temp_reset**

   If the temperature rises back above this threshold after becoming solid,
   it is regarded as having (partially) re-melted, and the previously computed
   G, L, and solidification time are erased. This is an optional variable;
   its value defaults to the value of end_temp, and if specified it must be
   greater than or equal to that value.

Namelist Variables for the Solid Fraction Based Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**begin_frac**

   Solidification is regarded as having started when the solid fraction first
   rises above this threshold. This is a required variable; its value must
   belong to the interval (0,1).

**end_frac**

   Solidification is regarded as complete when the solid fraction first rises
   above this threshold. This is a required variable; its value must belong
   to the interval (begin_frac, 1).

**gl_frac**

   When the solid fraction first rises above this threshold, G and L are
   set to the current value of the thermal gradient and cooling rate. The
   value must belong to the interval [begin_frac, end_frac). This is an
   optional variable; its value defaults to the value of begin_frac.

**begin_frac_reset**

   If the solid fraction falls back below this threshold while solidifying
   it is regarded as having returned to a liquid state and any G and L data
   that may have been computed is erased. This is an optional variable; its
   value defaults to the value of begin_frac, and if specified it must belong
   to the interval (0, begin_frac].

**end_frac_reset**

   If the solid fraction falls back below this threshold after becoming solid,
   it is regarded as having (partially) re-melted, and the previously computed
   G, L, and solidification time are erased. This is an optional variable; its
   value defaults to the value of end_frac, and if specified it must belong to
   the interval [begin_frac, end_frac].

Custom Microstructure Analysis Modules
--------------------------------------
The Truchas microstructure analysis framework is designed to readily accept
additional user-written analysis modules. Contact the Truchas developers for
details and refer to the template file located in the `src/truchas/ustruc`
directory.
