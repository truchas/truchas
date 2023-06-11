.. _PBF_MATERIAL_Namelist:

PBF_MATERIAL Namelist (Experimental)
====================================

In the powder bed fusion (PBF) additive manufacturing process, a layer of a
powdered form of a material is spread across the surface of the build volume
and selectively melted. When the melt solidifies it doesn't return to the
powdered form but rather to the solid material. Such history-dependent phase
changes cannot be described directly by the material model used in Truchas,
which assumes that material phase fractions are a function of temperature
alone. This namelist enables an indirect method for modeling such a
powder-liquid-solid history dependent material.

The method requires that two separate 2-phase materials are defined: a
*precursor material*(the powder) and an *end material*. The volume fraction
of solid precursor material represents the volume fraction of powder, the
volume fraction of solid end material represents the volume fraction of solid,
and the sum of the liquid volume fractions of the two materials represents the
liquid volume fraction. There are several requirements on the materials in
order for this to work properly:

* The two liquid phases must have identical properties.

* The two phase changes must define identical solid fraction functions; however
  the latent heats may differ.

* All phases must have the same density value. Truchas heat transfer and flow
  cannot directly model physical volume changes associated with phase change.

The user is responsible for ensuring these requirements are satisfied.

The history dependence algorithm follows two simple rules. On any given cell:

1. precursor material melts first during heating; and
2. the solid volume fraction of precursor material is constant during cooling.

.. note::

   :Required/Optional: Optional
   :Single/Multiple Instances: Single

Namelist Variables
------------------
The namelist consists of just two variables that specify the materials:

**material1**

      The name of the precursor material.

**material2**

      The name of the end material.
