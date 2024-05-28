.. _SPECIES_SOURCE_Namelist:

.. toctree::
   :maxdepth: 1

SPECIES_SOURCE Namelist
=======================

The SPECIES_SOURCE namelist is used to define external volumetric sources
for the species advection-diffusion model.

.. admonition:: Namelist Usage

   :Required/Optional: Optional
   :Single/Multiple Instances: Multiple

Namelist Variables
------------------

.. contents::
   :local:


comp_id
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The index of the species component to which this source applies.

:Type: integer
:Default: 1
:Note: The number of species components is defined by the PHYSICS namelist
       variable :ref:`number_of_species<PHYSICS_NOS>`.


cell_set_ids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A list of cell set IDs that define the subdomain where the source is applied.

:Type: a list of up to 32 integers
:Default: none
:Valid Values: any valid mesh cell set ID
:Note: Different instances of this namelist for a given component must apply
       to disjoint subdomains; overlapping of source functions is not allowed.

       ExodusII mesh element blocks are interpreted by Truchas as cell sets
       having the same IDs.


source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The constant value of the source. To specify a function use `source_func`_
instead.

:Type: real
:Default: none


source_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist that defines the
source function. The function is expected to be a function of :math:`(t,x,y,z)`.

:Type: string
:Default: none
