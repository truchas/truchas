.. _THERMAL_SOURCE_Namelist:

.. toctree::
   :maxdepth: 1

THERMAL_SOURCE Namelist
============================

Overview
------------
The THERMAL_SOURCE namelist is used to define external volumetric heat sources (power per unit volume). This source is in addition to any other sources coming from other physics, such as a Joule heat source. Each instance of this namelist defines a source, and the final source is the sum of all such sources (subject to some limitations).

Two forms of sources :math:`q` can be defined:

1. :math:`q(t,x) = f(t,x,T)\chi_S(x)`, where :math:`f` is a user-defined function and :math:`S` is a subdomain corresponding to one or more user-specified mesh cell sets. Here :math:`\chi_S` is the characteristic function on :math:`S: \chi_S = 1` for :math:`x \in S` and :math:`\chi_S = 0` for :math:`x \notin S`. Sources of this form can be summed as long as the interiors of their subdomains do not intersect.

2. :math:`q(t,x) = A(t)\Sigma_jq_j\chi_j(x)`, where :math:`A` is a user-defined time-dependent prefactor, :math:`q_j` is a constant source of mesh cell :math:`j`, and :math:`\chi_j` is the characteristic function on cell :math:`j`. The collection of values {:math:`q_j`} is read from a data file.

THERMAL_SOURCE Namelist Features
----------------------------------
| **Required/Optional        :** Optional
| **Single/Multiple Instances:** Multiple

Components
------------
* :ref:`name<TS_N>`
* :ref:`cell_set_ids<TS_CSI>`
* :ref:`data_file<TS_DF>`
* :ref:`prefactor<TS_P>`
* :ref:`prefactor_func<TS_PF>`
* :ref:`source<TS_S>`
* :ref:`source_func<TS_SF>`

.. _TS_N:

Equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A unique name used to identify a particular instance of this namelist.
| **Type**        : string (31 characters max)
| **Default**     : none

.. _TS_CSI:

cell_set_ids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : A list of cell set IDs that define the subdomain where the source is applied.
| **Type**        : a list of up to 32 integers
| **Default**     : none
| **Valid values**: any valid mesh cell set ID
| **Note**        : Different instances of this namelist must apply to disjoint subdomains; overlapping of source functions of this form is not supported. Exodus II mesh element blocks are interpreted by Truchas as cell sets having the same IDs.

.. _TS_S:

source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The constant value of the heat source. To specify a function, use :ref:`source_func<TS_SF>` instead.
| **Physical dimension**: :math:`E/T L^3`
| **Type**        : real
| **Default**     : none

.. _TS_SF:

source_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist that defines the source function. That function is expected to be a function of (t,x,y,z,T).
| **Type**        : string
| **Default**     : none

.. _TS_DF:

data_file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The path to the data file. It is expected to be a raw binary file consisting of a sequence of 8-byte floating point values, the number of which equals the number of mesh cells. The order of the values is assumed to correspond to the external ordering of the mesh cells. The file can be created with most any programming language. In Fortran use an unformatted stream access file.
| **Type**        : string
| **Default**     : none

.. _TS_P:

prefactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The constant value of the prefactor A. For a function use :ref:`prefactor_func<TS_PF>`.
| **Type**        : real
| **Default**     : none

.. _TS_PF:

prefactor_func
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Description** : The name of a :ref:`FUNCTION<FUNCTION_Namelist>` namelist that defines the function that computes the value of the time dependent prefactor :math:`A(t)`.
| **Type**        : string
| **Default**     : none
