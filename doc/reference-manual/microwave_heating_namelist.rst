.. _MICROWAVE_HEATING_Namelist:

.. toctree::
   :maxdepth: 1

MICROWAVE_HEATING Namelist
===========================
This namelist specifies the parameters that manage the computation of the
dielectric heat source for microwave heating simulations. It is required
whenever the PHYSICS namelist option :ref:`microwave_heating<physics-mwh>`
is enabled. The computation uses the frequency-domain solver for time-harmonic
Maxwell's equations, using the configuration specified by the
:ref:`FDME_SOLVER<FDME_SOLVER_Namelist>` namelist,
to solve an auxiliary electromagnetics (EM) problem.
The equations are solved on the tetrahedral mesh specified by the
:ref:`EM_MESH<EM_MESH_Namelist>` namelist, which is generally different
than the mesh used for heat transfer. EM boundary conditions are defined using
:ref:`ELECTROMAGNETIC_BC<ELECTROMAGNETIC_BC_Namelist>` namelists.

.. caution::

   The electromagnetic models assume SI units by default. In particular, the
   computed dielectric heat is a power density, W/m\ :sup:`3` in SI units. If
   you want to use a different system of units, you must assign appropriate
   values to the PHYSICAL_CONSTANTS namelist variables
   :ref:`vacuum_permittivity<PhyCo_VPm>` (:math:`\epsilon_0`) and
   :ref:`vacuum_permeability<PhyCo_VPt>` (:math:`\mu_0`).


Namelist Variables
------------------
.. contents::
   :local:

frequency
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The frequency (cycles per unit time).

:Type: real
:Default: none


wg_port_bc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A list of names of :ref:`"wg-port" type<em-bc-type>` ELECTROMAGNETIC_BC
namelists whose power is to be managed by this namelist.

:Type: string list
:Default: none


times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
An optional list of one or more heat transfer times :math:`t_1 < \ldots < t_n`
when the waveguide port feed input powers changes.

:Type: real list
:Default: none


powers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The waveguide input powers. These values are effectively assigned to the
EM boundary condition variable :ref:`em-bc-power` at the appropriate times.

:Type: real array
:Default: None
:Notes: This is formally a 2-dimensional array of values. When `times`_ is
   not specified, the powers are constant in time and are specified as

      |      ``powers(:,1) =`` :math:`p_1, \ldots, p_m`

   or more simply as

      |      ``powers =`` :math:`p_1, \ldots, p_m`

   with one value for each of the :math:`m` BC specified by `wg_port_bc`_.
   When `times`_ is specified, additional columns of powers to use at each
   of those times must also be specified, as in

      | ``powers(:,1) =`` :math:`p_1, \ldots, p_m`
      | ``powers(:,2) =`` :math:`p_1, \ldots, p_m`
      |     ...
      | ``powers(:,n+1) =`` :math:`p_1, \ldots, p_m`


prop_change_threshold
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This control mechanism determines, at each heat transfer time step, whether
the EM heat source is recomputed to account for temperature-dependent changes
in the EM material properties. The EM heat source is recomputed whenever the
maximum relative difference between the current properties and those from the
previous computation exceeds this threshold value. Otherwise the previously
computed source is retained. When the properties are independent of temperature,
the heat source is recomputed solely when there is a change in the waveguide
input powers.

:Type: real
:Default: 0.3
:Valid Values: > 0
:Notes: The relevant properties are the (real) magnetic permeability and the
        real and imaginary parts of the electric permittivity.


data_mapper_kind (experimental)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This specifies the tool used to map fields between the primary heat transfer
mesh and the EM mesh. Select "portage" to use an experimental data mapper based
on the Portage toolkit, https://laristra.github.io/portage. Otherwise, the
standard mapping tool will be used by default.

:Type: string
:Default: `"default"`
:Valid Values: `"default"`, `"portage"`
