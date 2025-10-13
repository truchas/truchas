.. _INDUCTION_HEATING_Namelist:

.. toctree::
   :maxdepth: 1

INDUCTION_HEATING Namelist
===========================
This namelist specifies the parameters that manage the computation of
Joule heat source for induction heating simulations. It is required whenever
the PHYSICS namelist option :ref:`induction_heating<physics-ih>` is enabled.
The computation solves an auxiliary electromagnetics (EM) problem that
simulates the eddy currents induced by an external low-frequency magnetic
field. The EM field equations are solved on the tetrahedral mesh specified
by the :ref:`EM_MESH<EM_MESH_Namelist>` namelist, which is generally different
than the mesh used for heat transfer. EM boundary conditions are defined using
:ref:`ELECTROMAGNETIC_BC<ELECTROMAGNETIC_BC_Namelist>` namelists. A properly
configured model will typically include a :ref:`magnetic induction
source<em-bc-ih-field>` boundary condition with an external field defined by
the :ref:`INDUCTION_SOURCE_FIELD<INDUCTION_SOURCE_FIELD_Namelist>` namelist.

The time scale of the forced EM equations is assumed to be very short compared
to the fundamental time scale of heat transfer (see the INDUCTION_SOURCE_FIELD
namelist parameter
:ref:`frequency<INDUCTION_SOURCE_FIELD_Namelist/index:frequency>`.) To avoid
introducing this short time scale into heat transfer, the rapid temporal
variation in the Joule heat is averaged over a period or cycle of the forcing.
In effect, the EM field equations evolve on a distinct inner "fast" time that
unfolds in an instant of the outer "slow" time of heat transfer.

.. caution::

   The electromagnetic models assume SI units by default. In particular, the
   computed Joule heat is a power density, W/m\ :sup:`3` in SI units. If you
   want to use a different system of units, you must assign appropriate values
   to the PHYSICAL_CONSTANTS namelist variables
   :ref:`vacuum_permittivity<PhyCo_VPm>` (:math:`\epsilon_0`) and
   :ref:`vacuum_permeability<PhyCo_VPt>` (:math:`\mu_0`).


Namelist Variables
------------------
.. contents::
   :local:


use_fd_solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Enabling this flag activates the frequency-domain solver, using the
configuration specified by the :ref:`FDME_SOLVER<FDME_SOLVER_Namelist>`
namelist. Otherwise the standard time-domain solver is used instead, using
the configuration specified by the
:ref:`TDME_JOULE_SOLVER<TDME_JOULE_SOLVER_Namelist>` namelist.
The frequency-domain solver should be considered experimental for this
application.

:Type: logical
:Default: ``.false.``


prop_change_threshold
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This control determines, at each heat transfer time step, whether the Joule
heat source is recomputed to account for temperature-dependent changes in the
EM material properties. The source is recomputed whenever the maximum relative
difference between the current properties and those from the previous
computation exceeds this threshold value. Otherwise, the previously computed
source is retained. When the properties are independent of temperature,
the Joule heat source is only recomputed when the forcing function changes.

At each heat transfer time step, this control governs the recomputation of
the Joule heat source based
:Type: real
:Default: 0.3
:Valid Values: > 0
:Notes: The electrical conductivity and magnetic susceptibility are the only
properties whose changes are monitored. The electric susceptibility
only enters the equations through the displacement current term,
which is an exceedingly small perturbation in the quasi-magnetostatic
regime, and has essentially no effect on the Joule heat.

For electric conductivity, only the conducting region (where the
value is positive) is considered when computing the difference.
An underlying assumption is that this region remains fixed throughout
the simulation.


data_mapper_kind (experimental)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This specifies the tool used to map fields between the primary heat transfer
mesh and the EM mesh. Select "portage" to use an experimental data mapper based
on the Portage toolkit, https://laristra.github.io/portage. Otherwise, the
standard mapping tool will be used by default.

:Type: string
:Default: "default"
:Valid Values: "default", "portage"
