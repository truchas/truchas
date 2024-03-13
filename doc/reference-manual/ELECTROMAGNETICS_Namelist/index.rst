.. _ELECTROMAGNETICS_Namelist:

.. toctree::
   :maxdepth: 1

ELECTROMAGNETICS Namelist
===========================

The ELECTROMAGNETICS namelist specifies the parameters for the induction
heating model and algorithms. This namelist is read whenever the
:ref:`PHYSICS<PHYSICS_Namelist>` option :ref:`electromagnetics<PHYSICS_EM>`
is enabled. The induction heating model solves the electromagnetic (EM) field
equations to compute the Joule heat source used by the heat transfer model
for induction heating simulations.
The EM field equations are solved on the tetrahedral mesh specified by the
:ref:`EM_MESH<EM_MESH_Namelist>` namelist, which is generally different
than the mesh used for heat transfer.
A properly configured model will include
a periodic boundary condition with a forcing function defined by the
:ref:`INDUCTION_SOURCE_FIELD<INDUCTION_SOURCE_FIELD_Namelist>` namelist.
EM boundary conditions are defined using
:ref:`ELECTROMAGNETIC_BC<ELECTROMAGNETIC_BC_Namelist>` namelists (but see
:ref:`legacy boundary conditions<legacy_em_bc>` below.)

The time scale of the forced EM equations is assumed to be very short compared
to the fundamental time scale of heat transfer (see the INDUCTION_SOURCE_FIELD
namelist parameter
:ref:`frequency<INDUCTION_SOURCE_FIELD_Namelist/index:frequency>`.) To avoid
introducing this short time scale into heat transfer, the rapid temporal
variation in the Joule heat is averaged over a period or cycle of the forcing.
In effect, the EM field equations evolve on a distinct inner "fast" time that
unfolds in an instant of the outer "slow" time of heat transfer.

.. caution::

   The induction heating model assumes SI units by default. In particular,
   the result of the Joule heat computation is a power density,
   W/m\ :sup:`3` in SI units. To use a different system of units,
   appropriate values must be assigned to the free-space constants
   :ref:`vacuum_permittivity<PhyCo_VPm>` and
   :ref:`vacuum_permeability<PhyCo_VPt>` in the
   :ref:`PHYSICAL_CONSTANTS<PHYSICAL_CONSTANTS_Namelist>` namelist.


.. admonition:: Namelist Usage

   :Required/Optional: Required when :ref:`electromagnetics<PHYSICS_EM>` is enabled.
   :Single/Multiple Instances: Single

Namelist Variables
------------------

.. contents::
   :local:

General Varibles
+++++++++++++++++++++++++++++++++

matl_change_threshold
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This controls, at each heat transfer time step, whether the Joule heat is
recomputed in response to temperature-induced changes in the EM material
properties. The Joule heat is recomputed whenever the maximum difference
between the current properties and those when the Joule heat was
last computed relative to the maximum of the current properties exceeds this
threshold value. Otherwise the previously computed Joule heat is used.
When the properties are independent of temperature, the Joule heat is only
recomputed when the forcing function changes.

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
This specifies the tool that will be used to map fields between the primary
mesh used for heat transfer and the EM mesh. If "portage" is selected,
an experimental data mapper based on the Portage toolkit,
https://laristra.github.io/portage, will be used. Otherwise, the normal
data mapping tool will be used by default.

:Type: string
:Default: `"default"`
:Valid Values: `"default"`, `"portage"`


Time-Domain Joule Heat Solver
+++++++++++++++++++++++++++++++++
The following parameters pertain to the time-domain Joule heat solver,
which integrates the periodically-forced EM field equations in time until
a periodic steady state solution is attained, and then the time-averaged
Joule heat is computed over one cycle of the forcing.


steps_per_cycle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The number of time steps per cycle of the periodic forcing used to integrate
the EM field equations.

:Type: integer
:Default: 20
:Valid Values: > 0
:Notes: Increasing the number of time steps per cycle increases the accuracy
   and stability of the Joule heat computation, while generally increasing
   the execution time. A reasonable range of values is [10,40]; anything
   less than 10 is severely discouraged.


steady_state_tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Convergence to the periodic steady state is measured by comparing the
computed Joule heat field averaged over the last source field cycle,
:math:`q_\text{last}`, with the result from the previous cycle,
:math:`q_\text{prev}`. The Joule heat computation is considered converged
when :math:`\lVert{q_\text{last}-q_\text{prev}}\rVert_\infty < \delta\,
\lVert{q_\text{last}}\rVert_\infty` where :math:`\delta` is the value of
this parameter.

:Type: real
:Default: :math:`10^{-2}`
:Valid Values: > 0
:Notes: Depending on the accuracy of the other physics, :math:`10^{-2}` or
        :math:`10^{-3}` are adequate values. The measured error is more
        properly the error in :math:`q_\text{prev}`; the actual error in
        :math:`q_\text{last}` will be significantly less.


max_source_cycles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Specifies the time limit, measured in cycles of the periodic forcing, allowed
for convergence to the periodic steady state. If convergence is not attained
within this limit, the last result is returned and a warning issued, but
execution continues.

:Type: integer
:Default: 5
:Valid Values: > 0
:Notes: To avoid ringing, the amplitude of the external source field is phased
        in and is not at full strength until after approximately two cycles
        have passed. Consequently this input variable should normally be at
        least 3. Convergence to a periodic steady state is commonly attained
        within 5 cycles; see `steady_state_tol`_.


cg_tol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The convergence tolerance :math:`\eta` for the conjugate gradient (GC) solution
of the linear time step system. The CG iteration is considered converged when
the residual :math:`r` satisfies :math:`\lVert{r}\rVert<\eta\lVert{r_0}\rVert`.

:Type: real
:Default: :math:`10^{-8}`
:Valid Values: :math:`(0, 0.1)`
:Notes: The numerical characteristics of the EM system require that the linear
        systems be solved to significantly greater accuracy than might
        otherwise be required. Too loose a tolerance will manifest itself in
        a significant build-up of noise in the solution of the electric field
        over the course of the simulation.


cg_max_iter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The maximum number of conjugate gradient (CG) iterations allowed. It is a
fatal error if convergence is not attained within this number of iterations.

:Type: integer
:Default: 500
:Valid Values: > 0


output_level
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Controls the verbosity of the time-domain Joule heat solver.

:Type: integer
:Default: 1
:Valid Values: 1, 2, 3 or 4
:Notes: At the default level, 1, a status message is output at the end of
  each source field cycle showing the progress toward steady state. Level 2
  adds a summary of the CG iteration for each time step. Level 3 adds the
  norm of the difference between the solution and extrapolated predictor
  for each time step. This gives an indication of the (time) truncation error,
  and if noise is accumulating in the system it will be seen here; see
  `cg_tol`_. Level 4 adds convergence info for each CG iterate. Levels 1
  and 2 are typical.


graphics_output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A flag that enables the output of the time-domain EM solution for graphics
visualization.

:Type: logical
:Default: `.false.`
:Notes: When enabled, the time-domain Joule heat solver will generate its
   own graphics data files in OpenDX format (see http://www.opendx.org).
   The files contain the material parameter fields, the averaged Joule heat
   field, and the time series of the EM fields, all on the EM mesh. The files
   are identified by the suffixes `-EM.dx` and `-EM.bin`. The value of this
   variable has no impact on the normal graphics output generated by Truchas,
   which is determined elsewhere, and the averaged Joule heat field mapped
   onto the heat transfer mesh will be output there in either case.

.. _legacy_em_bc:

Legacy Boundary Conditions (deprecated)
+++++++++++++++++++++++++++++++++++++++

When the induction heating model was originally introduced, it did not require
a specification of the EM boundary conditions, but inferred them for a small
number of cylinder-like EM mesh domain shapes. This feature is deprecated, but
still supported for the present. Users should migrate to defining the boundary
conditions explicitly using the :ref:`ELECTROMAGNETIC_BC
<ELECTROMAGNETIC_BC_Namelist>` namelist.


use_legacy_bc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
When this flag is enabled, the boundary conditions will be automatically
assigned based on hints given by the following parameters. Note that any
boundary condition specified by a :ref:`ELECTROMAGNETIC_BC
<ELECTROMAGNETIC_BC_Namelist>` namelist will be ignored in this case.

:Type: logical
:Default: `.false.`


symmetry_axis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Specifies which axis is the symmetry axis of the EM domain when
`use_legacy_bc`_ is enabled.

:Type: string
:Default: `"z"`
:Valid Values: `"x"`, `"y"`, or `"z"`
:Notes: The value assigned here should match the value assigned to the
   :ref:`orientation<INDUCTION_SOURCE_FIELD_Namelist/index:orientation>`
   parameter of the INDUCTION_SOURCE_FIELD namelist.


em_domain_type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Specifies the type of domain that is discretized by the EM mesh when
`use_legacy_bc`_ is enabled. The allowed values and corresponding domain when
`symmetry_axis`_ is `"z"` are given in the following table. The domains for
the other axis options are defined similarly and obtained by the appropriate
cyclic permutation of the coordinates.

.. csv-table::
   :header: "Value", "Domain"
   :class: tight-table
   :widths: auto

   `"full_cylinder"`,":math:`\Omega = \{(x,y,z) \,|\, x^2+y^2 \leq r^2, z_1\leq z\leq z_2 \}`"
   `"half_cylinder"`, ":math:`\Omega = \{(x,y,z) \,|\, x^2+y^2 \leq r^2, x\geq 0, z_1\leq z \leq z_2 \}`"
   `"quarter_cylinder"`, ":math:`\Omega = \{(x,y,z) \,|\, x^2+y^2 \leq r^2, x, y\geq 0, z_1\leq z \leq z_2 \}`"

:Type: string
:Default: none
:Notes:
   The values :math:`r>0` and :math:`z_1\gt z_2` are inferred from the mesh
   and are not specified directly.

   When `symmetry_axis`_ is `"z"`, the boundary condition
   :math:`\hat{n}\times\vec{H} = \hat{n}\times\vec{H}_\text{ext}` is imposed
   on the boundaries :math:`\{x^2+y^2=r^2\}` and :math:`\{z=z_1,z_2\}`, where
   :math:`\vec{H}_\text{ext}` is defined by the
   :ref:`INDUCTION_SOURCE_FIELD<INDUCTION_SOURCE_FIELD_Namelist>` namelist,
   and the boundary condition :math:`\hat{n}\times\vec{E}=0` is imposed on
   the symmetry planes :math:`\{x=0\}` and :math:`\{y=0\}` if present. And
   analogously for the other axis options.
