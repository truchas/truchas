# NEWS

## 2021-02-04 Version 21.02

## 2021-01-07 (77d200c2)

The `write-restart.py` utility for creating restart files has been updated
give the option of using Portage to compute the mapping of fields in "mapped
restart" use cases. Use the `-h` option for usage information. This also
applies to the underlying python tools that can be used to write custom
scripts.

## 2021-01-06 (fe855127)

Adds an option `metis_ptype` to METIS mesh partitioning to choose between
using multilevel recursive bisection, the new default, and multilevel k-way
partitioning, which was used originally. See the Reference Manual for details.

## 2021-01-05 Version 21.01

## 2021-01-05 (540d6109)

Significant performance improvements were made to the initialization of the
main mesh. Users can expect to see dramatically shorter mesh initialization
times for meshes larger than several million cells, with greater relative
improvements for larger meshes.

## 2020-12-24 (a0c0dd62)

Added an option to use METIS graph partitioning to compute the parallel
decomposition of the mesh. This appears to produce much higher quality
decompositions than the default Chaco partitioner, especially for large
meshes. Users are encouraged to try it. To use, set `partitioner = 'metis'`
in the `MESH` and `ALTMESH` namelists. See the Reference Manual for more
details and additional options for modifying the METIS graph partitioning
algorithm.

## 2020-12-23 (856efdcb)

Truchas can now be built for the first time using the GFortran compiler.
Versions 9.1, 9.2, 9.3, 10.1, and 10.2 have been tested. All but a couple
tests are passing.

## 2020-12-01 Version 20.12

## 2020-12-01 (1b1aff7b)

Added a new microstructure analysis module "gl", which collects the thermal
gradient vector G and the rate of cooling L at the onset of solidification,
and the total solidification time. This module uses temperature values for
specifying the thresholds rather than solid fraction values.

## 2020-11-23 Version 20.11.1

## 2020-11-23 (8dd75d28)

Added a new variable `event_lookahead` to the `SIMULATION_CONTROL` namelist,
which is the number of steps over which the time step size is gradually
adjusted in order to hit an event time -- typically an output time. This had
been hardwired to 5, which is now the default. See the Reference Manual.

## 2020-11-23 (2550bd1b)

Output of the temperature gradient and time derivative fields was restored.

## 2020-11-22 (42adcb44)

A minimum frequency with which the preconditioner is updated can now be
specified in the input file. See the Reference Manual description of the
`DIFFUSION_SOLVER` namelist variable `pc_freq` for details.

## 2020-11-09 Version 20.11

## 2020-11-08 (914f932c)

The recent configuration change to the Hypre BoomerAMG preconditioner used in
heat transport has be reverted to its original configuration. Further testing
on large problems showed little, if any, performance benefit and a loss in
robustness requiring more careful attention to the number of cycles.

## 2020-10-29 (78ce9948)

An option to apply rotations to the computational mesh has been added. See the
Reference Manual description of the `rotation_angles` variable in the `MESH`
and `ALTMESH` namelists.

## 2020-10-22 (8772728b)

The recently introduced Portage data mapper has been made an optional
compile-time component due to the difficulty in compiling the Portage library
on some platforms.

## 2020-10-21 (cbc562ea)

A Dirichlet temperature boundary condition can now be imposed on a surface
that participates in enclosure radiation. This greatly simplifies the modeling
of heating elements with a prescribed surface temperature.

## 2020-10-09 Version 20.10

## 2020-10-09 (122a475)

Induction heating simulations employ a tool for mapping data fields between
the main heat transfer mesh and the electromagnetic solver mesh. An alternative
data mapper based on the Portage toolkit, https://laristra.github.io/portage,
has been added. Unlike the existing data mapper, this new experimental data
mapper is capable of handling main meshes containing prism and pyramid cells.
It is *much* slower however, and not used by default. See the Reference Manual
chapter on the `ALTMESH` namelist on how to enable it.

## 2020-09-29 (48c2de2)

The configuration of the Hypre BoomerAMG preconditioner used in heat transport
has been updated to follow recommendations from the Hypre team. As a result the
default number of AMG cycles (`pc_amg_cycles`) has been increased from 2 to 4,
and users should expect to use significantly more cycles than previously.
Nevertheless users should hope to see a noticeable overall improvement in
performance.

## 2020-09-15 (1cf4670)

The method of computing solidification velocity for microstructure GV-analysis
has been greatly simplified. The new method is much less costly, but more
significantly it resolves a long-standing, difficult-to-locate, memory leak
that occurred when using the Intel compiler. Users will observe some differences
in the new solidification velocities.

## 2020-09-13 (a9d9c2a)

Fixed an error that occurred when using a polynomial-type specific heat
function for the high temperature phase of a 2-phase material. See
https://gitlab.com/truchas/truchas/-/issues/161

## 2020-09-11 (b7db4ef)

A new flow algorithm option has been added that helps suppress the formation
of tiny isolated fragments of fluid (wisps) during splashy filling simulations.
These wisps often cause significantly reduced time step sizes, and even
complete stalling on occasion. This algorithm option is enabled using the new
`FLOW` namelist flag `wisp_redistribution`. By default the option is disabled.
See the Reference Manual for a description of the associated `wisp_*` algorithm
parameters. A detailed description of the algorithm will be found in the
Physics and Algorithms manual.

## 2020-09-11 (1429272)

The `SOLID MECHANICS` namelist preconditioner variables `preconditioning_steps`
and `relaxation_parameter` have been removed. Internally they default to 1 and
1.0, respectively, which correspond to the diagonal preconditioning used in
practice.

## 2020-09-02 Version 20.09

## 2020-08-14 Version 20.08

## 2020-08-14 (3b61747)

The original legacy flow solver has been removed.

## 2020-08-14 (38145ea)

User-defined heat sources are now defined using the new `THERMAL_SOURCE`
namelist instead of the `DS_SOURCE` namelist. There is no change in the
variable names, except that the `equation` variable is no longer needed;
"temperature" is naturally assumed.

With the new namelist comes a new option to read the source from a data
file. See the Reference Manual description of `THERMAL_SOURCE` for details.

## 2020-08-13 (b6fa8e7)

The linear and nonlinear solver parameters for solid mechanics have been
moved into the `SOLID_MECHANICS` namelist from `LINEAR_SOLVER` and
`NONLINEAR_SOLVER`. The current solver uses an accelerated nonlinear Krylov
(NKA) method with simple diagonal preconditioning. For a linear problem this
effectively reduces to GMRES. Originally other options were available, but
this was the only one used in practice. See the Reference Manual description
of `SOLID_MECHANICS` for details.

## 2020-07-27 (5d41b69) Version 20.07

We are moving to a rolling release of Truchas and monthly date-based tagging
of the master branch using a YY.MM format.

## 2020-07-14 (8941aef)

Interface heat transfer boundary conditions may now be temperature dependent.
This can be used to model the change in heat conduction across a metal/mold
interface due to the shrinkage of the metal as it solidifies, for example.
A side-efect of this change is that time and spatially dependent Interface
heat transfer coefficients must be defined slightly differently. See the
Reference Manual description of interface heat transfer boundary conditions
in `THERMAL_BC` for details.

## 2020-02-08 (21c6bb7)

The input format for materials and properties was completely overhauled. The
new format is more sensible and far simpler, and adds new options for modeling
phase change, but it will require significant changes to existing input files.
Here is a brief summary of the input changes; see the Truchas Reference Manual
for details.

* The `MATERIAL_SYSTEM` namelist was renamed `MATERIAL` with the original
  `MATERIAL` namelist eliminated and its few remaining variables moved to
  other namelists.
* Materials and phases are no longer assigned numeric labels. Input variables
  that used material numbers now use the assigned material or phase name.
* The `immobile` flag is replaced by the `is_fluid` flag whose default value
  is false. Previously all materials were fluid by default.
* The requirement to appoint one material as the background material was
  eliminated. An optional background type `BODY` namelist serves the purpose.
* The `priority`, `sound_speed`, and `permeability_constant` variables are
  legacy flow algorithm parameters and have been moved to the `LEGACY_FLOW`
  namelist.
* A void material no longer needs to be defined (nor can be). Use the
  reserved name `"VOID"` to refer to the void material.
* Properties for single-phase materials are now defined directly in the
  `MATERIAL` namelist; an additional `PHASE` namelist is no longer required.
* Multi-phase material properties that apply to all phases can be defined
  directly in the `MATERIAL` namelist as well. The `PHASE` namelist continues
  to be used to define properties on a per-phase basis where needed.
* The generic `property_name`, `property_constant`, and `property_function`
  input arrays have been replaced by property-specific variables. For example,
  `viscosity` or `viscosity_func` for the fluid viscosity; the former taking
  a constant value, and the latter the name of a `FUNCTION` namelist similar
  to before. These are used in both the `MATERIAL` and `PHASE` namelists.
* The latent heat and low/high temperature phase change parameters were moved
  from the `MATERIAL` namelist to a new `PHASE_CHANGE` namelist which specifies
  the parameters associated with a single phase change that is identified by
  the name of its low and high temperature phases.

A significant change to the default phase change model was also made.
Previously the solid fraction varied linearly between the solidus and
liquidus temperatures with some optional smoothing at the corners. This
non-physical model was replaced by an equally non-physical, but simpler,
mathematical model that uses a smooth Hermite cubic polynomial to
interpolate the solid fraction between solidus and liquidus.

Complementing this mathematical phase change model is a new option to
explicitly define the temperature-dependent solid fraction as a data table.
This allows physics-based phase change models like lever and Scheil to be
defined, and data from CALPHAD tools to be used directly. Related to this
is the new capability to define a material's specific enthalpy (perhaps as
a table) as an alternative to defining its specific heat and latent heats
of its phase changes. See the `PHASE_CHANGE` namelist section in the latest
Truchas Reference manual for details.

## 2020-01-09 (056da49)

Improvements were made to the internal algorithm that initializes the material
volume fractions as specified by the input BODY namelists. Users may notice
small differences in the initial volume fractions for geometrically-defined
bodies. The new results are much more accurate. But otherwise there are just
a few minor changes that are unlikely to impact most users:
* Several rarely-used (if ever) geometric BODY types have been removed.
* While the order of BODY namelists has always been significant, the new
  behavior may be slightly different. BODY namelists are processed in the
  order they appear, and identify the specified part of the computational
  domain not claimed by any preceding BODY namelist. The "background" type
  BODY, if any, must come last.
* There is no longer a choice in volume fraction initialization method,
  and as a consequence the INTERFACES namelist has been removed.

See the BODY namelist section in the latest Truchas Reference Manual for
further details.

## 2019-09-24 (fd75405)

The current VOF flow algorithm effectively treats small fragments of void
entrained in fluid as incompressible, resulting in unphysical void "bubbles"
that persist in the flow. This is a particular problem for splashy filling
simulations. A new void collapse model has been implemented that drives
the collapse of these void fragments. While the model is still experimental
and subject to change, it appears to work remarkably well enough at this
point to recommend trying it. To activate the model set the new FLOW namelist
variable `void_collapse` to true (the default is false). Please send any
feedback about the model to truchas@lanl.gov, and see the latest Truchas
Reference Manual for further details.

## 2019-08-29 (3927011)

It is now possible to define a rectilinear hex mesh of a brick domain in
the MESH namelist as an alternative to reading an Exodus II mesh file.
Though not especially useful for most real applications, it does make
Truchas much more accessible for simple tests and demo problems by
avoiding the "meshing tool" obstacle. See the MESH namelist section
in the latest Truchas Reference Manual for details.

## 2019-07-25 (88211a5f)

Solution probe output was re-implemented. This should resolve several
long-standing bugs and memory leaks that had been reported. This involves
a number of user-visible changes:
* Each probe writes directly to its own text file in the output directory,
  and no longer writes into the Truchas HDF5 output file.
* Rather than write data for every possible solution quantity, a probe now
  outputs data for a specific quantity.
* There are several changes to the PROBE namelist:
  - The `probe_name` variable was removed.
  - `probe_coords` was changed to `coord`.
  - `probe_coords_scale` was changed to `coord_scale_factor`.
  - `probe_description` was changed to `description`.
  - The new variable `data` specifies the quantity to produce output for. The
    current choices are: `"temperature"`, `"pressure"`, and `"velocity"`.
    Additional choices will be added in the near future.
  - The new variable `data_file` specifies the name of the text output file.
* The `write-probes.py` utility will be retained for use with previous
  Truchas .h5 output files that contained the probe data.

See the latest Truchas Reference Manual section on the PROBE namelist for
further details.

## 2019-06-05 (ddec9290)

The configuration of boundary conditions for heat transfer and species
diffusion has undergone a major revision. The DS_BOUNDARY_CONDITION and
DS_INTERFACE_CONDITION namelists have been replaced by two new namelists:
THERMAL_BC for heat transfer boundary and interface conditions, and SPECIES_BC
for species boundary conditions. The content of the namelists is essentially
unchanged but some variables and their values are different:
* The `variable` variable is no longer needed and was removed.
* The `condition` variable was renamed to `type`. There is no change to the
  types of boundary and interface conditions available, but some keywords
  have changed.
* The THERMAL_BC namelist handles both boundary conditions and interface
  conditions. The mapping of type keywords is
  - DS_BOUNDARY_CONDITION: `"dirichlet"` --> `"temperature"`
  - DS_INTERFACE_CONDITION: `"htc"` --> `"interface-htc"`, and
    `"radiation"` --> `"gap-radiation"`
* The SPECIES_BC namelist includes a new variable `comp` to specify the
  species component, and the type keyword `"dirichlet"` --> `"concentration"`.
* For improved clarity, the generic `data_constant` and `data_function`
  arrays have been replaced by variables specific to each type of boundary
  condition; e.g., `temp` and `temp_func` for the temperature Dirichlet BC.

See the latest Truchas Reference Manual sections on the THERMAL_BC and
SPECIES_BC namelists for further details.

## 2019-05-09 (8814f471)

A very basic "mapped restart" capability has been restored to the new
`write-restart.py` utility. This generates a restart file for a new mesh,
interpolating field data from the output file on the original mesh. Use the
`--help` option to get the usage information.

## 2019-03-29 (d22d4a88)

The framework for post-processing the .h5 output file and for regression
testing has been entirely rewritten. Several consequences follow from
this change:

* Python 3.5 or later is required, together with the h5py and scipy Python
  packages. Previously 2.7 and numpy were required. The Python development
  files are also no longer required.

* The names of the post-processing utilities have changed, but they retain
  their original interface:

  - `xdmf-parser.py` --> `write-xdmf.py`
  - `truchas-gmv-parser.py` --> `write-gmv.py`
  - `write_restart` --> `write-restart.py`
  - `write_probe` --> `write-probe.py`

## 2019-03-26 (2f097cf5)

The inputs for flow and solid mechanics have been reorganized. **This change
will require the updating of all input files,** including those not using
either flow or solid mechanics.

* The PHYSICS namelist variable `fluid_flow` was renamed `legacy_flow` and
  its default value changed to `.false.`
* Input files that previously set `fluid_flow = .false.` can simply delete
  that line. That may be the only change required for such files.
* The PHYSICS namelist variable `body_force` was renamed `body_force_density`
  to more accurately reflect its meaning.
* All other flow-related PHYSICS namelist variables were moved into the new
  LEGACY_FLOW namelist.
* All flow-specific NUMERICS namelist variables were also moved into the
  LEGACY_FLOW namelist.
* The TURBULENCE namelist variables were renamed by dropping the redundant
  `turbulence_` prefix.
* The PHYSICS namelist variable `solid_mechanics_body_force` was moved to
  the new SOLID_MECHANICS namelist
* All solid mechanics-specific NUMERICS namelist variables were also moved
  into the SOLID_MECHANICS namelist.

See the documentation for the PHYSICS, NUMERICS, LEGACY_FLOW, SOLID_MECHANICS,
and TURBULENCE namelists in the current Reference Manual for details.
