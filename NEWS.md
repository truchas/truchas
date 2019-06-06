# NEWS

Important changes since the 3.0 release.

## 2019-06-05 (ddec9290)

The configuration of boundary conditions for heat transfer and species diffusion has undergone a major revision. The DS_BOUNDARY_CONDITION and DS_INTERFACE_CONDITION namelists have been replaced by two new namelists: THERMAL_BC for heat transfer boundary and interface conditions, and SPECIES_BC for species boundary conditions. The content of the namelists is essentially unchanged but some variables and their values are different:
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
