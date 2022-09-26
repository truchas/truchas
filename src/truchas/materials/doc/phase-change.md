# Phase Change
The `multiphase_matl` implementation of the `material` class describes a
multiphase material whose phases are linearly ordered from low to high
temperature phases, with separated temperature-dependent transformations
between consecutive pairs of phases. The abstract derived type `phase_change`
defines the interface that `multiphase_matl` type objects use to compute the
phase fractions as a function of temperature through one such transformation.

A transformation occurs over a given temperature interval. The lower limit
of the interval is referred to as the *solidus temperature*, and the upper
limit as the *liquidus temperature*. The temperature-dependent fraction of
the lower temperature phase is referred to as the *solid fraction*. While
this terminology reflects the usual solid-liquid transformation, it applies
equally well to solid-solid transformations with the obvious reinterpretations.

## The `phase_change` Abstract Class
> **NB:** This is an auxiliary class to the `multiphase_matl` type and is
not visible to application code.

### Initialization
Instances of the `phase_change` class are instantiated using the `alloc_phase_change`
procedure from the `phase_change_factory` module described below.

### Type Bound Procedures
The class has the following type bound procedures.

#### solidus_temp
```Fortran
class(phase_change) :: pc
pc%solidus_temp()
```
Returns the solidus temperature of the transformation `pc`.

#### liquidus_temp
```Fortran
class(phase_change) :: pc
pc%liquidus_temp()
```
Returns the liquidus temperature of the transformation `pc`.

#### solid_frac
```Fortran
class(phase_change) :: pc
pc%solid_frac(temp)
```
Returns the fraction of the solid phase at the given temperature `temp`.
The expectation is that this function is continuous and
* returns 1 for all temperatures less than or equal to the solidus
  temperature;
* returns 0 for all temperatures greater than or equal to the liquidus
  temperature; and
* returns a monotonically decreasing value between the solidus and liquidus
  temperatures of the transformation.

#### ref_liquidus_temp
```Fortran
class(phase_change) :: pc
pc%ref_liquidus_temp()
```
The latent heat is the difference in specific enthalpy between the liquid
and solid phases at the temperature value returned by this function. This
defaults to `liquidus_temp()`, but implementing classes may override it if
desired (none do currently). This is only used to define where to apply a
*given* latent heat in the construction of the specific enthalpy of the
liquid phase.

#### write_solid_frac_plotfile
```Fortran
class(phase_change) :: pc
character(*) :: filename
call pc%write_solid_frac_plotfile(filename, digits, npoints, iostat)
```
This subroutine writes a plot file of the solid fraction between the solidus
and liquidus temperatures of the phase change to the given text file
`filename`. Each line consists of a temperature value and the corresponding
solid fraction value. The number of equally-spaced data points to write is
given by `npoints` and the number of decimal digits to write for each value is
given by `digits`. The subroutine opens the specified file and the resulting
`iostat` value of that operation is returned by `iostat`.

### Concrete Implementations of the `phase_change` Class
There are two concrete extensions of the `phase_change` class.

* The [`smooth_phase_change`](../smooth_phase_change_type.F90) derived type
  implements the solid fraction as a simple smooth Hermite cubic polynomial
  of temperature between the solidus and liquidus temperatures.

* The [`tabular_phase_change`](../tabular_phase_change_type.F90) derived type
  uses a table of temperature-solid fraction values to define the solid
  fraction function.

  > **NOTE:**
  The function is smoothed using Akima interpolation which is an option for
  other tabular functions in Truchas. See the Reference Manual documentation
  for the
  [`tabular_interp`](https://www.truchas.org/docs/reference-manual/FUNCTION_Namelist/index.html#func-ti)
  parameter for more details on Akima interpolation. Note that the
  implementation of `tabular_phase_change` automatically takes care adding
  the extra data points outside the transformation interval that are needed
  to help avoid unphysical undulations in the interpolated solid fraction
  function near the endpoints of the transformation interval.

----

## The `phase_change_factory` Module
> **NB:** This is an auxiliary module to `material_factory` and is used in
the instantiation of multiphase material objects. Application code should
never be using this module directly.

The [`phase_change_factory`](../phase_change_factory.F90) module provides the
following subroutine for instantiating a `phase_change` class object. This is
the primary and recommended method for instantiating these objects.

```Fortran
class(phase_change), allocatable :: pc
type(parameter_list) :: params
call alloc_phase_change(pc, params, stat, errmsg)
```
This allocates the `phase_change` class variable `pc` according to the
specification given by the parameter list `params`. If an error is
encountered, the integer `stat` is assigned a non-zero value and the
allocatable deferred-length character `errmsg` is assigned an explanatory
message; otherwise `stat` is assigned the value 0.

There are two valid forms of the parameter list `params`, one for each of
the concrete extensions of the `phase_change` class. Their JSON-format
representations are:
* Simple smooth phase change:
  ```json
  {
    "solidus-temp": number,
    "liquidus-temp": number,
    "latent-heat": number
  }
  ```
* Tabular phase change:
  ```json
  {
    "solid-frac-table": [[number,number],[number,number], ...],
    "latent-heat": number
  }
  ```
The *latent-heat* parameter is required only when the specific enthalpy of the
liquid (i.e., high-temperature) phase is specified implicitly by its specific
heat. Otherwise, the *latent-heat* parameter is ignored.

The *solid-frac-table* parameter is a 2 x N array of (temperature, solid
fraction) pairs. The pairs may be given in either increasing or decreasing
temperature order, and the solid fraction values must be strictly monotone.
The table endpoints must specify the 0 and 1 solid fractions and these define
the liquidus and solidus temperatures.
