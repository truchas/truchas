# The `ustruc_analysis_factory` Module
The [`ustruc_analysis_factory`](../ustruc_analysis_factory.F90) module provides
the following procedure for instantiating a `ustruc_analysis` class object.

```Fortran
type(parameter_list) :: params
new_ustruc_analysis(n, params)
```
This returns a pointer to a new `class(ustruc_analysis)` object that will
perform microstructure analysis on an array of `n` independent points. The
parameter list `params` supplies the parameters that define the analysis
components to enable.

if the parameter *model-file* is defined then a custom component will be
enabled. Currently the only such component is the LDRD analysis component
described below. Otherwise the default GL analysis component described
next is enabled.

## The GL Analysis Component
The GL analysis component records the thermal gradient $G=\nabla T$ and cooling
rate $L=-\partial{T}/\partial{t}$ at the onset of solidification, and the total
time spent in the mushy zone at the completion of solidification. There are
several thresholds that determine when solidification is regarded as having
started and completed, and when the GL information is recorded. The value of
these thresholds are supplied by the `params` parameter list. There are two
mutually exclusive types of thresholds:

* Solid fraction based thresholds:
  ```json
  {
    "begin-frac": number,
    "end-frac": number,
    "gl-frac": number,          // optional
    "begin-frac-reset": number, // optional
    "end-frac-reset": number    // optional
  }
  ```
* Temperature based thresholds:
  ```json
  {
    "begin-temp": number,
    "end-temp": number,
    "gl-temp": number,          // optional
    "begin-temp-reset": number, // optional
    "end-temp-reset": number    // optional
  }
  ```
A full description of the these threshold parameters can be found in the
Reference Manual documentation of the
[`MICROSTRUCTURE` namelist](https://www.truchas.org/docs/reference-manual/MICROSTRUCTURE_Namelist/index.html).

The actual instantiation of this `class(ustruc_comp)` object is performed by
the `new_ustruc_gl` procedure defined by the
[`ustruc_gl_type`](../ustruc_gl_type.F90) module.

## Custom Analysis Components
If the *model-file* parameter is defined by `params`, the JSON-format file it
specifies is read into another `parameter_list` type variable. That parameter
list must define the *model-type* parameter, which is used to select between the
available custom analysis components, and its remaining parameters provide the
input to that component.

Currently the only recognized value of *model-type* is "ldrd" which selects the
analysis component described briefly in the
[ustruc_analysis.md](./ustruc_analysis.md#the-ustruc_ldrd-derived-type) page.
The actual instantiation of this `class(ustruc_comp)` object is performed by
the `new_ustruc_ldrd` procedure defined by the
[`ustruc_ldrd_type`](../ustruc_ldrd_type.F90) module.
