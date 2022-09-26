# The `material_utilities` Module

The `material_utilities` module provides a collect of high-level utility
procedures that operate on materials.

> **NOTE:** At present these are ad hoc procedures that supplied a convenient
bridge from the existing Truchas drivers and solvers to the new material
package implemented in the `material` and associated classes. Some of these
should probably be deleted and others perhaps absorbed by the `material_model`
type, while other Truchas-specific material procedures are moved here from
places they don't belong.

#### required_property_check
```Fortran
class(material_model) :: model
call required_property_check(model, name, stat, errmsg)
```
Checks whether the property `name` is defined for all non-void materials in
the material model. The integer stat is assigned the value 0 if so; otherwise
it is assigned a non-zero value and the deferred-length allocatable character
`errmsg` is assigned an explanatory message.

> **Temporary usage notes:**
> * init_module: density
> * htsd_model_factory: enthalpy, conductivity, diffusivity
> * fht_model_factory: enthalpy, conductivity

#### required_fluid_property_check
```Fortran
class(material_model) :: model
call required_fluid_property_check(model, name, stat, errmsg)
```
Checks whether the property `name` is defined for all non-void phases with
the "is-fluid" attribute in the material model. The integer stat is assigned
the value 0 if so; otherwise it is assigned a non-zero value and the
deferred-length allocatable character `errmsg` is assigned an explanatory
message.

> **Temporary usage notes:** NONE (flow???)

#### optional_property_check
```Fortran
class(material_model) :: model
call optional_property_check(model, name, stat, errmsg)
```
Examines all non-void phases of the material model for the existence of the
property `name`. If all have the property, `stat` is assigned the value 0.
If none have the property, `stat` is assigned the value 1. Otherwise some
have the property while others do not, and `stat` is assigned the value -1
and the deferred-length character `errmsg` is assigned an explanatory message.

> **Temporary usage notes:**
> * htsd_model_factory: soret-coef

#### constant_property_check
```Fortran
class(material_model) :: model
call constant_property_check(model, name, stat, errmsg)
```
Checks whether a constant-valued property `name` is defined for each non-void
material in the material model (the values may differ between materials).
The integer `stat` is assigned the value 0 if so; otherwise it is assigned a
non-zero value and the deferred-length allocatable character `errmsg` is
assigned an explanatory message.

> **Temporary usage notes:**
> * init_module: density
> * EM_properties: 3 EM properties

#### define_property_default
```Fortran
class(material_model) :: model
call define_property_default(model, name, default)
```
Ensure that property `name` is defined for all non-void phases in the
material model by assigning the real constant value `default` to those
phases without the property.

> **Temporary usage notes:**
> * EM: 3 EM properties
> * init_module: specific_enthalpy (for when HT is not enabled)

#### define_fluid_property_default
```Fortran
class(material_model) :: model
call define_fluid_property_default(model, name, default)
```
Ensure that property `name` is defined for all non-void phases in the
material model with the attribute "is-fluid" by assigning the real constant
value `default` to those phases without the property.

> **Temporary usage notes:** NONE (flow???)

#### add_enthalpy_prop
```Fortran
class(material_model) :: model
call add_enthalpy_prop(model, stat, errmsg)
```
Attempts to add the "enthalpy" property to each non-void material in the
material model. This uses the [`add_enthalpy`](./materials#add_enthalpy)
method of the `material` class to build the property from existing primitive
properties. If an error is encountered, the integer `stat` is assigned
a non-zero value and the deferred-length allocatable character `errmsg` is
assigned an explanatory message; otherwise `stat` is assigned the value 0.

> **Temporary usage notes:** init_module
