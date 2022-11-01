# The `material_factory` module
The [`material_factory` module](../material_factory.F90) provides the
following procedures for instantiating `material` class objects. Both use a
parameter list to describe the typically complex inputs that are required.
This aims to be as generic as possible, not imposing any pre-defined set of
valid material or property names, leaving that to the client application code.
The one exception is in the specification of phase changes which expose
parameters specific to the models used by Truchas.

#### alloc_material
```Fortran
class(material), allocatable :: matl
type(parameter_list) :: params
call alloc_material(matl, name, params, stat, errmsg)
```
This allocates the material `matl` specified by the parameter list `params`,
setting its name to `name`. If an error is encountered, the integer `stat`
is assigned a non-zero value and the allocatable deferred-length character
`errmsg` is assigned an explanatory message. The expected structure of
`params` is the *material* structure described in the next section.

#### load_material_database
This is the primary method application code will use to instantiate a
collection of `material` class objects that can then be referenced by name.

```Fortran
type(material_database) :: matl_db
type(parameter_list) :: params
call load_material_database(matl_db, params, stat, errmsg)
```
This instantiates the `material` class objects specified by the parameter
list `params` and adds them to the material database `matl_db`. If an error
is encountered, the integer `stat` is assigned a non-zero value and the
allocatable deferred-length character `errmsg` is assigned an explanatory
message. The expected structure of `params` is the *material-list* structure
described in the next section.

## Parameter List Structure

This section describes the expected structure of the `params` argument in
terms of its JSON-format representation. The `load_material_database`
subroutine expects a *material-list* structure and the `alloc_material`
subroutine a *material* structure.

* *prop-list* is:
  > { *prop-name*: *prop-func*, ... }

  *prop-name* is a user-choice **string** used to reference the associated
  *prop-func*. Each item in the list must have a distinct name.

* *prop-func* is:

  >  **true** | **false** | **number** | **string** | *func-specifier*

  If *prop-func* is a boolean value, then the associated *prop-name* is
  interpreted as an attribute name and retained as an attribute when the
  value is *true*; otherwise the associated *prop-list* item is ignored.
  (Ideally we would use the JSON **null** value instead of a boolean value
  to signify an attribute -- a valueless property is an attribute -- but
  the `parameter_list` type is not able to represent **null**.)

  In the remaining cases:
  - **number** is the constant value of the property function
  - **string** is the name of a `scalar_func` class object stored in the
    function table to use as the property function.
  - *func-specifier* is a parameter list that defines a `scalar_func` object.
    See the subroutine
    [`alloc_scalar_func`](../functions/scalar_func_factories.F90#L256) from
    the `scalar_func_factories` module for the supported forms of this
    parameter list. This is currently unused by Truchas; the `FUNCTION`
    namelists generate `scalar_func` objects that are added to the function
    table and referenced here by name.

  See [`get_scalar_func`](../functions/scalar_func_factories.F90#L327)
  from the `scalar_func_factories` module which implements the parsing of
  *prop-func*.

* *phase* is:
  > *prop-list*

* *material* is:
  >  *single-phase-material* | *multi-phase-material*

* *material-list* is:

  > { *material-name*: *material*, ... }

  *material-name* is a user-choice **string** used to reference the associated
  *material*. Each item in the list must have a distinct name.

* *single-phase-material* is:
  > { **"properties"**: *prop-list* }

* *multi-phase-material* is:
  > {\
      **"properties"**: *prop-list*,\
      **"phases"**: *phase-list*,\
      **"phase-changes"**: *phase-change-list*\
    }

* *phase-list* is:
  > { *phase-name*: *phase*, ... }

  *phase-name* is a user-choice **string** used to reference the associated
  *phase*. Each item in the list must have a distinct name.

* *phase-change-list* is:
  > { *name-pair* : *phase-change*, ... }

  *name-pair* is a **string** that identifies the ordered pair of phases to
  connect with the associated *phase-change*. It has the form "name1:name2"
  where name1 is the name of the low-temperature phase and name2 is the name
  of the high-temperature phase.

* *phase-change* is:
  > *phase-change-simple* | *phase-change-table*

* *phase-change-simple* is:

  > {\
      **"solidus-temp"**: **number**,\
      **"liquidus-temp"**: **number**,\
      **"latent-heat"**: **number**\
    }

  NB: The "latent-heat" parameter is not required (nor used) when the
  specific enthalpy of the "liquid" phase of the phase change is given
  explicitly, rather than being formed as a definite integral of the
  specific heat.

* *phase-change-table* is:

  > {\
      **"solid-frac-table"**: *2-column-table*,\
      **"latent-heat"**: **number**\
    }

  NB: The comment on the "latent-heat" parameter for *phase-change-simple*
  applies here as well.

* *2-column-table* is:

  > [ [**number**, **number**], [**number**, **number**], ...]

  This is a list of two or more (temperature, solid fraction) pairs. It may
  be given in either increasing or decreasing temperature order, and the
  solid fraction data must be strictly monotone. The table endpoints must
  specify the 0 and 1 solid fractions.

### Example

Here is a parameter list input example that defines the properties (*values
are bogus!*) for two materials: a single phase material "graphite", and a
2-phase material "uranium". Note the different forms of property values:
constant, reference into function table, and embedded specification. For the
2-phase material note the definition of properties common to all phases, and
the definition of phase-specific properties. The property names used here
reflect Truchas usage and are not imposed by the factory procedures provided
by this module.


```json
{
  "graphite": {
    "properties": {
      "density": 1750.0,           // constant value
      "specfic-heat": {            // embedded parameter list specification
        "type": "polynomial",
        "poly-coef": [0, 1],
        "poly-powers": [2.4e2, 2.5]
      },
      "conductivity": "graphite-k" // reference into function table
    }
  },
  "uranium":{
    "properties": {
      "density": 17.5e3,
      "conductivity": 46.0
    },
    "phases": {
      "gamma": {
        "specific-heat": 161.0
      },
      "liquid": {
        "specific-heat": 200.0,
        "fluid": true
      }
    },
    "phase-changes": {
      "gamma:liquid": {
        "solidus-temp": 1400.0,
        "liquidus-temp": 1410.0,
        "latent-heat": 35800.0
      }
    }
  }
}
```
