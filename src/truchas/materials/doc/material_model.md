# The `material_model` Derived Type

The `material_model` derived type is a container for a set of materials used
by a simulation. Materials are referenced by an integer *material index* in
a contiguous 1-based range. Likewise the phases that comprise the materials
are referenced by an integer *phase index* in a 1-based range. The phases are
numbered according to the order of materials, and those for each material by
temperature (low to high). This array-like indexing of materials (phases)
allows for corresponding arrays of material (phase) volume fractions that
describe mixed material (phase) computational cells. It is this feature,
and the methods that exploit it, that distinguishes this container from the
`material_database` container. The latter may contain a large collection of
materials, whereas the former is intended to include only those few used by
a particular simulation.

The derived type also provides for an optional virtual "void" material (and
corresponding phase) that is intended to represent the absence of material
as it pertains to volume fractions. There is no underlying `material` or
`phase` object associated with void; it serves merely as a placeholder in
the array of materials and phases. If void is included, it is always the
final material and phase. This makes it straightforward to limit loops to
only the real, non-void materials/phases when desired.

## Initialization
The `init` method is used to initialize a `material_model` type variable.
```Fortran
type(material_model) :: model
type(material_database) :: matl_db
call model%init(matl_names, matl_db, stat, errmsg)
```
The rank-1 character array `matl_names` is an array of unique material names
from the material database `matl_db` to include in the model, or the special
reserved name "VOID" that signifies the virtual void material. (Any material
with the name "VOID" in the material database is ignored.) The order of the
materials in the model is the order of the given names, except that "VOID",
if present, is moved to the last position. The model stores references to the
materials stored in the database, and thus `matl_db` must persist for the
lifetime of `model`. Once initialized, the set of materials contained in the
model is fixed; no materials may be added or removed (though the materials
themselves may be modified by the addition of properties or attributes). If an
error is encountered, the integer `stat` is assigned a non-zero value and the
allocatable deferred-length character `errmsg` is assigned an explanatory
message.

## Public Data Components
An initialized `material_model` variable has the following public data
components. These are **read-only** and must **never** be modified by client
code. In the future it will be possible to have the compiler enforce this by
use of the `protected` attribute being introduced in the upcoming F202X
standard.

* `nmatl`: The number of materials, including void. Valid material indices
  are in the range [1, `nmatl`]. If the model includes void, `nmatl` is the
  material index of void.

* `nmatl_real`: the number of materials, excluding void. If the model
  includes void, this will be 1 less than `nmatl`.

* `nphase`: the number of phases, including void. Valid phase indices are in
  the range [1, `nphase`]. If the model includes void, `nphase` is the phase
  index of void.

* `nphase_real`: the number of phases, excluding void. If the model includes
  void, this will be 1 less than `nphase`.

* `have_void`: true if the material model includes void, otherwise false.

* `void_index`: this is the phase index of void if the material model
  includes void, and will equal `nphase`; otherwise it is 0.

* `is_fluid`: an array indexed by phase index that is true if the phase
  has the "fluid" attribute, and is false otherwise. The void phase is
  regarded as having the "fluid" attribute.

## Type Bound Procedures

#### has_matl
```Fortran
type(material_model) :: model
model%has_material(name)
```
Returns true if material `name` exists in the model; otherwise false is
returned. This function is elemental.

#### matl_index
```Fortran
type(material_model) :: model
model%matl_index(name)
```
Returns the material index for material `name`, or 0 if no such material
exists. `name` may be the reserved name `"VOID"`. This function is elemental.

#### matl_name
```Fortran
type(material_model) :: model
model%matl_name(mid)
```
Returns the name of the material `mid`. It is an error if `mid` is not a valid
material index.

#### has_phase
```Fortran
type(material_model) :: model
model%has_phase(name)
```
Returns true if phase `name` exists in the model; otherwise false is
returned. This function is elemental.

#### phase_index
```Fortran
type(material_model) :: model
model%phase_index(name)
```
Returns the phase index for phase `name`, or 0 if no such phase exists.
`name` may be the reserved name `"VOID"`. This function is elemental.

#### phase_name
```Fortran
type(material_model) :: model
model%phase_name(pid)
```
Returns the name of phase `pid`. It is an error if `pid` is not a valid phase
index.

#### num_matl_phase
```Fortran
type(material_model) :: model
model%num_matl_phase(mid)
```
Returns the number of phases that comprise material `mid`. It is an error if
`mid` is not a valid material index.

#### get_matl_phase_index_range
```Fortran
type(material_model) :: model
call model%get_matl_phase_index_range(mid, first, last)
```
Returns the phase index range [`first`, `last`] of the phases that comprise
material `mid`. The phases are ordered from low-temperature to
high-temperature. It is an error if `mid` is not a valid material index.

#### get_matl_phase_frac
```Fortran
type(material_model) :: model
real(r8) :: temp, beta(:)
call model%get_matl_phase_frac(mid, temp, beta)
```
Returns the phase fractions `beta` for material `mid` at the temperature
`temp`. The order of values in the `beta` array correspond to the order of the
phases of the material. The size of `beta` may be larger than the number of
phases; unused trailing array elements are left unchanged. `mid` must be a
valid material index.

#### get_matl_ref
```Fortran
type(material_model) :: model
class(material), pointer :: matl
call model%get_matl_ref(mid, matl)
```
Returns a `material` class pointer to material `mid`; `mid` must be a valid
non-void material index.

#### get_phase_ref
```Fortran
type(material_model) :: model
type(phase), pointer :: phase
call model%get_phase_ref(pid, phase)
```
Returns a `phase` type pointer to phase `pid`; `pid` must be a valid non-void
phase index.

#### const_phase_prop(n, name)
```Fortran
type(material_model) :: model
const_phase_prop(pid, name)
```
Returns the constant value of the constant property `name` of phase `pid`. It
is an error if property does not exist or is not constant valued, or if `pid`
is not a valid non-void phase index.

#### get_phase_prop
```Fortran
type(material_model) :: model
class(scalar_func), allocatable :: prop
call model%get_phase_prop(pid, name, prop)
```
Allocates the scalar function `prop` that is a copy of property `name` of
phase `pid` in the model. `pid` must be a valid non-void phase index. `prop`
is returned unallocated if the property does not exist.

#### get_matl_prop
```Fortran
type(material_model) :: model
class(matl_prop), allocatable :: prop
call model%get_matl_prop(mid, name, prop, errmsg)
```
Allocates the material property `prop` that is a copy of the material property
`name` of material `mid` in the model. `mid` must be a valid non-void material
index. If an error is encountered, `prop` is returned unallocated and an
explanatory message is assigned to the deferred-length allocatable character
`errmsg`.
