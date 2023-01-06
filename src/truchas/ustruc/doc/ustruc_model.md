# The `ustruc_model` Derived Type
The [`ustruc_model`](../ustruc_model_type.F90) derived type encapsulates the
microstructure analysis functionality. Client code interacts with objects of
this type only through its type bound procedures; all data components are
private. The type handles the computation of derived state data, including
data requiring mesh-based computation, before handing off the actual point-wise
analyses to its `class(ustruc_analysis)` data component. That class is
described by the [ustruc_analysis.md](./ustruc_analysis.md) page.

## Initialization
The `init` method is used to initialize a `ustruc_model` type variable:
```Fortran
type(ustruc_model) :: model
type(unstr_mesh), target :: mesh
type(parameter_list) :: params
call model%init(mesh, params)
```
The variable `model` holds a reference to `mesh`. The procedure directly
accesses two parameters from the parameter list `params`:

* *cell-set-ids* is an array of cell set IDs that specify the part of the
domain that will be included in the microstructure analysis.

* *material-fraction-threshold* specifies the volume fraction of the subject
material below which no microstructure analysis will be done. If not specified
it defaults to 0.01. This only applies to cells included in the analysis and
such cells falling below this threshold, which can change from one time step
to the next, are tagged as invalid for the analysis components.

The parameter list `params` is also passed to the
[`new_ustruc_analysis`](./ustruc_analysis_factory.md)
factory procedure which consumes the remainder of its parameters.

## Type Bound Procedures

#### set_state
```Fortran
type(ustruc_model) :: model
real(r8) :: t, tcell(:), tface(:), liq_vf(:), sol_vf(:)
call model%set_state(t, tcell, tface, liq_vf, sol_vf)
```
Sets the internal state at time `t` to the specified values: `tcell` and
`tface` are the cell and face temperatures on the mesh, and `liq_vf` and
`sol_vf` are the liquid and solid volume fractions of the subject material
on the cells. The sum of `liq_vf` and `sol_vf` will not equal 1 for cells
that contain other materials, such as void. This also resets all analysis
data to their initial values.

#### update_state
```Fortran
type(ustruc_model) :: model
real(r8) :: t, tcell(:), tface(:), liq_vf(:), sol_vf(:)
call model%update_state(t, tcell, tface, liq_vf, sol_vf)
```
Updates the internal state at time `t` to the specified values. The interface
is identical to `set_state`. The analysis components are advanced from the
previous state to this new specified state.

#### has
```Fortran
type(ustruc_model) :: model
model%has(name)
```
Returns true if one of the analysis components in `model` produces the named
data item.

#### get
```Fortran
type(ustruc_model) :: model
call model%get(name, array)
```
Returns the value of the named analysis item in the provided `array`. This is a
generic procedure with `array` a rank-1 or 2 real cell-based array. Its size in
the last dimension must be at least the number of on-process cells. A dummy
value of 0 is returned for all cells not included in the analysis. The value
returned for cells included in the analysis but without valid data is
determined by the associated analysis component (but tends to be 0 as well).

> **NOTE:** Access to the analysis data is a weak aspect of the current
> implementation. The client must know *a priori* the names of the valid data
> items, and for rank-2 data the size in the first dimension (see
> [`ustruc_driver::ustruc_output`](../ustruc_driver.F90#L348)). Furthermore,
> the analysis components produce some data items (as reported by the `has`
> function) that are not accessible by the current generic `get` subroutine.

#### get_map
```Fortran
type(ustruc_model) :: model
integer, allocatable :: map(:)
call model%get_map(map)
```
Returns a copy of the internal `map` array. This array is a list of on-process
cell indices that are included in the analysis. It is used in the writing and
reading of restart/checkpoint data.

#### get_comp_list
```Fortran
type(ustruc_model) :: model
integer, allocatable :: list(:)
call model%get_comp_list(list)
```
Returns the list of analysis component IDs. Each component has a unique ID that
is used in the calls to the `serialize` and `deserialize` procedures.

#### serialize
```Fortran
type(ustruc_model) :: model
integer(int8), allocatable :: array(:,:)
call model%serialize(cid, array)
```
Returns the internal state of the analysis component with ID `cid` in the
allocatable byte array `array`. This serves as input to the companion
`deserialize` subroutine which restores the internal state.

#### deserialize
```Fortran
type(ustruc_model) :: model
integer(int8) :: array(:,:)
call model%deserialize(cid, array)
```
Restores the internal state of the analysis component with ID `cid` using the
data provided by the byte `array`. This is the companion to the `serialize`
subroutine.
