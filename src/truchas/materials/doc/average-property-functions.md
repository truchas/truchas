# Average Property Functions
A frequently encountered pattern is the evaluation of a linear combination
of the form

$$ f(s;w) = \sum_{i=1}^m w_i f_i(s) $$

for a fixed set of scalar valued functions $f_i$. In specific view here is
the case where the $f_i$ are material or phase property functions and the
weights $w_i$ are the volume fractions of the materials or phases on a
particular cell, so that the linear combination is volume fraction weighted
average of the individual property values. The independent state variable
$s$ typically consists of a single temperature component, but in general may
consist of multiple state components.

This linear combination function is encapsulated by two derived types: one
for linear combinations of functions of the `matl_prop` class used for
material properties, and the other for linear combinations of functions
of the `scalar_func` class used for phase properties.

> **NOTE:** It would be straightforward to extend these data types to support
additional averaging methods such as harmonic averaging.

## The `avg_matl_prop` Derived Type

The `avg_matl_prop` derived type evaluates linear combinations of `matl_prop`
class functions. It has no public data components, only the following type
bound procedures.

### Type Bound Procedures
#### init
The `init` method is used to initialize an `avg_matl_prop` type variable.
This is a generic subroutine that has two forms. The first is

```Fortran
type(avg_matl_prop) :: prop
character(*) :: name
integer :: mids(:)
type(material_model) :: model
call prop%init(name, mids, model, stat, errmsg)
```
This initializes the `avg_matl_prop` object `prop` to compute the linear
combination of the functions $f_i$ that are associated with property `name`
and material indices given by `mids` from the material `model`. The void
material index is not allowed. If an error is encountered, the integer `stat`
is assigned a nonzero value and the deferred-length allocatable character
`errmsg` is assigned an explanatory message; otherwise `stat` is assigned the
value 0.

The second form of this generic function omits the `mids` argument and behaves
as if all non-void materials were specified in material index order.

#### compute_value / compute_deriv
```Fortran
type(avg_matl_prop) :: prop

real(r8) :: w(:), state(:), value, deriv
call prop%compute_value(w, state, value)
call prop%compute_deriv(w, state, n, deriv)

real(r8) :: w(:,:), state(:,:), value(:), deriv(:)
call prop%compute_value(w, state, value)
call prop%compute_deriv(w, state, n, deriv)
```
These procedures compute the value of the average material property function
`prop` and its derivative with respect to the `n`th component of the state
vector. The value of the weights $w_i$ are given by the elements of the
array `w`, whose order corresponds to the order of the materials specified
by the `mids` array used in the allocation of `prop`. The value of the state
vector is given by the array `state`. The result is returned in `value` and
`deriv` respectively.

There are two forms of these generic procedures. The first, just described,
computes the value at a single input. The second, computes the value for a
one-dimensional array of inputs: `w(j,:)` are the weights and `state(:,j)` is
the state vector for the `j`th input.

> **NB:** Note carefully the different ordering of dimensions of the `w` and
`state` arrays. This reflects the organization of the actual arguments in
current Truchas usage.

---

## The `avg_phase_prop` Derived Type
The `avg_phase_prop` derived type evaluates linear combinations of
`scalar_func` class functions, which are used as phase property functions.
It has no public data components, only the following type bound procedures.

### Type Bound Procedures

#### init
The `init` method is used to initialize an `avg_phase_prop` type variable.
```Fortran
type(avg_phase_prop) :: prop
character(*) :: name
integer :: pids(:)
type(material_model) :: model
call prop%init(name, pids, model, stat, errmsg)
```
This initializes the `avg_phase_prop` object `prop` to compute the linear
combination of the functions $f_i$ that are associated with property `name`
and phase indices given by `pids` from the material `model`. The index for the
void phase *is allowed* for this procedure, with 0 taken for the value of its
property. If an error is encountered, the integer `stat` is assigned a nonzero
value and the deferred-length allocatable character `errmsg` is assigned an
explanatory message; otherwise `stat` is assigned the value 0.

#### compute_value
```Fortran
type(avg_phase_prop) :: prop
real(r8) :: w(:), state(:), value
call prop%compute_value(w, state, value)
```
This computes the value of the average phase property function `prop`,
returning the result in `value`. The value of the weights $w_i$ are given by
the elements of the array `w`, whose order corresponds to the order of the
phases specified by the `pids` array used in the allocation of `prop`. The
value of the state vector is given by the array `state`.
