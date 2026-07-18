# Boundary Function Evaluation Design Notes

This note captures a design discussion about species mass-transfer-coefficient
(MTC) boundary/interface functions and the broader shape of the boundary
function abstractions.

## Immediate Issue

The species-only SD model has no temperature field. Some species boundary and
interface condition helpers, however, currently use the same MTC function
classes as the coupled HTSD model. Those classes expect two state fields:

```fortran
call bc%compute(t, var1, var2)
```

For HTSD this is naturally:

```fortran
call bc%compute(t, Cface, Tface)
```

For SD, the only physical state field available is concentration:

```fortran
call bc%compute(t, Cface)
```

The recent code path avoided exposing a dummy temperature at the
`species_component` call site, but the lower-level MTC implementation was made
to optionally ignore its first state argument. That works mechanically, but it
is not a satisfying design because a logical flag changes the meaning and shape
of the user-function argument vector inside a low-level mesh-function class.

## Why Cell Material Properties Were Easier

The material property abstraction already has the right shape:

```fortran
class(cell_prop_func), allocatable :: diffusivity
```

The abstract `cell_prop_func` class declares generic evaluation methods for the
state combinations the models can provide:

```fortran
call prop%compute_value(temp, value)
call prop%compute_value(conc, value)
call prop%compute_value(temp, conc, value)
```

The concrete property evaluator can implement these as thin adapters that build
the anonymous material-function state vector cell-by-cell:

```fortran
state = [T]
state = [C1, C2, ...]
state = [T, C1, C2, ...]
```

The output value is returned through a call argument, so the object does not
need to store per-evaluation results.

## What Makes Boundary Functions Awkward

The boundary/interface function classes currently mix several responsibilities:

1. Static topology: the face or interface-link indices where the condition is
   defined.
2. Evaluation policy: how user functions are evaluated on those entities.
3. Mutable evaluation results: arrays such as `value`, `deriv`, and `deriv2`.

The current `bndry_func2` and `bndry_func3` classes differ mainly by the number
of state fields passed to `compute`:

```fortran
! bndry_func2
call bc%compute(t, var)

! bndry_func3
call bc%compute(t, var1, var2)
```

They also differ in stored derivative component names (`deriv` versus
`deriv2`). Since callers hold MTC conditions as `class(bndry_func3)` and
`class(intfc_multifield_func)`, Fortran only exposes the multifield bindings
declared by those abstract types. Adding a one-field generic binding to a
concrete MTC class would not help callers through the existing abstract type.

## Possible Direction

A cleaner design would keep static topology in the object, but return transient
evaluation data through call arguments.

For example, a boundary-state-function abstraction could expose the fixed index
set as a component:

```fortran
type, abstract :: bndry_state_func
  integer, allocatable :: index(:)
contains
  procedure(value_c_iface),  deferred :: compute_value_c
  procedure(value_tc_iface), deferred :: compute_value_tc
  procedure(deriv_c_iface),  deferred :: compute_deriv_c
  procedure(deriv_tc_iface), deferred :: compute_deriv_tc
  generic :: compute_value => compute_value_c, compute_value_tc
  generic :: compute_deriv => compute_deriv_c, compute_deriv_tc
end type
```

Callers would use the persistent index for assembly, but own the temporary
evaluated data:

```fortran
call bc%compute_value(t, Cface, value)

do j = 1, size(bc%index)
  n = bc%index(j)
  Fface(n) = Fface(n) + area(n)*value(j)
end do
```

For coupled models:

```fortran
call bc%compute_deriv1(t, Cface, Tface, deriv1)
call matrix%incr_face_diag(bc%index, deriv)
```

For interface functions, the same idea applies: keep `index(:,:)` in the object
and return `value(:)` or `deriv(:,:)` arrays ordered to match it.

## Benefits

- SD can call a genuine concentration-only interface.
- HTSD can call a genuine concentration/temperature interface.
- The user-function argument vector can match the model state instead of being
  padded with dummy values.
- Mutable result arrays are not hidden inside boundary-function objects.
- The derivative naming issue (`deriv` vs `deriv2`) goes away at the caller
  interface; the caller receives the derivative it requested.
- This design resembles the successful `cell_prop_func` pattern while retaining
  static boundary/interface topology in the function object.

## Open Questions

- Should this be introduced as a new abstract class and migrated gradually, or
  should the existing `bndry_func2/3`, `intfc_field_func`, and
  `intfc_multifield_func` classes be replaced?
- Should unsupported state signatures be deferred, or should the base class
  provide default implementations that fail with a clear programming error?
- How much caller-side workspace should be reused to avoid repeated allocation
  of returned value/derivative arrays?
- Are there boundary function types outside thermal/species transport that
  rely on the current cached `value`/`deriv` components in a way that would make
  migration disruptive?

## Near-Term Recommendation

Do not keep expanding the logical-flag approach in the MTC classes. If this is
worth pursuing, prototype a small species-MTC-specific abstraction first:

```fortran
call mtc%compute_value(t, Cface, value)
call mtc%compute_value(t, Cface, Tface, value)
call mtc%compute_deriv(t, Cface, deriv)
call mtc%compute_deriv1(t, Cface, Tface, deriv1)
```

with `index` retained as a component. That would test the design where the
current problem is most visible without forcing an immediate rewrite of all
mesh-function boundary abstractions.
