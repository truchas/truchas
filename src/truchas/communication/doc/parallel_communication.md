## The `parallel_communication` Module

The `parallel_communication` module is a simple abstraction layer over MPI
that presents a much friendlier and safer interface to some basic parallel
operations of MPI. It also provides some MPI-related data. Its features
are designed to meet the specific needs of Truchas, and aims to encapsulate
most access to MPI within its module procedures and variables. It is still
possible for Truchas application code to make direct MPI calls. However it
is preferable to use the facilities provided by this module if possible,
and if not, to instead consider extending them to meet the unmet need.
The goal is to encapsulate all direct MPI use within as few modules as
possible.

The module provides a few read-only module variables, and a few generic
collective procedures: broadcast, distribute, collate, and a selection
of global reductions. These are described in the following.

### Preliminary Notions
The module defines a public module variable `comm` that is used for
the communicator argument in all MPI calls. Currently it is set to
`MPI_COMM_WORLD`. Consistent use of this variable for all MPI calls
rather than hardwiring `MPI_COMM_WORLD` will make it simple in the
future to switch from `MPI_COMM_WORLD` to a different communicator.

The module defines a module variable `root` that is used for the root
argument in all MPI procedure calls that require one (like `MPI_Bcast`).
Currently it is hardwired to rank 0.

The module continues to use some of the same vocabulary and conventions
of PGSLib which this module in part replaces. The MPI notion of a rank is
replaced by the neutral term "processing element" or more simply "process",
and these are numbered starting at 1 instead of 0. Thus MPI rank *n* is
process *n+1*.

Many of the procedures involve a single distiguished process or rank.
In MPI this would be the root rank. Since input and output are frequently
performed in serial by this distinguished process, it has been called the
"I/O processing element", or "IO PE" or just "IOP" for short.

All procedures are collective and must be called by all processes in the
communicator.

The distribute and collate subroutines have an array argument that is only
referenced on the IO process with an expected size of `N` say. Other
processes, which still must pass an array to the argument, may pass a
0-sized array. The following idiom is a compact but still readable way
to allocate such an array:
```Fortran
real, allocatable :: array(:)
allocate(array(merge(N,0,is_iop)))
```
The logical variable `is_iop` is described below.

### Initialization, Finalization
Before any of the other module procedures are called or module variables
referenced, the module must be initialized with this call:

```fortran
call init_parallel_communication
```
The subroutine takes no arguments. If MPI has not been initialized then
the subroutine will call `MPI_init` to do so. It then initializes the
module variables. Temporarily, the subroutine also calls `PGSLib_init`
to make PGSLib procedures available to a few legacy Truchas components.

Before a graceful exit of the program, parallel communication should be
finalized with this call:
```fortran
call halt_parallel_communication
```
The subroutine takes no arguments and does nothing more than call
`MPI_finalize`, but only if MPI was initialized by the module. Otherwise
the client code, which initialized MPI, is responsible for calling
`MPI_finalize`. Temporarily, the subroutine also calls `PGSLib_finalize`.

For an ungraceful exit of the program, say due to a fatal error on a
single process, parallel communication can be aborted with this call:
```fortran
call abort_parallel_communication
```
The subroutine takes no arguments and merely calls `MPI_Abort`.

### Module variables
The module provides the following public read-only variables:

* `comm` -- The MPI communicator. Any application code that chooses to
  call MPI procedures directly should use this variable for the communicator
  argument.

* `npe` -- The number of processing elements or MPI ranks.

* `this_pe` -- The index of this processing element. This is 1 more than
  the MPI rank.

* `io_pe`  -- The index of the IO PE. It is currently hardwired to 1
  (consistent with the MPI root rank), but client code should always use
  this variable.

* `is_iop` -- A logical flag with the value `this_pe == io_pe`.

### Generic Broadcast Subroutine

```fortran
call broadcast(a)
```
copies the value of `a` on the IO process to the corresponding argument on all
other processes. `a` may be a scalar or an array of rank no greater than 3.
If an array, it must have the same shape on all processes. It must have the
same type and kind on all processes. Currently there are specific subroutines
for `int8`, `int32`, and `int64` integers; `real32` and `real64` reals;
default logical; and default character. A character argument must have the
same length on all processes.

### Generic Distribute and Collate Subroutines
These subroutines distribute data from the IO process to all processes
and the reverse operation of collating data from all processes onto the
IO process. Both arguments must have the same type and kind across all
processes. The following argument types and kinds are currently supported:
* `integer`: kinds `int8`, `int32`, and `int64`
* `real`: kinds `real32` and `real64`
* default `logical`
* default `character` (for `collate` only)

Character arguments must have the same length across all processes.

#### Scalar Data
```Fortran
call distribute(src, dest)
```
copies element *n* of the rank-1 array `src` on the IO process to the scalar
argument `dest` on process *n*. The `src` array is only accessed on the IO
process, where its size must be at least `npe`. `src` may be a 0-sized array
on other processes. This subroutine is analogous to MPI's `MPI_Scatter`.

```fortran
call collate(src, dest)
```
copies the value of the scalar argument `src` on process *n* to element *n*
of the rank-1 array `dest` on the IO process. `dest` is only referenced on
the IO process, where its size must be at least `npe` and only its first
`npe` elements are modified. `dest` may be a 0-sized array on other processes.
This subroutine is analogous to MPI's `MPI_Gather`.

#### Array Data
The generic `distribute` subroutine supports arrays of rank not greater
than 3, and the `collate` subroutine supports arrays of rank not greater
than 2.

```Fortran
call distribute(src, dest)
```
distributes elements of the rank-*n* array `src` on the IO process to the
`dest` array on all processes. Both arguments must have the same rank, and
must have the same extent in all but the last dimension. In the multi-rank
case, the copied "elements" are entire rank-(*n-1*) sections `src(:,...,j)`
for indices `j`. The number of elements copied to `dest` equals the extent
its last dimension. `src` is only referenced on the IO process, where its
extent in the last dimension must be at least the sum over all processes
of these `dest` extents. The elements of `src` on the IO process are copied
in order, first to `dest` on process 1, then to `dest` on process 2, and so
forth. `src` may be 0-sized on other processes.
This subroutine is analagous to MPI's `MPI_Scatterv`.


```Fortran
call collate(src, dest)
```
collates the elements of the rank-*n* array `src` on each process into the
array `dest` on the IO process. Both arguments must have the same rank, and
must have the same extent in all but the last dimension. In the multi-rank
case, the copied "elements" are entire rank-(*n-1*) sections `src(:,...,j)`
for indices `j`. The number of elements copied from `src` equals the extent
in its last dimension. `dest` is only referenced on the IO process, where
its extent in the last dimension must be at least the sum over all processes
of these `src` extents in order to receive all the elements of the `src`
arrays. The elements of `dest` on the IO process are assigned in order,
first from the elements of `src` on process 1, then those from `src` on
process 2, and so forth. `dest` may be 0-sized on other processes.
This subroutine is analogous to MPI's `MPI_Gatherv`

### Generic Global Reduction Functions
These reduction operations are the straightforward extensions of
similarly-named Fortran intrinsic reductions to distributed data.
Except for `global_count`, all functions return a result of the same type
and kind as their argument. The function result is returned on all processes.

#### Scalar Arguments
These functions are analagous to MPI's `MPI_allreduce`.

* ``global_all(a)``

  returns the value true if the logical scalar `a` has value true on all
  processes, and otherwise returns false.

* ``global_any(a)``

  returns the value true if the logical scalar `a` has value true on any
  process, and otherwise returns false.

* ``global_count(a)``

    returns the number of processes where the logical scalar `a` has
    value true.

* ``global_sum(a)``

  returns the sum over all processes of the scalar `a`.

* ``global_minval(a)``

  returns the minimum value over all processes of the scalar argument `a`.

* ``global_maxval(a)``

  returns the maximum value over all processes of the scalar argument `a`.

#### Array Arguments
The functions `global_all`, `global_any`, and `global_count` support
rank-1 and rank-2 arrays.

The functions `global_sum`, `global_minval`, and `global_maxval` support
rank-1 arrays only, of either integer (`int32` and `int64`) or real
(`real32` and `real64`) type.

The functions `global_sum`, `global_minval`, and `global_maxval` have
an optional logical argument `mask`. If present, it must be conformable
with the array argument.

* ``global_all(mask)``

  returns the value `global_all(all(mask))` for a logical array `mask`.

* ``global_any(mask)``

  returns the value `global_any(any(mask))` for a logical array `mask`.

* ``global_count(mask)``

  returns the value `global_sum(count(mask))` for a logical array `count`.

* ``global_sum(a [,mask])``

  returns the value `global_sum(sum(a,mask))` for an array `a`.

* ``global_minval(a [,mask])``

  returns the value `global_minval(minval(a,mask))` for an array `a`.

* ``global_maxval(a [,mask])``

  returns the value `global_maxval(maxval(a,mask))` for an array `a`.

* ``global_dot_product(a, b)``

  returns the value `global_sum(dot_product(a, b))` for rank-1 real arrays
  `a` and `b`. Kinds `real32` and `real64` are supported.

#### Ad Hoc Procedures
There are several ad hoc procedures for some specific one-off needs.

* `global_minloc(a)`

  returns the global index of the first element of a distributed rank-1
  `real(real64)` array that equals `global_minval(a)`. The global index
  is the index of that element in the array `b` were `b` the result of
  `call collate(a, b)`.

* `global_maxloc(a)`

  returns the global index of the first element of a distributed rank-1
  `real(real64)` array that equals `global_maxval(a)`. The global index
  is the index of that element in the array `b` were `b` the result of
  `call collate(a, b)`.

* `call global_maxloc_sub(array, pid, lindex [,mask])`

  finds the first element of the distributed rank-1 `real(real64)` array that
  equals `global_maxval(array, mask)`. The process that has this element and
  its index in `a` are returned in the integer arguments `pid` and `lindex`.
