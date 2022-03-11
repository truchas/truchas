## Module `index_map_type`

The elements of a Fortran array are accessed by an *index*, which is an
integer in a contiguous range that usually starts at 1. To describe the
distribution of array elements across one or more parallel processes, it
suffices to describe how its index set is distributed across, or mapped to,
processes; the type of the array data itself is irrelevant. As its name
suggests, the `index_map_type` module defines the `index_map` derived
type which describes the mapping of a contiguous index set to one or more
parallel processes. The mapping allows for overlap between processes, and
provides array data gather and scatter procedures that are associated with
that overlap. The derived type also provides procedures for scattering
and gathering array data based on the index set, and procedures for
localizing serial indirect indexing arrays.

### Concepts
A global index set is a 1-based contiguous range of integers, {1, 2, ..., N}.
Given *P* processes, an index is assigned to a unique process according to a
block partition of the index set: The first *N_1* indices are assigned to
process 1, the next *N_2* indices to process 2, and so forth. The sum of the
*block sizes* *N_1 + N_2 + ... + N_P = N*, and a block size of 0 is allowed.
An index so assigned to a process is said to be *owned* by that process. In
addition to these, a process may include indices owned by other processes.
The owned indices are said to be *on-process*, and these additional indices,
if any, are said to be *off-process*. It is important to note that a global
index is owned by a *unique* process; this is in contrast to other schemes
(e.g. Trilinos/Tpetra) where a global index may be owned by multiple processes.

**Process-local index sets.**
For the purposes of process-local computation, the collection of all global
indices known to a process are mapped to a 1-based, contiguous, *local index
set* as follows. The contiguous block of on-process global indices are assigned,
in order, consecutive local index numbers starting at 1, and the off-process
indices are assigned, in order, consecutive local index numbers immediately
following these.

### The `index_map` Derived Type
An `index_map` type variable is a distributed parallel object. Except for
the type-bound function `global_index`, all type-bound procedures must be
called collectively by all processes.

Some procedures involve a single distinguished process called the *root*
process. By default this is MPI rank 0. While there seems little reason to
choose a different rank, one can, by giving a different rank to the optional
`root` argument to the `init` procedure.

All MPI calls use the derived type component `%comm` for their communicator
argument. The value of this component is currently hardwired to MPI_COMM_WORLD,
but could in the future be easily passed as an argument to `init`.

The next sections describe the public components and type-bound procedures
of the type.

#### Initialization
A variable
```Fortran
type(index_map) :: imap
```
must be initialized before being used by calling its `init` method. There are
several possible versions. All are collective procedures that must be called
by all processes.

1. `call imap%init(bsize [,root])`

  Each process gives its own block size using the scalar integer `bsize`.
  A block size of 0 is allowed. The resulting `imap`includes no off-process
  indices. These can be added later using the `add_offp_index` method.

2. `call imap%init(bsizes [,root])`

  The root process gives the block sizes for all processes using the rank-1
  integer array `bsizes` whose size equals the number of processes. A block
  size of 0 is allowed. `bsize` is only accessed on the root process; other
  processes may pass a 0-sized array. The resulting `imap` includes no
  off-process indices. These can be added later using the `add_offp_index`
  method.

3. `call imap%init(bsize, offp_index [,root])`

  In this extension of version 1, each process also gives a list of its
  off-process global indices using the rank-1 integer array `offp_index`.
  This will be a 0-sized array if the process has no off-process indices.
  The global indices given in `offp_index` must belong to the global index
  set as defined by the collective values of `bsize`, and must not be
  on-process indices for the process.

4. `call imap%init(bsizes, offp_counts, offp_indices [,root])`

  In this extension of version 2, the root process also gives lists of
  off-process global indices to include for all process using the rank-1
  integer arrays `offp_counts` and `offp_indices`. The size of `offp_counts`
  must equal the number of processes. `offp_counts(n)` is the number of
  off-process indices for the nth process; a value of 0 is allowed. The
  off-process indices are stored in the `off-process-indices` array in
  order of process, so that its size equals the sum of the elements of
  `offp_counts`. The global off-process indices are subject to the same
  requirements as for version 3. The subroutine arguments are only accessed
  on the root process; other processes may pass 0-sized arrays.

5. `call imap%init(domain, g_count)`

  This is a specialized version associated with ragged indexing arrays.
  `domain` is an `index_map` type variable and `g_count` is a rank-1 integer
  array. Both are `intent(in)`. `g_count` is only referenced on the root
  process of `domain`, where its size must be at least `domain%global_size`;
  other processes can pass a 0-sized array. The global index set for `imap`
  will be derived from the `g_count` array: starting with 1, the first
  consecutive `g_count(1)` integers derive from `g_count(1)`, the next
  consecutive `g_count(2)` integers derive from `g_count(2)`, and so forth.
  The size of this global index set is thus the sum of the first
  `domain%global_size` elements of `g_count`. `imap` is initialized taking
  the sum of the on-process elements of `g_count` (were it distributed
  according to `domain`) as the block size for a process. If `domain`
  includes off-process indices, these give rise to off-process indices in
  `imap`: the global indices that derive from off-process `g_count` elements
  on this process are added as off-process indices to `imap`.

#### Data Components
An initialized `index_map` variable has the following public data components.
These are **read-only** and must **never** be modified by client code. It will
be possible to have the compiler enforce this in the future by use of the
`protected` attribute being introduced in the upcoming F202X standard.

* `%onp_size`: the number of on-process indices for this process
* `%offp_size`: the number of off-process indices for this process
* `%local_size`: the size of the local index set for this process. This is
  the sum of `%onp_size` and `%offp_size`
* `%global_size`: The size of the global index set.
* `%first_gid`: the global index of the first on-process index for this process
* `%last_gid`: the global index of the last on-process index for this process
* TODO: others?

#### Miscellaneous Procedures

```fortran
type(index_map) :: imap
imap%global_index(n)
```
Returns the global index corresponding to local index `n` with respect to
`imap`. This function is not collective and can be called from any single
process. This is also an elemental function.

```Fortran
type(index_map) :: imap
call imap%add_offp_index(offp_index)
```
`offp_index` is a rank-1 array of global off-process indices with respect to
`imap` to include on this process. This may be a 0-sized array. It is an
error if the given indices are not contained in the `imap` global index set,
or if any of the indices are on-process with respect to `imap`. This is a
collective procedure. It is not currently possible to add new off-process
indices to an `index_map` object that already includes some. Thus it is an
error to call this procedure if `imap` already includes off-process indices.

#### Scatter and Gather Subroutines
These type-bound generic subroutines scatter data from the root process
to all processes and the reverse operation of collecting data from all
processes onto the root process. These are collective procedures that must
be called from all processes. Both arguments must have the same type, kind
parameters, and rank. Currently there are specific procedures for `real`
arguments of `real32` and `real64` kinds, `integer` arguments of `int32` and
`int64` kinds, and default `logical` arguments, and array ranks up to 3.

```Fortran
type(index_map) :: imap
call imap%scatter(src, dest)
```
This scatters the elements of a global array `src` on the `imap` root
process to the local array `dest` on all processes according to `imap`.
Both arguments must have the same rank, and must have the same extent in all
but the last dimension. In the multi-rank case `imap` describes the indexing
of the last dimension, and the "elements" of `src` being scattered to
corresponding "elements" of `dest` are rank-(*n-1*) sections `src(:,...,j)`
for global `imap` indices `j`. `src` is only referenced on the `imap` root
process, where its extent in the last dimension must be at least
`imap%global_size`; other processes can pass a 0-sized array. The extent
of `dest` in its last dimension must be at least `imap%onp_size`. It is
`intent(inout)` so that the value of any trailing elements is unchanged.
This subroutine wraps `MPI_Scatterv`.

If `imap` includes off-process indices it often desired that the distributed
array `dest` also include elements for these indices. To accomplish this the
`dest` array needs to be properly sized (extent in its last dimension at
least `imap%local_size`) and then simply make a subsequent call to gather
these off-process values:
``call imap%gather_offp(dest)``

```Fortran
type(index_map) :: imap
call imap%gather(src, dest)
```
This is the reverse of `scatter`. This gathers the on-process elements of
the local `src` array on each process into the global array `dest` on the
`imap` root process according to `imap`. Both arguments must have the same
rank, and must have the same extent in all but the last dimension. In the
multi-rank case `imap` describes the indexing of the last dimension, and the
"elements" of `src` that are being gathered into corresponding "elements" of
`dest` are rank-(*n-1*) sections `src(:,...,j)` for global `imap` indices `j`.
The extent of `src` in its last dimension must be at least `imap%onp_size`;
any trailing elements are ignored. The `dest`array is only referenced on the
`imap` root process, where its extent in the last dimension must be at least
`imap%global_size`; other processes can pass a 0-sized array. It is
`intent(inout)` so that the value of any trailing elements is unchanged.
This subroutine wraps `MPI_Gatherv`.

#### Off-Process Gather Subroutines

Each global index is owned by exactly one process (where it is on-process)
but may exist on other processes as an off-process index. Applied to an array
distributed according to an `index_map` object, the generic `gather_offp`
subroutines replace the value of each off-process element with the value of
its corresponding on-process element.

These are collective procedures that must be called from all process. The
arguments must have the same rank, type, and kind on all processes. Currently
there are specific procedures for rank-1, 2, and 3 arrays of types: `real`,
kinds `real32` and `real64`; `integer`, kinds `int32` and `int64`; and
default `logical`.

There are two forms of the generic `gather_offp`. The most-used form is

```Fortran
type(index_map) :: imap
call imap%gather_offp(local_data)
```
`local_data` is a rank-*n* array distributed according to `imap`. Its extent
in all but the last dimension, if any, must be the same on all processes,
and its extent in the last dimension must be at least `imap%local_size`.
In the multi-rank case `imap` describes the indexing of the last dimension of
`local_data`, and the "elements" that are being gathered are rank-(*n-1*)
sections `local_data(:,...,j)` for local `imap` indices `j`.

The second form simply splits the on-process and off-process array elements
into separate array arguments:
```Fortran
type(index_map) :: imap
call imap%gather_offp(onp_data, offp_data)
```
This is useful in less-common cases where the on and off-process elements are
not stored in contiguously in memory. In this form `onp_data` is `intent(in)`
with extent in the last dimension at least `imap%onp_size`, while `offp_data`
is `intent(inout)`, with extent in the last dimension at least `imap%offp_size`.


#### Off-Process Scatter-Reduce Subroutines

These are in some sense the reverse of `gather_offp`. However since an
on-process global index may exist on multiple other processes as an
off-process index, the operation of copying an off-process array element
to its corresponding on-process array element is an ill-defined operation.
Instead, the off-process elements are reduced with the on-process element
using an associative, commutative binary operation.

These are collective procedures that must be called from all process. The
arguments must have the same rank, type, and kind on all processes. Currently
there are specific procedures for rank-1, 2, and 3 arrays of types: `real`,
kinds `real32` and `real64`; `integer`, kinds `int32` and `int64`; and
default `logical`.


There are two forms of the generic `scatter_offp_<op>` subroutines. The
most-used form is
```Fortran
type(index_map) :: imap
call imap%scatter_offp_<op>(local_data)
```
where `<op>` may be `sum`, `min`, or `max` when `local_data` is of integer
or real type, and `or` or `and` when the array is of logical type.
`local_data` is a rank-*n* array distributed according to `imap`. Its extent
in all but the last dimension, if any, must be the same on all processes,
and its extent in the last dimension must be at least `imap%local_size`.
In the multi-rank case `imap` describes the indexing of the last dimension of
`local_data`, and the "elements" that are being scattered are rank-(*n-1*)
sections `local_data(:,...,j)` for local `imap` indices `j`. In this case
the reduction operation is applied to corresponding section elements.

The second form simply splits the on-process and off-process array elements
into separate array arguments:
```Fortran
type(index_map) :: imap
call imap%scatter_offp_<op>(onp_data, offp_data)
```
This is useful in less-common cases where the on and off-process elements are
not stored in contiguously in memory. In this form `onp_data` is `intent(inout)`
with extent in the last dimension at least `imap%onp_size`, while `offp_data`
is `intent(in)`, with extent in the last dimension at least `imap%offp_size`.

#### Localization of Indirect Indexing Arrays

Indirect indexing arrays are a necessary feature of any algorithm that
deals with unstructured meshes. These are integer arrays that map one
index set into another; for example, a rank-2 array that gives the node
indices that are the vertices of each cell. When these refer to global
indices and perhaps also exist only as a global indexing array on one
process, as is usually the case, These need to "localized" for doing
process-local computation on distributed data; that is, they need to be
distributed, if necessary, according to a distributed mapping of their
domain index set, and have their values converted to local indices
according to a distributed mapping of their range index set. This
procedure results in process-local indirect indexing arrays. It is
usually the case that a distributed indirect indexing array will reference
global indices that are not local to the process. In order to obtain
closure for the local indexing array, these off-process references need
to be included in the local index set. This is one of the key aspects
of the localization procedure, and one of the main methods by which
off-process indices are added to `index_map` objects.

The generic `localize_index_array` subroutine has the following forms:

```fortran
type(index_map) :: range
call range%localize_index_array(index [,stat])
```
`index` is a distributed indirect indexing array that takes values in the
`range` global index set. The array may be rank-1 or rank-2. Each process
"localizes" its `intent(inout)` `index` array by replacing its global index
values with their corresponding local index values with respect to `range`.
When off-process global index values are encountered, which have no automatic
corresponding local index as do on-process global indices, they are added as
additional off-process indices to `range` if not already included as such.
Thus `range` may be modified by this subroutine. An `index` value of 0 is
treated specially and not modified by this subroutine. Otherwise all `index`
values must belong to the interval [1,`range%global_size`]. It is currently
not possible to add new off-process indices to an `index_map` object that
already includes some. Thus if `range` includes off-process indices before
calling this subroutine and it is found that additional off-process indices
need to be added in order to satisfy the values of `index`, this will result
in an error. If the optional integer `stat` argument is present, it will
return -1; otherwise the program will abort with an error message.

```Fortran
type(index_map) :: domain
call domain%localize_index_array(g_index, range, l_index [,stat])
```
`range` is an `index_map` type object and `g_index` is a serial indirect
indexing array indexed by global indices from the `domain` global index set
and taking values in the `range` global index set. The array may be rank-1
or rank-2. In the latter case it is the last dimension that is indexed by
`domain` and the section `g_index(:,j)` is an array of `range` global indices
associated with `domain` global index `j`. `g_index` is only referenced on
the root process of `domain`, where its extent in the last dimension must
be at least `domain%global_size` (any additional elements are ignored);
on other processes a 0-sized array can be passed.
`l_index` is an allocatable `intent(out)` array. Its rank must be the same
as `g_index`. In the first step each process allocates `l_index` to its
proper shape: the same extents as `g_index` in all but the last dimension,
and extent `domain%local_size` in the last dimension. `g_index` is then
distributed to `l_index` according to `domain`, including any off-process
elements that `domain` may include. The last step is to localize `l_array`
by calling the previous form of `localize_index_array`. As a result, this
subroutine may modify `range`, and the same limitation on pre-existing
off-process indices noted for the previous form applies here.

NB: It is possible for an indirect indexing array to have the same global
index set for both its domain and range. In such a case `domain` and `range`
are aliased arguments, one `intent(in)` and the other `intent(inout)`. This
should not be a source of error, however, because in the first step `range`
is not referenced, and in the last step `domain` is not referenced.

```Fortran
type(index_map) :: domain
call domain%localize_index_array(g_count, g_index, range, l_count, l_index [,stat])
```
This version handles the case of a ragged indirect indexing array that is
represented by a pair of rank-1 integer arrays. In the pair (`g_count`,
`g_index`) `g_count` is indexed by the `domain` global index set and
`g_index` takes values in the `range` global index set. The first
`g_count(1)` elements of `g_index` are associated with `domain` global
index 1, the next `g_count(2)` elements with global index 2, and so forth.
If all values of `g_count` were equal to `n`, say, then the indexing array
could be more easily be represented by a single rank-2 indexing array whose
first dimension had extent equal to `n`. A similar description hold for the
output localized pair (`l_count`, `l_index`) which is indexed by the `domain`
local index set and takes values in the `range` local index set. Apart from
the difference in representation of the indirect indexing arrays, this
subroutine operates in exactly the same way as the preceding form.
