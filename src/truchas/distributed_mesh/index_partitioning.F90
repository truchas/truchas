!!
!! INDEX_PARTITIONING
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Last revised 3 April 2007
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! The elements of a Fortran array are accessed or referenced by index, which
!! is an integer in a contiguous range that usually starts at 1.  To distribute
!! an array across two or more processes requires one to choose which elements
!! of the array will exist on each process.  For this choice the array data
!! itself is generally irrelevant.  All that matters is how the underlying
!! index set is distributed across the processes, and that choice determines
!! how the array, and any array based on that index set, is to be distributed.
!! This module implements the following scheme for partitioning an index set
!! across multiple processes.
!!
!! Block Partitioning With Selective Overlap.  An index set is a 1-based
!! contiguous range of integers, {1, 2, 3, ..., N}.  Given P processes, an
!! index is assigned to a unique process according to a block partition of
!! the index set: the first N1 indices are assigned to process 1, the next
!! N2 indices are assigned to process 2, etc.  The sum of the block sizes
!! N1+N2+...+NP = N, and a block size of 0 is allowed.  An index so assigned
!! to a process is said to be "owned" by that process.  In addition to these,
!! a process may include indices owned by other processes.  The owned indices
!! are said to be "on-process" and these additional indices (if any) are said
!! to be "off-process".  For the purposes of computation, the collection of
!! all indices known to a process are mapped to a process-local index set
!! (contiguous, 1-based) as follows.  The block of on-process indices are
!! assigned, in order, consecutive local index numbers starting at 1, and the
!! off-process indices are assigned, in no particular order, consecutive local
!! index numbers following the last on-process local index.
!!
!! Each global index is owned by exactly one process (where it is on-process)
!! but may exist on other processes as an off-process index.  The off-process
!! indices are construed as the overlap of the index set partition.  The
!! mapping from off-process indices to their corresponding on-process index
!! is a many-to-one mapping, and a fundamental operation associated with a
!! distributed array is to communicate data associated with these indices
!! between processes according this mapping.  Communicating data to an
!! off-process index from its corresponding on-process index is known as a
!! "gather" operation, and communicating data from an off-process index to
!! its corresponding on-process index is known as a "scatter" operation.  
!! Because the mapping is many-to-one, the scatter operation may have
!! multiple data elements arriving at an on-process index and so the scatter
!! operation must include some sort of reduction operation, like summation,
!! in order to form a single data element for the index.
!! 
!! The module provides the derived type IP_DESC with private components that
!! encapsulates all the info describing such an index set partition, including
!! the info needed to perform the inter-process communication.  The following
!! procedures operate on a variable of this type.
!!
!! 1. ACCESSORS
!!
!!  ONP_SIZE(THIS) returns the number of on-process indices for the process
!!    for the index partition THIS.  Thus local indices in the range 1 through
!!    ONP_SIZE(THIS) are on-process (or owned) indices.  The function is pure
!!    and can be used in a specification expression.
!!  
!!  OFFP_SIZE(THIS) returns the number of off-process indices for the process
!!    for the index partition THIS.  Thus local indices in the range
!!    ONP_SIZE(THIS) + 1 through ONP_SIZE(THIS) + OFFP_SIZE(THIS) correspond
!!    to the off-process indices.  The function is pure and can be used in a
!!    specification expression.
!!
!!  LOCAL_SIZE(THIS) returns the size of the local index set for the process
!!    for the index partition THIS.  The return value is just the sum of the
!!    on-process and off-process sizes.  The function is pure and can be used
!!    in a specification expression.
!!
!!  GLOBAL_SIZE(THIS) returns the size of the global index set corresponding
!!    to the index partition THIS.  All processes return the same value.
!!    The function is pure and can be used in a specification expression.
!!
!!  GLOBAL_INDEX(THIS, N) returns the global index corresponding to local
!!    index N of the process in the index partition THIS.  N must lie in the
!!    range [1, LOCAL_SIZE(THIS)], and the returned value is in the range
!!    [1, GLOBAL_SIZE(THIS)].  If N is out-of-range, -1 is returned.  The
!!    function is elemental.
!!
!!  FIRST_INDEX(THIS) returns the first index of the global range assigned to
!!    the process for the index partition THIS.  Note that FIRST_INDEX(THIS)=
!!    GLOBAL_INDEX(THIS, 1).
!!
!!  LAST_INDEX(THIS) returns the last index of the global range assigned to the
!!    process for the index partition THIS.  Note that LAST_INDEX(THIS)=
!!    GLOBAL_INDEX(THIS, ONP_SIZE(THIS)).
!!
!! 2. CONSTRUCTORS/DESTRUCTORS
!!
!!  CALL CREATE (THIS, BSIZE) creates the index partition THIS for the block
!!    partition defined by the block sizes given in BSIZE.  BSIZE is either
!!    an integer scalar or rank-1 array. In the scalar case the value of BSIZE
!!    is the block size assigned to the process.  In the array case BSIZE is
!!    ignored on all but the I/O process, where its size must equal the number
!!    of processes, and the value of its elements are the block sizes assigned
!!    to the corresponding processes.  A block size may be 0, and the sum of
!!    the block sizes, either across processes in the scalar case or of the
!!    I/O process array in the array case, is the size of the corresponding
!!    global index set.  Note that this creates an index partition without
!!    overlap; all local indices are on-process.  Additional off-process
!!    indices may be added later by calling ADD_OFFP_INDEX.
!!
!!  CALL CREATE (THIS, BSIZE, OFFP_INDEX) creates the index partition THIS
!!    where the integer BSIZE is the block size assigned to the process and
!!    the rank-1 integer array OFFP_INDEX contains the off-process indices
!!    to be included on the process.  BSIZE may be 0 and OFFP_INDEX may be
!!    a 0-sized array.  The size of the corresponding global index set is the
!!    sum of BSIZE across processes.  It is an error for any of the specified
!!    off-process indices to be owned by the process, or to not belong to the
!!    global index set.
!!
!!  CALL CREATE (THIS, BSIZE, OFFP_SIZE, OFFP_INDEX) creates the index
!!    partition THIS where BSIZE gives the block sizes and OFFP_SIZE/OFFP_INDEX
!!    give the off-process indices to include on the processes.  The arguments
!!    BSIZE, OFFP_SIZE and OFFP_INDEX are rank-1 integer arrays and are ignored
!!    on all but the I/O process, where the size of BSIZE and OFFP_SIZE must
!!    equal the number of processes and the size of OFFP_INDEX must equal the
!!    sum of the values of OFFP_SIZE.  On the I/O process, the elements of BSIZE
!!    are the block sizes of the corresponding processes, and OFFP_INDEX gives
!!    the off-process indices in packed format: the first OFFP_SIZE(1) elements
!!    are the off-process indices for the first process, the next OFFP_SIZE(2)
!!    elements are the off-process indices for the second processes, and so on.
!!    This version of CREATE is just like the preceding version except that all
!!    the partition info is provided on the I/O process rather being distributed
!!    across the processes.
!!    
!!  CALL ADD_OFFP_INDEX (THIS, OFFP_INDEX) adds the global indices given in the
!!    rank-1 integer array OFFP_INDEX to the existing off-process indices
!!    for the process in the index partition THIS.  THIS must have been defined
!!    by one of the previous CREATE calls, and it is an error for any of the
!!    off-process indices specified to be owned by the process or not belong to
!!    the global index set.  Currently, use is restricted to adding off-process
!!    indices to a non-overlapping index partition created by the first version
!!    of the CREATE call; once off-process indices have been added no further
!!    ones can be added.  This restriction may be relaxed in the future.
!!
!!  CALL DESTROY (THIS) deallocates any storage associated with the index
!!    partition THIS and returns it to its default initialization state.
!!
!! The method for deciding which indices should be included as off-process
!! is often simply to include just those indices that need to be referenced
!! but aren't owned by the process.  Consider an indirect indexing array that
!! maps cell indices into node indices. When the array is distributed according
!! to a partition of the cell index set, the portion of the array that exists
!! on a process may refer to nodes that are not owned by the process.  These
!! node indices should be included as off-process indices in the partition
!! of the node index set in order to obtain closure of the local portion of
!! the indirect indexing array.  In addition, the global values of the local
!! indexing array need to be mapped to the corresponding local node indices.
!! The following routine performs both these tasks.
!!
!!  CALL LOCALIZE_INDEXING_ARRAY (G_INDEX, DOMAIN, RANGE, L_INDEX, OFFP_INDEX)
!!    takes the global indexing array G_INDEX presented on the I/O process and
!!    distributes it across processes according to the index partition DOMAIN.
!!    The local array values, which are global indices in an index set with
!!    partition RANGE, are then mapped to local indices, and the resulting
!!    localized indexing array returned in L_INDEX.  Any unknown index that
!!    is referenced is identified, and a list of all such indices returned
!!    in OFFP_INDEX.  In the localization process these unknown indices are
!!    mapped, in order, to new local indices starting from the last known
!!    local index.  It is necessary to subsequently call ADD_OFFP_INDEX,
!!    passing the returned OFFP_INDEX array, in order to make the necessary
!!    updates to the RANGE index partition.  The arrays G_INDEX and L_INDEX
!!    must be of the same rank, and may be either rank-1 or rank-2.  In the
!!    multidimensional case the last dimension is taken to be the distributed
!!    axis.  L_INDEX is a pointer and is allocated by the subroutine.
!!    OFFP_INDEX is a rank-1 integer pointer array which is allocated by the
!!    subroutine.  A 0-sized array is returned if no unknown local indices
!!    were found.  Only when a 0-sized array is returned on all processes is
!!    it unnecessary to make a subsequent call to ADD_OFFP_INDEX.  The caller
!!    is responsible for deallocating the OFFP_INDEX array when it is no
!!    longer needed.
!!
!!    Added 5/11/2009: The special index array value of 0 is now allowed with
!!    the interpretation that it is to be ignored.  The 0 value is propogated
!!    to the localized index array unchanged, and since it is not regarded as
!!    belonging to the range index set, it is never included in the returned
!!    OFFP_INDEX array.
!!
!! INTERPROCESS COMMUNICATION
!!
!!  CALL GATHER_BOUNDARY (THIS, ONP_DATA, OFFP_DATA) gathers off-process data
!!    from the corresponding on-process data on other processes.  This is a
!!    global communication procedure.  The process' on-process data is passed
!!    in the intent-in array ONP_DATA, and its off-process data is returned in
!!    the intent-out array OFFP_DATA.  THIS is the index partition for the
!!    distributed array.  The arrays ONP_DATA and OFFP_DATA must have the same
!!    rank and type.  The last dimension is the distributed axis and its size
!!    must equal the on-process and off-process sizes, respectively, of the
!!    index partition THIS.  ONP_DATA and OFFP_DATA must have the same extents
!!    in each of the remaining dimensions (if any). Currently, the data arrays
!!    may be rank 1, 2 or 3, and of default logical, integer, real, or double
!!    precision kinds.
!!
!!  CALL GATHER_BOUNDARY (THIS, LOCAL_DATA).  This is a variation of the
!!    preceding procedure in which the on and off-process data is packed into
!!    a single intent-inout data array LOCAL_DATA.  The on-process data forms
!!    the first part of the array and is followed by the off-process data.
!!    This packing of data corresponds to the organization of the on and
!!    off-process indices within the local index set.
!!
!!  CALL SCATTER_BOUNDARY_<op> (THIS, ONP_DATA, OFFP_DATA) scatters off-process
!!    data into the corresponding on-process data on other processes.  This is
!!    a global communications procedure. The process' on-process data is passed
!!    in the intent-inout array ONP_DATA, and its off-process data is passed in
!!    the intent-in array OFFP_DATA.  THIS is the index partition for the
!!    distributed array.  The arrays ONP_DATA and OFFP_DATA must have the same
!!    rank and type.  The last dimension is the distributed axis and its size
!!    must equal the on-process and off-process sizes, respectively, of the
!!    index partition THIS.  ONP_DATA and OFFP_DATA must have the same extents
!!    in the remaining dimensions (if any).  Currently the data arrays must be
!!    rank-1 and of default logical, integer, real, or double precision type.
!!    The on-process data is combined with its corresponding off-process data
!!    from other processes using a commutative, associative, reduction operator
!!    <op>.  Currently <op> may be SUM, MIN, or MAX for the arithmetic data
!!    types, and OR or AND for the logical data type.
!!
!!    N.B.: The off-process data is left inconsistent with its corresponding
!!    on-process data after this call.  If needed, make a subsequent call to
!!    GATHER_BOUNDARY.
!!
!!  CALL SCATTER_BOUNDARY_<op> (THIS, LOCAL_DATA).  This is a variation of the
!!    preceding procedure in which the on and off-process data is packed into
!!    a single intent-inout data array LOCAL_DATA.  The on-process data forms
!!    the first part of the array and is followed by the off-process data.
!!    This packing of data corresponds to the organization of the on and
!!    off-process indices within the local index set.
!!
!!    N.B.: The off-process data is left inconsistent with its corresponding
!!    on-process data after this call.  If needed, make a subsequent call to
!!    GATHER_BOUNDARY.
!!

#include "f90_assert.fpp"

module index_partitioning

  use pgslib_module
  use parallel_communication

  implicit none
  private

  public :: localize_index_array, localize_index_struct
  public :: gather_boundary, boundary_is_current
  public :: scatter_boundary_sum, scatter_boundary_min, scatter_boundary_max
  public :: scatter_boundary_and, scatter_boundary_or

  type, public :: ip_desc
    private
    integer :: onP_size_ = 0    ! number of indices assigned to this process (on-process)
    integer :: offP_size_ = 0   ! number of off-process indices referenced from this process
    integer :: local_size_ = 0  ! number of local indices (on and off-process)
    integer :: global_size_ = 0 ! size of the global index set
    integer :: first = 0        ! first global index of the range assigned to this process
    integer :: last  = 0        ! last global index of the range assigned to this process
    integer, pointer :: offP_index(:) => null() ! off-process indices referenced from this process
    integer, pointer :: dup_index(:) => null()  ! index from duplicate buffer to on-process data
    integer, pointer :: sup_index(:) => null()  ! index from off-process data to supplement buffer
    type(PGSLib_GS_Trace), pointer :: trace => null() ! PGSLib communication trace
  contains
    procedure, private :: create_1
    procedure, private :: create_2
    procedure, private :: create_3
    procedure, private :: create_4
    generic   :: init => create_1, create_2, create_3, create_4
    procedure :: add_offP_index
    procedure :: onP_size
    procedure :: offP_size
    procedure :: local_size
    procedure :: global_size
    procedure :: global_index
    procedure :: first_index
    procedure :: last_index
    procedure :: defined => defined_ip_desc
    final :: ip_desc_delete
  end type ip_desc

  interface localize_index_array
    module procedure localize_index_array_1, localize_index_array_2, localize_index_array_3
  end interface

  interface gather_boundary
    module procedure gather_boundary_L1, gather_boundary_XL1
    module procedure gather_boundary_L2, gather_boundary_XL2
    module procedure gather_boundary_L3, gather_boundary_XL3
    module procedure gather_boundary_I1, gather_boundary_XI1
    module procedure gather_boundary_I2, gather_boundary_XI2
    module procedure gather_boundary_I3, gather_boundary_XI3
    module procedure gather_boundary_S1, gather_boundary_XS1
    module procedure gather_boundary_S2, gather_boundary_XS2
    module procedure gather_boundary_S3, gather_boundary_XS3
    module procedure gather_boundary_D1, gather_boundary_XD1
    module procedure gather_boundary_D2, gather_boundary_XD2
    module procedure gather_boundary_D3, gather_boundary_XD3
  end interface

  interface scatter_boundary_sum
    module procedure scatter_boundary_sum_I1, scatter_boundary_sum_XI1
    module procedure scatter_boundary_sum_S1, scatter_boundary_sum_XS1
    module procedure scatter_boundary_sum_D1, scatter_boundary_sum_XD1
  end interface

  interface scatter_boundary_min
    module procedure scatter_boundary_min_I1, scatter_boundary_min_XI1
    module procedure scatter_boundary_min_S1, scatter_boundary_min_XS1
    module procedure scatter_boundary_min_D1, scatter_boundary_min_XD1
  end interface

  interface scatter_boundary_max
    module procedure scatter_boundary_max_I1, scatter_boundary_max_XI1
    module procedure scatter_boundary_max_S1, scatter_boundary_max_XS1
    module procedure scatter_boundary_max_D1, scatter_boundary_max_XD1
  end interface

  interface scatter_boundary_and
    module procedure scatter_boundary_AND_L1, scatter_boundary_AND_XL1
  end interface

  interface scatter_boundary_or
    module procedure scatter_boundary_OR_L1, scatter_boundary_OR_XL1
  end interface

  interface boundary_is_current
    module procedure boundary_is_current_L1, boundary_is_current_L2!, boundary_is_current_L3
    module procedure boundary_is_current_S1, boundary_is_current_S2!, boundary_is_current_S3
    module procedure boundary_is_current_D1, boundary_is_current_D2!, boundary_is_current_D3
  end interface

contains

  !! Final subroutine for IP_DESC objects.
  subroutine ip_desc_delete (this)
    type(ip_desc), intent(inout) :: this
    if (associated(this%offP_index)) deallocate(this%offP_index)
    if (associated(this%dup_index))  deallocate(this%dup_index)
    if (associated(this%sup_index))  deallocate(this%sup_index)
    if (associated(this%trace)) call PGSLib_Deallocate_Trace (this%trace)
  end subroutine ip_desc_delete

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! CREATE specific procedures
 !!

  subroutine create_1 (this, bsize)

    class(ip_desc), intent(out) :: this
    integer, intent(in) :: bsize(:)

    integer :: sizes(nPE)

    if (is_IOP) then
      ASSERT( size(bsize) == nPE )
      ASSERT( all(bsize >= 0) )
      sizes = bsize
    end if

    call broadcast (sizes)

    this%onP_size_ = sizes(this_PE)
    this%last  = sum(sizes(1:this_PE))
    this%first = this%last - this%onP_size_ + 1
    this%local_size_ = this%onP_size_
    this%global_size_ = sum(sizes)

  end subroutine create_1


  subroutine create_2 (this, bsize)

    class(ip_desc), intent(out) :: this
    integer, intent(in) :: bsize

    integer :: sizes(nPE)

    call collate (sizes, bsize)
    call create_1 (this, sizes)

  end subroutine create_2


  subroutine create_3 (this, bsize, offP_size, offP_index)

    class(ip_desc), intent(out) :: this
    integer, intent(in) :: bsize(:)
    integer, intent(in) :: offP_size(:)
    integer, intent(in) :: offP_index(:)

    integer :: n
    integer, allocatable :: index(:)

    call create_1 (this, bsize)

    if (is_IOP) then
      ASSERT( size(offP_size) == nPE )
      ASSERT( all(offP_size >= 0) )
      ASSERT( sum(offP_size) == size(offP_index) )
    end if

    call distribute (n, offP_size)
    allocate(index(n))
    call distribute (index, offP_index)
    call add_offP_index (this, index)
    deallocate(index)

  end subroutine create_3


  subroutine create_4 (this, bsize, offP_index)

    class(ip_desc), intent(out) :: this
    integer, intent(in) :: bsize
    integer, intent(in) :: offP_index(:)

    call create_2 (this, bsize)
    call add_offP_index (this, offP_index)

  end subroutine create_4

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Accessors to the public read-only components of the partition THIS.
 !!

  pure integer function onP_size (this)
    class(ip_desc), intent(in) :: this
    onP_size = this%onP_size_
  end function onP_size

  pure integer function offP_size (this)
    class(ip_desc), intent(in) :: this
    offP_size = this%offP_size_
  end function offP_size

  pure integer function local_size (this)
    class(ip_desc), intent(in) :: this
    local_size = this%local_size_
  end function local_size

  pure integer function global_size (this)
    class(ip_desc), intent(in) :: this
    global_size = this%global_size_
  end function global_size

  elemental integer function global_index (this, n)
    class(ip_desc), intent(in) :: this
    integer, intent(in) :: n
    global_index = -1
    if (n < 1) return
    if (n <= this%onP_size_) then
      global_index = this%first + n - 1
    else if (n <= this%local_size_) then
      global_index = this%offP_index(n-this%onP_size_)
    end if
  end function global_index
  
  pure integer function first_index (this)
    class(ip_desc), intent(in) :: this
    first_index = this%first
  end function first_index

  pure integer function last_index (this)
    class(ip_desc), intent(in) :: this
    last_index = this%last
  end function last_index

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ADD_OFFP_INDEX
 !!
 !! Should we require distinct indices in the list?
 !! We require no previous off-process indices.
 !!
 !! NOTES
 !!
 !! (1) PGSLib rewrites the index array.  Off-process references are mapped to
 !! negative values intended to point into a special 'supplement' buffer used
 !! to communicate data between processes.  Since by design all references
 !! are off-process, -SUP_INDEX is the desired index array.
 !!
 !! (2) If the original off-process index values were distinct, as is normally
 !! the case, then SUP_INDEX should just be a permutation.  If in addition the
 !! values were suitably ordered, SUP_INDEX would be the identity permutation,
 !! and we can dispense entirely with the usual indirect copy between the
 !! supplement buffer and off-process data that occurs during gathers and
 !! scatters.  This special case is signaled by an unassociated SUP_INDEX
 !! component.
 !!
 !! (3) With regard to (2), experimentation has revealed that the off-process
 !! indices must be in ascending order of their assigned process, and in
 !! descending order of their values among indices assigned to the same
 !! process in order to obtain an identity mapping.  For now we have not gone
 !! to the effort of reordering the off-process index vector accordingly.
 !!
 !! (4) To address the issue in (3), it might be appropriate to modify PGSLib
 !! directly.  It appears that index values are inserted into some sort of
 !! data table as they are encountered in the indexing array, and the
 !! PGSLib-internal order of values is generated by reading the values out
 !! of the table in a LIFO fashion (thus reversing their order) instead of
 !! a FIFO fashion.  If the latter were done, then an ordered off-process
 !! index array, as generated by LOCALIZE_INDEX_ARRAY for example, would be
 !! exactly what is required to obtain an identity permutation.
 !!

  subroutine add_offP_index (this, offP_index)

    use permutations, only: is_identity_perm

    class(ip_desc), intent(inout) :: this
    integer, intent(in) :: offP_index(:)
    
    integer, pointer :: old_offP_index(:)

    ASSERT( this%defined() )
    ASSERT( minval(offP_index) >= 1 )
    ASSERT( maxval(offP_index) <= this%global_size_ )
    ASSERT( all((offP_index < this%first) .or. (offP_index > this%last)) )

    !! Record the off-process indices to be referenced from this process.
    if (associated(this%offP_index)) then ! append to the existing list
      old_offP_index => this%offP_index
      allocate(this%offP_index(this%offP_size_+size(offP_index)))
      this%offP_index(:this%offP_size_) = old_offP_index
      this%offP_index(this%offP_size_+1:) = offP_index
      deallocate(old_offP_index)
    else  ! create a new list
      allocate(this%offP_index(size(offP_index)))
      this%offP_index = offP_index
    end if
    this%offP_size_ = size(offP_index)
    this%local_size_ = this%onP_size_ + this%offP_size_

    !! Generate the PGSLib communication trace for the off-process indices.
    if (associated(this%sup_index)) deallocate(this%sup_index)
    allocate(this%sup_index(size(this%offP_index)))
    this%sup_index = this%offP_index
    if (associated(this%trace)) call PGSLib_Deallocate_Trace (this%trace)
    this%trace => pgslib_setup_trace(this%sup_index, this%onP_size_)

    if (associated(this%dup_index)) deallocate(this%dup_index)
    allocate(this%dup_index(pgslib_size_of_dup(this%trace)))
    this%dup_index = pgslib_dup_index(this%trace)

    !! Normalize the rewritten off-process index array (see Note 1).
    this%sup_index = -this%sup_index
    INSIST( all(this%sup_index > 0) )

    !! Eliminate SUP_INDEX if possible (see Note 2).
    if (is_identity_perm(this%sup_index)) then
      INSIST( pgslib_size_of_sup(this%trace) == size(this%sup_index) )
      deallocate(this%sup_index)
    end if

    ASSERT( gather_boundary_verified(this) )

  end subroutine add_offP_index

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! LOCALIZE_INDEX_ARRAY specific procedures
 !!
 !! This procedure takes a global indexing array G_INDEX presented on
 !! the IO processor and distributes it across the processors according to the
 !! partition, DOMAIN, of its last dimension (the distributed axis).  The local
 !! indexing array values, which are global indices in an index set with
 !! partition RANGE, are then mapped to local indices, and the resulting local
 !! indexing arrays are returned in L_INDEX.  Any unknown off-process index that
 !! is referenced is identified, and a list of all such indices is returned in
 !! OFFP_INDEX.  These unknown off-process indices are mapped to local indices
 !! starting from the last known local index in sequence.
 !!
 !! IMPLEMENTATION NOTES
 !!
 !! (1) The identification and manipulation of the off-process indices is
 !! greatly simplified by the use of a temporary integer vector whose size is
 !! the global size of the range index set, even though only a relatively small
 !! number of the elements are typically used.  Memory requirements could be
 !! reduced, if necessary, by employing a more sophisticated data structure.
 !! Regardless, the ultimate list of off-process indices is stored compactly.
 !!

  subroutine localize_index_array_1 (g_index, domain, range, l_index, offP_index)

    integer, intent(in) :: g_index(:)
    type(ip_desc), intent(in) :: domain, range
    integer, allocatable :: l_index(:), offP_index(:)

    integer :: j
    integer, allocatable :: map(:)

    ASSERT( domain%defined() )

    if (is_IOP) then
      ASSERT( size(g_index) == domain%global_size_ )
    end if

    !! Distribute the global indexing array according to the domain partition.
    allocate(l_index(domain%local_size_))
    call distribute (l_index(:domain%onP_size_), g_index)
    if(associated(domain%offP_index)) call gather_boundary (domain, l_index)

    call localize_index_array_2 (range, l_index, offP_index)

  end subroutine localize_index_array_1


  subroutine localize_index_array_2 (range, index, offP_index)

    type(ip_desc), intent(in) :: range
    integer, intent(inout) :: index(:)
    integer, allocatable :: offP_index(:)

    integer :: j
    integer, allocatable :: map(:)

    ASSERT( range%defined() )
    ASSERT( minval(index) >= 0 )
    ASSERT( maxval(index) <= range%global_size_ )

    !! Identify all unknown off-process index references (map>0).
    allocate(map(0:range%global_size_))
    map = 0
    do j = 1, size(index)
      map(index(j)) = index(j)
    end do
    map(range%first:range%last) = 0   ! on-process indices are known
    if (associated(range%offP_index)) then  ! known off-process indices
      do j = 1, size(range%offP_index)
        map(range%offP_index(j)) = 0
      end do
    end if

    !! Extract the list of unknown referenced indices.
    !allocate(offP_index(count(map>0)))
    offP_index = pack(map, mask=(map>0))

    !! Local numbering of the known off-process indices; by convention these
    !! were numbered sequentially following the on-process indices.
    if (associated(range%offP_index)) then
      do j = 1, size(range%offP_index)
        map(range%offP_index(j)) = range%onP_size_ + j
      end do
    end if

    !! Generate a local numbering of the unknown off-process indices; these
    !! will be numered sequentially following the 'known' local indices.
    do j = 1, size(offP_index)
      map(offP_index(j)) = range%local_size_ + j
    end do

    !! Remap the local index array values to the local numbering.
    do j = 1, size(index)
      if ((index(j) < range%first) .or. (index(j) > range%last)) then
        index(j) = map(index(j))
      else
        !! 1-based numbering of the on-process indices.
        index(j) = index(j) - range%first + 1
      end if
    end do

    deallocate(map)

  end subroutine localize_index_array_2


  subroutine localize_index_array_3 (g_index, domain, range, l_index, offP_index)

    integer, intent(in) :: g_index(:,:)
    type(ip_desc), intent(in) :: domain, range
    integer, allocatable :: l_index(:,:), offP_index(:)

    integer :: i, j
    integer, allocatable :: map(:)

    ASSERT( domain%defined() )
    ASSERT( range%defined() )

    if (is_IOP) then
      ASSERT( size(g_index,2) == domain%global_size_ )
      ASSERT( minval(g_index) >= 0 )
      ASSERT( maxval(g_index) <= range%global_size_ )
    end if

    !! Distribute the global indexing array according to the domain partition.
    allocate(l_index(size(g_index,1),domain%local_size_))
    call distribute (l_index(:,:domain%onP_size_), g_index)
    if(associated(domain%offP_index)) call gather_boundary (domain, l_index)

    !! Identify all unknown off-process index references (map>0).
    allocate(map(0:range%global_size_))
    map = 0
    do j = 1, size(l_index, dim=2)
      do i = 1, size(l_index, dim=1)
        map(l_index(i,j)) = l_index(i,j)
      end do
    end do
    map(range%first:range%last) = 0   ! on-process indices are known
    if (associated(range%offP_index)) then  ! known off-process indices
      do j = 1, size(range%offP_index)
        map(range%offP_index(j)) = 0
      end do
    end if

    !! Extract the list of unknown referenced indices.
    !allocate(offP_index(count(map>0)))
    offP_index = pack(map, mask=(map>0))

    !! Local numbering of the known off-process indices; by convention these
    !! were numbered sequentially following the on-process indices.
    if (associated(range%offP_index)) then
      do j = 1, size(range%offP_index)
        map(range%offP_index(j)) = range%onP_size_ + j
      end do
    end if

    !! Generate a local numbering of the unknown off-process indices; these
    !! will be numered sequentially following the 'known' local indices.
    do j = 1, size(offP_index)
      map(offP_index(j)) = range%local_size_ + j
    end do

    !! Remap the local index array values to the local numbering.
    do j = 1, size(l_index, dim=2)
      do i = 1, size(l_index, dim=1)
        if ((l_index(i,j) < range%first) .or. (l_index(i,j) > range%last)) then
          l_index(i,j) = map(l_index(i,j))
        else
          !! 1-based numbering of the on-process indices.
          l_index(i,j) = l_index(i,j) - range%first + 1
        end if
      end do
    end do

    deallocate(map)

  end subroutine localize_index_array_3


  subroutine localize_index_struct (g_count, g_index, domain, range, l_count, l_index, offP_index)

    integer, intent(in) :: g_index(:), g_count(:)
    type(ip_desc), intent(in) :: domain, range
    integer, allocatable, intent(out) :: l_index(:), l_count(:), offP_index(:)

    integer :: j, n, offset
    integer, allocatable :: map(:), bsize(:), domain_offP_index(:), g_ghosts(:), x_index(:)

    ASSERT( domain%defined() )
    ASSERT( range%defined() )

    if (is_IOP) then
      ASSERT( minval(g_count) >= 0 )
      ASSERT( sum(g_count) == size(g_index) )
      ASSERT( size(g_count) == domain%global_size_ )
      ASSERT( minval(g_index) >= 0 )
      ASSERT( maxval(g_index) <= range%global_size_ )
    end if
    
    !! Distribute the global count array according to the domain partition.
    allocate(l_count(domain%local_size_))
    call distribute (l_count(:domain%onP_size_), g_count)
    if (associated(domain%offP_index)) call gather_boundary (domain, l_count)
    
    !! Distribute the global indexing array according to the domain partition.
    allocate(l_index(sum(l_count)))
    call distribute (l_index(:sum(l_count(:domain%onP_size_))), g_index)
    ! if (associated(domain%offP_index)) "gather off-process part of l_index"
    !! Here we need to gather the off-process indexing data, ideally via a
    !! call to gather_boundary as we do in elsewhere (the commented line above).
    !! But no such version of gather_boundary exists; it would require the
    !! ability to communicate variable amounts of data per off-process index
    !! and that capability doesn't exist through PGSLib.  We would be stuck
    !! were it not for having the global indexing data G_INDEX at hand.  So
    !! instead of the processes gathering their off-process indexing data
    !! from the distributed indexing data, we collect the off-process
    !! indexing data on the IO process from G_INDEX, and then distribute it.
    !! The somewhat complicated block of code below does this.
    if (associated(domain%offP_index)) then
      !! Collate the off-process index sets.
      allocate(bsize(merge(nPE,0,is_IOP)))
      call collate (bsize, size(domain%offP_index))
      allocate(domain_offP_index(merge(sum(bsize),0,is_IOP)))
      call collate (domain_offP_index, domain%offP_index)
      !! Allocate space to hold the collated ghost indexing data.
      call collate (bsize, size(l_index) - sum(l_count(:domain%onP_size_)))
      allocate(g_ghosts(merge(sum(bsize),0,is_IOP)))
      deallocate(bsize)
      if (is_IOP) then
        !! Generate the indexing array into g_index.
        allocate(x_index(size(g_count)+1))
        x_index(1) = 1
        do j = 1, size(g_count)
          x_index(j+1) = x_index(j) + g_count(j)
        end do
        !! Collect the ghost indexing data.
        offset = 0
        do j = 1, size(domain_offP_index)
          n = domain_offP_index(j)
          associate (list => g_index(x_index(n):x_index(n+1)-1))
            g_ghosts(offset+1:offset+size(list)) = list
            offset = offset + size(list)
          end associate
        end do
        deallocate(x_index)
      end if
      !! Distribute the collected ghost indexing data.
      associate (l_index_offP => l_index(1+sum(l_count(:domain%onP_size_)):))
        call distribute (l_index_offP, g_ghosts)
      end associate
      deallocate(domain_offP_index, g_ghosts)
    end if

    !! Identify all unknown off-process index references (map>0).
    allocate(map(0:range%global_size_))
    map = 0
    do j = 1, size(l_index)
      map(l_index(j)) = l_index(j)
    end do
    map(range%first:range%last) = 0   ! on-process indices are known
    if (associated(range%offP_index)) then  ! known off-process indices
      do j = 1, size(range%offP_index)
        map(range%offP_index(j)) = 0
      end do
    end if

    !! Extract the list of unknown referenced indices.
    !allocate(offP_index(count(map>0)))
    offP_index = pack(map, mask=(map>0))

    !! Local numbering of the known off-process indices; by convention these
    !! were numbered sequentially following the on-process indices.
    if (associated(range%offP_index)) then
      do j = 1, size(range%offP_index)
        map(range%offP_index(j)) = range%onP_size_ + j
      end do
    end if

    !! Generate a local numbering of the unknown off-process indices; these
    !! will be numered sequentially following the 'known' local indices.
    do j = 1, size(offP_index)
      map(offP_index(j)) = range%local_size_ + j
    end do

    !! Remap the local index array values to the local numbering.
    do j = 1, size(l_index)
      if ((l_index(j) < range%first) .or. (l_index(j) > range%last)) then
        l_index(j) = map(l_index(j))
      else
        !! 1-based numbering of the on-process indices.
        l_index(j) = l_index(j) - range%first + 1
      end if
    end do

    deallocate(map)

  end subroutine localize_index_struct

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DEFINED
 !!
 !! Return the value true if the partition descriptor THIS is defined;
 !! otherwise return the value false.  Defined means that the variable has
 !! data and that the data is consistent.  This is a global procedure, and
 !! returns a global result.
 !!
 !! This procedure is mainly intended to be used in assertion checks.
 !!
 !! We do the best that we can to verify that the component data is consistent
 !! with the pgslib_gs_trace component, but have no way of checking whether
 !! that component is internally well-defined.
 !!

  logical function defined_ip_desc (this)

    class(ip_desc), intent(in) :: this

    integer :: sizes(nPE)

    call collate (sizes, this%onP_size_)
    call broadcast (sizes)

    defined_ip_desc = .false.
    CHECKLIST_1: do ! the cheap bits
      if (this%onP_size_ < 0) exit
      if (this%offP_size_ < 0) exit
      if (this%local_size_ /= this%onP_size_ + this%offP_size_) exit
      if (this%global_size_ /= sum(sizes)) exit
      if (this%last /= sum(sizes(:this_PE))) exit
      if (this%first /= this%last - this%onP_size_ + 1) exit
      if (associated(this%offP_index)) then
        if (this%offP_size_ /= size(this%offP_index)) exit
        if (.not.associated(this%trace)) exit
        if (.not.associated(this%dup_index)) exit
      else
        if (this%offP_size_ /= 0) exit
        if (associated(this%trace)) exit
        if (associated(this%dup_index)) exit
        if (associated(this%sup_index)) exit
      end if
      defined_ip_desc = .true.
      exit
    end do CHECKLIST_1
    defined_ip_desc = global_all(defined_ip_desc)

    if (.not.defined_ip_desc) return
    if (global_all(.not.associated(this%offP_index))) return  ! We're good; nothing else to check.

    defined_ip_desc = .false.
    !! It would be nice to check the trace component but PGSlib doesn't provide a function.
    !! if (.not.defined(this%trace)) return
    CHECKLIST_2: do ! the more expensive bits
      if (.not.associated(this%offP_index)) exit
      if (minval(this%offP_index) < 1 .or. maxval(this%offP_index) > this%global_size_) exit
      if (any(this%offP_index >= this%first .and. this%offP_index <= this%last)) exit
      if (size(this%dup_index) /= pgslib_size_of_dup(this%trace)) exit
      if (minval(this%dup_index) < 1 .or. maxval(this%dup_index) > this%onP_size_) exit
      if (.not.associated(this%sup_index)) then
        if (this%offP_size_ /= pgslib_size_of_sup(this%trace)) exit
      else
        if (size(this%sup_index) /= this%offP_size_) exit
        if (minval(this%sup_index) < 1 .or. maxval(this%sup_index) > this%offP_size_) exit
      end if
      defined_ip_desc = .true.
      exit
    end do CHECKLIST_2
    defined_ip_desc = global_all(defined_ip_desc)

  end function defined_ip_desc


  logical function gather_boundary_verified (this)

    type(ip_desc), intent(in) :: this

    integer :: j
    integer, allocatable :: onP_data(:), offP_data(:)

    !ASSERT( defined(this) )

    allocate(onP_data(this%onP_size_), offP_data(this%offP_size_))
    do j = 1, this%onP_size_
      onP_data(j) = global_index(this, j)
    end do
    call gather_boundary (this, onP_data, offP_data)
    gather_boundary_verified = global_all(offP_data == this%offP_index)
    deallocate(onP_data, offP_data)

  end function gather_boundary_verified

!  subroutine verify_boundary_scatter (this, verified)
!
!    type(ip_desc), intent(in) :: this
!    logical, intent(out) :: verified
!
!    integer :: j
!    integer, allocatable :: onP_data(:), offP_data(:)
!
!    ASSERT( defined(this) )
!
!    allocate(onP_data(this%onP_size_), offP_data(this%offP_size_))
!    do j = 1, this%onP_size_
!      onP_data(j) = this%first + j - 1
!    end do
!    call boundary_gather (this, onP_data, offP_data)
!    verified = global_all(offP_data == offP_data)
!    deallocate(onP_data, offP_data)
!
!  end subroutine verify_boundary_scatter

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! GATHER_BOUNDARY specific procedures.
 !! BOUNDARY_IS_CURRENT specific procedures.
 !!
 !! NOTES
 !!
 !! (1) The multidimensional versions are just do-loop-wrapped calls to the
 !! scalar versions.  This is ugly.  The core issue is that PGSLib gathers
 !! and scatters ought to be extended to handle general blocks of data per
 !! index, instead of a scalar.  This would be far more efficient.
 !!

#define _INTEGER_DATA_
#include "gather_boundary.fpp"

#define _SINGLE_DATA_
#include "gather_boundary.fpp"

#define _DOUBLE_DATA_
#include "gather_boundary.fpp"

#define _LOGICAL_DATA_
#include "gather_boundary.fpp"


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SCATTER_BOUNDARY_<op> specific procedures.
 !!

#define _SUM_REDUCTION_
#define _INTEGER_DATA_
#include "scatter_boundary.fpp"

#define _SUM_REDUCTION_
#define _SINGLE_DATA_
#include "scatter_boundary.fpp"

#define _SUM_REDUCTION_
#define _DOUBLE_DATA_
#include "scatter_boundary.fpp"

#define _MIN_REDUCTION_
#define _INTEGER_DATA_
#include "scatter_boundary.fpp"

#define _MIN_REDUCTION_
#define _SINGLE_DATA_
#include "scatter_boundary.fpp"

#define _MIN_REDUCTION_
#define _DOUBLE_DATA_
#include "scatter_boundary.fpp"

#define _MAX_REDUCTION_
#define _INTEGER_DATA_
#include "scatter_boundary.fpp"

#define _MAX_REDUCTION_
#define _SINGLE_DATA_
#include "scatter_boundary.fpp"

#define _MAX_REDUCTION_
#define _DOUBLE_DATA_
#include "scatter_boundary.fpp"

#define _AND_REDUCTION_
#include "scatter_boundary.fpp"

#define _OR_REDUCTION_
#include "scatter_boundary.fpp"

end module index_partitioning
