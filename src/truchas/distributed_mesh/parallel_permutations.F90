!!
!! PARALLEL_PERMUTATIONS
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Last revised 10 Apr 2007
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL CREATE_PAR_PERM (THIS, PERM, DEST_BSIZE, SRC_BSIZE) takes the global
!!    permutation vector PERM presented on the IO processor and creates the
!!    PAR_PERM object THIS that describes the parallel form of the permutation.
!!    The permutation is regarded as a mapping from destination index space
!!    with partition block sizes given by DEST_BSIZE into a source index space
!!    with partition block sizes given by SRC_BSIZE.  The sizes of the source
!!    and destination index spaces must of course be the same; it is their
!!    partitioning that may differ.  PERM, DEST_BSIZE and SRC_BSIZE are all
!!    intent-in.  The intent-out result THIS is used as input to REORDER to
!!    perform the reordering a source-distributed array into a destination-
!!    distributed array.  This is a global procedure.
!!
!!  CALL CREATE_PAR_PERM (P1, P2, P12 [, P21])
!!
!!    In this form of the call, the intent-in integer arrays P1 and P2 are
!!    distributed permutation vectors that map onto the same reference index
!!    set, and they can be viewed as describing two different labelings (of
!!    cells, for example).  This routine returns the PAR_PERM object P12 that
!!    describes the parallel form of the permutation from the P1-labeling to
!!    the P2-labeling.  The parallel from of the reverse permutation from
!!    P2-labeling to P1-labeling is also returned in P21 if it is present.
!!    This is a global procedure.
!!
!!  CALL REORDER (THIS, DEST, SRC) reorders the values of the distributed
!!    array SRC according to the parallel permutation THIS that describes
!!    the mapping from the DEST index space into the SRC index space,
!!    returning the result in the distributed array DEST.  This is a global
!!    procedure.
!!
!!  CALL DESTROY (THIS) deallocates all allocated storage associated with the
!!    PAR_PERM object THIS.  This is a global procedure.
!!
!!  DEFINED(THIS) returns the value true if the PAR_PERM object THIS is
!!    defined; otherwise it returns the value false.  Defined means that
!!    the components of the variable have been consistently defined.
!!    This is a global procedure.
!!
!! The above routines facilitate the 'complete' rearrangement of distributed
!! data from one numbering to another, hence the mappings are permutations.
!! A need has arisen to handle 'partial' rearrangements of data.  A case in
!! point is rearranging cell-based data from the old Truchas mesh which may
!! contain gap elements to the new DS mesh which has dropped these elements.
!! In this case there is no complete 1-1 correspondence between the two
!! numberings, but we can still perform the rearranngement of data where
!! there is a correspondence.  The following routines provide this capability.
!!
!!  CALL CREATE_PAR_PERM (P1, P2, P12, LIST1 [, P21, LIST2])
!!
!!    In this variation of call, the intent-in integer arrays P1 and P2 are
!!    the distributed mapping vectors (not necessarily permutations) into the
!!    same reference index set.  The mappings are required to be 1-1 in
!!    keeping with the notion that these represent renumberings.  LIST1 and
!!    LIST2 are integer array pointers that will be allocated by the routine
!!    and assigned the indices for which the 1-to-2 mapping and 2-to-1 mapping,
!!    respectively are not defined (see REARRANGE below).  The PAR_PERM objects
!!    P12 and P21 returned by the call should not be passed to REORDER -- use
!!    REARRANGE instead.
!!
!!  CALL REARRANGE (THIS, DEST, SRC [,DEFAULT]) rearranges values from the
!!    distributed array SRC to the distributed array DEST according to the
!!    parallel permutation THIS that describes the mapping from the DEST index
!!    space into the SRC index space.  This is analogous to REORDER (and in fact
!!    equivalent if THIS was created with either of the first two forms of
!!    CREATE_PAR_PERM) except that DEST is intent-inout and not all elements
!!    of DEST may be assigned a value by the call.  The elements of DEST not
!!    assigned to will be those in the index list corresponding to THIS that
!!    was returned by the third form of CREATE_PAR_PERM.  If the optional
!!    scalar value DEFAULT is specified, those elements in DEST not assigned
!!    a value from SRC will be assigned this value.
!!

#include "f90_assert.fpp"

module parallel_permutations

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use permutations
  use index_map_type
  use parallel_communication
  implicit none
  private

  public :: create_par_perm, defined, reorder, rearrange

  type, public :: par_perm
    private
    integer, allocatable :: perm(:) ! dest-to-src mapping (local indices)
    type(index_map) :: src_imap ! partition descriptor for the src index set
  end type par_perm

  interface defined
    module procedure defined_par_perm
  end interface

  interface reorder
    module procedure reorder_r0, reorder_i0
  end interface

  interface create_par_perm
    module procedure create_par_perm_1, create_par_perm_2, create_par_perm_3
  end interface

  interface rearrange
    module procedure rearrange_r0, rearrange_r1, rearrange_i0, rearrange_i1
    module procedure rearrange_l0, rearrange_l1
  end interface

contains

  subroutine create_par_perm_1 (this, perm, dest_bsize, src_bsize)

    type(par_perm), intent(out) :: this
    integer, intent(in) :: perm(:)  ! global dest-to-src permutation mapping
    integer, intent(in) :: dest_bsize(:)  ! dest index partition block sizes
    integer, intent(in) :: src_bsize(:)   ! src index partition block sizes

    type(index_map) :: dest_imap

    if (is_IOP) then
      ASSERT( size(dest_bsize) == nPE )
      ASSERT( all(dest_bsize >= 0) )
      ASSERT( size(src_bsize) == nPE )
      ASSERT( all(src_bsize >= 0) )
      ASSERT( sum(src_bsize) == sum(dest_bsize) )
      ASSERT( size(perm) == sum(dest_bsize) )
      ASSERT( is_perm(perm) )
    end if

    call dest_imap%init (dest_bsize)
    call this%src_imap%init (src_bsize)

    call dest_imap%localize_index_array (perm, this%src_imap, this%perm)

    ASSERT( minval(this%perm) > 0 )
    ASSERT( maxval(this%perm) <= this%src_imap%local_size )

  end subroutine create_par_perm_1

  subroutine create_par_perm_2 (dp1, dp2, pp12, pp21)

    integer, intent(in) :: dp1(:), dp2(:)
    type(par_perm), intent(out) :: pp12
    type(par_perm), intent(out), optional :: pp21

    integer :: n, bsize1(nPE), bsize2(nPE)
    integer, pointer :: p1(:) => null(), p2(:) => null()

    n = global_sum(size(dp1))
    ASSERT( n == global_sum(size(dp2)) )

    !! Assemble the 1-labeling-to-external permutation vector.
    allocate(p1(merge(n,0,is_iop)))
    call gather (dp1, p1)

    !! Assemble the 2-labeling-to-external permutation vector.
    allocate(p2(merge(n,0,is_iop)))
    call gather (dp2, p2)

    !! Construct the permutation from the 1-to-2 labeling.
    if (is_IOP) then
      ASSERT( is_perm(p1) )
      ASSERT( is_perm(p2) )
      call invert_perm (p2)
      call reorder (p2, perm=p1)  ! pull-back the inverse of p2 by p1
      ASSERT( is_perm(p2) )
    end if

    !! Form the partition block-size arrays inferred from the two labelings.
    call gather (size(dp1), bsize1)
    call gather (size(dp2), bsize2)

    !! Create the 1-to-2 parallel permutation structure.
    call create_par_perm (pp12, p2, bsize1, bsize2)

    !! Create the 2-to-1 parallel permutation structure.
    if (present(pp21)) then
      if (is_IOP) call invert_perm (p2, invp=p1)
      call create_par_perm (pp21, p1, bsize2, bsize1)
    end if

    !! Clean-up.
    deallocate(p1, p2)

  end subroutine create_par_perm_2

  subroutine create_par_perm_3 (dm1, dm2, pp12, list1, pp21, list2)

    integer, intent(in) :: dm1(:), dm2(:)
    type(par_perm), intent(out) :: pp12
    integer, pointer :: list1(:)
    type(par_perm), intent(out), optional :: pp21
    integer, pointer, optional :: list2(:)

    integer :: j, n, bsize1(nPE), bsize2(nPE), nmin, nmax
    integer, pointer :: m1(:) => null(), m2(:) => null()
    integer, allocatable :: m2inv(:)

    !! Assemble the 1-labeling-to-external mapping vector M1.
    n = global_sum(size(dm1))
    allocate(m1(merge(n,0,is_iop)))
    call gather (dm1, m1)

    !! Assemble the 2-labeling-to-external mapping vector M2.
    n = global_sum(size(dm2))
    allocate(m2(merge(n,0,is_iop)))
    call gather (dm2, m2)

    !! Overwrite M1 with the 1-to-2-labeling mapping.
    !! Zero values in M1 mark 1-labels without a corresponding 2-label.
    if (is_IOP) then
      nmin = min(minval(m1), minval(m2))
      nmax = max(maxval(m1), maxval(m2))
      ASSERT(nmin > 0)
      allocate(m2inv(nmin:nmax))
      ASSERT(is_one_to_one(m1).or. .not.present(pp21))
      ASSERT(is_one_to_one(m2))
      m2inv = 0
      do j = 1, size(m2)
        m2inv(m2(j)) = j
      end do
      do j = 1, size(m1)
        m1(j) = m2inv(m1(j))
      end do
      deallocate(m2inv)
    end if

    !! Form the partition block-size arrays inferred from the two labelings.
    call gather (size(dm1), bsize1)
    call gather (size(dm2), bsize2)

    !! Create the 1-to-2 parallel permutation structure.
    call create_par_perm_var (pp12, m1, bsize1, bsize2)
    call generate_undefined_map_list (pp12%perm, list1)

    !! Create the 2-to-1 parallel permutation structure.
    if (present(pp21)) then
      INSIST(present(list2))
      if (is_IOP) then
        !! Overwrite M2 with the inverse of M1.
        !! Zero values in M2 mark 2-labels without a corresponding 1-label.
        m2 = 0
        do j = 1, size(m1)
          if (m1(j) /= 0) m2(m1(j)) = j
        end do
      end if
      call create_par_perm_var (pp21, m2, bsize2, bsize1)
      call generate_undefined_map_list (pp21%perm, list2)
    end if

    !! Clean-up.
    deallocate(m1, m2)

  contains

    subroutine generate_undefined_map_list (map, list)
      integer, intent(in) :: map(:)
      integer, pointer :: list(:)
      integer :: j, n
      n = count(map == 0)
      allocate(list(n))
      if (n == 0) return
      n = 0
      do j = 1, size(map)
        if (map(j) == 0) then
          n = n + 1
          list(n) = j
        end if
      end do
    end subroutine

    !! We need to encapsulate this bit of code in a function
    !! in order to use its result in an assertion statement.

    logical function is_one_to_one (m)
      integer, intent(in) :: m(:)
      integer :: j
      is_one_to_one = .false.
      m2inv = 0
      do j = 1, size(m)
        if (m2inv(m(j)) /= 0) return
        m2inv(m(j)) = j
      end do
      is_one_to_one = .true.
    end function

  end subroutine create_par_perm_3

  subroutine create_par_perm_var (this, perm, dest_bsize, src_bsize)

    type(par_perm), intent(out) :: this
    integer, intent(in) :: perm(:)  ! global dest-to-src mapping
    integer, intent(in) :: dest_bsize(:)  ! dest index partition block sizes
    integer, intent(in) :: src_bsize(:)   ! src index partition block sizes

    type(index_map) :: dest_imap

    call dest_imap%init (dest_bsize)
    call this%src_imap%init (src_bsize)

    call dest_imap%localize_index_array (perm, this%src_imap, this%perm)

    ASSERT( minval(this%perm) >= 0 )
    ASSERT( maxval(this%perm) <= this%src_imap%local_size )

  end subroutine create_par_perm_var

  logical function defined_par_perm (this)
    type(par_perm), intent(in) :: this
    defined_par_perm = .false.
    CHECKLIST: block
      if (.not.allocated(this%perm)) exit CHECKLIST
!      if (.not.this%src_imap%defined()) exit CHECKLIST
      defined_par_perm = .true.
    end block CHECKLIST
    defined_par_perm = global_all(defined_par_perm)
  end function defined_par_perm

  subroutine reorder_r0 (this, dest, src)

    type(par_perm), intent(in) :: this
    real(r8), intent(out) :: dest(:)
    real(r8), intent(in)  :: src(:)

    integer :: j, k
    real(r8) :: src_offP(this%src_imap%onp_size+1:this%src_imap%local_size)

    ASSERT( defined(this) )
    ASSERT( size(dest) == size(this%perm) )
    ASSERT( size(src) == this%src_imap%onp_size )

    call this%src_imap%gather_offp(src, src_offP)

    do j = 1, size(dest)
      k = this%perm(j)
      if (k > this%src_imap%onp_size) then
        dest(j) = src_offP(k)
      else
        dest(j) = src(k)
      end if
    end do

  end subroutine reorder_r0

  subroutine reorder_i0 (this, dest, src)

    type(par_perm), intent(in) :: this
    integer, intent(out) :: dest(:)
    integer, intent(in)  :: src(:)

    integer :: j, k
    integer :: src_offP(this%src_imap%onp_size+1:this%src_imap%local_size)

    ASSERT( defined(this) )
    ASSERT( size(dest) == size(this%perm) )
    ASSERT( size(src) == this%src_imap%onp_size )

    call this%src_imap%gather_offp(src, src_offP)

    do j = 1, size(dest)
      k = this%perm(j)
      if (k > this%src_imap%onp_size) then
        dest(j) = src_offP(k)
      else
        dest(j) = src(k)
      end if
    end do

  end subroutine reorder_i0

  subroutine rearrange_r0 (this, dest, src, default)

    type(par_perm), intent(in) :: this
    real(r8), intent(inout) :: dest(:)
    real(r8), intent(in) :: src(:)
    real(r8), intent(in), optional :: default

    integer :: j, k
    real(r8) :: src_offP(this%src_imap%onp_size+1:this%src_imap%local_size)

    ASSERT( defined(this) )
    ASSERT( size(dest) == size(this%perm) )
    ASSERT( size(src) == this%src_imap%onp_size )

    call this%src_imap%gather_offp(src, src_offP)

    do j = 1, size(dest)
      k = this%perm(j)
      if (k > this%src_imap%onp_size) then
        dest(j) = src_offP(k)
      else if (k > 0) then
        dest(j) = src(k)
      else if (present(default)) then
        dest(j) = default
      end if
    end do

  end subroutine rearrange_r0

  subroutine rearrange_r1 (this, dest, src, default)

    type(par_perm), intent(in) :: this
    real(r8), intent(inout) :: dest(:,:)
    real(r8), intent(in) :: src(:,:)
    real(r8), intent(in), optional :: default

    integer :: j, k
    real(r8) :: src_offP(size(src,1),this%src_imap%onp_size+1:this%src_imap%local_size)

    ASSERT( defined(this) )
    ASSERT( size(dest,1) == size(src,1) )
    ASSERT( size(dest,2) == size(this%perm) )
    ASSERT( size(src,2) == this%src_imap%onp_size )

    call this%src_imap%gather_offp(src, src_offP)

    do j = 1, size(dest,2)
      k = this%perm(j)
      if (k > this%src_imap%onp_size) then
        dest(:,j) = src_offP(:,k)
      else if (k > 0) then
        dest(:,j) = src(:,k)
      else if (present(default)) then
        dest(:,j) = default
      end if
    end do

  end subroutine rearrange_r1

  subroutine rearrange_i0 (this, dest, src, default)

    type(par_perm), intent(in) :: this
    integer, intent(inout) :: dest(:)
    integer, intent(in) :: src(:)
    integer, intent(in), optional :: default

    integer :: j, k
    integer :: src_offP(this%src_imap%onp_size+1:this%src_imap%local_size)

    ASSERT( defined(this) )
    ASSERT( size(dest) == size(this%perm) )
    ASSERT( size(src) == this%src_imap%onp_size )

    call this%src_imap%gather_offp(src, src_offP)

    do j = 1, size(dest)
      k = this%perm(j)
      if (k > this%src_imap%onp_size) then
        dest(j) = src_offP(k)
      else if (k > 0) then
        dest(j) = src(k)
      else if (present(default)) then
        dest(j) = default
      end if
    end do

  end subroutine rearrange_i0

  subroutine rearrange_i1 (this, dest, src, default)

    type(par_perm), intent(in) :: this
    integer, intent(inout) :: dest(:,:)
    integer, intent(in) :: src(:,:)
    integer, intent(in), optional :: default

    integer :: j, k
    integer :: src_offP(size(src,1),this%src_imap%onp_size+1:this%src_imap%local_size)

    ASSERT( defined(this) )
    ASSERT( size(dest,1) == size(src,1) )
    ASSERT( size(dest,2) == size(this%perm) )
    ASSERT( size(src,2) == this%src_imap%onp_size )

    call this%src_imap%gather_offp(src, src_offP)

    do j = 1, size(dest,2)
      k = this%perm(j)
      if (k > this%src_imap%onp_size) then
        dest(:,j) = src_offP(:,k)
      else if (k > 0) then
        dest(:,j) = src(:,k)
      else if (present(default)) then
        dest(:,j) = default
      end if
    end do

  end subroutine rearrange_i1

  subroutine rearrange_l0 (this, dest, src, default)

    type(par_perm), intent(in) :: this
    logical, intent(inout) :: dest(:)
    logical, intent(in) :: src(:)
    logical, intent(in), optional :: default

    integer :: j, k
    logical :: src_offP(this%src_imap%onp_size+1:this%src_imap%local_size)

    ASSERT( defined(this) )
    ASSERT( size(dest) == size(this%perm) )
    ASSERT( size(src) == this%src_imap%onp_size )

    call this%src_imap%gather_offp(src, src_offP)

    do j = 1, size(dest)
      k = this%perm(j)
      if (k > this%src_imap%onp_size) then
        dest(j) = src_offP(k)
      else if (k > 0) then
        dest(j) = src(k)
      else if (present(default)) then
        dest(j) = default
      end if
    end do

  end subroutine rearrange_l0

  subroutine rearrange_l1 (this, dest, src, default)

    type(par_perm), intent(in) :: this
    logical, intent(inout) :: dest(:,:)
    logical, intent(in) :: src(:,:)
    logical, intent(in), optional :: default

    integer :: j, k
    logical :: src_offP(size(src,1),this%src_imap%onp_size+1:this%src_imap%local_size)

    ASSERT( defined(this) )
    ASSERT( size(dest,1) == size(src,1) )
    ASSERT( size(dest,2) == size(this%perm) )
    ASSERT( size(src,2) == this%src_imap%onp_size )

    call this%src_imap%gather_offp(src, src_offP)

    do j = 1, size(dest,2)
      k = this%perm(j)
      if (k > this%src_imap%onp_size) then
        dest(:,j) = src_offP(:,k)
      else if (k > 0) then
        dest(:,j) = src(:,k)
      else if (present(default)) then
        dest(:,j) = default
      end if
    end do

  end subroutine rearrange_l1

end module parallel_permutations
