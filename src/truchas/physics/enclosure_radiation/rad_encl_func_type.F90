!!
!! RAD_ENCL_FUNC_TYPE
!!
!! Neil N. Carlson <nnc@lanl.gov>, 5 Apr 2010
!!
!! This module provides a derived type for describing scalar, face-based
!! surface data that depends on time.  An instance of this type should be
!! viewed as the discretization of a time-space function over the faces
!! of a distributed enclosure surface mesh.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! The derived type ENCL_FUNC is an opaque structure with private
!! components.  No defined assignment is provided; do not use instances of
!! this type in an assignment statement unless you really know what default
!! assignment does (it's almost certainly not what you want in this case).
!! The enclosure surface data is defined by an initial call to EF_PREP
!! followed by one or more calls to EF_ADD_FUNCTION and a final call to
!! EF_DONE.  This definition phase is mandatory.  At this point the object
!! may be evaluated using EF_EVAL.
!!
!!  CALL EF_PREP (THIS, ENCL) prepares the ENCL_FUNC object THIS to begin
!!    receiving the data functions to be associated with groups of enclosure
!!    surface faces.  ENCL is a pointer to the distributed enclosure mesh
!!    object over which the boundary condition applies.  The implementation
!!    assumes that the enclosure surface is never modified.
!!
!!  CALL EF_ADD_FUNCTION (THIS, F, BLOCKID, STAT, ERRMSG) associates the
!!    SCAFUN function object F to the group of surface faces that belong to
!!    a block with an ID from the list given by rank-1 integer array BLOCKID.
!!    The possible blocks are those defined in ENCL.  The function object
!!    should expect (/ t, x, y, z /) as the array argument.  An independent
!!    copy of F is stored so that the actual argument may be destroyed
!!    afterwards if desired.  A nonzero value is assigned to the integer STAT
!!    if an error condition occurs, and an explanatory message is assigned to
!!    the character string ERRMSG.  Possible errors are referring to an invalid
!!    block ID or attempting to associate a face to more than one function.
!!    This routine may be called multiple times as long as the specified face
!!    blocks do not overlap from one call to the next.
!!
!!  CALL EF_DONE (THIS, STAT, ERRMSG) performs the final configuration of the
!!    ENCL_FUNC object THIS after all the desired calls to EF_ADD_FUNCTION
!!    have been made.  A nonzero value is assigned to the integer STAT if an
!!    error condition occurs, and an explanatory error message is assigned to
!!    the character string ERRMSG.  It is an error for a surface face to not be
!!    associated with a function.  After this call the object THIS may be
!!    evaluated using EF_EVAL.
!!
!!  CALL EF_EVAL (THIS, T, VALUE) evaluates the ENCL_FUNC object THIS at
!!    the time T and returns the computed data values in the face-based array
!!    VALUE.  A simple 1-point (face centroid) quadrature is currently used for
!!    functions that are spatially dependent.
!!
!!  CALL EF_DESTROY (THIS) frees the resources belonging to the ENCL_FUNC
!!    object THIS, returning it to its default initialization state.
!!
!! IMPLEMENTATION NOTE
!!
!! The evaluation procedure is written to short-circuit some potentially costly
!! calculations as directed by some optimization "hinting".  However, this
!! hinting is not currently exposed in the API for the user to control.  But
!! since we can directly interrogate SCAFUN objects to see if they are of the
!! constant flavor, we do set the optimization hinting for those functions at
!! least.
!!

#include "f90_assert.fpp"

module rad_encl_func_type

  use kinds
  use scalar_func_containers
  use rad_encl_type
  implicit none
  private

  type, public :: rad_encl_func
    private
    real(r8), allocatable, public :: values(:)
    !! the rest are private
    type(rad_encl), pointer :: encl => null() ! reference only -- do not own
    real(r8) :: tlast = -huge(1.0_r8)
    integer :: ngroup = -1
    integer, allocatable :: group(:)
    integer, allocatable :: faces(:)
    type(scalar_func_box), allocatable :: farray(:)
    integer, allocatable :: hint(:)
    !! temporary components used during initialization
    integer, allocatable :: tag(:)
    type(scalar_func_list) :: flist
  contains
    procedure :: prep
    procedure :: add_function
    procedure :: done
    procedure :: eval
  end type rad_encl_func

  !! Optimization hint values.
  integer, parameter, public :: EF_HINT_NONE    = 0
  integer, parameter, public :: EF_HINT_CONST   = 1
  integer, parameter, public :: EF_HINT_T_INDEP = 2
  integer, parameter, public :: EF_HINT_X_INDEP = 3

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! EF_PREP
 !!

  subroutine prep (this, encl)

    class(rad_encl_func), intent(out) :: this
    type(rad_encl), target, intent(in) :: encl

    this%encl => encl
    this%ngroup = 0
    allocate(this%tag(encl%nface))
    this%tag = 0

  end subroutine prep

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! EF_ADD_FUNCTION
 !!

  subroutine add_function (this, f, blockID, stat, errmsg)

    class(rad_encl_func), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: f
    integer, intent(in) :: blockID(:)
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    !! Verify that THIS is in the correct state.
    INSIST(this%ngroup>=0 .and. allocated(this%tag))

    !! Tag the faces associated with this function.
    call set_tag_array (this, blockID, stat, errmsg)
    if (stat /= 0) return

    !! Append the function to the temporary list.
    call this%flist%append (f)

  end subroutine add_function

  !!
  !! This auxillary subroutine tags the faces identified by the list
  !! of face block IDs, checking for error conditions in the process.
  !!

  subroutine set_tag_array (this, blockID, stat, errmsg)

    type(rad_encl_func), intent(inout) :: this
    integer, intent(in) :: blockID(:)
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: i, j
    logical :: mask(this%encl%nface), block_mask(size(this%encl%face_block_id))

    !! Create the mask corresponding to blockID.
    block_mask = .false.
    do i = 1, size(blockID)
      do j = size(this%encl%face_block_id), 1, -1
        if (blockID(i) == this%encl%face_block_id(j)) exit
      end do
      if (j == 0) then
        stat = 1
        write(errmsg,'(a,i0)') 'unknown enclosure face block ID: ', blockID(i)
        return
      end if
      block_mask(j) = .true.
    end do

    !! Identify the faces specified by blockID.
    mask = block_mask(this%encl%face_block)

    !! Check that these faces don't overlap those from preceding calls.
    if (any(mask .and. this%tag /= 0)) then
      stat = 1
      errmsg = 'enclosure face already associated with a function'
      return
    end if

    !! Set the tag array.
    this%ngroup = 1 + this%ngroup
    this%tag = merge(this%ngroup, this%tag, mask)

    stat = 0
    errmsg = ''

  end subroutine set_tag_array

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! EF_DONE
 !!

  subroutine done (this, stat, errmsg)
  
    use const_scalar_func_type

    class(rad_encl_func), intent(inout) :: this
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: n, j

    !! Verify that THIS is in the correct state.
    INSIST(this%ngroup>=0 .and. allocated(this%tag))

    ASSERT(minval(this%tag)>=0 .and. maxval(this%tag)<=this%ngroup)

    !! Verify that every face has been associated with a function.
    n = count(this%tag > 0)
    if (n /= this%encl%nface) then
      stat = 1
      errmsg = 'some enclosure faces not associated with a function'
      return
    end if

    allocate(this%faces(n), this%group(this%ngroup+1))

    !! Prepare GROUP: faces of group N will be FACES(GROUP(N):GROUP(N+1)-1).
    this%group(1) = 1
    do n = 1, this%ngroup
      this%group(n+1) = this%group(n) + count(this%tag == n)
    end do

    !! Fill the FACES array; GROUP(N) stores the next free location for group N.
    do j = 1, size(this%tag)
      n = this%tag(j)
      if (n == 0) cycle
      this%faces(this%group(n)) = j
      this%group(n) = 1 + this%group(n)
    end do
    deallocate(this%tag)

    !! Restore GROUP; GROUP(N) is now the start of condition N+1 instead of N
    do n = this%ngroup, 1, -1
      this%group(n+1) = this%group(n)
    end do
    this%group(1) = 1

    !! Convert the function list into the final function array.
    call scalar_func_list_to_box_array (this%flist, this%farray)

    !! For now we don't expose optimization hinting to the user,
    !! but we can determine directly which functions are constant.
    allocate(this%hint(this%ngroup))
    do n = 1, this%ngroup
      select type (f => this%farray(n)%f)
      type is (const_scalar_func)
        this%hint(n) = EF_HINT_CONST
      class default
        this%hint(n) = EF_HINT_NONE
      end select
    end do

    stat = 0
    errmsg = ''

  end subroutine done

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! EF_EVAL
 !!

  subroutine eval (this, t)

    class(rad_encl_func), intent(inout) :: this
    real(r8), intent(in)  :: t

    integer :: j, n
    logical :: unevaluated
    real(r8) :: args(0:size(this%encl%coord,dim=1))

    !! Verify that THIS is in the correct state.
    INSIST(allocated(this%faces) .and. .not.allocated(this%tag))

    unevaluated = .not.allocated(this%values)
    if (unevaluated) allocate(this%values(this%encl%nface))

    if (unevaluated .or. t /= this%tlast) then
      args(0) = t
      do n = 1, this%ngroup
        associate (faces => this%faces(this%group(n):this%group(n+1)-1))
          select case (this%hint(n))
          case (EF_HINT_CONST)
            if (unevaluated) this%values(faces) = this%farray(n)%f%eval(args)
          case (EF_HINT_X_INDEP)
            this%values(faces) = this%farray(n)%f%eval(args)
          case (EF_HINT_T_INDEP)
            if (unevaluated) then
              do j = 1, size(faces)
                associate (fnode => this%encl%fnode(this%encl%xface(faces(j)):this%encl%xface(faces(j)+1)-1))
                  args(1:) = sum(this%encl%coord(:,fnode),dim=2) / size(fnode)
                  this%values(faces(j)) = this%farray(n)%f%eval(args)
                end associate
              end do
            end if
          case default
            do j = 1, size(faces)
              associate (fnode => this%encl%fnode(this%encl%xface(faces(j)):this%encl%xface(faces(j)+1)-1))
                args(1:) = sum(this%encl%coord(:,fnode),dim=2) / size(fnode)
                this%values(faces(j)) = this%farray(n)%f%eval(args)
              end associate
            end do
          end select
        end associate
      end do
      this%tlast = t
    end if

  end subroutine eval

end module rad_encl_func_type
