!!
!! INTERFACE_DATA
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Last revised 22 June 2009
!!
!! This module provides a derived type for describing scalar, face-based
!! interface data that depends on time.  An instance of this type should be
!! viewed as the discretization of a time-space function over the faces
!! on a portion of the mesh interfaces.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! An instance of the derived type IF_DATA is intended to encapsulate the data
!! associated with a particular interface condition.  Once initialized, an
!! instance of this type has two public array components, FACES and VALUES:
!! for each j, %FACES(1:2,j) are the indices of an opposing pair of faces on
!! the interface and %VALUES(:,j) are the corresponding condition parameter
!! values.  FACES has shape [2,n] and VALUES has shape [p,n], where n is the
!! number of interface face pairs associated with the instance and p is the
!! number of parameters associated with the condition.  Each type of condition
!! will have 0 or more parameters associated with it.  For example, a heat
!! transfer condition, in which the heat flux is proportional to the temperature
!! difference across the interface, would have the proportionality constant as
!! the only parameter.  Other conditions may have more or less parameters.
!! The parameter values are assumed to be functions of time and space.
!! Note that the interpretation of the parameter values is solely determined by
!! the code that uses them.  All other components of an IF_DATA object are
!! intended to be private.  No defined assignment is provided; do not use
!! instances of this type in an assignment statement unless you really
!! understand what the default assignment does.
!! 
!!  CALL IF_DATA_PREP (THIS, MESH, NPAR) prepares the IF_DATA object THIS to
!!    begin receiving the condition data.  MESH is a pointer to the distributed
!!    mesh object over which the interface condition applies.  The object THIS
!!    maintains an internal pointer to MESH (read-only).  The number of
!!    condition parameters is specified by NPAR.
!!
!!  CALL IF_DATA_ADD (THIS, F, SETID, STAT, ERRMSG) associates the rank-1 array
!!    of parameter functions F with the face pairs identified by the link set
!!    IDs given in the rank-1 integer array SETID.  The possible link sets
!!    are those defined in MESH.  The size of the array F must equal the number
!!    of condition parameters specified in the call to IF_DATA_PREP.  The
!!    functions are SCAFUN-type objects that should expect (/ t, x, y, z /) as
!!    the vector argument.  This subroutine must be called after IF_DATA_PREP,
!!    and can be called multiple times, or even not at all.  The function
!!    object F can be destroyed afterwards; THIS maintains an independent copy
!!    of that object.  A nonzero value is assigned to the integer STAT if an
!!    error condition occurs, and an explanatory message is assigned to the
!!    string ERRMSG. It is an error to refer to an unknown link set ID or
!!    attempt to associate a face pair with more than one function array.
!!
!!  CALL IF_DATA_DONE (THIS) performs the final configuration of the IF_DATA
!!    object THIS after all the desired calls to IF_DATA_ADD have been made
!!    (if any).  After this call,  THIS%FACES is the list of face pairs that
!!    were associated with some parameter functions in one of the calls to
!!    IF_DATA_ADD.  The array THIS%VALUES will hold the corresponding parameter
!!    values and are set by the subroutine IF_DATA_EVAL.
!!
!!  CALL IF_DATA_EVAL (THIS, T) evaluates the VALUES component of the IF_DATA
!!    object THIS at the time T.  THIS%VALUES(:,j) is assigned the result of
!!    evaluating the function vector F(:) associated with the pair of faces
!!    THIS%FACES(:,j).  In order to reduce the value of a function over a face
!!    to a single value, the function is evaluated at the average of the
!!    midpoints of the pair of faces where necessary.
!!
!!  CALL IF_DATA_DESTROY (THIS) frees the resources belonging to the IF_DATA
!!    object THIS, returning it to its default initialization state.
!!

#include "f90_assert.fpp"

module interface_data

  use kinds
  use scalar_func_containers
  use unstr_mesh_type
  use string_utilities, only: i_to_c
  implicit none
  private

  public :: if_data_prep, if_data_add, if_data_done, if_data_eval, if_data_destroy

  type, public :: if_data
    integer,  pointer :: faces(:,:)  => null()
    real(r8), pointer :: values(:,:) => null()
    !! rest are private
    integer :: npar
    type(unstr_mesh), pointer :: mesh => null()
    logical :: evaluated = .false.
    real(r8) :: tlast = -huge(1.0_r8)
    integer :: ngroup = -1
    integer, pointer :: xgroup(:) => null()
    type(scalar_func_box), allocatable :: farray(:,:)
    integer, pointer :: hint(:,:) => null()
    integer, pointer :: tag(:) => null()
    type(scalar_func_vlist) :: flist
  end type if_data

  !! Optimization hint values.
  integer, parameter, public :: IF_DATA_HINT_NONE    = 0
  integer, parameter, public :: IF_DATA_HINT_CONST   = 1
  integer, parameter, public :: IF_DATA_HINT_T_INDEP = 2
  integer, parameter, public :: IF_DATA_HINT_X_INDEP = 3

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! IF_DATA_PREP
 !!

  subroutine if_data_prep (this, mesh, npar)

    type(if_data), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    integer, intent(in) :: npar

    INSIST(npar >= 0)

    this%npar = npar
    this%mesh => mesh

    this%ngroup = 0
    allocate(this%tag(size(mesh%link_set_mask)))
    this%tag = 0

  end subroutine if_data_prep

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! IF_DATA_ADD
 !!

  subroutine if_data_add (this, f, setID, stat, errmsg)

    type(if_data), intent(inout) :: this
    type(scalar_func_box), allocatable, intent(inout) :: f(:)
    integer, intent(in) :: setID(:)
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    !! Verify that THIS is in the correct state.
    INSIST(this%ngroup>=0 .and. associated(this%tag))

    INSIST(size(f) == this%npar)

    call set_tag_array (this, setID, stat, errmsg)
    if (stat /= 0) return

    !! Append the function array F to the temporary list.
    call this%flist%append (f)

  end subroutine if_data_add

  !!
  !! This auxillary subroutine tags the faces identified by the list of
  !! face set IDs, checking for error conditions in the process.
  !!

  subroutine set_tag_array (this, setID, stat, errmsg)

    use bitfield_type
  
    type(if_data), intent(inout) :: this
    integer, intent(in) :: setID(:)
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: i, j
    logical :: mask(this%mesh%nlink)
    type(bitfield) :: bitmask

    !! Create the bitmask corresponding to setID.
    bitmask = ZERO_BITFIELD
    do i = 1, size(setID)
      do j = size(this%mesh%link_set_ID), 1, -1
        if (setID(i) == this%mesh%link_set_ID(j)) exit
      end do
      if (j == 0) then
        stat = 1
        errmsg = 'unknown link set ID: ' // i_to_c(setID(i))
        return
      end if
      bitmask = ibset(bitmask, j)
    end do

    !! Identify the link specified by setID.
    mask = (iand(bitmask, this%mesh%link_set_mask) /= ZERO_BITFIELD)

    !! Check that these links don't overlap those from preceding calls.
    if (any(mask .and. this%tag /= 0)) then
      stat = 1
      errmsg = 'link already associated with a function'
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
 !! IF_DATA_DONE
 !!

  subroutine if_data_done (this)
  
    use const_scalar_func_type

    type(if_data), intent(inout) :: this

    integer :: n, j

    !! Verify that THIS is in the correct state.
    INSIST(this%ngroup>=0 .and. associated(this%tag))

    ASSERT(minval(this%tag)>=0 .and. maxval(this%tag)<=this%ngroup)

    n = count(this%tag > 0)
    allocate(this%faces(2,n), this%values(this%npar,n), this%xgroup(this%ngroup+1))

    !! Prepare XGROUP: faces of group N will be FACES(XGROUP(N):XGROUP(N+1)-1).
    this%xgroup(1) = 1
    do n = 1, this%ngroup
      this%xgroup(n+1) = this%xgroup(n) + count(this%tag == n)
    end do

    !! Fill the FACES array; XGROUP(N) stores the next free location for condition N.
    do j = 1, size(this%tag)
      n = this%tag(j)
      if (n == 0) cycle
      this%faces(:,this%xgroup(n)) = this%mesh%lface(:,j)
      this%xgroup(n) = 1 + this%xgroup(n)
    end do

    !! Restore XGROUP; XGROUP(N) is now the start of condition N+1 instead of N
    do n = this%ngroup, 1, -1
      this%xgroup(n+1) = this%xgroup(n)
    end do
    this%xgroup(1) = 1

    !! Convert the function list into the final function array.
    call scalar_func_vlist_to_box_array (this%flist, this%farray)

    !! For now we don't expose optimization hinting to the user,
    !! but we can determine directly which functions are constant.
    allocate(this%hint(this%npar,this%ngroup))
    do n = 1, this%ngroup
      do j = 1, this%npar
        select type (f => this%farray(j,n)%f)
        type is (const_scalar_func)
          this%hint(j,n) = IF_DATA_HINT_CONST
        class default
          this%hint(j,n) = IF_DATA_HINT_NONE
        end select
      end do
    end do
    
    deallocate(this%tag)

  end subroutine if_data_done

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! IF_DATA_EVAL
 !!

  subroutine if_data_eval (this, t)
  
    type(if_data), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: n, i, j
    integer, pointer :: faces(:,:)
    real(r8), pointer :: values(:,:)
    real(r8) :: args(0:size(this%mesh%x,dim=1))

    !! Verify that THIS is in the correct state.
    ASSERT(associated(this%faces) .and. .not.associated(this%tag))

    if (this%evaluated .and. t == this%tlast) return  ! values already set for this T

    args(0) = t
    do n = 1, this%ngroup
      faces => this%faces(:,this%xgroup(n):this%xgroup(n+1)-1)
      values => this%values(:,this%xgroup(n):this%xgroup(n+1)-1)
      do i = 1, this%npar
        select case (this%hint(i,n))
        case (IF_DATA_HINT_CONST)
          if (.not.this%evaluated) values(i,:) = this%farray(i,n)%f%eval(args)
        case (IF_DATA_HINT_X_INDEP)
          values(i,:) = this%farray(i,n)%f%eval(args)
        case (IF_DATA_HINT_T_INDEP)
          if (.not.this%evaluated) then
            do j = 1, size(faces,dim=2)
              associate (fnode1 => this%mesh%fnode(this%mesh%xfnode(faces(1,j)):this%mesh%xfnode(faces(1,j)+1)-1), &
                         fnode2 => this%mesh%fnode(this%mesh%xfnode(faces(2,j)):this%mesh%xfnode(faces(2,j)+1)-1))
                args(1:) = (sum(this%mesh%x(:,fnode1),dim=2) + sum(this%mesh%x(:,fnode2),dim=2)) / &
                            (2*size(fnode1))
                values(i,j) = this%farray(i,n)%f%eval(args)
              end associate
            end do
          end if
        case default
          do j = 1, size(faces,dim=2)
            associate (fnode1 => this%mesh%fnode(this%mesh%xfnode(faces(1,j)):this%mesh%xfnode(faces(1,j)+1)-1), &
                       fnode2 => this%mesh%fnode(this%mesh%xfnode(faces(2,j)):this%mesh%xfnode(faces(2,j)+1)-1))
              args(1:) = (sum(this%mesh%x(:,fnode1),dim=2) + sum(this%mesh%x(:,fnode2),dim=2)) / &
                          (2*size(fnode1))
              values(i,j) = this%farray(i,n)%f%eval(args)
            end associate
          end do
        end select
      end do
    end do

    this%tlast = t
    this%evaluated = .true.

  end subroutine if_data_eval

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! IF_DATA_DESTROY
 !!

  subroutine if_data_destroy (this)

    type(if_data), intent(inout) :: this

    type(if_data) :: default

    if (associated(this%faces))  deallocate(this%faces)
    if (associated(this%values)) deallocate(this%values)
    if (associated(this%xgroup))  deallocate(this%xgroup)
    if (associated(this%hint)) deallocate(this%hint)
    if (associated(this%tag)) deallocate(this%tag)

    this = default  ! assign default initialization values

  end subroutine if_data_destroy

end module interface_data
