!!
!! BNDRY_FACE_FUNC_TYPE
!!
!! This module defines an implementation of the base class BNDRY_FUNC that
!! describes a time-dependent function on a subset of the boundary faces of
!! a mesh of type UNSTR_MESH.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for Fortran 2008, February 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! An instance of the derived type BNDRY_FACE_FUNC is designed to encapsulate
!! the time dependent data associated with a boundary condition, along with the
!! means of evaluating that data.  Once initialized, an instance of this type
!! has two public array components, INDEX and VALUE: for each j, %INDEX(j) is
!! the index of a boundary face and %VALUE(j) is the datum associated with
!! that face.  Note that the interpretation of the data is up to the client
!! code; Dirichlet BC data or flux BC data, for example.
!!
!! A BNDRY_FACE_FUNC object is defined incrementally using the following type
!! bound subroutines:
!!
!!  INIT(MESH) initializes the object to begin receiving the specification of
!!    the boundary face data.  MESH is the UNSTR_MESH object on which the data
!!    is defined. The object maintains an internal read-only reference to MESH.
!!
!!  ADD(F, SETIDS, STAT, ERRMSG) specifies that the SCALAR_FUNC class object F
!!    should be used to compute the data for the mesh faces belonging to the
!!    face sets with IDs given by the rank-1 array SETIDS.  The EVAL method of
!!    F will be passed [t, x, y, z] as its vector argument. This method must be
!!    called after INIT, and can be called multiple times, or even not at all.
!!    F is an allocatable variable and its allocation is moved into the object
!!    and returned unallocated. It is an error to reference an unknown face set
!!    ID, associate a face with more than one function, or specify a face that
!!    is not a boundary face.  A nonzero value is assigned to the integer STAT
!!    if an error occurs, and an explanatory message assigned to the deferred
!!    length allocatable character argument ERRMSG.  The returned STAT and
!!    ERRMSG values are collective across all processes.
!!
!!  ADD_COMPLETE() performs the final configuration of the object after all the
!!    desired calls to ADD have been made (if any).
!!
!! Once defined, the %VALUE component of the object is computed at time T by
!! calling the COMPUTE(T) type bound subroutine.
!!

#include "f90_assert.fpp"

module bndry_face_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use bndry_func_class
  use unstr_mesh_type
  use scalar_func_class
  use scalar_func_containers
  use bndry_face_group_builder_type
  implicit none
  private

  type, extends(bndry_func), public :: bndry_face_func
    private
    type(unstr_mesh), pointer :: mesh => null() ! reference only - do not own
    logical :: evaluated = .false.
    real(r8) :: tlast = -huge(1.0_r8)
    integer :: ngroup
    integer, allocatable :: xgroup(:)
    type(scalar_func_box), allocatable :: farray(:)
    integer, allocatable :: hint(:)
    ! temporaries used during construction
    type(bndry_face_group_builder), allocatable :: builder
    type(scalar_func_list) :: flist
  contains
    procedure :: init
    procedure :: add
    procedure :: add_complete
    procedure :: compute
  end type bndry_face_func

  !! Optimization hint values.
  integer, parameter :: DATA_HINT_NONE    = 0
  integer, parameter :: DATA_HINT_CONST   = 1
  integer, parameter :: DATA_HINT_T_INDEP = 2
  integer, parameter :: DATA_HINT_X_INDEP = 3

contains

  subroutine init(this, mesh)
    class(bndry_face_func), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    this%mesh => mesh
    allocate(this%builder)
    call this%builder%init(mesh)
  end subroutine init

  subroutine add(this, f, setids, stat, errmsg)
    class(bndry_face_func), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: f
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_face_group(setids, stat, errmsg)
    if (stat /= 0) return
    call this%flist%append(f)
  end subroutine add

  subroutine add_complete(this)
    use const_scalar_func_type
    class(bndry_face_func), intent(inout) :: this
    integer :: n
    ASSERT(allocated(this%builder))
    call this%builder%get_face_groups(this%ngroup, this%xgroup, this%index)
    deallocate(this%builder)
    allocate(this%value(size(this%index)))
    call scalar_func_list_to_box_array(this%flist, this%farray)
    !! For now we don't expose optimization hinting to the user,
    !! but we can determine directly which functions are constant.
    allocate(this%hint(this%ngroup))
    do n = 1, this%ngroup
      select type (f => this%farray(n)%f)
      type is (const_scalar_func)
        this%hint(n) = DATA_HINT_CONST
      class default
        this%hint(n) = DATA_HINT_NONE
      end select
    end do
  end subroutine add_complete


  subroutine compute(this, t)

    class(bndry_face_func), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: n, j
    real(r8) :: args(0:size(this%mesh%x,dim=1))

    !! Verify that THIS is in the correct state.
    ASSERT(allocated(this%index))

    if (this%evaluated .and. t == this%tlast) return  ! values already set for this T

    args(0) = t
    do n = 1, this%ngroup
      associate(faces => this%index(this%xgroup(n):this%xgroup(n+1)-1), &
                value => this%value(this%xgroup(n):this%xgroup(n+1)-1))
        select case (this%hint(n))
        case (DATA_HINT_CONST)
          if (.not.this%evaluated) value = this%farray(n)%f%eval(args)
        case (DATA_HINT_X_INDEP)
          value = this%farray(n)%f%eval(args)
        case (DATA_HINT_T_INDEP)
          if (.not.this%evaluated) then
            do j = 1, size(faces)
              associate (fnode => this%mesh%fnode(this%mesh%xfnode(faces(j)):this%mesh%xfnode(faces(j)+1)-1))
                args(1:) = sum(this%mesh%x(:,fnode),dim=2) / size(fnode)
                value(j) = this%farray(n)%f%eval(args)
              end associate
            end do
          end if
        case default
          do j = 1, size(faces)
            associate (fnode => this%mesh%fnode(this%mesh%xfnode(faces(j)):this%mesh%xfnode(faces(j)+1)-1))
              args(1:) = sum(this%mesh%x(:,fnode),dim=2) / size(fnode)
              value(j) = this%farray(n)%f%eval(args)
            end associate
          end do
        end select
      end associate
    end do

    this%tlast = t
    this%evaluated = .true.

  end subroutine compute

end module bndry_face_func_type
