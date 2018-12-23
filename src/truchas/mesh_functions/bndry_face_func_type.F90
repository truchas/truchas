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
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
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
!!  INIT(MESH [,BNDRY_ONLY]) initializes the object to begin receiving the
!!    specification of the boundary face data. MESH is the UNSTR_MESH object
!!    on which the data s defined. The object maintains an internal read-only
!!    reference to MESH. Faces will be required to be boundary faces unless
!!    the option BNDRY_ONLY is present with value .false.
!!
!!  ADD(F, SETIDS, STAT, ERRMSG) specifies that the SCALAR_FUNC class object F
!!    should be used to compute the data for the mesh faces belonging to the
!!    face sets with IDs given by the rank-1 array SETIDS.  The EVAL method of
!!    F will be passed [t, x, y, z] as its vector argument. This method must be
!!    called after INIT, and can be called multiple times, or even not at all.
!!    F is an allocatable variable and its allocation is moved into the object
!!    and returned unallocated. It is an error to reference an unknown face set
!!    ID, associate a face with more than one function, or specify a face that
!!    is not a boundary face unless BNDRY_ONLY was specified with value .false.
!!    A nonzero value is assigned to the integer STAT if an error occurs, and
!!    an explanatory message assigned to the deferred length allocatable
!!    character argument ERRMSG.  The returned STAT and ERRMSG values are
!!    collective across all processes.
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
    procedure :: add_face_list
    procedure :: compute
  end type bndry_face_func

  !! Optimization hint values.
  integer, parameter :: DATA_HINT_NONE    = 0
  integer, parameter :: DATA_HINT_CONST   = 1
  integer, parameter :: DATA_HINT_T_INDEP = 2
  integer, parameter :: DATA_HINT_X_INDEP = 3

contains

  subroutine init(this, mesh, bndry_only)
    class(bndry_face_func), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    logical, intent(in), optional :: bndry_only
    this%mesh => mesh
    allocate(this%builder)
    call this%builder%init(mesh, bndry_only)
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
        ! these hints don't play well with the dp pressure boundary conditions
!!$     type is (const_scalar_func)
!!$        this%hint(n) = DATA_HINT_CONST
      class default
        this%hint(n) = DATA_HINT_NONE
      end select
    end do
  end subroutine add_complete

  !! Optionally called after init, add, and add_complete
  !! to append a list of faces to the object.
  subroutine add_face_list(this, f, faceids)

    class(bndry_face_func), intent(inout) :: this
    class(scalar_func), allocatable, intent(in) :: f
    integer, intent(in) :: faceids(:)

    integer, allocatable :: tmp(:)
    type(scalar_func_box), allocatable :: ftmp(:)

    ASSERT(allocated(this%value) .and. allocated(this%index))
    ASSERT(.not.allocated(this%builder))

    !! Update the index array
    tmp = this%index
    call move_alloc(this%index, tmp)
    allocate(this%index(size(tmp) + size(faceids)))
    this%index(:size(tmp)) = tmp
    this%index(size(tmp)+1:) = faceids
    deallocate(tmp)

    !! Update the value array
    deallocate(this%value)
    allocate(this%value(size(this%index)))

    !! Update the hint array
    call move_alloc(this%hint, tmp)
    allocate(this%hint(size(tmp)+1))
    this%hint(:size(tmp)) = tmp
    this%hint(size(tmp)+1) = DATA_HINT_NONE
    deallocate(tmp)

    !! Update the group array
    this%ngroup = this%ngroup + 1
    call move_alloc(this%xgroup, tmp)
    allocate(this%xgroup(size(tmp)+1))
    this%xgroup(:size(tmp)) = tmp
    this%xgroup(size(tmp)+1) = size(this%index) + 1
    deallocate(tmp)

    !! Update the function list
    call move_alloc(this%farray, ftmp)
    allocate(this%farray(size(ftmp)+1))
    this%farray(:size(ftmp)) = ftmp
    this%farray(size(ftmp)+1)%f = f
    deallocate(ftmp)

  end subroutine add_face_list


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
