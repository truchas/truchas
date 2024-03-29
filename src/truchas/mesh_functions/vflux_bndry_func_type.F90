!!
!! VFLUX_BNDRY_FUNC_TYPE
!!
!! This module defines an extension of the base class BNDRY_FUNC2 that
!! implements the VFLUX boundary condition flux function on a subset of the
!! boundary faces of a mesh type that extends the UNSTR_BASE_MESH class.
!!
!! Narendran Raghavan <naren@lanl.gov>
!! June 2022
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module vflux_bndry_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use bndry_func2_class
  use unstr_base_mesh_class
  use scalar_func_containers
  use bndry_face_group_builder_type
  use vector_func_class
  use vector_func_containers
  implicit none
  private

  type, extends(bndry_func2), public :: vflux_bndry_func
    private
    class(unstr_base_mesh), pointer :: mesh => null() ! reference only - do not own
    integer :: ngroup
    integer, allocatable :: xgroup(:)
    type(scalar_func_box), allocatable :: farray(:)
    type(vector_func_box), allocatable :: garray(:)
    ! temporaries used during construction
    type(bndry_face_group_builder), allocatable :: builder
    type(scalar_func_list) :: flist
    type(vector_func_list) :: glist
  contains
    procedure :: init
    procedure :: add
    procedure :: add_complete
    procedure :: compute
    procedure :: compute_value => compute
    procedure :: compute_deriv => compute
  end type vflux_bndry_func

contains

  subroutine init(this, mesh)
    class(vflux_bndry_func), intent(out) :: this
    class(unstr_base_mesh), intent(in), target :: mesh
    this%mesh => mesh
    allocate(this%builder)
    call this%builder%init(mesh, no_overlap=.false.)
  end subroutine init

  subroutine add(this, f, g, setids, stat, errmsg)
    class(vflux_bndry_func), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: f
    class(vector_func), allocatable, intent(inout) :: g
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_face_group(setids, stat, errmsg)
    if (stat /= 0) return
    call this%flist%append(f)
    call this%glist%append(g)
  end subroutine add

  subroutine add_complete(this)
    class(vflux_bndry_func), intent(inout) :: this
    ASSERT(allocated(this%builder))
    call this%builder%get_face_groups(this%ngroup, this%xgroup, this%index)
    deallocate(this%builder)
    allocate(this%value(size(this%index)), this%deriv(size(this%index)))
    call scalar_func_list_to_box_array(this%flist, this%farray)
    call vector_func_list_to_box_array(this%glist, this%garray)
  end subroutine add_complete

  subroutine compute(this, t, var)
    class(vflux_bndry_func), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: var(:)
    integer :: n, j
    real(r8) :: args(0:size(this%mesh%x,dim=1))
    real(r8) :: absorptivity, irrad(3)
    ASSERT(allocated(this%index))
    args(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(this%xgroup(n):this%xgroup(n+1)-1), &
                value => this%value(this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index)
          associate (fnode => this%mesh%face_node_list_view(index(j)))
            args(1:3) = sum(this%mesh%x(:,fnode),dim=2) / size(fnode)
          end associate
          absorptivity = this%farray(n)%eval(var)
          irrad = this%garray(n)%f%eval(args)
          value(j) = absorptivity*dot_product(this%mesh%normal(:,index(j)), irrad)
        end do
      end associate
    end do
  end subroutine compute

end module vflux_bndry_func_type
