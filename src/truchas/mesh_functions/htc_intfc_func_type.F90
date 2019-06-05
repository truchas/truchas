!!
!! HTC_INTFC_FUNC_TYPE
!!
!! This module defines an extension of the base class INTFC_FUNC2 that
!! implements the HTC interface condition flux function on a subset of the
!! interface faces of a mesh of type UNSTR_MESH.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module htc_intfc_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use intfc_func2_class
  use unstr_mesh_type
  use scalar_func_containers
  use intfc_link_group_builder_type
  implicit none
  private

  type, extends(intfc_func2), public :: htc_intfc_func
    private
    type(unstr_mesh), pointer :: mesh => null() ! reference only - do not own
    integer :: ngroup
    integer, allocatable :: xgroup(:)
    type(scalar_func_box), allocatable :: f(:)
    ! temporaries used during construction
    type(intfc_link_group_builder), allocatable :: builder
    type(scalar_func_list) :: flist
  contains
    procedure :: init
    procedure :: add
    procedure :: add_complete
    procedure :: compute
    procedure :: compute_value => compute
    procedure :: compute_deriv => compute
  end type htc_intfc_func

contains

  subroutine init(this, mesh)
    class(htc_intfc_func), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    this%mesh => mesh
    allocate(this%builder)
    call this%builder%init(mesh)
  end subroutine init

  subroutine add(this, htc, setids, stat, errmsg)
    class(htc_intfc_func), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: htc
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_link_group(setids, stat, errmsg)
    if (stat /= 0) return
    call this%flist%append(htc)
  end subroutine add

  subroutine add_complete(this)
    class(htc_intfc_func), intent(inout) :: this
    ASSERT(allocated(this%builder))
    call this%builder%get_link_groups(this%ngroup, this%xgroup, this%index)
    deallocate(this%builder)
    allocate(this%value(size(this%index,2)), this%deriv(2,size(this%index,2)))
    call scalar_func_list_to_box_array(this%flist, this%f)
  end subroutine add_complete

  subroutine compute(this, t, var)
    class(htc_intfc_func), intent(inout) :: this
    real(r8), intent(in) :: t, var(:)
    integer :: n, j
    real(r8) :: args(0:size(this%mesh%x,dim=1)), c
    ASSERT(allocated(this%index))
    args(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(:,this%xgroup(n):this%xgroup(n+1)-1), &
                value => this%value(this%xgroup(n):this%xgroup(n+1)-1), &
                deriv => this%deriv(:,this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index,dim=2)
          associate(fnode => this%mesh%fnode(this%mesh%xfnode(index(1,j)):this%mesh%xfnode(index(1,j)+1)-1))
            args(1:) = sum(this%mesh%x(:,fnode),dim=2)/size(fnode)
          end associate
          c = this%f(n)%eval(args)*this%mesh%area(index(1,j))
          value(j) = c*(var(index(1,j)) - var(index(2,j)))
          deriv(1,j) = c
          deriv(2,j) = -c
        end do
      end associate
    end do
  end subroutine compute

end module htc_intfc_func_type
