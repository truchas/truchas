!!
!! MTC_INTFC_FUNC_TYPE
!!
!! This module defines an extension of the base class INTFC_FUNC3 that
!! implements the MTC interface condition flux function on a subset of the
!! interface faces of a mesh type that extends the UNSTR_BASE_MESH class.
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

module mtc_intfc_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use intfc_func3_class
  use unstr_base_mesh_class
  use scalar_func_containers
  use intfc_link_group_builder_type
  implicit none
  private

  type, extends(intfc_func3), public :: mtc_intfc_func
    private
    class(unstr_base_mesh), pointer :: mesh => null() ! reference only - do not own
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
    procedure :: compute_deriv2 => compute
  end type mtc_intfc_func

contains

  subroutine init(this, mesh)
    class(mtc_intfc_func), intent(out) :: this
    class(unstr_base_mesh), intent(in), target :: mesh
    this%mesh => mesh
    allocate(this%builder)
    call this%builder%init(mesh)
  end subroutine init

  subroutine add(this, mtc, setids, stat, errmsg)
    class(mtc_intfc_func), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: mtc
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_link_group(setids, stat, errmsg)
    if (stat /= 0) return
    call this%flist%append(mtc)
  end subroutine add

  subroutine add_complete(this)
    class(mtc_intfc_func), intent(inout) :: this
    ASSERT(allocated(this%builder))
    call this%builder%get_link_groups(this%ngroup, this%xgroup, this%index)
    deallocate(this%builder)
    allocate(this%value(size(this%index,2)), this%deriv2(2,size(this%index,2)))
    call scalar_func_list_to_box_array(this%flist, this%f)
  end subroutine add_complete

  subroutine compute(this, t, var1, var2)
    class(mtc_intfc_func), intent(inout) :: this
    real(r8), intent(in) :: t, var1(:), var2(:)
    integer :: n, j
    real(r8) :: args(-2:size(this%mesh%x,dim=1)), c1, c2, fp, fm, df, fdinc
    ASSERT(allocated(this%index))
    args(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(:,this%xgroup(n):this%xgroup(n+1)-1), &
                value => this%value(this%xgroup(n):this%xgroup(n+1)-1), &
                deriv2 => this%deriv2(:,this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index,dim=2)
          associate(fnode => this%mesh%face_node_list_view(index(1,j)))
            args(1:) = sum(this%mesh%x(:,fnode),dim=2)/size(fnode)
          end associate
          associate (v1 => var2(index(1,j)), v2 => var2(index(2,j)))
            args(-2) = max(var1(index(1,j)), var1(index(2,j)))
            args(-1) = max(v1, v2)
            c1 = this%f(n)%eval(args) * this%mesh%area(index(1,j))
            value(j) = c1*(v1 - v2)
            fdinc = max(1.0_r8, abs(args(-1))) * sqrt(epsilon(1.0_r8))
            fdinc = scale(1.0_r8,exponent(fdinc))
            args(-1) = max(v1, v2) + fdinc
            fp = this%f(n)%eval(args)
            args(-1) = max(v1, v2) - fdinc
            fm = this%f(n)%eval(args)
            df = (fp - fm) / (2*fdinc)
            c2 = df * (v1 - v2) * this%mesh%area(index(1,j))
            deriv2(1,j) =  c1 + merge(c2, 0.0_r8, v1 > v2)
            deriv2(2,j) = -c1 + merge(c2, 0.0_r8, v2 > v1)
          end associate
        end do
      end associate
    end do
  end subroutine compute

end module mtc_intfc_func_type
