!!
!! RAD_BNDRY_FUNC_TYPE
!!
!! This module defines an extension of the base class BNDRY_FUNC2 that
!! implements the radiation boundary condition flux function on a subset of
!! the boundary faces of a mesh type that extends the UNSTR_BASE_MESH class.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2017
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module rad_bndry_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use bndry_func2_class
  use unstr_base_mesh_class
  use scalar_func_containers
  use bndry_face_group_builder_type
  implicit none
  private

  type, extends(bndry_func2), public :: rad_bndry_func
    private
    class(unstr_base_mesh), pointer :: mesh => null() ! reference only - do not own
    real(r8) :: sigma, abszero
    integer :: ngroup
    integer, allocatable :: xgroup(:)
    type(scalar_func_box), allocatable :: f(:), g(:)
    ! temporaries used during construction
    type(bndry_face_group_builder), allocatable :: builder
    type(scalar_func_list) :: flist, glist
  contains
    procedure :: init
    procedure :: add
    procedure :: add_complete
    procedure :: compute
    procedure :: compute_value => compute
    procedure :: compute_deriv => compute
  end type rad_bndry_func

contains

  subroutine init(this, mesh, sigma, abszero)
    class(rad_bndry_func), intent(out) :: this
    class(unstr_base_mesh), intent(in), target :: mesh
    real(r8), intent(in) :: sigma, abszero
    this%mesh => mesh
    this%sigma = sigma
    this%abszero = abszero
    allocate(this%builder)
    call this%builder%init(mesh)
  end subroutine init

  subroutine add(this, eps, tamb, setids, stat, errmsg)
    class(rad_bndry_func), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: eps, tamb
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_face_group(setids, stat, errmsg)
    if (stat /= 0) return
    call this%flist%append(eps)
    call this%glist%append(tamb)
  end subroutine add

  subroutine add_complete(this)
    class(rad_bndry_func), intent(inout) :: this
    ASSERT(allocated(this%builder))
    call this%builder%get_face_groups(this%ngroup, this%xgroup, this%index)
    deallocate(this%builder)
    allocate(this%value(size(this%index)), this%deriv(size(this%index)))
    call scalar_func_list_to_box_array(this%flist, this%f)
    call scalar_func_list_to_box_array(this%glist, this%g)
  end subroutine add_complete

  subroutine compute(this, t, var)
    class(rad_bndry_func), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: var(:)
    integer :: n, j
    real(r8) :: args(-1:size(this%mesh%x,dim=1)), temp, eps, tamb, c, tmp, fdinc, eps_p, eps_m, deps
    ASSERT(allocated(this%index))
    args(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(this%xgroup(n):this%xgroup(n+1)-1), &
                value => this%value(this%xgroup(n):this%xgroup(n+1)-1), &
                deriv => this%deriv(this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index)
          associate (fnode => this%mesh%face_node_list_view(index(j)))
            args(1:3) = sum(this%mesh%x(:,fnode),dim=2) / size(fnode)
          end associate
          temp = var(index(j))
          args(-1) = temp
          eps  = this%f(n)%eval(args)
          tamb = this%g(n)%eval(args(0:))
          c = this%mesh%area(index(j)) * this%sigma
          tmp = (temp-this%abszero)**4 - (tamb-this%abszero)**4
          value(j) = c*eps*tmp
          fdinc = max(1.0_r8, abs(temp)) * sqrt(epsilon(1.0_r8))
          fdinc = scale(1.0_r8,exponent(fdinc))
          args(-1) = temp + fdinc
          eps_p = this%f(n)%eval(args)
          args(-1) = temp - fdinc
          eps_m = this%f(n)%eval(args)
          deps = (eps_p - eps_m) / (2*fdinc)
          deriv(j) = 4*c*eps*(temp-this%abszero)**3 + c*deps*tmp
        end do
      end associate
    end do
  end subroutine compute

end module rad_bndry_func_type
