!!
!! RAD_INTFC_FUNC_TYPE
!!
!! This module defines an extension of the base class INTFC_FIELD_FUNC that
!! implements the gap radiation interface condition flux function on a subset of
!! the interface faces of a mesh type that extends the UNSTR_BASE_MESH class.
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

module rad_intfc_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use intfc_field_func_class
  use unstr_base_mesh_class
  use scalar_func_containers
  use intfc_link_group_builder_type
  implicit none
  private

  type, extends(intfc_field_func), public :: rad_intfc_func
    private
    class(unstr_base_mesh), pointer :: mesh => null() ! reference only - do not own
    real(r8) :: sigma, abszero
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
    procedure :: compute_value
    procedure :: compute_deriv
  end type rad_intfc_func

contains

  subroutine init(this, mesh, sigma, abszero)
    class(rad_intfc_func), intent(out) :: this
    class(unstr_base_mesh), intent(in), target :: mesh
    real(r8), intent(in) :: sigma, abszero
    this%mesh => mesh
    this%sigma = sigma
    this%abszero = abszero
    allocate(this%builder)
    call this%builder%init(mesh)
  end subroutine init

  subroutine add(this, eps, setids, stat, errmsg)
    class(rad_intfc_func), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: eps
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_link_group(setids, stat, errmsg)
    if (stat /= 0) return
    call this%flist%append(eps)
  end subroutine add

  subroutine add_complete(this)
    class(rad_intfc_func), intent(inout) :: this
    ASSERT(allocated(this%builder))
    call this%builder%get_link_groups(this%ngroup, this%xgroup, this%index)
    deallocate(this%builder)
    call scalar_func_list_to_box_array(this%flist, this%f)
  end subroutine add_complete

  subroutine compute_value(this, t, u, value)
    class(rad_intfc_func), intent(inout) :: this
    real(r8), intent(in) :: t, u(:)
    real(r8), allocatable, intent(out) :: value(:)
    integer :: n, j
    real(r8) :: args(-1:size(this%mesh%x,dim=1)), c, f, g
    ASSERT(allocated(this%index))
    allocate(value(size(this%index,2)))
    args(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(:,this%xgroup(n):this%xgroup(n+1)-1), &
                value => value(this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index,dim=2)
          associate(fnode => this%mesh%face_node_list_view(index(1,j)))
            args(1:) = sum(this%mesh%x(:,fnode),dim=2)/size(fnode)
          end associate
          associate (v1 => u(index(1,j)), v2 => u(index(2,j)))
            args(-1) = max(v1, v2)
            c = this%sigma * this%mesh%area(index(1,j))
            f = this%f(n)%eval(args)
            g = (v1-this%abszero)**4 - (v2-this%abszero)**4
            value(j) = c*f*g
          end associate
        end do
      end associate
    end do
  end subroutine compute_value

  subroutine compute_deriv(this, t, u, deriv)
    class(rad_intfc_func), intent(inout) :: this
    real(r8), intent(in) :: t, u(:)
    real(r8), allocatable, intent(out) :: deriv(:,:)
    integer :: n, j
    real(r8) :: args(-1:size(this%mesh%x,dim=1))
    real(r8) :: c, f, g, fp, fm, df, fdinc, tmp, umax
    ASSERT(allocated(this%index))
    allocate(deriv(2,size(this%index,2)))
    args(0) = t
    do n = 1, this%ngroup
      associate(index => this%index(:,this%xgroup(n):this%xgroup(n+1)-1), &
                deriv => deriv(:,this%xgroup(n):this%xgroup(n+1)-1))
        do j = 1, size(index,dim=2)
          associate(fnode => this%mesh%face_node_list_view(index(1,j)))
            args(1:) = sum(this%mesh%x(:,fnode),dim=2)/size(fnode)
          end associate
          associate (v1 => u(index(1,j)), v2 => u(index(2,j)))
            umax = max(v1, v2)
            args(-1) = umax
            c = this%sigma * this%mesh%area(index(1,j))
            f = this%f(n)%eval(args)
            g = (v1-this%abszero)**4 - (v2-this%abszero)**4
            fdinc = max(1.0_r8, abs(umax)) * sqrt(epsilon(1.0_r8))
            fdinc = scale(1.0_r8,exponent(fdinc))
            args(-1) = umax + fdinc
            fp = this%f(n)%eval(args)
            args(-1) = umax - fdinc
            fm = this%f(n)%eval(args)
            df = (fp - fm) / (2*fdinc)
            tmp = c*df*g
            deriv(1,j) =  4*c*f*(v1-this%abszero)**3 + merge(tmp, 0.0_r8, v1 > v2)
            deriv(2,j) = -4*c*f*(v2-this%abszero)**3 + merge(tmp, 0.0_r8, v2 > v1)
          end associate
        end do
      end associate
    end do
  end subroutine compute_deriv

end module rad_intfc_func_type
