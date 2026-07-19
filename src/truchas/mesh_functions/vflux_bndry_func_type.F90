!!
!! VFLUX_BNDRY_FUNC_TYPE
!!
!! This module defines an extension of the base class BNDRY_FIELD_FUNC that
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
  use bndry_field_func_class
  use unstr_base_mesh_class
  use scalar_func_containers
  use bndry_face_group_builder_type
  use vector_func_class
  use vector_func_containers
  implicit none
  private

  type, extends(bndry_field_func), public :: vflux_bndry_func
    private
    class(unstr_base_mesh), pointer :: mesh => null() ! reference only - do not own
    integer :: ngroup
    integer, allocatable :: xgroup(:), map(:) ! XGROUP partitions MAP
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
    procedure :: compute_value
    procedure :: compute_deriv
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
    call this%builder%get_unique_face_groups(this%ngroup, this%xgroup, this%map, this%index)
    deallocate(this%builder)
    call scalar_func_list_to_box_array(this%flist, this%farray)
    call vector_func_list_to_box_array(this%glist, this%garray)
  end subroutine add_complete

  subroutine compute_value(this, t, u, value)
    class(vflux_bndry_func), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: u(:)
    real(r8), allocatable, intent(out) :: value(:)
    integer :: n, i, j, face
    real(r8) :: args(0:size(this%mesh%x,dim=1)), abs_args(1)
    real(r8) :: absorptivity, irrad(3)
    ASSERT(allocated(this%index))
    allocate(value(size(this%index)), source=0.0_r8)
    args(0) = t
    do n = 1, this%ngroup
      do i = this%xgroup(n), this%xgroup(n+1)-1
        j = this%map(i)
        face = this%index(j)
        associate (fnode => this%mesh%face_node_list_view(face))
          args(1:3) = sum(this%mesh%x(:,fnode),dim=2) / size(fnode)
        end associate
        abs_args(1) = u(face)
        absorptivity = this%farray(n)%eval(abs_args)
        irrad = this%garray(n)%f%eval(args)
        value(j) = value(j) + absorptivity*dot_product(this%mesh%normal(:,face), irrad)
      end do
    end do
  end subroutine compute_value

  subroutine compute_deriv(this, t, u, deriv)
    class(vflux_bndry_func), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: u(:)
    real(r8), allocatable, intent(out) :: deriv(:)
    integer :: n, i, j, face
    real(r8) :: args(0:size(this%mesh%x,dim=1)), abs_args(1)
    real(r8) :: dabs, irrad(3)
    ASSERT(allocated(this%index))
    allocate(deriv(size(this%index)), source=0.0_r8)
    args(0) = t
    do n = 1, this%ngroup
      do i = this%xgroup(n), this%xgroup(n+1)-1
        j = this%map(i)
        face = this%index(j)
        associate (fnode => this%mesh%face_node_list_view(face))
          args(1:) = sum(this%mesh%x(:,fnode),dim=2) / size(fnode)
        end associate
        abs_args(1) = u(face)
        call compute_first_arg_deriv(this%farray(n), abs_args, dabs)
        irrad = this%garray(n)%f%eval(args)
        deriv(j) = deriv(j) + dabs*dot_product(this%mesh%normal(:,face), irrad)
      end do
    end do
  end subroutine compute_deriv

  subroutine compute_first_arg_deriv(f, args, deriv)
    type(scalar_func_box), intent(in) :: f
    real(r8), intent(inout) :: args(:)
    real(r8), intent(out) :: deriv
    real(r8) :: fdinc, fm, fp, u

    ! Differentiate with respect to the first function argument.
    u = args(1)
    fdinc = max(1.0_r8, abs(u)) * sqrt(epsilon(1.0_r8))
    fdinc = scale(1.0_r8, exponent(fdinc))
    args(1) = u + fdinc
    fp = f%eval(args)
    args(1) = u - fdinc
    fm = f%eval(args)
    args(1) = u
    deriv = (fp - fm) / (2*fdinc)
  end subroutine compute_first_arg_deriv

end module vflux_bndry_func_type
