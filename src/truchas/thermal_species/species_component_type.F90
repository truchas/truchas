!!
!! SPECIES_COMPONENT_TYPE
!!
!! This module defines the derived type that collects the data common to the
!! species transport part of the SD and HTSD models: species material
!! properties, source terms, and boundary/interface conditions. Model and
!! factory code are responsible for defining the components and interpreting
!! them in residual and preconditioner assembly; the type provides focused
!! helper operations for behavior that belongs to its own boundary/interface
!! condition objects.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module species_component_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use cell_prop_func_class
  use scalar_mesh_func_class
  use bndry_func1_class
  use bndry_func3_class
  use intfc_multifield_func_class
  use mfd_diff_matrix_type, only: mfd_diff_matrix
  implicit none
  private

  type, public :: species_component
    class(cell_prop_func), allocatable :: diffusivity
    real(r8), allocatable :: ext_rate(:)
    class(scalar_mesh_func), allocatable :: src
    class(cell_prop_func), allocatable :: soret
    class(bndry_func1), allocatable :: bc_dir
    class(bndry_func1), allocatable :: bc_flux
    class(bndry_func3), allocatable :: bc_mtc
    class(intfc_multifield_func), allocatable :: ic_mtc
  contains
    procedure, private :: add_flux_bc_residual_c
    procedure, private :: add_flux_bc_residual_tc
    generic :: add_flux_bc_residual => add_flux_bc_residual_c, add_flux_bc_residual_tc
    procedure, private :: add_flux_bc_jacobian_c
    procedure, private :: add_flux_bc_jacobian_tc
    generic :: add_flux_bc_jacobian => add_flux_bc_jacobian_c, add_flux_bc_jacobian_tc
  end type species_component

contains

  subroutine add_flux_bc_residual_c(this, t, Cface, Fface, area, void_face)
    class(species_component), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: Cface(:), area(:)
    real(r8), intent(inout) :: Fface(:)
    logical, intent(in), optional :: void_face(:)

    integer :: j, n

    !! Simple flux BC contribution.
    if (allocated(this%bc_flux)) then
      call this%bc_flux%compute(t)
      do j = 1, size(this%bc_flux%index)
        n = this%bc_flux%index(j)
        Fface(n) = Fface(n) + area(n)*this%bc_flux%value(j)
      end do
    end if

    !! Mass-transfer-coefficient (MTC) BC flux contribution.
    if (allocated(this%bc_mtc)) then
      block
        real(r8), allocatable :: value(:)
        call this%bc_mtc%compute_value(t, Cface, value)
        do j = 1, size(this%bc_mtc%index)
          n = this%bc_mtc%index(j)
          Fface(n) = Fface(n) + value(j)
        end do
      end block
    end if

    !! Internal MTC flux contribution.
    if (allocated(this%ic_mtc)) then
      block
        real(r8), allocatable :: value(:)
        call this%ic_mtc%compute_value(t, Cface, value)
        do j = 1, size(this%ic_mtc%index,2)
          if (present(void_face)) then
            if (any(void_face(this%ic_mtc%index(:,j)))) cycle
          end if
          call add_interface_flux(this%ic_mtc%index(:,j), value(j))
        end do
      end block
    end if

  contains

    subroutine add_interface_flux(face, value)
      integer, intent(in) :: face(2)
      real(r8), intent(in) :: value
      Fface(face(1)) = Fface(face(1)) + value
      Fface(face(2)) = Fface(face(2)) - value
    end subroutine

  end subroutine

  subroutine add_flux_bc_residual_tc(this, t, Tface, Cface, Fface, area, void_face)
    class(species_component), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: Tface(:), Cface(:), area(:)
    real(r8), intent(inout) :: Fface(:)
    logical, intent(in), optional :: void_face(:)

    integer :: j, n

    !! Simple flux BC contribution.
    if (allocated(this%bc_flux)) then
      call this%bc_flux%compute(t)
      do j = 1, size(this%bc_flux%index)
        n = this%bc_flux%index(j)
        Fface(n) = Fface(n) + area(n)*this%bc_flux%value(j)
      end do
    end if

    !! Mass-transfer-coefficient (MTC) BC flux contribution.
    if (allocated(this%bc_mtc)) then
      block
        real(r8), allocatable :: value(:)
        call this%bc_mtc%compute_value(t, Cface, Tface, value)
        do j = 1, size(this%bc_mtc%index)
          n = this%bc_mtc%index(j)
          Fface(n) = Fface(n) + value(j)
        end do
      end block
    end if

    !! Internal MTC flux contribution.
    if (allocated(this%ic_mtc)) then
      block
        real(r8), allocatable :: value(:)
        call this%ic_mtc%compute_value(t, Cface, Tface, value)
        do j = 1, size(this%ic_mtc%index,2)
          if (present(void_face)) then
            if (any(void_face(this%ic_mtc%index(:,j)))) cycle
          end if
          call add_interface_flux(this%ic_mtc%index(:,j), value(j))
        end do
      end block
    end if

  contains

    subroutine add_interface_flux(face, value)
      integer, intent(in) :: face(2)
      real(r8), intent(in) :: value
      Fface(face(1)) = Fface(face(1)) + value
      Fface(face(2)) = Fface(face(2)) - value
    end subroutine

  end subroutine

  subroutine add_flux_bc_jacobian_c(this, t, Cface, matrix, void_face)
    class(species_component), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: Cface(:)
    type(mfd_diff_matrix), intent(inout) :: matrix
    logical, intent(in), optional :: void_face(:)

    !! Mass-transfer-coefficient (MTC) BC flux contribution.
    if (allocated(this%bc_mtc)) then
      block
        real(r8), allocatable :: deriv(:)
        call this%bc_mtc%compute_deriv(t, Cface, deriv)
        call matrix%incr_face_diag(this%bc_mtc%index, deriv)
      end block
    end if

    !! Internal MTC flux contribution.
    if (allocated(this%ic_mtc)) then
      block
        real(r8), allocatable :: deriv(:,:)
        call this%ic_mtc%compute_deriv(t, Cface, deriv)
        call add_interface_jacobian(this%ic_mtc%index, deriv)
      end block
    end if

  contains

    subroutine add_interface_jacobian(index, deriv)
      integer, intent(in) :: index(:,:)
      real(r8), intent(inout) :: deriv(:,:)
      integer :: j

      if (present(void_face)) then
        do j = 1, size(index,2)
          if (any(void_face(index(:,j)))) deriv(:,j) = 0.0_r8
        end do
      end if
      call matrix%incr_interface_flux3(index, deriv)
    end subroutine

  end subroutine

  subroutine add_flux_bc_jacobian_tc(this, t, Tface, Cface, matrix, void_face)
    class(species_component), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: Tface(:), Cface(:)
    type(mfd_diff_matrix), intent(inout) :: matrix
    logical, intent(in), optional :: void_face(:)

    !! Mass-transfer-coefficient (MTC) BC flux contribution.
    if (allocated(this%bc_mtc)) then
      block
        real(r8), allocatable :: deriv1(:)
        call this%bc_mtc%compute_deriv1(t, Cface, Tface, deriv1)
        call matrix%incr_face_diag(this%bc_mtc%index, deriv1)
      end block
    end if

    !! Internal MTC flux contribution.
    if (allocated(this%ic_mtc)) then
      block
        real(r8), allocatable :: deriv1(:,:)
        call this%ic_mtc%compute_deriv1(t, Cface, Tface, deriv1)
        call add_interface_jacobian(this%ic_mtc%index, deriv1)
      end block
    end if

  contains

    subroutine add_interface_jacobian(index, deriv)
      integer, intent(in) :: index(:,:)
      real(r8), intent(inout) :: deriv(:,:)
      integer :: j

      if (present(void_face)) then
        do j = 1, size(index,2)
          if (any(void_face(index(:,j)))) deriv(:,j) = 0.0_r8
        end do
      end if
      call matrix%incr_interface_flux3(index, deriv)
    end subroutine

  end subroutine

end module species_component_type
