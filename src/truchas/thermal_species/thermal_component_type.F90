!!
!! THERMAL_COMPONENT_TYPE
!!
!! This module defines the derived type that collects the data common to the
!! enthalpy transport part of the HT, HTSD, and FHT models: thermal material
!! properties, source terms, and boundary/interface conditions. Model and
!! factory code are responsible for defining the components and interpreting
!! them in residual and preconditioner assembly; the type provides focused
!! helper operations for behavior that belongs to its own boundary/interface
!! condition objects.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module thermal_component_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use cell_prop_func_class
  use scalar_mesh_multifunc_type
  use bndry_func1_class
  use bndry_func2_class
  use intfc_field_func_class
  use mfd_diff_matrix_type, only: mfd_diff_matrix
  implicit none
  private

  type, public :: thermal_component
    class(cell_prop_func), allocatable :: conductivity
    class(cell_prop_func), allocatable :: H_of_T
    real(r8), allocatable :: ext_rate(:)
    type(scalar_mesh_multifunc), allocatable :: src
    class(bndry_func1), allocatable :: bc_dir
    class(bndry_func1), allocatable :: bc_flux
    class(bndry_func2), allocatable :: bc_vflux
    class(bndry_func2), allocatable :: bc_htc
    class(bndry_func2), allocatable :: bc_rad
    class(intfc_field_func), allocatable :: ic_htc
    class(intfc_field_func), allocatable :: ic_rad
    class(bndry_func2), allocatable :: evap_flux
  contains
    procedure :: add_flux_bc_residual
    procedure :: add_flux_bc_jacobian
  end type

contains

  subroutine add_flux_bc_residual(this, t, Tface, Fface, area, void_face)
    class(thermal_component), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: Tface(:), area(:)
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

    !! Oriented flux BC contribution.
    if (allocated(this%bc_vflux)) then
      call this%bc_vflux%compute(t, Tface)
      do j = 1, size(this%bc_vflux%index)
        n = this%bc_vflux%index(j)
        Fface(n) = Fface(n) + this%bc_vflux%value(j)
      end do
    end if

    !! External HTC flux contribution.
    if (allocated(this%bc_htc)) then
      call this%bc_htc%compute(t, Tface)
      do j = 1, size(this%bc_htc%index)
        n = this%bc_htc%index(j)
        Fface(n) = Fface(n) + this%bc_htc%value(j)
      end do
    end if

    !! Ambient radiation BC flux contribution.
    if (allocated(this%bc_rad)) then
      call this%bc_rad%compute(t, Tface)
      do j = 1, size(this%bc_rad%index)
        n = this%bc_rad%index(j)
        Fface(n) = Fface(n) + this%bc_rad%value(j)
      end do
    end if

    !! Experimental evaporation heat flux contribution.
    if (allocated(this%evap_flux)) then
      call this%evap_flux%compute_value(t, Tface)
      do j = 1, size(this%evap_flux%index)
        n = this%evap_flux%index(j)
        Fface(n) = Fface(n) + area(n)*this%evap_flux%value(j)
      end do
    end if

    !! Internal HTC flux contribution.
    if (allocated(this%ic_htc)) call add_interface_flux(this%ic_htc)

    !! Internal gap radiation contribution.
    if (allocated(this%ic_rad)) call add_interface_flux(this%ic_rad)

  contains

    subroutine add_interface_flux(bc)
      class(intfc_field_func), intent(in) :: bc
      integer :: j, n1, n2
      real(r8) :: value(size(bc%index,2))
      call bc%compute_value(t, Tface, value)
      do j = 1, size(bc%index,2)
        if (present(void_face)) then
          if (any(void_face(bc%index(:,j)))) cycle
        end if
        n1 = bc%index(1,j)
        n2 = bc%index(2,j)
        Fface(n1) = Fface(n1) + value(j)
        Fface(n2) = Fface(n2) - value(j)
      end do
    end subroutine

  end subroutine

  subroutine add_flux_bc_jacobian(this, t, Tface, area, matrix, void_face)
    class(thermal_component), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: Tface(:), area(:)
    type(mfd_diff_matrix), intent(inout) :: matrix
    logical, intent(in), optional :: void_face(:)

    !! TODO: Oriented flux currently has no Jacobian contribution. Preserve
    !! that existing behavior until the appropriate approximation is clarified.

    !! External HTC boundary condition contribution.
    if (allocated(this%bc_htc)) then
      call this%bc_htc%compute_deriv(t, Tface)
      call matrix%incr_face_diag(this%bc_htc%index, this%bc_htc%deriv)
    end if

    !! Simple radiation boundary condition contribution.
    if (allocated(this%bc_rad)) then
      call this%bc_rad%compute_deriv(t, Tface)
      call matrix%incr_face_diag(this%bc_rad%index, this%bc_rad%deriv)
    end if

    !! Experimental evaporation heat flux contribution.
    if (allocated(this%evap_flux)) then
      call this%evap_flux%compute_deriv(t, Tface)
      associate (index => this%evap_flux%index, deriv => this%evap_flux%deriv)
        call matrix%incr_face_diag(index, area(index)*deriv)
      end associate
    end if

    !! Internal HTC interface condition contribution.
    if (allocated(this%ic_htc)) call add_interface_jacobian(this%ic_htc)

    !! Internal gap radiation condition contribution.
    if (allocated(this%ic_rad)) call add_interface_jacobian(this%ic_rad)

  contains

    subroutine add_interface_jacobian(bc)
      class(intfc_field_func), intent(in) :: bc
      integer :: j
      real(r8) :: deriv(2,size(bc%index,2))

      call bc%compute_deriv(t, Tface, deriv)
      if (present(void_face)) then
        do j = 1, size(bc%index,2)
          if (any(void_face(bc%index(:,j)))) deriv(:,j) = 0.0_r8
        end do
      end if
      call matrix%incr_interface_flux3(bc%index, deriv)
    end subroutine

  end subroutine

end module thermal_component_type
