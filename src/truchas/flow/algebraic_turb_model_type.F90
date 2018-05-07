!!
!! ALGEBRAIC_TURB_MODEL_TYPE
!!
!! This module provides an implementation of the abstract turbulence_model_class
!!
!! Peter Brady <ptb@lanl.gov>
!! May 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  Calculates cell-centered turbulent diffusivity based on the
!!  algebraic turbulence model:
!!
!!  Nu_Turb = cmu * turbulence_length * turbulence_velocity
!!
!!   where turbulence_length is specified via input, and
!!   turbulence_velocity is determined from algebraic relations with the
!!   local velocities approximating the turbulence intensities.  No
!!   solution of transport equations is required.
!!
!!   Effective (Apparent) Diffusivity = Nu_Turb + Molecular Diffusivity
!!
!!   In the case of momentum transport, the diffusivity is also called
!!   kinematic viscosity (nu), which is given by mu / rho.
!!
module algebraic_turb_model_type
  use kinds
  use turbulence_model_class
  use parameter_list_type
  use truchas_logging_services
  implicit none
  private

  public :: algebraic_turb_model, alloc_algebraic_turb_moel

  type, extends(turbulence_model) :: algebraic_turb_model
    real(r8) :: length
    real(r8) :: cmu
    real(r8) :: ke_fraction
    real(r8), allocatable :: nu_turb(:)
  end type algebraic_turb_model

contains

  subroutine alloc_algebraic_turb_moel(t)
    class(turbulence_model), allocatable, intent(out) :: t
    type(algebraic_turb_model), allocatable :: m

    allocate(m)
    call move_alloc(m, t)
  end subroutine alloc_algebraic_turb_moel


  subroutine read_params(this, params)
    class(algebraic_turb_model), intent(inout) :: this
    type(parameter_list), pointer, intent(in) :: params

    call params%get('length', this%length)
    call params%get('cmu', this%cmu, 0.05_r8)
    call params%get('ke-fraction', this%ke_fraction, 0.1_r8)

    if (this%length <= 0.0_r8) call TLS_fatal("turbulence length must be > 0")
    if (this%cmu <= 0.0_r8) call TLS_fatal("turbulence cmu must > 0")
    if (this%ke_fraction <= 0.0_r8 .or. this%ke_fraction >= 1.0_r8) &
        call TLS_fatal("turbulence ke-fraction must be in (0,1)")

  end subroutine read_params

  subroutine init(this, mesh)
    class(algebraic_turb_model), intent(inout) :: this
    type(flow_mesh), pointer, intent(in) :: mesh

    this%mesh => mesh
    allocate(this%nu_turb(this%mesh%mesh%ncell))
  end subroutine init

  subroutine setup(this, vel_cc)
    class(algebraic_turb_model), intent(inout) :: this
    real(r8), intent(in) :: vel_cc(:,:)
    !-
    integer :: i

    associate(cmu => this%cmu, l => this%length, f => this%ke_fraction)
      do i = 1, this%mesh%mesh%ncell
        nu_turb(i) = cmu*l*sqrt(0.5_r8*f*sum(vel_cc(:,i)**2))
      end do
    end associate

  end subroutine setup

  subroutine apply
    class(algebraic_turb_model), intent(inout) :: this
    real(r8), intent(inout) :: visc_cc(:)
  end subroutine apply

end module algebraic_turb_model_type
