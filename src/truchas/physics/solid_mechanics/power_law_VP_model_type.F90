!!
!! POWER_LAW_VP_MODEL_TYPE
!!
!! The implementation of the power law viscoplastic model.  This wraps
!! existing code by Dave Korzekwa to compute the strain rate.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module power_law_VP_model_type

  use kinds, only: r8
  use VP_model_class
  implicit none
  private
  
  type, extends(VP_model), public :: power_law_VP_model
    private
    real(r8) :: A, N, Q, R
  contains
    procedure :: strain_rate
  end type power_law_VP_model
  
  public :: new_power_law_VP_model

contains
  
  function new_power_law_VP_model (A, N, Q, R) result (model)
    real(r8), intent(in) :: A, N, Q, R
    class(VP_model), pointer :: model
    allocate(model, source=power_law_VP_model(A, N, Q, R))
  end function new_power_law_VP_model


  function strain_rate (this, stress, temp, dt)

    class(power_law_VP_model), intent(in) :: this
    real(r8), intent(in) :: stress  ! von Mises stress
    real(r8), intent(in) :: temp    ! temperature [K]
    real(r8), intent(in) :: dt      ! time step size
    real(r8) :: strain_rate

    real(r8) :: epsdot_min
    
    strain_rate = this%A * stress**this%N * exp(- this%Q / this%R / temp)

    strain_rate = min(1.0d7, strain_rate)

    epsdot_min = 1.0e-20 / dt

    if (strain_rate < epsdot_min) strain_rate = 0.0

  end function strain_rate

end module power_law_VP_model_type
