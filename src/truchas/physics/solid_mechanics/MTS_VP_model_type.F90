!!
!! MTS_VP_MODEL_TYPE
!!
!! The implementation of the MTS viscoplastic model.  This wraps existing
!! code by Dave Korzekwa to compute the strain rate.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module MTS_VP_model_type

  use kinds, only: r8
  use VP_model_class
  use truchas_logging_services
  implicit none
  private
  
  type, extends(VP_model), public :: MTS_VP_model
    private
    real(r8) :: kb, mu_0, sig_a, d, temp_0, b, edot_0i, g_0i, q_i, p_i, sig_i
  contains
    procedure :: strain_rate
  end type MTS_VP_model
  
  public :: new_MTS_VP_model
  
contains
  
  function new_MTS_VP_model (kb, mu_0, sig_a, d, temp_0, b, edot_0i, g_0i, q_i, p_i, sig_i) result (model)
    real(r8), intent(in) :: kb, mu_0, sig_a, d, temp_0, b, edot_0i, g_0i, q_i, p_i, sig_i
    class(VP_model), pointer :: model
    allocate(model, source=MTS_VP_model(kb, mu_0, sig_a, d, temp_0, b, edot_0i, g_0i, q_i, p_i, sig_i))
  end function new_MTS_VP_model


  function strain_rate (this, stress, temp, dt)

    class(MTS_VP_model), intent(in) :: this
    real(r8), intent(in) :: stress  ! von Mises stress
    real(r8), intent(in) :: temp    ! temperature [K]
    real(r8), intent(in) :: dt      ! time step size
    real(r8) :: strain_rate
    
    real(r8) :: a, c, mu, epsdot_min, eps_sig_a, k, tmp

    !! Temperature dependent shear modulus
    mu = this%mu_0 - this%d/(exp(this%temp_0/temp) - 1.0)

    if (stress > (this%sig_a + this%sig_i*mu/this%mu_0)) then
      strain_rate = this%edot_0i
      return
    end if

    ! We need these terms in any case
    c = this%kb * temp / (mu * this%b**3 * this%g_0i)
    epsdot_min = 1.0e-20 / dt

    if (stress < this%sig_a) then
      !! Get strain rate for the athermal stress sig_a
      eps_sig_a = this%edot_0i / exp(1.0 /c)
      if (eps_sig_a < epsdot_min) then
        strain_rate = 0.0
        return
      end if
      !! Assume a sigma^5 relationship for strain rate
      k = log(eps_sig_a) - 5 * log(this%sig_a)
      k = exp(k)
      strain_rate = k * stress**5
    else

      !! Rate dependence of the flow (yield) strength
      a = this%mu_0/mu/this%sig_i
      tmp = 1.0-(a*(stress - this%sig_a))**this%p_i
      if (tmp < 0.0_r8) then
        call TLS_panic ('MTS_VP_MODEL:STRAIN_RATE: invalid exponentiation')
      else
        strain_rate = this%edot_0i / (exp(tmp**this%q_i /c))
      end if
    end if

    if (strain_rate < epsdot_min) strain_rate = 0.0_r8
    
  end function strain_rate

end module MTS_VP_model_type
