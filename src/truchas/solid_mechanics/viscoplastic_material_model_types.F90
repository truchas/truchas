!!
!! Zach Jibben <zjibben@lanl.gov>
!! October 2021
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module viscoplastic_material_model_types

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use truchas_logging_services
  !use ieee_arithmetic, only: ieee_is_normal
  implicit none
  private

  public :: alloc_viscoplastic_material_model

  type, abstract, public :: viscoplastic_material_model
    logical, public :: is_elastic = .false.
  contains
    procedure(init), deferred :: init
    procedure(strain_rate_st), deferred :: strain_rate_st
    procedure(strain_rate_s), deferred :: strain_rate_s
    procedure(strain_rate_st), deferred :: strain_rate_precon
    procedure(set_temperature), deferred :: set_temperature
    generic :: strain_rate => strain_rate_st, strain_rate_s
  end type viscoplastic_material_model

  abstract interface
    subroutine init(this, params)
      import viscoplastic_material_model, parameter_list
      class(viscoplastic_material_model), intent(out) :: this
      type(parameter_list), intent(inout) :: params
    end subroutine
    pure function strain_rate_st(this, stress, temperature)
      import viscoplastic_material_model, r8
      class(viscoplastic_material_model), intent(in) :: this
      real(r8), intent(in) :: stress, temperature
      real(r8) :: strain_rate_st
    end function
    pure subroutine set_temperature(this, temperature)
      import viscoplastic_material_model, r8
      class(viscoplastic_material_model), intent(inout) :: this
      real(r8), intent(in) :: temperature
    end subroutine
    pure function strain_rate_s(this, stress)
      import viscoplastic_material_model, r8
      class(viscoplastic_material_model), intent(in) :: this
      real(r8), intent(in) :: stress
      real(r8) :: strain_rate_s
    end function
  end interface


  type, extends(viscoplastic_material_model), public :: elastic_model
  contains
    procedure :: init => init_el
    procedure :: strain_rate_st => strain_rate_el_st
    procedure :: strain_rate_s => strain_rate_el_s
    procedure :: strain_rate_precon => strain_rate_precon_el
    procedure :: set_temperature => set_temperature_el
  end type elastic_model

  type, extends(viscoplastic_material_model), public :: viscoplastic_power_law_model
    private
    real(r8) :: A, n, Q, R
    real(r8) :: A2
  contains
    procedure :: init => init_pl
    procedure :: strain_rate_st => strain_rate_pl_st
    procedure :: strain_rate_s => strain_rate_pl_s
    procedure :: strain_rate_precon => strain_rate_precon_pl
    procedure :: set_temperature => set_temperature_pl
  end type viscoplastic_power_law_model

  type, extends(viscoplastic_material_model), public :: viscoplastic_mts_model
    private
    real(r8) :: mu0, D, b, k, T0, sigma_i, sigma_a, epsdot0i, g0i, pi, qi

    real(r8) :: mu, tmp0, tmp1, K1
    real(r8) :: b3_g0i_k, sigma_i_mu0, sigma_a_5
  contains
    procedure :: init => init_mts
    procedure :: strain_rate_st => strain_rate_mts_st
    procedure :: strain_rate_s => strain_rate_mts_s
    procedure :: strain_rate_precon => strain_rate_precon_mts
    procedure :: set_temperature => set_temperature_mts
  end type viscoplastic_mts_model

  type, public :: viscoplastic_material_model_box
    class(viscoplastic_material_model), allocatable :: m
  end type viscoplastic_material_model_box

contains

  subroutine alloc_viscoplastic_material_model(plist, m)

    type(parameter_list), intent(inout) :: plist
    class(viscoplastic_material_model), intent(out), allocatable :: m

    character(:), allocatable :: model

    call plist%get('model', model, default='elastic')
    select case (model)
    case ('power law')
      allocate(viscoplastic_power_law_model :: m)
    case ('MTS')
      allocate(viscoplastic_mts_model :: m)
    case ('elastic')
      allocate(elastic_model :: m)
    case default
      INSIST(.false.) ! TODO: proper error message
    end select
    if (allocated(m)) call m%init(plist)

  end subroutine alloc_viscoplastic_material_model

  ! Power-Law strain rate model
  subroutine init_pl(this, params)
    class(viscoplastic_power_law_model), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    call params%get('A', this%A)
    call params%get('n', this%n)
    call params%get('Q', this%Q)
    call params%get('R', this%R)

    if (this%a < 0) call tls_fatal("pwr_law_a must be >= 0")
    if (this%n < 0) call tls_fatal("pwr_law_n must be >= 0")
    if (this%r == 0) call tls_fatal("pwr_law_r must be /= 0")
  end subroutine init_pl

  pure function strain_rate_pl_st(this, stress, temperature) result(strain_rate)
    class(viscoplastic_power_law_model), intent(in) :: this
    real(r8), intent(in) :: stress, temperature
    real(r8) :: strain_rate
    strain_rate = this%A * stress**this%n * exp(-this%Q / (this%R * temperature))
  end function strain_rate_pl_st


  pure function strain_rate_precon_pl(this, stress, temperature) result(strain_rate_precon)
    class(viscoplastic_power_law_model), intent(in) :: this
    real(r8), intent(in) :: stress, temperature
    real(r8) :: strain_rate_precon
    if (this%n == 0) then
      strain_rate_precon = 0
    else
      strain_rate_precon = this%A * exp(-this%Q / (this%R * temperature)) * (this%n * stress**(this%n-1))
    end if
  end function strain_rate_precon_pl


  pure subroutine set_temperature_pl(this, temperature)
    class(viscoplastic_power_law_model), intent(inout) :: this
    real(r8), intent(in) :: temperature
    this%A2 = this%A * exp(-this%Q / (this%R * temperature))
  end subroutine set_temperature_pl

  pure function strain_rate_pl_s(this, stress) result(strain_rate)
    class(viscoplastic_power_law_model), intent(in) :: this
    real(r8), intent(in) :: stress
    real(r8) :: strain_rate
    strain_rate = this%A2 * stress**this%n
  end function strain_rate_pl_s


  !!!!! Elastic "clear" model !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_el(this, params)
    class(elastic_model), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    this%is_elastic = .true.
  end subroutine
  pure function strain_rate_el_st(this, stress, temperature) result(strain_rate)
    class(elastic_model), intent(in) :: this
    real(r8), intent(in) :: stress, temperature
    real(r8) :: strain_rate
    strain_rate = 0
  end function
  pure function strain_rate_precon_el(this, stress, temperature) result(strain_rate_precon)
    class(elastic_model), intent(in) :: this
    real(r8), intent(in) :: stress, temperature
    real(r8) :: strain_rate_precon
    strain_rate_precon = 0
  end function
  pure subroutine set_temperature_el(this, temperature)
    class(elastic_model), intent(inout) :: this
    real(r8), intent(in) :: temperature
  end subroutine
  pure function strain_rate_el_s(this, stress) result(strain_rate)
    class(elastic_model), intent(in) :: this
    real(r8), intent(in) :: stress
    real(r8) :: strain_rate
    strain_rate = 0
  end function


  ! MTS (mechanical threshold stress) strain rate model
  subroutine init_mts(this, params)
    class(viscoplastic_mts_model), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    call params%get('mu0', this%mu0)
    call params%get('D', this%D)
    call params%get('b', this%b)
    call params%get('k', this%k)
    call params%get('T0', this%T0)
    call params%get('sigma_i', this%sigma_i)
    call params%get('sigma_a', this%sigma_a)
    call params%get('epsdot0i', this%epsdot0i)
    call params%get('g0i', this%g0i)
    call params%get('pi', this%pi)
    call params%get('qi', this%qi)

    if (this%b <= 0) call tls_fatal("MTS_B must be > 0")
    if (this%epsdot0i < 0) call tls_fatal("MTS_epsdot0i must be >= 0")
    if (this%g0i <= 0) call tls_fatal("MTS_g_0i must be > 0")
    if (this%k <= 0) call tls_fatal("MTS_k must be > 0")
    if (this%mu0 <= 0) call tls_fatal("MTS_mu0 must be > 0")
    if (this%pi <= 0) call tls_fatal("MTS_pi must be > 0")
    if (this%qi <= 0) call tls_fatal("MTS_qi must be > 0")
    if (this%sigma_a < 0) call tls_fatal("MTS_sigma_a must be >= 0")
    if (this%sigma_i <= 0) call tls_fatal("MTS_sigma_i must be > 0")
    if (this%T0 <= 0) call tls_fatal("MTS_T0 must be > 0")

    this%b3_g0i_k = this%b**3 * this%g0i / this%k
    this%sigma_a_5 = this%sigma_a**5
    this%sigma_i_mu0 = this%sigma_i / this%mu0
  end subroutine init_mts

  pure function strain_rate_mts_st(this, stress, temperature) result(strain_rate)

    class(viscoplastic_mts_model), intent(in) :: this
    real(r8), intent(in) :: stress, temperature
    real(r8) :: strain_rate

    real(r8) :: mu, K1, tmp0, tmp1, tmp2

    mu = this%mu0 - this%D / (exp(this%T0/temperature) - 1)
    tmp0 = mu * this%sigma_i_mu0
    tmp1 = mu * this%b3_g0i_k / temperature

    if (stress < this%sigma_a) then
      ! K1 chosen for continuity with the bottom eq @ stress = sigma_a
      K1 = this%epsdot0i * exp(-tmp1) / this%sigma_a_5
      strain_rate = K1 * stress**5
    else if (stress > tmp0 + this%sigma_a) then
      strain_rate = this%epsdot0i
    else
      tmp2 = (stress - this%sigma_a) / tmp0
      strain_rate = this%epsdot0i * exp(-tmp1 * (1 - tmp2**this%pi)**this%qi)
    end if

    !ASSERT(ieee_is_normal(strain_rate))

  end function strain_rate_mts_st


  pure subroutine set_temperature_mts(this, temperature)
    class(viscoplastic_mts_model), intent(inout) :: this
    real(r8), intent(in) :: temperature
    ! K1 chosen for continuity with the bottom eq @ stress = sigma_a
    this%mu = this%mu0 - this%D / (exp(this%T0/temperature) - 1)
    this%tmp0 = this%mu * this%sigma_i_mu0
    this%tmp1 = this%mu * this%b3_g0i_k / temperature
    !this%eps_sig_a = this%epsdot0i * exp(-this%tmp1)
    this%K1 = this%epsdot0i * exp(-this%tmp1) / this%sigma_a_5
  end subroutine set_temperature_mts

  pure function strain_rate_mts_s(this, stress) result(strain_rate)

    class(viscoplastic_mts_model), intent(in) :: this
    real(r8), intent(in) :: stress
    real(r8) :: strain_rate

    real(r8) :: tmp2

    if (stress < this%sigma_a) then
      strain_rate = this%K1 * stress**5
    else if (stress > this%tmp0 + this%sigma_a) then
      strain_rate = this%epsdot0i
    else
      tmp2 = (stress - this%sigma_a) / this%tmp0
      strain_rate = this%epsdot0i * exp(-this%tmp1 * (1 - tmp2**this%pi)**this%qi)
    end if

    !ASSERT(ieee_is_normal(strain_rate))

  end function strain_rate_mts_s


  pure function strain_rate_precon_mts(this, stress, temperature) result(strain_rate_precon)

    class(viscoplastic_mts_model), intent(in) :: this
    real(r8), intent(in) :: stress, temperature
    real(r8) :: strain_rate_precon

    real(r8) :: mu, K1, tmp0, tmp1, tmp2

    mu = this%mu0 - this%D / (exp(this%T0/temperature) - 1)
    tmp0 = mu * this%sigma_i / this%mu0
    tmp1 = mu * this%b**3 * this%g0i / (this%k*temperature)

    if (stress < this%sigma_a) then
      ! K1 chosen for continuity with the bottom eq @ stress = sigma_a
      K1 = this%epsdot0i * exp(-tmp1) / this%sigma_a**5
      strain_rate_precon = K1 * 5*stress**4
    else if (stress - this%sigma_a > tmp0) then
      strain_rate_precon = 0
    else
      tmp2 = (stress - this%sigma_a) / tmp0
      strain_rate_precon = this%epsdot0i * exp(-tmp1 * (1 - tmp2**this%pi)**this%qi) &
          * (-tmp1 * this%qi * (1 - tmp2**this%pi)**(this%qi-1) * (-this%pi * tmp2**(this%pi-1)) / tmp0)
    end if

    !ASSERT(ieee_is_normal(strain_rate_precon))

  end function strain_rate_precon_mts

end module viscoplastic_material_model_types
