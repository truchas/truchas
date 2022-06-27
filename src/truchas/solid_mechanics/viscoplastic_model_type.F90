!!
!! Zach Jibben <zjibben@lanl.gov>
!! May 2022
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module viscoplastic_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use sm_material_model_type
  implicit none
  private

  type, public :: viscoplastic_model
    private
    ! unowned references used to avoid copies
    real(r8), dimension(:), pointer :: vof, strain_total_old, strain_total_new, &
        strain_thermal_old, strain_thermal_new, strain_pc
    ! real(r8), allocatable :: vof(:)
    ! real(r8), dimension(6) :: strain_total_old, strain_total_new, &
    !     strain_thermal_old, strain_thermal_new, strain_pc
    real(r8) :: temperature
    real(r8), public :: lame1, lame2 ! expose to preconditioner
    real(r8) :: dt

    type(sm_material_model), pointer :: matl_model => null() ! unowned reference
  contains
    procedure :: init
    procedure :: set_state
    procedure :: compute_udot
    procedure :: strain_rate
    procedure :: strain_rate2
    procedure :: is_elastic
    procedure :: compute_stresses
    procedure :: g
    procedure, private :: g_precon
    procedure :: compute_precon
  end type viscoplastic_model

contains

  subroutine init(this, matl_model)
    class(viscoplastic_model), intent(out) :: this
    type(sm_material_model), intent(in), target :: matl_model
    this%matl_model => matl_model
  end subroutine init


  subroutine set_state(this, dt, temperature, lame1, lame2, vof, strain_pc, &
      strain_total_old, strain_thermal_old, strain_total_new, strain_thermal_new)

    class(viscoplastic_model), intent(inout) :: this
    real(r8), intent(in) :: dt, temperature, lame1, lame2
    real(r8), intent(in), target :: vof(:), strain_pc(:), &
        strain_total_old(:), strain_thermal_old(:), strain_total_new(:), strain_thermal_new(:)

    integer :: m

    ASSERT(6 == size(strain_total_old))
    ASSERT(6 == size(strain_total_new))
    ASSERT(6 == size(strain_thermal_old))
    ASSERT(6 == size(strain_thermal_new))
    ASSERT(6 == size(strain_pc))
    ASSERT(all(vof >= 0))

    this%dt = dt
    this%temperature = temperature
    this%lame1 = lame1
    this%lame2 = lame2

    ! this%vof = vof
    ! this%strain_total_old = strain_total_old
    ! this%strain_total_new = strain_total_new
    ! this%strain_thermal_old = strain_thermal_old
    ! this%strain_thermal_new = strain_thermal_new
    ! this%strain_pc = strain_pc

    this%vof => vof
    this%strain_pc => strain_pc
    this%strain_total_old => strain_total_old
    this%strain_total_new => strain_total_new
    this%strain_thermal_old => strain_thermal_old
    this%strain_thermal_new => strain_thermal_new

    do m = 1, this%matl_model%nmat
      !if (.not.allocated(this%matl_model%vp(m)%m) .or. this%vof(m) == 0) cycle
      if (this%vof(m) == 0) cycle
      call this%matl_model%vp(m)%m%set_temperature(this%temperature)
    end do

  end subroutine set_state


  ! Returns whether the initialized state is elastic
  pure logical function is_elastic(this)
    class(viscoplastic_model), intent(in) :: this
    integer :: m
    is_elastic = .true.
    do m = 1, this%matl_model%nmat
      !if (.not.allocated(this%matl_model%vp(m)%m)) cycle
      is_elastic = is_elastic .and. (this%vof(m) == 0 .or. this%matl_model%vp(m)%m%is_elastic)
    end do
  end function is_elastic


  !! Expose a method to compute the time derivative of the plastic strain.
  !! This is used to save and set the "initial" udot at the next time step.
  !!
  !! du_i / dt = F_i = 3/2 * G * deviatoric_stress_i / vmstress
  !! G = sum_m { vof_m * f_m(vmstress, temperature) }
  pure subroutine compute_udot(this, t, u, udot)
    class(viscoplastic_model), intent(in) :: this
    real(r8), intent(in) :: t, u(:)
    real(r8), intent(out) :: udot(:)
    real(r8) :: stress(6), vmstress, dstress(6)
    call this%compute_stresses(t, u, stress, vmstress, dstress)
    !if (vmstress /= 0) then
    udot = (1.5_r8 * this%g(vmstress) / vmstress) * dstress
    ! else
    !   ! If the Von Mises stress is 0, set udot to 0.
    !   udot = 0
    ! end if
  end subroutine compute_udot


  real(r8) pure function strain_rate(this, t, u)
    class(viscoplastic_model), intent(in) :: this
    real(r8), intent(in) :: t, u(:)
    real(r8) :: stress(6), vmstress, dstress(6)
    call this%compute_stresses(t, u, stress, vmstress, dstress)
    strain_rate = this%g(vmstress)
  end function strain_rate

  !! A routine for computing strain rate from a fresh state for visualization
  function strain_rate2(this, temperature, lame1, lame2, vof, stress) result(strain_rate)
    use sm_bc_utilities, only: von_mises_stress
    class(viscoplastic_model), intent(inout) :: this
    real(r8), intent(in) :: temperature, lame1, lame2, stress(:)
    real(r8), intent(in), target :: vof(:)
    real(r8) :: strain_rate
    real(r8), parameter :: zeros(6) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8]
    call this%set_state(0.0_r8, temperature, lame1, lame2, vof, zeros, zeros, zeros, zeros, zeros)
    strain_rate = this%g(von_mises_stress(stress))
  end function strain_rate2


  !! G = sum_m { vof_m * f_m(vmstress, temperature) }
  real(r8) pure function g(this, vmstress)
    class(viscoplastic_model), intent(in) :: this
    real(r8), intent(in) :: vmstress
    integer :: m
    ! ASSERT(associated(this%matl_model))
    ! ASSERT(allocated(this%matl_model%vp))
    g = 0
    do m = 1, this%matl_model%nmat
      !if (.not.allocated(this%matl_model%vp(m)%m) .or. this%vof(m) == 0) cycle
      if (this%vof(m) == 0) cycle
      g = g + this%vof(m) * this%matl_model%vp(m)%m%strain_rate(vmstress)
    end do
  end function g

  !! G = sum_m { vof_m * f_m(vmstress, temperature) }
  !! compute dG / dvmstress
  real(r8) pure function g_precon(this, vmstress)
    class(viscoplastic_model), intent(in) :: this
    real(r8), intent(in) :: vmstress
    integer :: m
    g_precon = 0
    do m = 1, this%matl_model%nmat
      !if (.not.allocated(this%matl_model%vp(m)%m) .or. this%vof(m) == 0) cycle
      if (this%vof(m) == 0) cycle
      g_precon = g_precon + this%vof(m) * this%matl_model%vp(m)%m%strain_rate_precon(vmstress, this%temperature)
    end do
  end function g_precon


  !! Compute the stress, von mises stress (effective stress), and deviatoric stress.
  pure subroutine compute_stresses(this, t, strain_plastic, stress, vmstress, dstress)

    use sm_bc_utilities, only: compute_stress, von_mises_stress, compute_deviatoric_stress

    class(viscoplastic_model), intent(in) :: this
    real(r8), intent(in) :: t, strain_plastic(:)
    real(r8), intent(out) :: stress(:), vmstress, dstress(:)

    real(r8) :: a
    real(r8), dimension(6) :: strain_total, strain_thermal, strain_elastic

    a = merge(t / this%dt, 0.0_r8, this%dt > 0) ! dt might be zero during initialization
    strain_total = a * this%strain_total_new + (1 - a) * this%strain_total_old
    strain_thermal = a * this%strain_thermal_new + (1 - a) * this%strain_thermal_old
    strain_elastic = strain_total - strain_plastic - strain_thermal - this%strain_pc

    call compute_stress(this%lame1, this%lame2, strain_elastic, stress)

    vmstress = von_mises_stress(stress)
    call compute_deviatoric_stress(stress, dstress)

  end subroutine compute_stresses


  !! Compute the Jacobian (derivative of plastic strain rate wrt plastic strain).
  pure subroutine compute_precon(this, t, u, precon)

    class(viscoplastic_model), intent(in) :: this
    real(r8), intent(in) :: t, u(:)
    real(r8), intent(out) :: precon(:,:)

    integer :: i, j
    real(r8) :: stress(6), vmstress, dstress(6), strain_rate
    real(r8) :: stress_precon(6,6), deviatoric_precon(6,6)
    real(r8) :: von_mises_precon(6), strain_rate_precon(6)

    call this%compute_stresses(t, u, stress, vmstress, dstress)
    strain_rate = this%g(vmstress)

    ! dstress_i / dstrain_j
    stress_precon = 0
    do i = 1, 6
      stress_precon(i,i) = 2 * this%lame2
    end do
    stress_precon(:3,:3) = stress_precon(:3,:3) + this%lame1

    ! dvmstress / dstrain_i = (dvmstress / dstress_k) * (dstress_k / dstrain_i)
    von_mises_precon(1) = (stress(3) - stress(2)) / (2*vmstress)
    von_mises_precon(2) = (stress(1) - stress(3)) / (2*vmstress)
    von_mises_precon(3) = (stress(2) - stress(1)) / (2*vmstress)
    von_mises_precon(4:) = 3 * stress(4:) / vmstress
    von_mises_precon = matmul(transpose(stress_precon), von_mises_precon)

    ! ddeviatoric_stress_i / dstrain_j
    deviatoric_precon = stress_precon
    do i = 1, 3
      deviatoric_precon(i,:) = deviatoric_precon(i,:) - sum(stress_precon(:3,:), dim=1) / 3
    end do

    ! dG / dstrain_i = (dG / dvmstress) * (dvmstress / dstrain_i)
    strain_rate_precon = this%g_precon(vmstress) * von_mises_precon

    ! u: plastic strain
    ! du_i / dt = F_i = 3/2 * G * deviatoric_stress_i / vmstress
    ! G = sum_m { vof_m * f_m(vmstress, temperature) }
    ! residual: R_i = du_i/dt - udot_i = F_i - (u_i - uprev) / dt
    ! (NOTE: the udot part of the residual is left to the jacobian idaesol model)
    !
    ! dR_i / du_j = dF_i / du_j - dudot_i / du_j
    ! udot = (u - uprev) / dt  -->  dudot_i / du_j = delta_ij / dt
    ! dF_i/du_j = 3/2 * (dG/du_j * deviatoric_stress_i / vmstress
    !                   + G * ddeviatoric_stress_i/du_j / vmstress
    !                   - G * deviatoric_stress_i * dvmstress/du_j / vmstress**2)
    do j = 1, 6
      do i = 1, 6
        precon(i,j) = 1.5_r8 * (strain_rate_precon(j) * dstress(i) / vmstress &
            + strain_rate * deviatoric_precon(i,j) / vmstress &
            - strain_rate * dstress(i) * von_mises_precon(j) / vmstress**2)
      end do
    end do

  end subroutine compute_precon

end module viscoplastic_model_type
