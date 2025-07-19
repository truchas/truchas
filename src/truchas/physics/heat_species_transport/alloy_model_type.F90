!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module alloy_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use alloy_vector_type
  use multicomp_lever_type
  use mfd_disc_type
  use prop_mesh_func_type
  use scalar_mesh_multifunc_type
  use bndry_func1_class
  use bndry_func2_class
  use intfc_func2_class
!  use index_partitioning
  use TofH_type
  use truchas_timers
  use parameter_list_type
  implicit none
  private

  type, public :: alloy_model
    type(unstr_mesh), pointer :: mesh => null()
    type(mfd_disc),   pointer :: disc => null()
    !! Equation parameters
    type(multicomp_lever) :: alloy
    type(prop_mesh_func) :: conductivity ! thermal conductivity
    type(prop_mesh_func) :: H_of_T       ! enthalpy as a function of temperature
    type(TofH) :: T_of_H
    real(r8), allocatable :: q_adv(:) ! advective source
    type(scalar_mesh_multifunc), allocatable :: src ! external heat source
    !! Boundary condition data
    class(bndry_func1), allocatable :: bc_dir  ! Dirichlet
    class(bndry_func1), allocatable :: bc_flux ! simple flux
    class(bndry_func2), allocatable :: bc_vflux ! oriented flux
    class(bndry_func2), allocatable :: bc_htc  ! external HTC
    class(bndry_func2), allocatable :: bc_rad  ! simple radiation
    class(intfc_func2), allocatable :: ic_htc  ! internal HTC
    class(intfc_func2), allocatable :: ic_rad  ! internal gap radiation
    class(bndry_func2), allocatable :: evap_flux
  contains
    procedure :: init
    procedure :: init_vector
    procedure :: compute_f
    procedure :: set_alloy_adv_source
  end type alloy_model

contains

  subroutine init_vector(this, vec)
    class(alloy_model), intent(in) :: this
    type(alloy_vector), intent(out) :: vec
    call vec%init(this%mesh)
  end subroutine

  subroutine init(this, disc)
    class(alloy_model), intent(out), target :: this
    type(mfd_disc), intent(in), target :: disc
    this%disc => disc
    this%mesh => disc%mesh
    block
      !TODO: pass parameters
      !TODO: redesign the TofH type after redesign of prop_mesh_func type
      real(r8) :: eps, delta
      integer  :: max_try
      type(parameter_list) :: params
      call params%get('tofh-tol', eps, default=0.0_r8)
      call params%get('tofh-max-try', max_try, default=50)
      call params%get('tofh-delta', delta, default=1.0_r8)
      call this%T_of_H%init(this%H_of_T, eps=eps, max_try=max_try, delta=delta)
    end block
  end subroutine

  subroutine set_alloy_adv_source(this, q_adv)
    class(alloy_model), intent(inout) :: this
    real(r8), intent(in) :: q_adv(:)
    ASSERT(size(q_adv) == this%mesh%ncell)
    if (.not.allocated(this%q_adv)) allocate(this%q_adv(this%mesh%ncell))
    this%q_adv = q_adv
  end subroutine

  subroutine compute_f(this, t, u, udot, f)

    use mfd_disc_type

    class(alloy_model), intent(inout) :: this
    real(r8), intent(in) :: t
    type(alloy_vector), intent(inout) :: u, udot ! data is intent(in)
    type(alloy_vector), intent(inout) :: f       ! data is intent(out)
    target :: u

    integer :: j, n, n1, n2
    real(r8), pointer :: qrad(:)
    real(r8), dimension(this%mesh%ncell) :: value
    real(r8), allocatable :: Tdir(:)
    real(r8), pointer :: state(:,:)

    call start_timer('ht-function')

    !call this%mesh%cell_imap%gather_offp(u%lf)
    !call this%mesh%cell_imap%gather_offp(u%hc)
    !call this%mesh%cell_imap%gather_offp(u%tc)
    !call this%mesh%face_imap%gather_offp(u%tf)
    call u%gather_offp

    !TODO: The existing prop_mesh_func%compute_value function expects a rank-2
    !      state array. This is a workaround until prop_mesh_func is redesigned.
    state(1:this%mesh%ncell,1:1) => u%tc

    !! Residual of the liquid fraction relation
    call this%alloy%compute_r(u%hc, u%lf, f%lf)

    ! Residual of the algebraic enthalpy-temperature relation
    call this%alloy%compute_H(u%tc, u%lf, f%hc)
    f%hc = u%hc - f%hc

  !!!! RESIDUAL OF THE HEAT EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Overwrite the temperature on Dirichlet faces with the boundary
    !! data, saving the original values to restore them later.
    if (allocated(this%bc_dir)) then
      call this%bc_dir%compute(t)
      allocate(Tdir(size(this%bc_dir%index)))
      do j = 1, size(this%bc_dir%index)
        n = this%bc_dir%index(j)
        Tdir(j) = u%tf(n)
        u%tf(n) = this%bc_dir%value(j)
      end do
    end if

    !! Compute the generic heat equation residual.
    call this%conductivity%compute_value(state, value)
    call this%disc%apply_diff(value, u%tc, u%tf, f%tc, f%tf)
    call this%mesh%cell_imap%gather_offp(udot%hc)
    if (allocated(this%q_adv)) then
      f%tc = f%tc + this%mesh%volume*(udot%hc - this%q_adv)
    else
      f%tc = f%tc + this%mesh%volume*udot%hc
    end if

    !! Additional heat source
    if (allocated(this%src)) then
      call this%src%compute(t, u%tc)
      f%tc = f%tc - this%mesh%volume*this%src%value
    end if

    !! Overwrite face residuals with the Dirichlet BC residual and
    !! restore the face temperatures to their original input values.
    if (allocated(this%bc_dir)) then
      do j = 1, size(this%bc_dir%index)
        n = this%bc_dir%index(j)
        f%tf(n) = Tdir(j) - this%bc_dir%value(j)
        u%tf(n) = Tdir(j)
      end do
      deallocate(Tdir)
    end if

    !! Simple flux BC contribution.
    if (allocated(this%bc_flux)) then
      call this%bc_flux%compute(t)
      do j = 1, size(this%bc_flux%index)
        n = this%bc_flux%index(j)
        f%tf(n) = f%tf(n) + this%mesh%area(n) * this%bc_flux%value(j)
      end do
    end if

    !! Oriented flux BC contribution.
    if (allocated(this%bc_vflux)) then
      call this%bc_vflux%compute(t, u%tf)
      do j = 1, size(this%bc_vflux%index)
        n = this%bc_vflux%index(j)
        f%tf(n) = f%tf(n) + this%bc_vflux%value(j)
      end do
    end if

    !! External HTC flux contribution.
    if (allocated(this%bc_htc)) then
      call this%bc_htc%compute(t, u%tf)
      do j = 1, size(this%bc_htc%index)
        n = this%bc_htc%index(j)
        f%tf(n) = f%tf(n) + this%bc_htc%value(j)
      end do
    end if

    !! Ambient radiation BC flux contribution.
    if (allocated(this%bc_rad)) then
      call this%bc_rad%compute(t, u%tf)
      do j = 1, size(this%bc_rad%index)
        n = this%bc_rad%index(j)
        f%tf(n) = f%tf(n) + this%bc_rad%value(j)
      end do
    end if

    !! Experimental evaporation heat flux
    if (allocated(this%evap_flux)) then
      call this%evap_flux%compute_value(t, u%tf)
      do j = 1, size(this%evap_flux%index)
        n = this%evap_flux%index(j)
        f%tf(n) = f%tf(n) + this%mesh%area(n)*this%evap_flux%value(j)
      end do
    end if

    !! Internal HTC flux contribution.
    if (allocated(this%ic_htc)) then
      call this%ic_htc%compute(t, u%tf)
      do j = 1, size(this%ic_htc%index,2)
        n1 = this%ic_htc%index(1,j)
        n2 = this%ic_htc%index(2,j)
        f%tf(n1) = f%tf(n1) + this%ic_htc%value(j)
        f%tf(n2) = f%tf(n2) - this%ic_htc%value(j)
      end do
    end if

    !! Internal gap radiation contribution.
    if (allocated(this%ic_rad)) then
      call this%ic_rad%compute(t, u%tf)
      do j = 1, size(this%ic_rad%index,2)
        n1 = this%ic_rad%index(1,j)
        n2 = this%ic_rad%index(2,j)
        f%tf(n1) = f%tf(n1) + this%ic_rad%value(j)
        f%tf(n2) = f%tf(n2) - this%ic_rad%value(j)
      end do
    end if

    !TODO: is this necessary? Off-process values are not needed, but may be
    !      used in dummy vector operations, and we don't want fp exceptions.
    call f%gather_offp

    call stop_timer('ht-function')

  end subroutine compute_f

end module alloy_model_type
