!!
!! FHT_MODEL_TYPE
!!
!! This module defines the discrete thermal transport model used by the
!! fluid/heat transport solver. It owns the thermal component, mimetic
!! diffusion discretization, optional view-factor radiation coupling, and the
!! helper used to invert enthalpy density to temperature. It refers to the
!! dynamic void masks owned and updated by fht_solver.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! Residual assembly gathers the off-process temperature state needed by the
!! mimetic discretization, imposes thermal Dirichlet face values while
!! evaluating temperature-dependent face terms, and restores the input vector
!! before returning.
!!

#include "f90_assert.fpp"

module fht_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfd_disc_type
  use unstr_mesh_type
  use thermal_component_type
  use thermal_component_factory, only: define_thermal_component
  use thermal_view_factor_coupling_type, only: thermal_view_factor_coupling
  use truchas_timers
  use matl_mesh_func_type
  use fht_vector_type
  use TofH_type
  implicit none
  private

  type, public :: fht_model
    type(mfd_disc) :: disc
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    !! Current void masks owned by fht_solver; valid after solver void-state update.
    logical, pointer :: void_cell(:) => null(), void_face(:) => null()
    real(r8) :: void_temp = 0.0_r8
    type(thermal_component) :: thermal
    type(thermal_view_factor_coupling) :: vf_rad
    type(TofH) :: T_of_H
  contains
    procedure :: init
    procedure :: residual
    procedure :: init_vector
    procedure :: set_initial_time
    procedure :: update_moving_vf
    procedure :: add_moving_vf_events
    procedure :: set_ext_enthalpy_rate
  end type fht_model

contains

  subroutine init(this, mesh, mmf, params, stat, errmsg)

    use parameter_list_type
    use thermal_bc_factory_type
    use thermal_source_factory_type

    !! TARGET is needed because T_of_H stores a persistent pointer to
    !! this%thermal%H_of_T.
    class(fht_model), intent(out), target :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(matl_mesh_func), intent(in), target :: mmf
    type(parameter_list), intent(inout), target :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(thermal_bc_factory)    :: tbc_fac
    type(thermal_source_factory) :: tsrc_fac
    type(parameter_list), pointer :: bc_params, source_params, radiation_params
    real(r8) :: stefan_boltzmann, absolute_zero, eps, delta
    integer :: max_try

    call params%get('stefan-boltzmann', stefan_boltzmann, stat, errmsg, default=5.67e-8_r8)
    if (stat /= 0) return
    call params%get('absolute-zero', absolute_zero, stat, errmsg, default=0.0_r8)
    if (stat /= 0) return
    call params%get('void-temperature', this%void_temp, stat, errmsg, default=0.0_r8)
    if (stat /= 0) return
    call params%get('tofh-tol', eps, stat, errmsg, default=0.0_r8)
    if (stat /= 0) return
    call params%get('tofh-max-try', max_try, stat, errmsg, default=50)
    if (stat /= 0) return
    call params%get('tofh-delta', delta, stat, errmsg, default=1.0e-3_r8)
    if (stat /= 0) return

    bc_params => params%sublist('bc')
    source_params => params%sublist('sources')
    radiation_params => params%sublist('radiation')
    call tbc_fac%init(mesh, stefan_boltzmann, absolute_zero, bc_params)
    call tsrc_fac%init(mesh, source_params)

    call define_thermal_component(mesh, mmf, tbc_fac, tsrc_fac, this%thermal, stat, errmsg)
    if (stat /= 0) return

    call this%vf_rad%init(mesh, radiation_params, stat, errmsg)
    if (stat /= 0) return

    call this%vf_rad%validate_bc(mesh, this%thermal, stat, errmsg)
    if (stat /= 0) return

    call this%T_of_H%init(this%thermal%H_of_T, eps=eps, max_try=max_try, delta=delta)

    call this%disc%init(mesh)
    this%mesh => mesh

    stat = 0

  end subroutine init

  subroutine set_ext_enthalpy_rate(this, enthalpy_rate)
    class(fht_model), intent(inout) :: this
    real(r8), intent(in) :: enthalpy_rate(:)
    ASSERT(size(enthalpy_rate) == this%mesh%ncell_onP)
    this%thermal%ext_rate = enthalpy_rate
  end subroutine

  subroutine init_vector(this, vec)
    class(fht_model), intent(in) :: this
    type(fht_vector), intent(out) :: vec

    if (this%vf_rad%is_active()) then
      call vec%init(this%mesh, this%vf_rad%size())
    else
      call vec%init(this%mesh)
    end if
  end subroutine init_vector

  subroutine residual(this, t, u, hdot, r)

    class(fht_model), intent(inout) :: this
    real(r8), intent(in) :: t, hdot(:)
    type(fht_vector), intent(inout) :: u ! data is intent(in)
    type(fht_vector), intent(inout) :: r ! data is intent(out)

    integer :: ncell_onP, nface_onP
    real(r8), dimension(this%mesh%ncell) :: Rcell, value
    real(r8), dimension(this%mesh%nface) :: Rface
    real(r8), allocatable :: Tdir(:)

    ASSERT(associated(this%void_cell))
    ASSERT(associated(this%void_face))

    call start_timer('residual')

    ncell_onP = this%mesh%ncell_onP
    nface_onP = this%mesh%nface_onP

    !! IDAESOL callbacks require only on-process entries to be valid on entry.
    !! Gather all temperature fields used by the residual here; HDOT is only
    !! used on-process below.
    call this%mesh%cell_imap%gather_offp(u%tc)
    call this%mesh%face_imap%gather_offp(u%tf)

    call r%setval(0.0_r8)

    !! Thermal Dirichlet face values are imposed for thermal residual assembly
    !! because flux-type boundary conditions and view-factor radiation may
    !! depend on face temperature.
    call impose_thermal_dirichlet

    call compute_thermal_residual

    !! Restore the input vector before returning; vector data is conceptually
    !! input here despite the intent(inout) required by the vector type.
    call restore_thermal_dirichlet

    r%tc(:ncell_onP) = Rcell(:ncell_onP)
    r%tf(:nface_onP) = Rface(:nface_onP)

    call stop_timer('residual')

  contains

    subroutine impose_thermal_dirichlet
      integer :: j, n
      if (allocated(this%thermal%bc_dir)) then
        call this%thermal%bc_dir%compute(t)
        allocate(Tdir(size(this%thermal%bc_dir%index)))
        do j = 1, size(this%thermal%bc_dir%index)
          n = this%thermal%bc_dir%index(j)
          Tdir(j) = u%tf(n)
          u%tf(n) = this%thermal%bc_dir%value(j)
        end do
      end if
    end subroutine

    subroutine restore_thermal_dirichlet
      integer :: j, n
      if (allocated(this%thermal%bc_dir)) then
        do j = 1, size(this%thermal%bc_dir%index)
          n = this%thermal%bc_dir%index(j)
          u%tf(n) = Tdir(j)
        end do
        deallocate(Tdir)
      end if
    end subroutine

    subroutine compute_thermal_residual
      integer :: j, n

      !! Compute the generic heat equation residual.
      call this%thermal%conductivity%compute_value(u%tc, value)
      where (this%void_cell) value = 0.0_r8
      call this%disc%apply_diff(value, u%tc, u%tf, Rcell, Rface)

      !! Time derivative and external enthalpy rate contributions.
      Rcell(:ncell_onP) = Rcell(:ncell_onP) &
                        + this%mesh%volume(:ncell_onP)*hdot(:ncell_onP)

      !! External rate contribution; typically from advection.
      if (allocated(this%thermal%ext_rate)) then
        Rcell(:ncell_onP) = Rcell(:ncell_onP) - this%thermal%ext_rate(:ncell_onP)
      end if

      !! User-defined heat source contribution.
      if (allocated(this%thermal%src)) then
        call this%thermal%src%compute(t, u%tc)
        Rcell(:ncell_onP) = Rcell(:ncell_onP) &
                          - this%mesh%volume(:ncell_onP)*this%thermal%src%value(:ncell_onP)
      end if

      !! Overwrite face residuals with the Dirichlet BC residual.
      if (allocated(this%thermal%bc_dir)) then
        do j = 1, size(this%thermal%bc_dir%index)
          n = this%thermal%bc_dir%index(j)
          Rface(n) = Tdir(j) - this%thermal%bc_dir%value(j)
        end do
      end if

      !! Flux-type thermal BC contributions.
      call this%thermal%add_flux_bc_residual(t, u%tf, Rface, this%mesh%area, this%void_face)

      !! Overwrite residual on void cells and faces with dummy equation T=T_void.
      where (this%void_cell(:ncell_onP)) Rcell(:ncell_onP) = u%tc(:ncell_onP) - this%void_temp
      where (this%void_face(:nface_onP)) Rface(:nface_onP) = u%tf(:nface_onP) - this%void_temp

      if (this%vf_rad%is_active()) then
        !! Radiative heat flux contribution to the heat conduction face residual.
        call this%vf_rad%add_heat_flux_residual(t, u%tf, u%qrad, this%mesh%area, Rface)
        !! Residual of the algebraic radiosity system.
        call this%vf_rad%compute_radiosity_residual(t, u%tf, u%qrad, r%qrad)
        r%qrad(:) = -r%qrad(:)
      end if

    end subroutine compute_thermal_residual

  end subroutine residual

  subroutine set_initial_time(this, t)
    class(fht_model), intent(inout) :: this
    real(r8), intent(in) :: t
    call this%vf_rad%set_initial_time(t)
  end subroutine set_initial_time

  subroutine update_moving_vf(this)
    class(fht_model), intent(inout) :: this
    call this%vf_rad%update_moving_vf
  end subroutine

  subroutine add_moving_vf_events(this, eventq, rank)
    use sim_event_queue_type
    class(fht_model), intent(in) :: this
    type(sim_event_queue), intent(inout) :: eventq
    integer, intent(in), optional :: rank
    call this%vf_rad%add_moving_vf_events(eventq, rank)
  end subroutine

end module fht_model_type
