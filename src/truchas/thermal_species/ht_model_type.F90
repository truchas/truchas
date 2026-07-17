!!
!! HT_MODEL_TYPE
!!
!! This module defines the discrete thermal transport model used by the
!! implicit integrator. It owns the thermal component, mimetic diffusion
!! discretization, static void masks, optional view-factor radiation coupling
!! data, and the helper used to invert enthalpy density to temperature during
!! initial-condition setup.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! Residual assembly gathers the off-process state needed by the mimetic
!! discretization, imposes thermal Dirichlet face values while evaluating
!! temperature-dependent face terms, and restores the input vector before
!! returning.
!!

#include "f90_assert.fpp"

module ht_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use ht_vector_type
  use mfd_disc_type
  use thermal_component_type
  use thermal_component_factory, only: define_thermal_component
  use thermal_view_factor_coupling_type, only: thermal_view_factor_coupling
  use TofH_type
  use parameter_list_type
  use matl_mesh_func_type
  use thermal_bc_factory_type
  use thermal_source_factory_type
  implicit none
  private

  type, public :: ht_model
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    type(mfd_disc) :: disc
    logical, allocatable :: void_cell(:), void_face(:)
    integer, allocatable :: void_dir_faces(:)
    real(r8) :: void_temp = 0.0_r8
    type(thermal_component) :: thermal
    type(thermal_view_factor_coupling) :: vf_rad
    type(TofH) :: T_of_H
  contains
    procedure :: init
    procedure :: init_vector
    procedure :: define_void_masks
    procedure :: set_initial_time
    procedure :: residual
    procedure :: update_moving_vf
    procedure :: add_moving_vf_events
    procedure :: set_ext_enthalpy_rate
  end type ht_model

contains

  subroutine init(this, mesh, mmf, params, stat, errmsg)

    !! TARGET is needed because T_of_H stores a persistent pointer to
    !! this%thermal%H_of_T.
    class(ht_model), intent(out), target :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(matl_mesh_func), intent(in), target :: mmf
    type(parameter_list), intent(inout), target :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(thermal_bc_factory) :: bc_fac
    type(thermal_source_factory) :: src_fac
    type(parameter_list), pointer :: sublist
    real(r8) :: stefan_boltzmann, absolute_zero, eps, delta
    integer :: max_try

    call this%disc%init(mesh)
    this%mesh => mesh

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
    call params%get('tofh-delta', delta, stat, errmsg, default=1.0_r8)
    if (stat /= 0) return

    sublist => params%sublist('bc')
    call bc_fac%init(mesh, stefan_boltzmann, absolute_zero, sublist)

    sublist => params%sublist('sources')
    call src_fac%init(mesh, sublist)

    call define_thermal_component(mesh, mmf, bc_fac, src_fac, this%thermal, stat, errmsg)
    if (stat /= 0) return

    sublist => params%sublist('radiation')
    call this%vf_rad%init(mesh, sublist, stat, errmsg)
    if (stat /= 0) return

    call this%vf_rad%validate_bc(mesh, this%thermal, stat, errmsg)
    if (stat /= 0) return

    call this%T_of_H%init(this%thermal%H_of_T, eps=eps, max_try=max_try, delta=delta)

    stat = 0

  end subroutine

  subroutine init_vector(this, vec)
    class(ht_model), intent(in) :: this
    type(ht_vector), intent(out) :: vec
    if (this%vf_rad%is_active()) then
      call vec%init(this%mesh, this%vf_rad%size())
    else
      call vec%init(this%mesh)
    end if
  end subroutine

  subroutine define_void_masks(this, mmf)

    class(ht_model), intent(inout) :: this
    type(matl_mesh_func), intent(in) :: mmf

    integer :: j, n

    if (allocated(this%void_dir_faces)) deallocate(this%void_dir_faces)
    call mmf%get_void_masks(this%mesh, this%void_cell, this%void_face)
    if (.not. allocated(this%void_cell)) return

    n = count(this%void_face)
    if (n > 0) then
      allocate(this%void_dir_faces(n))
      n = 0
      do j = 1, this%mesh%nface
        if (this%void_face(j)) then
          n = n + 1
          this%void_dir_faces(n) = j
        end if
      end do
    end if

  end subroutine define_void_masks

  subroutine set_initial_time(this, t)
    class(ht_model), intent(inout) :: this
    real(r8), intent(in) :: t
    call this%vf_rad%set_initial_time(t)
  end subroutine

  subroutine set_ext_enthalpy_rate(this, enthalpy_rate)
    class(ht_model), intent(inout) :: this
    real(r8), intent(in) :: enthalpy_rate(:)
    ASSERT(size(enthalpy_rate) == this%mesh%ncell_onP)
    this%thermal%ext_rate = enthalpy_rate
  end subroutine

  subroutine residual(this, t, u, udot, r)

    class(ht_model), intent(inout) :: this
    real(r8), intent(in) :: t
    type(ht_vector), intent(inout) :: u, udot ! data is intent(in)
    type(ht_vector), intent(inout) :: r       ! data is intent(out)

    real(r8), allocatable :: Tdir(:)

    !! IDAESOL callbacks require only on-process entries to be valid on entry.
    !! Gather all temperature fields used by the thermal residual here; UDOT
    !! is only used on-process below.
    call this%mesh%cell_imap%gather_offp(u%tc)
    call this%mesh%face_imap%gather_offp(u%tf)

    !! Thermal Dirichlet face values are imposed for thermal residual assembly
    !! because flux-type boundary conditions and view-factor radiation may
    !! depend on face temperature.
    call impose_thermal_dirichlet

    call compute_thermal_residual

    !! Restore the input vector before returning; vector data is conceptually
    !! input here despite the intent(inout) required by the vector type.
    call restore_thermal_dirichlet
    call overwrite_thermal_void_residual

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

    subroutine overwrite_thermal_void_residual
      associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
        if (allocated(this%void_cell)) &
            where (this%void_cell(:ncell_onP)) r%tc(:ncell_onP) = u%tc(:ncell_onP) - this%void_temp
        if (allocated(this%void_face)) &
            where (this%void_face(:nface_onP)) r%tf(:nface_onP) = u%tf(:nface_onP) - this%void_temp
      end associate
    end subroutine

    subroutine compute_thermal_residual

      integer :: j, n, ncell_onP
      real(r8), dimension(this%mesh%ncell) :: value

      ncell_onP = this%mesh%ncell_onP

      !! Residual of the algebraic enthalpy-temperature relation
      call this%thermal%H_of_T%compute_value(u%tc(:ncell_onP), value(:ncell_onP))
      r%hc(:ncell_onP) = u%hc(:ncell_onP) - value(:ncell_onP)

      !! Overwrite residual on void cells with dummy equation H=0.
      if (allocated(this%void_cell)) &
          where (this%void_cell(:ncell_onP)) r%hc(:ncell_onP) = u%hc(:ncell_onP)

      !! Heat equation residual

      !! Compute the generic heat equation residual.
      call this%thermal%conductivity%compute_value(u%tc, value)
      if (allocated(this%void_cell)) where (this%void_cell) value = 0.0_r8
      call this%disc%apply_diff(value, u%tc, u%tf, r%tc, r%tf)
      r%tc(:ncell_onP) = r%tc(:ncell_onP) &
          + this%mesh%volume(:ncell_onP)*udot%hc(:ncell_onP)

      !! External rate contribution; typically from advection or electromagnetics.
      if (allocated(this%thermal%ext_rate)) then
        r%tc(:ncell_onP) = r%tc(:ncell_onP) - this%thermal%ext_rate(:ncell_onP)
      end if

      !! User-defined heat source contribution.
      if (allocated(this%thermal%src)) then
        call this%thermal%src%compute(t, u%tc)
        r%tc(:ncell_onP) = r%tc(:ncell_onP) &
            - this%mesh%volume(:ncell_onP)*this%thermal%src%value(:ncell_onP)
      end if

      !! Overwrite face residuals with the Dirichlet BC residual.
      if (allocated(this%thermal%bc_dir)) then
        do j = 1, size(this%thermal%bc_dir%index)
          n = this%thermal%bc_dir%index(j)
          r%tf(n) = Tdir(j) - this%thermal%bc_dir%value(j)
        end do
      end if

      !! Flux-type thermal BC contributions.
      call this%thermal%add_flux_bc_residual(t, u%tf, r%tf, this%mesh%area, this%void_face)

      !! Residual of the algebraic enclosure radiation system.
      if (this%vf_rad%is_active()) then
        !! Radiative heat flux contribution to the heat conduction face residual.
        call this%vf_rad%add_heat_flux_residual(t, u%tf, u%qrad, this%mesh%area, r%tf)
        !! Residual of the algebraic radiosity system.
        call this%vf_rad%compute_radiosity_residual(t, u%tf, u%qrad, r%qrad)
        r%qrad(:) = -r%qrad(:)
      end if

    end subroutine compute_thermal_residual

  end subroutine residual

  subroutine update_moving_vf(this)
    class(ht_model), intent(inout) :: this
    call this%vf_rad%update_moving_vf
  end subroutine

  subroutine add_moving_vf_events(this, eventq, rank)
    use sim_event_queue_type
    class(ht_model), intent(in) :: this
    type(sim_event_queue), intent(inout) :: eventq
    integer, intent(in), optional :: rank
    call this%vf_rad%add_moving_vf_events(eventq, rank)
  end subroutine

end module ht_model_type
