!!
!! HTSD_MODEL_TYPE
!!
!! This module defines the discrete coupled thermal/species transport model
!! used by the implicit integrator. It owns the shared thermal component, the
!! species components, the mimetic diffusion discretization, static void masks,
!! and optional view-factor radiation coupling data needed to assemble the
!! coupled residual.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! Residual assembly gathers the off-process state needed by the mimetic
!! discretization, imposes thermal Dirichlet face values while evaluating the
!! coupled thermal/species terms that may depend on face temperature, and
!! restores the input vector before returning.
!!

#include "f90_assert.fpp"

module htsd_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use mfd_disc_type
  use thermal_component_type
  use thermal_component_factory, only: define_thermal_component
  use thermal_view_factor_coupling_type, only: thermal_view_factor_coupling
  use TofH_type
  use species_component_type
  use matl_mesh_func_type
  use thermal_bc_factory_type
  use species_bc_factory_type
  use thermal_source_factory_type
  use species_source_factory_type
  use species_component_factory, only: define_species_component
  use htsd_vector_type
  implicit none
  private

  type, public :: htsd_model
    integer :: num_comp = 0
    type(thermal_component) :: thermal
    type(TofH) :: T_of_H
    type(thermal_view_factor_coupling) :: vf_rad
    type(species_component), allocatable :: species(:)
    type(mfd_disc) :: disc
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    logical, allocatable :: void_cell(:), void_face(:)
    integer, allocatable :: void_dir_faces(:)
    real(r8) :: void_temp = 0.0_r8
  contains
    procedure :: init
    procedure :: init_vector
    procedure :: define_void_masks
    procedure :: residual
    procedure :: set_initial_time
    procedure :: update_moving_vf
    procedure :: add_moving_vf_events
    procedure :: set_ext_enthalpy_rate
    procedure :: set_ext_species_rate
  end type

contains

  subroutine init(this, mesh, mmf, params, stat, errmsg)

    use parameter_list_type
    use species_bc_factory_type
    use thermal_bc_factory_type

    !! TARGET is needed because T_of_H stores a persistent pointer to
    !! this%thermal%H_of_T.
    class(htsd_model), intent(out), target :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(matl_mesh_func), intent(in), target :: mmf
    type(parameter_list), intent(inout), target :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(species_bc_factory) :: sbc_fac
    type(species_source_factory) :: ssrc_fac
    type(thermal_bc_factory) :: tbc_fac
    type(thermal_source_factory) :: tsrc_fac
    type(parameter_list), pointer :: heat_params, species_params, sublist
    real(r8) :: stefan_boltzmann, absolute_zero, eps, delta
    integer :: max_try
    integer :: n

    call this%disc%init(mesh)
    this%mesh => mesh

    heat_params => params%sublist('heat')
    call heat_params%get('stefan-boltzmann', stefan_boltzmann, stat, errmsg, default=5.67e-8_r8)
    if (stat /= 0) return
    call heat_params%get('absolute-zero', absolute_zero, stat, errmsg, default=0.0_r8)
    if (stat /= 0) return
    call heat_params%get('void-temperature', this%void_temp, stat, errmsg, default=0.0_r8)
    if (stat /= 0) return
    call heat_params%get('tofh-tol', eps, stat, errmsg, default=0.0_r8)
    if (stat /= 0) return
    call heat_params%get('tofh-max-try', max_try, stat, errmsg, default=50)
    if (stat /= 0) return
    call heat_params%get('tofh-delta', delta, stat, errmsg, default=1.0_r8)
    if (stat /= 0) return

    sublist => heat_params%sublist('bc')
    call tbc_fac%init(mesh, stefan_boltzmann, absolute_zero, sublist)

    sublist => heat_params%sublist('sources')
    call tsrc_fac%init(mesh, sublist)

    call define_thermal_component(mesh, mmf, tbc_fac, tsrc_fac, this%thermal, stat, errmsg)
    if (stat /= 0) return

    sublist => heat_params%sublist('radiation')
    call this%vf_rad%init(mesh, sublist, stat, errmsg)
    if (stat /= 0) return

    call this%vf_rad%validate_bc(mesh, this%thermal, stat, errmsg)
    if (stat /= 0) return

    call this%T_of_H%init(this%thermal%H_of_T, eps=eps, max_try=max_try, delta=delta)

    species_params => params%sublist('species')
    call species_params%get('num-species', this%num_comp, stat, errmsg)
    if (stat /= 0) return
    if (this%num_comp <= 0) then
      stat = -1
      errmsg = 'species model requires at least one species component'
      return
    end if

    sublist => species_params%sublist('bc')
    call sbc_fac%init(mesh, sublist)

    sublist => species_params%sublist('sources')
    call ssrc_fac%init(mesh, sublist)

    allocate(this%species(this%num_comp))
    do n = 1, this%num_comp
      call define_species_component(mesh, mmf, sbc_fac, ssrc_fac, n, this%species(n), &
          stat, errmsg, include_soret=.true.)
      if (stat /= 0) return
    end do

    stat = 0

  end subroutine init

  subroutine init_vector(this, vec)
    class(htsd_model), intent(in) :: this
    type(htsd_vector), intent(out) :: vec
    if (this%vf_rad%is_active()) then
      call vec%init(this%mesh, this%num_comp, this%vf_rad%size())
    else
      call vec%init(this%mesh, this%num_comp)
    end if
  end subroutine

  subroutine define_void_masks(this, mmf)

    class(htsd_model), intent(inout) :: this
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

  subroutine set_ext_enthalpy_rate(this, enthalpy_rate)
    class(htsd_model), intent(inout) :: this
    real(r8), intent(in) :: enthalpy_rate(:)
    ASSERT(size(enthalpy_rate) == this%mesh%ncell_onP)
    this%thermal%ext_rate = enthalpy_rate
  end subroutine

  subroutine set_ext_species_rate(this, n, species_rate)
    class(htsd_model), intent(inout) :: this
    integer, intent(in) :: n
    real(r8), intent(in) :: species_rate(:)
    ASSERT(n > 0 .and. n <= this%num_comp)
    ASSERT(size(species_rate) == this%mesh%ncell_onP)
    this%species(n)%ext_rate = species_rate
  end subroutine

  subroutine residual(this, t, u, udot, r)

    class(htsd_model), intent(inout) :: this
    real(r8), intent(in) :: t
    type(htsd_vector), intent(inout) :: u, udot ! data is intent(in)
    type(htsd_vector), intent(inout) :: r       ! data is intent(out)

    integer :: n, k
    real(r8), allocatable :: Tdir(:)

    !! IDAESOL callbacks require only on-process entries to be valid on entry.
    !! Gather all state fields used by the coupled thermal/species residuals
    !! here; UDOT is only used on-process below.
    call this%mesh%cell_imap%gather_offp(u%tc)
    call this%mesh%face_imap%gather_offp(u%tf)
    do k = 1, this%num_comp
      call this%mesh%cell_imap%gather_offp(u%cc(:,k))
      call this%mesh%face_imap%gather_offp(u%cf(:,k))
    end do

    !! Thermal Dirichlet face values are imposed for the whole coupled residual
    !! assembly because species fluxes and Soret terms may depend on face
    !! temperature.
    call impose_thermal_dirichlet

    call compute_thermal_residual

    do n = 1, this%num_comp
      call compute_species_residual(n)
    end do

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
      call this%thermal%H_of_T%compute_value(u%tc(:ncell_onP), u%cc(:ncell_onP,:), value(:ncell_onP))
      r%hc(:ncell_onP) = u%hc(:ncell_onP) - value(:ncell_onP)

      !! Overwrite residual on void cells with dummy equation H=0.
      if (allocated(this%void_cell)) &
          where (this%void_cell(:ncell_onP)) r%hc(:ncell_onP) = u%hc(:ncell_onP)

      !! Heat equation residual

      !! Compute the generic heat equation residual.
      call this%thermal%conductivity%compute_value(u%tc, u%cc, value)
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

    subroutine compute_species_residual(index)

      integer, intent(in) :: index

      integer :: j, n, ncell_onP, nface_onP
      real(r8), dimension(this%mesh%ncell) :: D, value
      real(r8), allocatable :: Cdir(:)

      ncell_onP = this%mesh%ncell_onP
      nface_onP = this%mesh%nface_onP

      associate (comp => this%species(index))

      !! Overwrite the concentration on Dirichlet faces with the boundary
      !! data, saving the original values to restore them later.
      if (allocated(comp%bc_dir)) then
        call comp%bc_dir%compute(t)
        allocate(Cdir(size(comp%bc_dir%index)))
        do j = 1, size(comp%bc_dir%index)
          n = comp%bc_dir%index(j)
          Cdir(j) = u%cf(n,index)
          u%cf(n,index) = comp%bc_dir%value(j)
        end do
      end if

      !! Compute the residual of the basic species diffusion equation.
      call comp%diffusivity%compute_value(u%tc, u%cc, D)
      if (allocated(this%void_cell)) where (this%void_cell) D = 0.0_r8
      call this%disc%apply_diff(D, u%cc(:,index), u%cf(:,index), r%cc(:,index), r%cf(:,index))

      !! Time derivative and external species rate contributions.
      r%cc(:ncell_onP,index) = r%cc(:ncell_onP,index) &
          + this%mesh%volume(:ncell_onP)*udot%cc(:ncell_onP,index)

      !! External rate contribution; typically from advection
      if (allocated(comp%ext_rate)) then
        r%cc(:ncell_onP,index) = r%cc(:ncell_onP,index) - comp%ext_rate(:ncell_onP)
      end if

      !! User-defined species source contribution.
      if (allocated(comp%src)) then
        call comp%src%compute(t)
        r%cc(:ncell_onP,index) = r%cc(:ncell_onP,index) &
            - this%mesh%volume(:ncell_onP)*comp%src%value(:ncell_onP)
      end if

      !! Overwrite face residuals with the Dirichlet BC residual and
      !! restore the face concentrations to their original input values.
      if (allocated(comp%bc_dir)) then
        do j = 1, size(comp%bc_dir%index)
          n = comp%bc_dir%index(j)
          r%cf(n,index) = Cdir(j) - comp%bc_dir%value(j)
          u%cf(n,index) = Cdir(j)
        end do
      end if

      !! Flux-type species BC contributions.
      call comp%add_flux_bc_residual(t, u%tf, u%cf(:,index), r%cf(:,index), this%mesh%area, this%void_face)

      !! Soret thermodiffusion contribution.
      if (allocated(comp%soret)) then
        call comp%soret%compute_value(u%tc, u%cc, value)
        value = D*value
        call this%disc%add_diff(value, u%tc, u%tf, r%cc(:,index), r%cf(:,index))
        if (allocated(comp%bc_dir)) then
          call comp%bc_dir%compute(t)
          do j = 1, size(comp%bc_dir%index)
            n = comp%bc_dir%index(j)
            r%cf(n,index) = Cdir(j) - comp%bc_dir%value(j)
          end do
        end if
      end if

      !! Overwrite residual on void cells and faces with dummy equation C=0.
      if (allocated(this%void_cell)) &
          where (this%void_cell(:ncell_onP)) r%cc(:ncell_onP,index) = u%cc(:ncell_onP,index)
      if (allocated(this%void_face)) &
          where (this%void_face(:nface_onP)) r%cf(:nface_onP,index) = u%cf(:nface_onP,index)

      end associate

    end subroutine compute_species_residual

  end subroutine residual

  subroutine set_initial_time(this, t)
    class(htsd_model), intent(inout) :: this
    real(r8), intent(in) :: t
    call this%vf_rad%set_initial_time(t)
  end subroutine

  subroutine update_moving_vf(this)
    class(htsd_model), intent(inout) :: this
    call this%vf_rad%update_moving_vf
  end subroutine

  subroutine add_moving_vf_events(this, eventq, rank)
    use sim_event_queue_type
    class(htsd_model), intent(in) :: this
    type(sim_event_queue), intent(inout) :: eventq
    integer, intent(in), optional :: rank
    call this%vf_rad%add_moving_vf_events(eventq, rank)
  end subroutine

end module htsd_model_type
