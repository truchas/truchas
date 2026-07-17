!!
!! SD_MODEL_TYPE
!!
!! This module defines the discrete species transport model used by the
!! implicit integrator. It owns the species components, the mimetic diffusion
!! discretization, and static void masks needed to assemble the concentration
!! residuals.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! Residual assembly gathers the off-process concentration fields needed by
!! the mimetic discretization. Species Dirichlet face values are imposed while
!! each component residual is evaluated and then restored before returning.
!!

#include "f90_assert.fpp"

module sd_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use mfd_disc_type
  use sd_vector_type
  use species_component_type
  use species_component_factory, only: define_species_component
  use matl_mesh_func_type
  use species_bc_factory_type
  use species_source_factory_type
  implicit none
  private

  type, public :: sd_model
    integer :: num_comp = 0
    type(species_component), allocatable :: species(:)
    type(mfd_disc) :: disc
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    logical, allocatable :: void_cell(:), void_face(:)
    integer, allocatable :: void_dir_faces(:)
  contains
    procedure :: init
    procedure :: init_vector
    procedure :: define_void_masks
    procedure :: residual
    procedure :: set_ext_species_rate
  end type

contains

  subroutine init(this, mesh, mmf, params, stat, errmsg)

    use parameter_list_type

    class(sd_model), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(matl_mesh_func), intent(in), target :: mmf
    type(parameter_list), intent(inout), target :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(species_bc_factory) :: bc_fac
    type(species_source_factory) :: src_fac
    type(parameter_list), pointer :: sublist
    integer :: n

    call this%disc%init(mesh)
    this%mesh => mesh

    call params%get('num-species', this%num_comp, stat, errmsg)
    if (stat /= 0) return
    if (this%num_comp <= 0) then
      stat = -1
      errmsg = 'species model requires at least one species component'
      return
    end if

    sublist => params%sublist('bc')
    call bc_fac%init(mesh, sublist, use_temperature=.false.)

    sublist => params%sublist('sources')
    call src_fac%init(mesh, sublist)

    allocate(this%species(this%num_comp))
    do n = 1, this%num_comp
      call define_species_component(mesh, mmf, bc_fac, src_fac, n, this%species(n), stat, errmsg)
      if (stat /= 0) return
    end do
    stat = 0

  end subroutine init

  subroutine init_vector(this, vec)
    class(sd_model), intent(in) :: this
    type(sd_vector), intent(out) :: vec
    call vec%init(this%mesh, this%num_comp)
  end subroutine

  subroutine define_void_masks(this, mmf)

    class(sd_model), intent(inout) :: this
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

  subroutine set_ext_species_rate(this, n, species_rate)
    class(sd_model), intent(inout) :: this
    integer, intent(in) :: n
    real(r8), intent(in) :: species_rate(:)
    ASSERT(n > 0 .and. n <= this%num_comp)
    ASSERT(size(species_rate) == this%mesh%ncell_onP)
    this%species(n)%ext_rate = species_rate
  end subroutine

  subroutine residual(this, t, u, udot, r)

    class(sd_model), intent(inout) :: this
    real(r8), intent(in) :: t
    type(sd_vector), intent(inout) :: u, udot ! data is intent(in)
    type(sd_vector), intent(inout) :: r       ! data is intent(out)

    integer :: n

    !! IDAESOL callbacks require only on-process entries to be valid on entry.
    !! Gather all cell and face concentration fields used by the species
    !! residuals here; UDOT is only used on-process below.
    call u%gather_offp

    do n = 1, this%num_comp
      call compute_species_residual(n)
    end do

  contains

    subroutine compute_species_residual(index)
      integer, intent(in) :: index

      integer :: j, n, ncell_onP, nface_onP
      real(r8), dimension(this%mesh%ncell) :: D
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
      call comp%diffusivity%compute_value(u%cc, D)
      if (allocated(this%void_cell)) where (this%void_cell) D = 0.0_r8
      call this%disc%apply_diff(D, u%cc(:,index), u%cf(:,index), &
          r%cc(:,index), r%cf(:,index))

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
      call comp%add_flux_bc_residual(t, u%cf(:,index), &
          r%cf(:,index), this%mesh%area, this%void_face)

      !! Overwrite residual on void cells and faces with dummy equation C=0.
      if (allocated(this%void_cell)) &
          where (this%void_cell(:ncell_onP)) r%cc(:ncell_onP,index) = u%cc(:ncell_onP,index)
      if (allocated(this%void_face)) &
          where (this%void_face(:nface_onP)) r%cf(:nface_onP,index) = u%cf(:nface_onP,index)

      end associate

    end subroutine compute_species_residual

  end subroutine residual

end module sd_model_type
