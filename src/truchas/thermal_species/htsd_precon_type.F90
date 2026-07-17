!!
!! HTSD_PRECON_TYPE
!!
!! This module defines the preconditioner used by the implicit integrator for
!! the nonlinear coupled thermal/species transport time-step system. It
!! assembles and applies approximate Jacobians for the thermal and species
!! transport residuals, including thermal and species boundary contributions,
!! optional Soret thermodiffusion coupling, and optional view-factor radiation
!! coupling.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module htsd_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use htsd_model_type
  use htsd_vector_type
  use unstr_mesh_type
  use mfd_diff_precon_type
  use mfd_diff_matrix_type
  use thermal_view_factor_coupling_type, only: VFR_JAC, VFR_FGS, VFR_BGS, VFR_FAC
  implicit none
  private

  type, public :: htsd_precon
    type(htsd_model), pointer :: model => null() ! unowned reference
    type(unstr_mesh), pointer :: mesh  => null() ! unowned reference
    ! enthalpy transport preconditioning data
    real(r8) :: dt
    real(r8), allocatable :: dHdT(:)
    type(mfd_diff_precon) :: ht_pc
    integer :: vfr_precon_coupling = VFR_BGS
    ! species transport preconditioning data
    type(mfd_diff_precon), allocatable :: sd_pc(:)
    logical :: have_soret_coupling = .false.
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
  end type

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type
    use thermal_view_factor_coupling_type, only: vfr_precon_coupling_from_string

    class(htsd_precon), intent(out) :: this
    type(htsd_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    type(mfd_diff_matrix), allocatable, target :: matrix
    type(mfd_diff_matrix), pointer :: mold_matrix => null()
    character(:), allocatable :: method

    this%model => model
    this%mesh  => model%mesh

    !! Enthalpy density/temperature relation derivative.
    allocate(this%dHdT(this%mesh%ncell))

    !! Create the heat equation preconditioner; compute defines its values.
    allocate(matrix)
    call matrix%init(model%disc)
    call this%ht_pc%init(matrix, params, stat, errmsg)
    if (stat /= 0) return
    mold_matrix => this%ht_pc%matrix_ref()

    !! Heat equation/view factor radiation preconditioner coupling method.
    if (model%vf_rad%is_active()) then
      call params%get('vfr-precon-coupling', method, default='backward-gs')
      this%vfr_precon_coupling = vfr_precon_coupling_from_string(method)
    end if

    !! Create the species transport preconditioner; compute defines its values.
    allocate(this%sd_pc(model%num_comp))
    do n = 1, this%model%num_comp
      allocate(matrix)
      call matrix%init(mold=mold_matrix)
      call this%sd_pc(n)%init(matrix, params, stat, errmsg)
      if (stat /= 0) return
      if (allocated(model%species(n)%soret)) this%have_soret_coupling = .true.
    end do

  end subroutine init

  subroutine compute (this, t, u, dt)

    class(htsd_precon), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    type(htsd_vector), intent(inout) :: u

    integer :: n
    real(r8) :: D(this%mesh%ncell), A(this%mesh%ncell)

    ASSERT(dt > 0.0_r8)

    !! Preconditioner assembly evaluates material properties over the
    !! off-process extended temperature and concentration fields.
    call u%gather_offp

    !! Compute the heat equation preconditioner.
    call compute_heat_precon

    !! Compute the species transport preconditioner.
    do n = 1, size(this%sd_pc)
      call compute_species_precon(n)
    end do

  contains

    subroutine compute_heat_precon

      type(mfd_diff_matrix), pointer :: matrix

      !! The time step size.
      this%dt = dt

      !! Derivative with respect to temperature, the first material-function argument.
      call this%model%thermal%H_of_T%compute_deriv(u%tc, u%cc, 1, this%dHdT)
      call this%model%thermal%conductivity%compute_value(u%tc, u%cc, D)
      A = this%mesh%volume * this%dHdT / dt

      !! Correct data on void cells.
      if (allocated(this%model%void_cell)) then
        where (this%model%void_cell)
          this%dHdT = 0.0_r8
          D = 0.0_r8
          A = 1.0_r8
        end where
      end if

      matrix => this%ht_pc%matrix_ref()

      !! Jacobian of the basic heat equation that ignores nonlinearities
      !! in the conductivity.  This has the algebraic H/T relation eliminated.
      call matrix%compute(D)
      call matrix%incr_cell_diag (A)

      !! Dirichlet boundary condition fixups.
      if (allocated(this%model%thermal%bc_dir)) then
        call this%model%thermal%bc_dir%compute(t)
        call matrix%set_dir_faces(this%model%thermal%bc_dir%index)
      end if

      !! Flux-type thermal boundary/interface condition Jacobian contributions.
      call this%model%thermal%add_flux_bc_jacobian(t, u%tf, this%mesh%area, matrix, &
          this%model%void_face)

      !! Dirichlet fix-ups for void faces.
      if (allocated(this%model%void_dir_faces)) call matrix%set_dir_faces(this%model%void_dir_faces)

      !! Enclosure radiation contributions to the preconditioner.
      !! TODO: what about factorization coupling?  Is this still correct?
      if (this%model%vf_rad%is_active()) then
        call this%model%vf_rad%add_heat_precon_deriv(t, u%tf, this%mesh%area, matrix)
      end if

      !! The matrix is now complete; re-compute the preconditioner.
      call this%ht_pc%compute

    end subroutine compute_heat_precon

    subroutine compute_species_precon(index)

      integer, intent(in) :: index

      type(mfd_diff_matrix), pointer :: matrix

      matrix => this%sd_pc(index)%matrix_ref()

      !! Jacobian of the diffusion operator that ignores nonlinearities.
      call this%model%species(index)%diffusivity%compute_value(u%tc, u%cc, D)
      if (allocated(this%model%void_cell)) where (this%model%void_cell) D = 0.0_r8
      call matrix%compute(D)

      !! Time derivative contribution to the diffusion equation Jacobian.
      D = this%mesh%volume / dt
      if (allocated(this%model%void_cell)) where (this%model%void_cell) D = 1.0_r8
      call matrix%incr_cell_diag (D)

      !! Dirichlet BC fixups.
      if (allocated(this%model%species(index)%bc_dir)) then
        call this%model%species(index)%bc_dir%compute(t)
        call matrix%set_dir_faces(this%model%species(index)%bc_dir%index)
      end if
      if (allocated(this%model%void_dir_faces)) call matrix%set_dir_faces(this%model%void_dir_faces)

      !! Flux-type species boundary/interface condition Jacobian contributions.
      call this%model%species(index)%add_flux_bc_jacobian(t, u%tf, &
          u%cf(:,index), matrix, this%model%void_face)

      !! The matrix is now complete; re-compute the preconditioner.
      call this%sd_pc(index)%compute

    end subroutine compute_species_precon

  end subroutine compute

  subroutine apply(this, t, u, f)

    class(htsd_precon), intent(inout) :: this
    real(r8), intent(in) :: t
    type(htsd_vector), intent(inout) :: u
    type(htsd_vector), intent(inout) :: f

    integer :: index, ncell_onP, nface_onP
    real(r8), allocatable :: FTdir(:)

    ncell_onP = this%mesh%ncell_onP
    nface_onP = this%mesh%nface_onP

    !! Precondition the HT component.
    call apply_heat_precon

    !! Precondition the SD components.
    !! Initialize extra data needed to handle Soret coupling.
    if (this%have_soret_coupling) then
      call f%gather_offp
      call u%gather_offp
      !! Off-process-extended preconditioned HT components. Temperature
      !! Dirichlet projection is applied directly and restored below.
      if (allocated(this%model%thermal%bc_dir)) then
        allocate(FTdir(size(this%model%thermal%bc_dir%index)))
        FTdir = f%tf(this%model%thermal%bc_dir%index)
        f%tf(this%model%thermal%bc_dir%index) = 0.0_r8 ! temperature Dirichlet projection
      end if
      !TODO! void face dirichlet projection?
    end if
    do index = 1, this%model%num_comp
      call apply_species_precon (index)
    end do
    if (allocated(FTdir)) then
      f%tf(this%model%thermal%bc_dir%index) = FTdir
      deallocate(FTdir)
    end if

  contains

    subroutine apply_heat_precon ()

      real(r8), allocatable :: z(:)

      !! Precondition the radiosity components:
      !! Factorization and forward Gauss-Seidel coupling methods.
      if (this%model%vf_rad%is_active()) then
        if (this%vfr_precon_coupling == VFR_FGS .or. &
            this%vfr_precon_coupling == VFR_FAC) then
          z = f%qrad
          call this%model%vf_rad%apply_rad_precon(t, z)
          if (this%vfr_precon_coupling == VFR_FGS) f%qrad(:) = z
          !! Update the heat equation face residual.
          call this%model%vf_rad%apply_rad_precon_matvec1(t, z)
          call this%model%vf_rad%add_heat_flux_to_residual(this%mesh%area, z, f%tf)
          deallocate(z)
        end if
      end if

      !! Heat equation cell residual with the H/T relation residual eliminated.
      if (allocated(this%model%void_cell)) then
        where (.not.this%model%void_cell(:ncell_onP))
          f%tc(:ncell_onP) = f%tc(:ncell_onP) &
              - (this%mesh%volume(:ncell_onP)/this%dt)*f%hc(:ncell_onP)
        end where
      else
        f%tc(:ncell_onP) = f%tc(:ncell_onP) &
            - (this%mesh%volume(:ncell_onP)/this%dt)*f%hc(:ncell_onP)
      end if

      !! Precondition the heat equation. The diffusion preconditioner owns any
      !! off-process gathers it needs; only on-process results are significant
      !! after the call.
      call this%ht_pc%apply(f%tc, f%tf)

      !! Backsubstitute to obtain the preconditioned H/T-relation residual.
      f%hc(:ncell_onP) = f%hc(:ncell_onP) + this%dHdT(:ncell_onP)*f%tc(:ncell_onP)

      !! Precondition the radiosity components:
      !! Factorization, backward Gauss-Seidel and Jacobi coupling methods.
      if (this%model%vf_rad%is_active()) then
        if (this%vfr_precon_coupling == VFR_JAC .or. &
            this%vfr_precon_coupling == VFR_BGS .or. &
            this%vfr_precon_coupling == VFR_FAC) then
          !! Update the radiosity residual components.
          if (this%vfr_precon_coupling /= VFR_JAC) then
            call this%model%vf_rad%add_rhs_deriv_times_face_vector(t, u%tf, f%tf, f%qrad)
          end if
          call this%model%vf_rad%apply_rad_precon(t, f%qrad)
        end if
      end if

    end subroutine apply_heat_precon

    subroutine apply_species_precon(index)

      integer, intent(in) :: index

      real(r8) :: Fcell(this%mesh%ncell), Fface(this%mesh%nface)
      real(r8), allocatable :: value(:)

      !! Lower triangle coupling of HT and SD; forward elimination.
      if (allocated(this%model%species(index)%soret)) then
        !! Compute the update.
        allocate(value(this%mesh%ncell))
        call this%model%species(index)%soret%compute_value(u%tc, u%cc, value)
        call this%model%species(index)%diffusivity%compute_value(u%tc, u%cc, Fcell) ! Fcell used as temporary
        value = value * Fcell
        call this%model%disc%apply_diff (value, f%tc, f%tf, Fcell, Fface)
        if (allocated(this%model%species(index)%bc_dir)) then
          Fface(this%model%species(index)%bc_dir%index) = 0.0_r8 ! concentration Dirichlet projection
        end if
        !TODO! void face dirichlet projection?
        !! Apply the update.
        f%cc(:ncell_onP,index) = f%cc(:ncell_onP,index) - Fcell(:ncell_onP)
        f%cf(:nface_onP,index) = f%cf(:nface_onP,index) - Fface(:nface_onP)
      end if

      !TODO! HT and SD coupling due to temperature-dependence of the MTC coefficient.

      !! Precondition the diffusion equation for this component. The component
      !! preconditioner owns any off-process gathers it needs; only on-process
      !! results are significant after the call.
      call this%sd_pc(index)%apply(f%cc(:,index), f%cf(:,index))

    end subroutine apply_species_precon

  end subroutine apply

end module htsd_precon_type
