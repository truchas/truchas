!!
!! FDME_SOLVER_TYPE
!!
!! This module provides a type for solving the frequency-domain Maxwell
!! Equations.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! April 2023
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module fdme_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simpl_mesh_type
  use fdme_model_type
  use fdme_minres_solver_type
  use fdme_mixed_minres_solver_type
#ifdef USE_MUMPS
  use fdme_mumps_solver_type
#endif
  implicit none
  private

  type, public :: fdme_solver
    private
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(fdme_model), pointer :: model => null()
    type(fdme_minres_solver), allocatable :: minres
    type(fdme_mixed_minres_solver), allocatable :: mixed_minres
#ifdef USE_MUMPS
    type(fdme_mumps_solver), allocatable :: mumps
#endif
    integer :: print_level
    complex(r8), allocatable, public :: efield(:), bfield(:)
  contains
    procedure :: init
    procedure :: solve
    procedure :: get_heat_source
    procedure :: get_cell_efield, get_cell_hfield, get_div_dfield
  end type fdme_solver

contains

  subroutine init(this, mesh, omega, epsr, epsi, mu, sigma, params, stat, errmsg)

    use em_bc_factory_type
    use parameter_list_type

    class(fdme_solver), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    real(r8), intent(in) :: omega, epsr(:), epsi(:), mu(:), sigma(:)
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(em_bc_factory) :: bc_fac
    type(parameter_list), pointer :: plist
    character(:), allocatable :: solver_type

    this%mesh => mesh

    plist => params%sublist('bc')
    call bc_fac%init(this%mesh, omega, plist)

    plist => params%sublist('fd-solver')

    allocate(this%model)
    call this%model%init(mesh, bc_fac, plist, stat, errmsg)
    if (stat /= 0) return

    call this%model%setup(omega, epsr, epsi, mu, sigma)

    call plist%get('print-level', this%print_level, stat, errmsg, default=0)
    if (stat /= 0) return
    call plist%set('print-level', max(0,this%print_level-1))

    call plist%get('solver-type', solver_type, stat, errmsg)
    if (stat /= 0) return

    select case (solver_type)
    case ('minres')
      if (this%model%use_mixed_form) then
        allocate(this%mixed_minres)
      else
        allocate(this%minres)
      end if
    case ('mumps')
#ifdef USE_MUMPS
      allocate(this%mumps)
#else
      stat = 1
      errmsg = 'executable compiled without support for solver-type "mumps"'
      return
#endif
    case default
      stat = 1
      errmsg = 'invalid solver-type value: ' // solver_type
      return
    end select

    if (allocated(this%minres)) then
      call this%minres%init(this%model, plist, stat, errmsg)
    else if (allocated(this%mixed_minres)) then
      call this%mixed_minres%init(this%model, plist, stat, errmsg)
#ifdef USE_MUMPS
    else if (allocated(this%mumps)) then
      call this%mumps%init(this%model, plist, stat, errmsg)
#endif
    else
      INSIST(.false.)
    end if

    allocate(this%efield(this%mesh%nedge), this%bfield(this%mesh%nface))

  end subroutine init


  subroutine solve(this, stat, errmsg)

    use string_utilities, only: i_to_c

    class(fdme_solver), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    if (allocated(this%minres)) then
      call this%minres%solve(this%efield, stat, errmsg)
    else if (allocated(this%mixed_minres)) then
      call this%mixed_minres%solve(this%efield, stat, errmsg)
#ifdef USE_MUMPS
    else if (allocated(this%mumps)) then
      call this%mumps%solve(this%efield, stat)
      if (stat /= 0) errmsg = 'MUMPS solve failed: stat=' // i_to_c(stat)
#endif
    else
      INSIST(.false.)
    end if

    block
      use parallel_communication, only: global_norm2
      use truchas_logging_services
      complex(r8) :: r(this%mesh%nedge), div_efield(this%mesh%nnode)
      real(r8) :: rnorm, bnorm, dnorm
      character(128) :: msg
      call this%model%residual(this%efield, r)
      call this%model%compute_div(this%efield, div_efield)
      rnorm = global_norm2(r(:this%mesh%nedge_onP))
      bnorm = global_norm2(this%model%rhs(:this%mesh%nedge_onP))
      dnorm = global_norm2(div_efield(:this%mesh%nnode_onP))
      write(msg,'(*(a,es8.2))') 'FDME_SOLVE: |r|=', rnorm, &
          ', |b|=', bnorm, ', |r|/|b|=', rnorm/bnorm, ', |Div E|=', dnorm
      call TLS_info(msg)
    end block

    call this%mesh%edge_imap%gather_offp(this%efield)
    call this%model%compute_bfield(this%efield, this%bfield)
    call this%mesh%face_imap%gather_offp(this%bfield)

  end subroutine

  subroutine get_heat_source(this, q)
    class(fdme_solver), intent(in) :: this
    real(r8), intent(out) :: q(:)
    ASSERT(size(q) >= this%mesh%ncell_onp)
    call this%model%compute_heat_source(this%efield, q)
  end subroutine

  subroutine get_cell_efield(this, efield)
    class(fdme_solver), intent(in) :: this
    complex(r8), intent(out) :: efield(:,:)
    ASSERT(size(efield,1) == 3)
    ASSERT(size(efield,2) >= this%mesh%ncell_onp)
    call this%model%cell_avg_efield(this%efield, efield)
  end subroutine

  subroutine get_cell_hfield(this, hfield)
    class(fdme_solver), intent(in) :: this
    complex(r8), intent(out) :: hfield(:,:)
    ASSERT(size(hfield,1) == 3)
    ASSERT(size(hfield,2) >= this%mesh%ncell_onp)
    call this%model%cell_avg_hfield(this%bfield, hfield)
  end subroutine

  subroutine get_div_dfield(this, div_dfield)
    class(fdme_solver), intent(in) :: this
    complex(r8), intent(out) :: div_dfield(:)
    ASSERT(size(div_dfield) ==  this%mesh%nnode)
    call this%model%compute_div(this%efield, div_dfield)
  end subroutine

end module fdme_solver_type
