!TODO: finish documentation
!!
!! HT_2D_PRECON_TYPE
!!
!! This module defines a derived type that encapsulates the preconditioner for
!! the 2D heat conduction model.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! April 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! The module defines the derived type HT_2D_PRECON_TYPE.  It has the following
!! type bound procedures.
!!
!!  INIT(MODEL, PARAMS) initializes the object.  MODEL is of type PC_MODEL.
!!    The object will hold a reference to the model, and so the actual argument
!!    must have the target attribute and persist for the lifetime of the object.
!!    The PARAMETER_LIST type argument PARAMS gives the parameters for the
!!    preconditioner.  For this model there is only the single heat equation
!!    and the expected parameters are those described for MFD_2D_DIFF_PRECON.
!!
!!  COMPUTE(T, U, DT) computes the preconditioner for the model at time T,
!!    unknown vector U, and time step DT.  It must be called before calling
!!    the APPLY procedure.
!!
!!  APPLY(T, U, F) applies the preconditioner at time T and state U to the
!!    vector F, which is overwritten with the result.
!!

#include "f90_assert.fpp"

module ht_2d_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ht_2d_model_type
  use ht_2d_vector_type
  use unstr_2d_mesh_type
  use mfd_2d_diff_precon_type
  use mfd_2d_diff_matrix_type
  implicit none
  private

  type, public :: ht_2d_precon
    type(ht_2d_model),   pointer :: model => null()  ! reference only -- do not own
    type(unstr_2d_mesh), pointer :: mesh  => null()  ! reference only -- do not own
    real(r8) :: dt  ! time step
    real(r8), allocatable :: dHdT(:)  ! derivative of the enthalpy/temperature relation
    type(mfd_2d_diff_precon) :: hcprecon ! heat equation preconditioner
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
  end type ht_2d_precon

contains

  subroutine init(this, model, params)

    use parameter_list_type
    use truchas_logging_services

    class(ht_2d_precon), intent(out) :: this
    type(ht_2d_model), intent(in), target :: model
    type(parameter_list) :: params

    integer :: stat
    type(mfd_2d_diff_matrix), allocatable :: dm
    character(:), allocatable :: errmsg

    this%model => model
    this%mesh  => model%mesh

    !! Heat density/temperature relation derivative.
    allocate(this%dHdT(this%mesh%ncell))

    !! Create the preconditioner for the heat equation.
    !! The preconditioner assumes ownership of the matrix.
    allocate(dm)
    call dm%init(model%disc)
    call this%hcprecon%init(dm, params, stat, errmsg)
    if (stat /= 0) call TLS_fatal('ht_2d_PRECON%INIT: ' // errmsg)

  end subroutine init


  subroutine compute(this, t, u, dt)

    class(ht_2d_precon), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    type(ht_2d_vector), intent(in) :: u

    real(r8) :: coef(this%mesh%ncell)
    real(r8), allocatable :: state(:,:)
    type(mfd_2d_diff_matrix), pointer :: dm

    ASSERT(dt > 0.0_r8)

    allocate(state(this%mesh%ncell,0:0))
    state(:,0) = u%tc
    call this%mesh%cell_imap%gather_offp(state(:,0))

    this%dt = dt
    dm => this%hcprecon%matrix()
    call this%model%conductivity%compute_value(state, coef)
    call dm%compute(coef)
    call this%model%H_of_T%compute_deriv(state, 1, this%dHdT)
    call dm%incr_cell_diag(this%mesh%volume*this%dHdT/dt)

    if (allocated(this%model%bc_dir)) then
      call this%model%bc_dir%compute(t)
      call dm%set_dir_faces(this%model%bc_dir%index)
    end if
    call this%hcprecon%compute

  end subroutine compute


  subroutine apply(this, t, u, f)

    class(ht_2d_precon), intent(in) :: this
    real(r8), intent(in) :: t
    type(ht_2d_vector), intent(inout) :: u
    type(ht_2d_vector), intent(inout) :: f

    associate (mesh => this%mesh)
      !! Eliminate the enthalpy-temperature residual from the heat equation.
      f%tc(:mesh%ncell_onP) = f%tc(:mesh%ncell_onP) &
                             - (mesh%volume(:mesh%ncell_onP)/this%dt)*f%hc(:mesh%ncell_onP)
      call mesh%cell_imap%gather_offp(f%tc)
      call mesh%face_imap%gather_offp(f%tf)

      call this%hcprecon%apply(f%tc, f%tf)

      !! Backsubstitute the enthalpy-temperature relation residual.
      f%hc(:mesh%ncell_onP) = f%hc(:mesh%ncell_onP) &
                             + this%dHdT(:mesh%ncell_onP)*f%tc(:mesh%ncell_onP)
    end associate

  end subroutine apply

end module ht_2d_precon_type
