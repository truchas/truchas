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
!!  APPLY(F) applies the preconditioner for the model to the vector F, which is
!!    overwritten with the result.
!!

#include "f90_assert.fpp"

module HT_2d_precon_type

  use kinds, only: r8
  use HT_2d_model_type
  use unstr_2d_mesh_type
  use mfd_2d_diff_precon_type
  use mfd_2d_diff_matrix_type
  use parameter_list_type
  use index_partitioning
  implicit none
  private

  type, public :: HT_2d_precon
    type(HT_2d_model),   pointer :: model => null()  ! reference only -- do not own
    type(unstr_2d_mesh), pointer :: mesh  => null()  ! reference only -- do not own
    real(r8) :: dt  ! time step
    real(r8), allocatable :: dHdT(:)  ! derivative of the enthalpy/temperature relation
    type(mfd_2d_diff_precon) :: hcprecon ! heat equation preconditioner
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
  end type HT_2d_precon

contains

  subroutine init(this, model, params)

    class(HT_2d_precon), intent(out) :: this
    type(HT_2d_model), intent(in), target :: model
    type(parameter_list) :: params

    type(mfd_2d_diff_matrix), allocatable :: dm

    this%model => model
    this%mesh  => model%mesh

    !! Heat density/temperature relation derivative.
    allocate(this%dHdT(this%mesh%ncell))

    !! Create the preconditioner for the heat equation.
    !! The preconditioner assumes ownership of the matrix.
    allocate(dm)
    call dm%init(model%disc)
    call this%hcprecon%init(dm, params)

  end subroutine init


  !TODO: handle void cells?
  subroutine compute(this, t, u, dt)

    class(HT_2d_precon), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), target :: u(:)

    real(r8) :: coef(this%mesh%ncell), Tface(this%mesh%nface)
    real(r8), allocatable :: state(:,:)
    type(mfd_2d_diff_matrix), pointer :: dm

    ASSERT(size(u) == this%model%num_dof())
    ASSERT(dt > 0.0_r8)

    call this%model%get_face_temp_copy(u, Tface)
    call gather_boundary(this%mesh%face_ip, Tface)

    call this%model%new_state_array(u, state)

    this%dt = dt  ! the time step size

    !! Grab the matrix; we will update its values.
    dm => this%hcprecon%matrix()

    !! Jacobian of the heat equation diffusion operator ignoring nonlinearities
    !! in the conductivity.
    call this%model%conductivity%compute_value(state, coef)
    call dm%compute(coef)

    !! Contribution from the time derivative (H/T relation eliminated).
    !TODO: fix hardwired assumption that T is the first component of state
    call this%model%H_of_T%compute_deriv(state, 1, this%dHdT)
    call dm%incr_cell_diag(this%mesh%volume * this%dHdT / dt)

    !! Dirichlet boundary condition fixups.
    if (allocated(this%model%bc_dir)) then
      call this%model%bc_dir%compute(t)
      call dm%set_dir_faces(this%model%bc_dir%index)
    end if

    !! The matrix is now complete; re-compute the preconditioner.
    call this%hcprecon%compute

  end subroutine compute


  subroutine apply(this, f)

    class(HT_2d_precon), intent(in) :: this
    real(r8), intent(inout), target :: f(:)

    real(r8), pointer :: f0(:), f1(:)  ! on-process segments
    real(r8) :: f1x(this%mesh%ncell), f2x(this%mesh%nface)  ! off-process extended copies

    !! Heat equation cell residual with the H/T relation residual eliminated.
    call this%model%get_cell_heat_view(f, f0)
    call this%model%get_cell_temp_view(f, f1)
    f1x(:this%mesh%ncell_onP) = f1 - (this%mesh%volume(:this%mesh%ncell_onP)/this%dt)*f0
    call gather_boundary(this%mesh%cell_ip, f1x)

    !! Heat equation face residual.
    call this%model%get_face_temp_copy(f, f2x)
    call gather_boundary(this%mesh%face_ip, f2x)

    !! Precondition the heat equation.
    call this%hcprecon%apply(f1x, f2x)
    call this%model%set_cell_temp(f1x, f)
    call this%model%set_face_temp(f2x, f)

    !! Backsubstitute to obtain the preconditioned H/T-relation residual.
    f0 = f0 + this%dHdT(:this%mesh%ncell_onP)*f1

  end subroutine apply

end module HT_2d_precon_type
