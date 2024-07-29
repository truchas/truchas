!!
!! Abstract solid mechanics preconditioner class
!!
!! Zach Jibben <zjibben@lanl.gov>
!! April 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module sm_precon_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use sm_model_type
  use parameter_list_type
  implicit none
  private

  type, abstract, public :: sm_precon
  contains
    procedure(init), deferred :: init
    procedure(compute), deferred :: compute
    procedure(apply), deferred :: apply
  end type

  abstract interface
    subroutine init(this, model, params)
      import sm_precon, sm_model, parameter_list
      class(sm_precon), intent(out) :: this
      type(sm_model), intent(in), target :: model
      type(parameter_list), intent(inout) :: params
    end subroutine

    subroutine compute(this, t, dt, displ)
      import sm_precon, r8
      class(sm_precon), intent(inout) :: this
      real(r8), intent(in) :: t, dt
      real(r8), intent(inout) :: displ(:,:) ! need to update halo
    end subroutine

    subroutine apply(this, u, f)
      import sm_precon, r8
      class(sm_precon), intent(in), target :: this
      real(r8), intent(in), contiguous, target :: u(:,:) ! current displacement guess
      real(r8), intent(inout), contiguous, target :: f(:,:) ! in residual, out next displacement guess
    end subroutine
  end interface

end module sm_precon_class
