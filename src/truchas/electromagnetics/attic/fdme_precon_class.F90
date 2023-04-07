module fdme_precon_class

  use fdme_model_type
  use simpl_mesh_type
  implicit none
  private

  type, abstract, public :: fdme_precon
  contains
    procedure(init),  deferred :: init
    procedure(setup), deferred :: setup
    procedure(apply), deferred :: apply
  end type

  abstract interface
    subroutine init(this, model, params) !TODO: add stat, errmsg
      use parameter_list_type
      import
      class(fdme_precon), intent(out) :: this
      type(fdme_model), intent(in), target :: model
      type(parameter_list), intent(inout) :: params
    end subroutine
    subroutine setup(this)
      import
      class(fdme_precon), intent(inout) :: this
    end subroutine
    subroutine apply(this, x)
      use,intrinsic :: iso_fortran_env, only: r8 => real64
      import
      class(fdme_precon), intent(inout) :: this
      real(r8), intent(inout) :: x(:,:)
    end subroutine
  end interface

end module fdme_precon_class
