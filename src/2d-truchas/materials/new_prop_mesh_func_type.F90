#include "f90_assert.fpp"

module new_prop_mesh_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use matl_prop_class
  private

  type, public :: prop_mesh_func
    !type(matl_mesh_func), pointer :: mmf => null() ! reference only -- not owned
    !type(avg_matl_prop)
    class(matl_prop), allocatable :: prop
  contains
    procedure :: init
    generic :: compute_value => compute, compute_cell
    procedure, private :: compute, compute_cell
    procedure :: compute_deriv
  end type

contains

  subroutine init(this, matl_model, name, stat, errmsg)

    use material_model_type

    class(prop_mesh_func), intent(out) :: this
    type(material_model), intent(in) :: matl_model
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    stat = 1
    if (matl_model%nmatl == 1 .and. matl_model%nmatl_real == 1) then
      call matl_model%get_matl_prop(1, name, this%prop, errmsg)
      if (.not.allocated(this%prop)) return
    else
      errmsg = 'indeterminate material'
      return
    end if
    stat = 0

  end subroutine init

  subroutine compute(this, state, value)

    class(prop_mesh_func), intent(in) :: this
    real(r8), intent(in) :: state(:,:)
    real(r8), intent(out) :: value(:)

    integer :: j

    ASSERT(size(state,dim=1) == size(value))

    !TODO: version of compute_value that applies to a vector of states
    do j = 1, size(value)
      call this%prop%compute_value(state(j,:), value(j))
    end do

  end subroutine compute

  subroutine compute_cell(this, n, state, value)

    class(prop_mesh_func), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(in) :: state(:)
    real(r8), intent(out) :: value

    call this%prop%compute_value(state, value)

  end subroutine compute_cell

  subroutine compute_deriv(this, state, index, value)

    class(prop_mesh_func), intent(in) :: this
    real(r8), intent(in)  :: state(:,:)
    integer,  intent(in)  :: index
    real(r8), intent(out) :: value(:)

    integer :: j

    ASSERT(size(state,dim=1) == size(value))

    do j = 1, size(value)
      call this%prop%compute_deriv(state(j,:), index, value(j))
    end do

  end subroutine compute_deriv

end module new_prop_mesh_func_type
