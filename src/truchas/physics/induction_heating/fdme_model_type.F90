module fdme_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simpl_mesh_type
  use pcsr_matrix_type
  implicit none
  private

  type, public :: fdme_model
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(pcsr_matrix) :: A(2,2)
  contains
    procedure :: init
  end type

contains

  subroutine init(this, mesh, params, stat, errmsg)
  
    use parameter_list_type

    class(fdme_model), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
  
    this%mesh => mesh
    
  end subroutine init

end module fdme_model_type
