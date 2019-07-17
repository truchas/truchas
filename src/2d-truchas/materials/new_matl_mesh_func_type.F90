module new_matl_mesh_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  private

  type, public :: matl_mesh_func
    integer, allocatable :: matl_id(:)
    real(r8), allocatable :: vfrac(:,:)
  contains
    procedure :: init
  end type

contains

  subroutine init(this, mesh, matl_id)
    use base_mesh_class
    type(matl_mesh_func), intent(out) :: this
    class(base_mesh), intent(in), target :: mesh
    integer, intent(in) :: matl_id(:)
    this%matl_id = matl_id
    allocate(this%vfrac(size(matl_id),mesh%ncell_onP)
  end subroutine

end module new_matl_mesh_func_type
