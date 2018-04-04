module flow_bc_type

  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  use flow_mesh_type
  use flow_bc_factory_type
  implicit none
  private

  public :: flow_bc

  type :: flow_bc
    class(bndry_func), allocatable :: p_dirichlet
    class(bndry_vfunc), allocatable :: v_dirichlet

  contains
    procedure :: init
  end type flow_bc

contains

  subroutine init(this, mesh, p)
    class(flow_bc), intent(out) :: this
    type(flow_mesh), pointer, intent(in) :: mesh
    type(parameter_list), intent(inout) :: p
    type(flow_bc_factory) :: f

    call f%init(mesh, p)
    call f%alloc_bc("velocity-dirichlet", v_dirichlet)
    call f%alloc_bc("pressure-dirichlet", p_dirichlet)

  end subroutine init


end module flow_bc_type
