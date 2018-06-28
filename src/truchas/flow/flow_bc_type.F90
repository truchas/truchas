module flow_bc_type

  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  use parameter_list_type
  use bndry_func_class
  use bndry_vfunc_class
  use flow_mesh_type
  use flow_bc_factory_type
  use parallel_communication
  implicit none
  private

  public :: flow_bc

  type :: flow_bc
    class(bndry_func), allocatable :: p_dirichlet, p_neumann, v_zero_normal
    class(bndry_vfunc), allocatable :: v_dirichlet
    logical :: pressure_d
  contains
    procedure :: init
    procedure :: compute
  end type flow_bc

contains

  subroutine init(this, mesh, p)
    class(flow_bc), intent(out) :: this
    type(flow_mesh), pointer, intent(in) :: mesh
    type(parameter_list), pointer, intent(in) :: p
!!$    type(parameter_list), pointer :: pn, pp
    type(flow_bc_factory) :: f

    ! need to generalize this to allow user input:
    ! - no slip walls (=> velocity-dirichlet + pressure_neumann)
    ! - free slip walls (=> v_zero_normal + pressure_neumann)
    ! - pressure inlet/outlet (=> pressure_dirichlet + v_neumann)
    ! - velocity inlet (=> v_dirichlet + pressure_neumann)


    call f%init(mesh, p)
    call f%alloc_vector_bc("velocity-dirichlet", this%v_dirichlet)
    call f%alloc_scalar_bc("pressure-dirichlet", this%p_dirichlet)
    ! need to have a p_neumann here... to complement velocity-dirichlet

!!$    allocate(pn)
!!$    pp => pn%sublist("")
    call f%alloc_scalar_bc("pressure-neumann", this%p_neumann)

    this%pressure_d = global_sum(size(this%p_dirichlet%index)) > 0

  end subroutine init

  subroutine compute(this, t)
    class(flow_bc), intent(inout) :: this
    real(r8), intent(in) :: t

    call this%p_dirichlet%compute(t)
    call this%p_neumann%compute(t)
    call this%v_dirichlet%compute(t)
    ! skip compute call for v_zero_normal for now.  perhaps there is a better abstraction
    ! for this, some sort of cell_face_boundary_condition which takes the current value
    ! at a cell and returns the appropriate value at a face...
  end subroutine compute


end module flow_bc_type
