#include "f90_assert.fpp"
module flow_bc_type

  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  use parameter_list_type
  use bndry_func_class
  use bndry_vfunc_class
  use unstr_mesh_type
  use flow_bc_factory_type
  use parallel_communication
  implicit none
  private

  type, public :: flow_bc
    class(bndry_func), allocatable :: p_dirichlet, dp_dirichlet, p_neumann, v_zero_normal
    class(bndry_vfunc), allocatable :: v_dirichlet
    logical :: pressure_d
    logical :: fix_neumann
    logical :: is_p_neumann_fix_PE
  contains
    procedure :: init
    procedure :: compute
    procedure :: compute_initial
  end type flow_bc

contains

  subroutine init(this, mesh, params)

    use parameter_list_type

    class(flow_bc), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout) :: params

    type(flow_bc_factory) :: f
    integer :: nc, flag, has_p_neumann(nPE)
    type(parameter_list), pointer :: plist

    ! need to generalize this to allow user input:
    ! - no slip walls (=> velocity-dirichlet + pressure_neumann)
    ! - free slip walls (=> v_zero_normal + pressure_neumann)
    ! - pressure inlet/outlet (=> pressure_dirichlet + v_neumann)
    ! - velocity inlet (=> v_dirichlet + pressure_neumann)

    this%fix_neumann = .false.

    plist => params%sublist('bc')
    call f%init(mesh, plist)
    call f%alloc_vector_bc( &
        [character(len=32) :: "velocity dirichlet", "no slip"], &
        this%v_dirichlet, &
        default=0.0_r8)
    call f%alloc_scalar_bc(["pressure dirichlet"], this%p_dirichlet)
    call f%alloc_scalar_bc(["pressure dirichlet"], this%dp_dirichlet)
    ! need to have a p_neumann here... to complement velocity-dirichlet

    call f%alloc_scalar_bc(&
        [character(len=32) :: "pressure neumann", "no slip", "slip"], this%p_neumann, default=0.0_r8)
    call f%alloc_scalar_bc(["slip"], this%v_zero_normal, default=0.0_r8)

    this%pressure_d = global_sum(size(this%p_dirichlet%index)) > 0

    if (.not.this%pressure_d) then
      ! find the lowest PE with a pressure neumann boundary
      flag = merge(1, 0, size(this%p_neumann%index) > 0)
      call collate(has_p_neumann, flag)
      if (is_IOP) then
        ! note: can replace this with Fortran 2008's findloc intrinsic once compilers support it
        do flag = 1,nPE
          if (has_p_neumann(flag) == 1) exit
        end do
        INSIST(flag <= nPE)
      end if
      call broadcast(flag)
      this%is_p_neumann_fix_PE = (flag == this_PE)
    end if
#ifndef NDEBUG
    print *, "size of p dirichlet: ", size(this%p_dirichlet%index)
    print *, "size of p neumann: ", size(this%p_neumann%index)
    print *, "size of v dirichlet: ", size(this%v_dirichlet%index)
#endif
  end subroutine init

  subroutine compute(this, t, dt)
    class(flow_bc), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    call this%p_dirichlet%compute(t)
    call this%p_neumann%compute(t)
    call this%v_dirichlet%compute(t)
    call this%dp_dirichlet%compute(t+dt)
    this%dp_dirichlet%value = this%dp_dirichlet%value - this%p_dirichlet%value
  end subroutine compute

  !! Compute BC data for calculating the initial flow state. Only the Dirichlet
  !! data for dp is different. Normally the data is the difference in pressure
  !! data over the time step, but in this case we are not taking an actual time
  !! step, so dp is 0 on pressure Dirichlet boundaries.

  subroutine compute_initial(this, t)
    class(flow_bc), intent(inout) :: this
    real(r8), intent(in) :: t
    call this%p_dirichlet%compute(t)
    call this%p_neumann%compute(t)
    call this%v_dirichlet%compute(t)
    this%dp_dirichlet%value = 0.0_r8
  end subroutine compute_initial

end module flow_bc_type
