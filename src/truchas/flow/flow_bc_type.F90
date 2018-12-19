#include "f90_assert.fpp"
module flow_bc_type

  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  use parameter_list_type
  use bndry_func_class
  use bndry_vfunc_class
  use flow_surface_tension_bc_type
  use unstr_mesh_type
  use flow_bc_factory_type
  implicit none
  private

  type, public :: flow_bc
    class(bndry_func), allocatable :: p_dirichlet, dp_dirichlet, p_neumann, v_zero_normal
    class(bndry_vfunc), allocatable :: v_dirichlet
    type(surface_tension_bc) :: surface_tension
    logical :: pressure_d
  contains
    procedure :: init
    procedure :: compute
    procedure :: compute_initial
    procedure :: is_p_neumann_fix_pe
    procedure, private :: apply_default
  end type flow_bc

contains

  subroutine init(this, mesh, params, vof, temperature_fc)

    use parameter_list_type
    use flow_props_type
    use parallel_communication, only: global_sum

    class(flow_bc), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout) :: params
    real(r8), intent(in), target :: vof(:), temperature_fc(:)

    type(flow_bc_factory) :: f
    type(parameter_list), pointer :: plist

    ! need to generalize this to allow user input:
    ! - no slip walls (=> velocity-dirichlet + pressure_neumann)
    ! - free slip walls (=> v_zero_normal + pressure_neumann)
    ! - pressure inlet/outlet (=> pressure_dirichlet + v_neumann)
    ! - velocity inlet (=> v_dirichlet + pressure_neumann)

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
        [character(len=32) :: "pressure neumann", "no slip", "slip", "surface tension"], &
        this%p_neumann, default=0.0_r8)
    call f%alloc_scalar_bc([character(len=32) :: "slip", "surface tension"], &
        this%v_zero_normal, default=0.0_r8)

    call this%surface_tension%init(plist, mesh, vof, temperature_fc)

    ! apply the default boundary condition on faces with no user-specified BC
    call this%apply_default(mesh)

    this%pressure_d = global_sum(size(this%p_dirichlet%index)) > 0

#ifndef NDEBUG
    print *, "size of p dirichlet: ", size(this%p_dirichlet%index)
    print *, "size of p neumann: ", size(this%p_neumann%index)
    print *, "size of v dirichlet: ", size(this%v_dirichlet%index)
    print *, "size of surface tension: ", size(this%surface_tension%index)
#endif
  end subroutine init

  ! Assigns an MPI rank to do the pressure-neumann null space fixup.
  ! The rank must have fluid in the current time step, and otherwise
  ! neumann pressure boundaries.
  logical function is_p_neumann_fix_pe(this, any_real_fluid_onP)

    use parallel_communication
    use f08_intrinsics, only: findloc

    class(flow_bc), intent(inout) :: this
    logical, intent(in) :: any_real_fluid_onP

    integer :: fix_pe
    logical :: is_valid_pe, has_p_neumann(nPE)

    is_p_neumann_fix_PE = .false.

    if (.not.this%pressure_d) then
      ! find the lowest PE with a pressure neumann boundary
      is_valid_pe = size(this%p_neumann%index) > 0 .and. any_real_fluid_onP
      call collate(has_p_neumann, is_valid_pe)
      if (is_IOP) then
        fix_pe = findloc(has_p_neumann, .true.)
        INSIST(fix_pe > 0 .and. fix_pe <= nPE)
      end if
      call broadcast(is_valid_pe)
      is_p_neumann_fix_PE = (fix_pe == this_PE)
    end if

  end function is_p_neumann_fix_pe

  ! apply the default boundary condition (slip) to
  ! boundary faces with unspecified conditions
  subroutine apply_default(this, mesh)

    use bndry_face_func_type
    use scalar_func_class
    use scalar_func_factories, only: alloc_const_scalar_func

    class(flow_bc), intent(inout) :: this
    type(unstr_mesh), intent(in) :: mesh

    integer :: f, nf
    integer, allocatable :: faces(:)
    class(scalar_func), allocatable :: func

    ! get list of boundary faces with no user-applied BC
    allocate(faces(mesh%nface))

    nf = 0
    do f = 1, mesh%nface_onP
      if (mesh%fcell(2,f) /= 0) cycle ! skip internal faces

      ! skip faces with user-given BC
      if (any(this%p_dirichlet%index == f)) cycle
      if (any(this%p_neumann%index == f)) cycle ! catches v_dirichlet & surface_tension

      nf = nf + 1
      faces(nf) = f
    end do

    if (nf == 0) return

#ifndef NDEBUG
    print '(a,i8,a)', "applying default BC to ", nf, " faces"
#endif

    ! apply slip bc (p_neumann + v_zero_normal) to faces in faces(:nf)
    call alloc_const_scalar_func(func, 0.0_r8)

    select type (bc => this%p_neumann)
    type is (bndry_face_func)
      call bc%add_face_list(func, faces(:nf))
    class default
      ASSERT(.false.)
    end select

    select type (bc => this%v_zero_normal)
    type is (bndry_face_func)
      call bc%add_face_list(func, faces(:nf))
    class default
      ASSERT(.false.)
    end select

  end subroutine apply_default

  subroutine compute(this, t, dt)
    class(flow_bc), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    call this%p_dirichlet%compute(t)
    call this%p_neumann%compute(t)
    call this%v_dirichlet%compute(t)
    call this%surface_tension%compute(t)
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
    call this%surface_tension%compute(t)
    this%dp_dirichlet%value = 0.0_r8
  end subroutine compute_initial

end module flow_bc_type
