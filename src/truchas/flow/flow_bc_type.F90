#include "f90_assert.fpp"
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

  public :: flow_bc, read_flow_bc_namelist

  type :: flow_bc
    class(bndry_func), allocatable :: p_dirichlet, dp_dirichlet, p_neumann, v_zero_normal
    class(bndry_vfunc), allocatable :: v_dirichlet
    logical :: pressure_d
    logical :: fix_neumann
    type(parameter_list), pointer :: plist => null()
  contains
    procedure :: read_params
    procedure :: init
    procedure :: compute
  end type flow_bc

contains

  subroutine read_flow_bc_namelist(lun, p)
    use string_utilities, only: i_to_c
    use parallel_communication, only: is_IOP, broadcast
    use flow_input_utils

    integer, intent(in) :: lun
    type(parameter_list), pointer, intent(inout) :: p
    type(parameter_list), pointer :: pp, bc
    integer :: ios
    logical :: found
    character(128) :: iom
    character(128) :: condition
    character(3) :: bc_label
    real(r8) :: data(3)
    integer :: face_sets(10), i

    namelist /flow_bc/ condition, data, face_sets

    bc => p%sublist("bc")

    if (is_IOP) then
      rewind lun
    end if

    i = 1
    do while(.true.)
      ! advance to next flow bc_namelist or exit
      if (is_IOP) then
        call seek_to_namelist(lun, 'flow_bc', found, iostat=ios)
      end if
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('Error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      call TLS_info("")
      call TLS_info("Reading FLOW_BC namelist...")

      !! Reset defaults and read the namelist.
      condition = null_c
      data = null_r
      face_sets = null_i

      if (is_IOP) then
        read(lun,nml=flow_bc,iostat=ios,iomsg=iom)
      end if
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading FLOW_BC namelist: ' // trim(iom))

      call broadcast(condition)
      call broadcast(face_sets)
      call broadcast(data)

      write(bc_label,'("bc",i1)') i
      pp => bc%sublist(bc_label)
      call plist_set_if(pp, 'condition', condition)
      call plist_set_if(pp, 'data', data)
      call plist_set_if(pp, 'face sets', face_sets)

      i = i + 1
    end do

  end subroutine read_flow_bc_namelist

  subroutine read_params(this, p)
    class(flow_bc), intent(inout) :: this
    type(parameter_list), pointer, intent(in) :: p

    this%plist => p%sublist("bc")
#ifndef NDEBUG
    print *, "size of bc plist", this%plist%count()
#endif
  end subroutine read_params

  subroutine init(this, mesh)
    class(flow_bc), intent(inout) :: this
    type(flow_mesh), pointer, intent(in) :: mesh
    type(flow_bc_factory) :: f
    integer :: nc
    integer, allocatable :: neumann_count(:)

    ! need to generalize this to allow user input:
    ! - no slip walls (=> velocity-dirichlet + pressure_neumann)
    ! - free slip walls (=> v_zero_normal + pressure_neumann)
    ! - pressure inlet/outlet (=> pressure_dirichlet + v_neumann)
    ! - velocity inlet (=> v_dirichlet + pressure_neumann)

    ASSERT(associated(this%plist))

    this%fix_neumann = .false.

    call f%init(mesh, this%plist)
    call f%alloc_vector_bc( &
        [character(len=32) :: "velocity dirichlet", "no slip"], &
        this%v_dirichlet, &
        default=0.0_r8)
    call f%alloc_scalar_bc(["pressure dirichlet"], this%p_dirichlet)
    call f%alloc_scalar_bc(["pressure dirichlet"], this%dp_dirichlet)
    ! need to have a p_neumann here... to complement velocity-dirichlet

    call f%alloc_scalar_bc(&
        [character(len=32) :: "pressure neumann", "no slip"], this%p_neumann, default=0.0_r8)
    call f%alloc_scalar_bc(["slip"], this%v_zero_normal, default=0.0_r8)

    this%pressure_d = global_sum(size(this%p_dirichlet%index)) > 0

    if (.not.this%pressure_d) then
      ! this should be replaced with something smarter
      if (is_IOP) then
        INSIST(size(this%p_neumann%index) > 0)
      end if
    end if
#ifndef NDEBUG
    print *, "size of p dirichlet: ", size(this%p_dirichlet%index)
    print *, "size of p neumann: ", size(this%p_neumann%index)
    print *, "size of v dirichlet: ", size(this%v_dirichlet%index)
#endif
  end subroutine init

  subroutine compute(this, t, dt, initial)
    class(flow_bc), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    logical, optional, intent(in) :: initial
    !-
    logical :: init
    integer :: i

    if (present(initial)) then
      init = initial
    else
      init = .false.
    end if

    call this%p_dirichlet%compute(t)
    call this%p_neumann%compute(t)
    call this%v_dirichlet%compute(t)


    if (init) then
      call this%dp_dirichlet%compute(t)
    else
      call this%dp_dirichlet%compute(t+dt)
      this%dp_dirichlet%value(:) = this%dp_dirichlet%value(:) - this%p_dirichlet%value(:)
    end if

    ! skip compute call for v_zero_normal for now.  perhaps there is a better abstraction
    ! for this, some sort of cell_face_boundary_condition which takes the current value
    ! at a cell and returns the appropriate value at a face...
  end subroutine compute


end module flow_bc_type
