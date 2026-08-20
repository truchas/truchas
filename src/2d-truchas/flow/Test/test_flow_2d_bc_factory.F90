program test_flow_2d_bc_factory

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use parameter_list_type
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use bndry_func1_class
  use bndry_vfunc_class
  use flow_2d_bc_factory_type
  implicit none

  integer :: status

  call init_parallel_communication
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NORMAL)

  status = 0
  call test_factory

  call TLS_finalize
  call halt_parallel_communication
  stop status

contains

  subroutine test_factory
    type(unstr_2d_mesh), pointer :: mesh
    type(parameter_list), target :: velocity_params, pressure_params, slip_params
    type(parameter_list), pointer :: plist
    type(flow_2d_bc_factory) :: factory
    class(bndry_vfunc), allocatable :: velocity
    class(bndry_func1), allocatable :: pressure, zero_normal, pressure_neumann
    character(:), allocatable :: errmsg
    integer :: stat

    mesh => new_unstr_2d_mesh([0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [4, 4], 0.0_r8, 0.0_r8)

    plist => velocity_params%sublist('inlet')
    call plist%set('type', 'velocity')
    call plist%set('face-set-ids', [1])
    call plist%set('velocity', [1.5_r8, -0.75_r8])
    plist => pressure_params%sublist('outlet')
    call plist%set('type', 'pressure')
    call plist%set('face-set-ids', [1])
    call plist%set('pressure', 3.0_r8)
    plist => slip_params%sublist('symmetry')
    call plist%set('type', 'free-slip')
    call plist%set('face-set-ids', [1])

    call factory%init(mesh, velocity_params)
    call factory%alloc_dir_vel_bc(velocity, stat, errmsg)
    call require(stat == 0 .and. allocated(velocity), 'velocity BC was not allocated')
    call velocity%compute(0.0_r8)
    call require(global_any(size(velocity%index) > 0), 'velocity BC has no faces')
    call require(all(abs(velocity%value - spread([1.5_r8, -0.75_r8], 2, size(velocity%index))) < 1.0e-12_r8), &
        'velocity BC values are incorrect')

    call factory%init(mesh, pressure_params)
    call factory%alloc_dir_prs_bc(pressure, stat, errmsg)
    call require(stat == 0 .and. allocated(pressure), 'pressure BC was not allocated')
    call pressure%compute(0.0_r8)
    call require(all(abs(pressure%value - 3.0_r8) < 1.0e-12_r8), 'pressure BC values are incorrect')

    call factory%init(mesh, slip_params)
    call factory%alloc_zero_vn_bc(zero_normal, stat, errmsg)
    call require(stat == 0 .and. allocated(zero_normal), 'free-slip BC was not allocated')
    call zero_normal%compute(0.0_r8)
    call require(all(abs(zero_normal%value) < 1.0e-12_r8), 'free-slip BC values are incorrect')

    call factory%init(mesh, velocity_params)
    call factory%alloc_neu_prs_bc(pressure_neumann, stat, errmsg)
    call require(stat == 0 .and. allocated(pressure_neumann), 'pressure Neumann BC was not allocated')
    call pressure_neumann%compute(0.0_r8)
    call require(all(abs(pressure_neumann%value) < 1.0e-12_r8), 'pressure Neumann BC values are incorrect')
  end subroutine


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine

end program test_flow_2d_bc_factory
