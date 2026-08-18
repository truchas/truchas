program test_flow_2d_bc

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use parameter_list_type
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use flow_2d_bc_type
  implicit none

  integer :: status

  call init_parallel_communication
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NORMAL)

  status = 0
  call test_pressure_pin

  call halt_parallel_communication
  stop status

contains

  subroutine test_pressure_pin
    type(unstr_2d_mesh), pointer :: mesh
    type(parameter_list), target :: velocity_params, pressure_params
    type(parameter_list), pointer :: plist
    type(flow_2d_bc) :: bc
    character(:), allocatable :: errmsg
    integer :: stat, pin_face, f
    logical :: defaults_complete

    mesh => new_unstr_2d_mesh([0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [4, 4], 0.0_r8, 0.0_r8)
    plist => velocity_params%sublist('wall')
    call plist%set('type', 'no-slip')
    call plist%set('face-set-ids', [1])

    call bc%init(mesh, velocity_params, stat, errmsg)
    call require(stat == 0, 'velocity boundary condition initialization failed')
    call bc%compute(0.0_r8, 1.0_r8)
    defaults_complete = .true.
    do f = 1, mesh%nface_onP
      if (mesh%fcell(2,f) /= 0) cycle
      defaults_complete = defaults_complete .and. any(bc%pressure_neumann%index == f) .and. &
          (any(bc%velocity_dirichlet%index == f) .or. any(bc%velocity_zero_normal%index == f))
    end do
    call require(defaults_complete, 'default boundary conditions do not cover every boundary face')
    if (size(bc%pressure_correction_dirichlet%value) > 0) then
      call require(maxval(abs(bc%pressure_correction_dirichlet%value)) < 1.0e-12_r8, &
          'static pressure data produced a nonzero pressure correction')
    end if
    pin_face = bc%pressure_pin_face()
    call require(global_sum(merge(1, 0, pin_face > 0)) == 1, &
        'all-Neumann pressure conditions did not select exactly one pin face')

    plist => pressure_params%sublist('outlet')
    call plist%set('type', 'pressure')
    call plist%set('face-set-ids', [1])
    call plist%set('pressure', 0.0_r8)
    call bc%init(mesh, pressure_params, stat, errmsg)
    call require(stat == 0, 'pressure boundary condition initialization failed')
    pin_face = bc%pressure_pin_face()
    call require(pin_face == 0, 'pressure Dirichlet condition should suppress pressure pinning')
  end subroutine


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine

end program test_flow_2d_bc
