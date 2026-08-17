program test_cell_matl_prop_func

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use parameter_list_json
  use unstr_2d_mesh_factory
  use material_database_type
  use material_model_type
  use material_factory, only: load_material_database
  use material_utilities, only: add_enthalpy_prop
  use material_composition_type
  use cell_matl_prop_func_type
  use parallel_communication
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  implicit none

  type(unstr_2d_mesh), pointer :: mesh
  type(material_database) :: matl_db
  type(material_model) :: matl_model
  type(material_composition), pointer :: composition
  type(cell_matl_prop_func) :: conductivity, enthalpy
  type(parameter_list), pointer :: params
  real(r8), allocatable :: state(:), value(:), deriv(:)
  character(:), allocatable :: errmsg
  integer :: stat, ncell

  call init_parallel_communication
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_SILENT)

  mesh => new_unstr_2d_mesh([0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [4,4])
  call parameter_list_from_json_string( &
    '{"mat-a":{"properties":{"conductivity":2.0,"density":3.0,"specific-heat":4.0}},&
     &"mat-b":{"properties":{"conductivity":10.0,"density":5.0,"specific-heat":6.0}}}', params, errmsg)
  if (.not.associated(params)) call fail('could not parse materials: ' // errmsg)
  call load_material_database(matl_db, params, stat, errmsg)
  if (stat /= 0) call fail('could not load materials: ' // errmsg)
  call matl_model%init(['mat-a', 'mat-b'], matl_db, stat, errmsg)
  if (stat /= 0) call fail('could not initialize material model: ' // errmsg)
  call add_enthalpy_prop(matl_model, stat, errmsg)
  if (stat /= 0) call fail('could not add enthalpy properties: ' // errmsg)

  allocate(composition)
  call composition%init_uniform(mesh, matl_model, 1, stat, errmsg)
  if (stat /= 0) call fail('could not initialize composition: ' // errmsg)
  composition%vfrac(1,:) = 0.25_r8
  composition%vfrac(2,:) = 0.75_r8

  call conductivity%init(matl_model, composition, 'conductivity', stat, errmsg)
  if (stat /= 0) call fail('could not initialize conductivity: ' // errmsg)
  call enthalpy%init(matl_model, composition, 'enthalpy', stat, errmsg)
  if (stat /= 0) call fail('could not initialize enthalpy: ' // errmsg)

  ncell = mesh%ncell_onP
  allocate(state(ncell+1), value(ncell), deriv(ncell))
  state = 2.0_r8

  call conductivity%compute_value(state, value)
  if (global_any(abs(value - 8.0_r8) > 16.0_r8 * epsilon(1.0_r8))) call fail('wrong mixed conductivity')
  call enthalpy%compute_value(state, value)
  if (global_any(abs(value - 51.0_r8) > 16.0_r8 * epsilon(1.0_r8))) call fail('wrong mixed enthalpy')
  call enthalpy%compute_deriv(state, 1, deriv)
  if (global_any(abs(deriv - 25.5_r8) > 16.0_r8 * epsilon(1.0_r8))) call fail('wrong mixed enthalpy derivative')

  call halt_parallel_communication

contains

  subroutine fail(message)
    character(*), intent(in) :: message
    if (is_IOP) write(*,'(a)') 'test_cell_matl_prop_func: ' // message
    error stop 1
  end subroutine fail

end program test_cell_matl_prop_func
