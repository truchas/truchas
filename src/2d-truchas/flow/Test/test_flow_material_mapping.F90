program test_flow_material_mapping

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use parameter_list_json
  use material_database_type
  use material_model_type
  use material_distribution_type
  use material_factory, only: load_material_database
  use flow_material_mapping_type
  implicit none

  type(material_database) :: database
  type(material_model) :: matl_model
  type(material_distribution) :: matl_dist
  type(flow_material_mapping) :: mapping
  type(material_database) :: phase_database
  type(material_model) :: phase_model
  type(material_distribution) :: phase_dist
  type(flow_material_mapping) :: phase_mapping
  type(parameter_list), pointer :: matl_params, tracking_params
  type(parameter_list), pointer :: phase_params
  character(:), allocatable :: errmsg
  real(r8) :: vfrac(4,1)
  real(r8) :: phase_vfrac(3,1), temperature(1)
  integer :: stat, priority(4)

  call parameter_list_from_json_string( &
      '{"water":{"properties":{"fluid":true}},"wall":{"properties":{}},"oil":{"properties":{"fluid":true}}}', &
      matl_params, errmsg)
  if (.not.associated(matl_params)) error stop 'parsing material database: ' // errmsg
  call load_material_database(database, matl_params, stat, errmsg)
  if (stat /= 0) error stop 'loading material database: ' // errmsg
  call matl_model%init([character(len=5) :: 'water', 'wall ', 'oil  ', 'VOID '], database, stat, errmsg)
  if (stat /= 0) error stop 'initializing material model: ' // errmsg

  call parameter_list_from_json_string( &
      '{"material-priority":["SOLID","VOID","oil","water"]}', tracking_params, errmsg)
  if (.not.associated(tracking_params)) error stop 'parsing tracking input: ' // errmsg
  call mapping%init(matl_model, stat, errmsg)
  if (stat /= 0) error stop 'initializing material mapping: ' // errmsg
  call mapping%set_priority(tracking_params, stat, errmsg)
  if (stat /= 0) error stop 'setting material priority: ' // errmsg
  call require(mapping%num_real_fluid() == 2, 'wrong real fluid count')
  call require(mapping%num_fluid() == 3, 'wrong moving material count')
  call require(mapping%num_material() == 4, 'wrong reduced material count')
  call mapping%get_priority(priority)
  call require(all(priority == [4,3,2,1]), 'wrong full material priority')

  allocate(matl_dist%vfrac(4,1))
  matl_dist%vfrac(:,1) = [0.1_r8, 0.3_r8, 0.4_r8, 0.2_r8]
  call mapping%get_reduced_volume_fractions(matl_dist, vfrac)
  call require(maxval(abs(vfrac(:,1) - [0.1_r8,0.4_r8,0.2_r8,0.3_r8])) < 1.0e-14_r8, &
      'wrong reduced volume fractions')
  vfrac(:,1) = [0.15_r8, 0.25_r8, 0.30_r8, 0.30_r8]
  call mapping%put_reduced_volume_fractions(vfrac, matl_dist)
  call require(maxval(abs(matl_dist%vfrac(:,1) - [0.15_r8,0.3_r8,0.25_r8,0.30_r8])) < 1.0e-14_r8, &
      'wrong accepted volume fractions')
  call require(abs(sum(matl_dist%vfrac(:,1)) - 1.0_r8) < 16.0_r8 * epsilon(1.0_r8), &
      'accepted volume fractions do not sum to one')

  call parameter_list_from_json_string( &
      '{"pcm":{"phases":{"solid":{},"liquid":{"fluid":true}},' // &
      '"phase-changes":{"solid:liquid":{"solidus-temp":0.0,"liquidus-temp":1.0,"latent-heat":1.0}}},' // &
      '"water":{"properties":{"fluid":true}}}', phase_params, errmsg)
  if (.not.associated(phase_params)) error stop 'parsing phase material database: ' // errmsg
  call load_material_database(phase_database, phase_params, stat, errmsg)
  if (stat /= 0) error stop 'loading phase material database: ' // errmsg
  call phase_model%init(['pcm  ', 'water'], phase_database, stat, errmsg)
  if (stat /= 0) error stop 'initializing phase material model: ' // errmsg
  call phase_mapping%init(phase_model, stat, errmsg)
  if (stat /= 0) error stop 'initializing phase material mapping: ' // errmsg
  call require(phase_mapping%num_real_fluid() == 2, 'wrong phase fluid count')
  call require(phase_mapping%num_material() == 3, 'wrong phase reduced count')

  allocate(phase_dist%vfrac(2,1))
  phase_dist%vfrac(:,1) = [0.6_r8, 0.4_r8]
  temperature = 0.5_r8
  call phase_mapping%get_phase_volume_fractions(phase_model, phase_dist, temperature, phase_vfrac)
  call require(maxval(abs(phase_vfrac(:,1) - [0.3_r8,0.4_r8,0.3_r8])) < 1.0e-14_r8, &
      'wrong phase-aware volume fractions')
contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (.not.condition) error stop 'FAIL: ' // message
  end subroutine

end program test_flow_material_mapping
