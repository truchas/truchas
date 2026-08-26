program test_flow_2d_material_layout

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use parameter_list_json
  use material_database_type
  use material_model_type
  use material_composition_type
  use material_factory, only: load_material_database
  use flow_2d_material_layout_type
  implicit none

  type(material_database) :: database
  type(material_model) :: matl_model
  type(material_composition) :: composition
  type(flow_2d_material_layout) :: layout
  type(parameter_list), pointer :: matl_params, tracking_params
  character(:), allocatable :: errmsg
  real(r8) :: vfrac(4,1)
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
  call layout%init(matl_model, tracking_params, stat, errmsg)
  if (stat /= 0) error stop 'initializing material layout: ' // errmsg
  call require(layout%num_real_fluid() == 2, 'wrong real fluid count')
  call require(layout%num_fluid() == 3, 'wrong moving material count')
  call require(layout%num_material() == 4, 'wrong reduced material count')
  call layout%get_priority(priority)
  call require(all(priority == [4,3,2,1]), 'wrong full material priority')

  allocate(composition%vfrac(4,1))
  composition%vfrac(:,1) = [0.1_r8, 0.3_r8, 0.4_r8, 0.2_r8]
  call layout%get_reduced_volume_fractions(composition, vfrac)
  call require(maxval(abs(vfrac(:,1) - [0.1_r8,0.4_r8,0.2_r8,0.3_r8])) < 1.0e-14_r8, &
      'wrong reduced volume fractions')
  vfrac(:,1) = [0.15_r8, 0.35_r8, 0.25_r8, 0.25_r8]
  call layout%put_reduced_volume_fractions(vfrac, composition)
  call require(maxval(abs(composition%vfrac(:,1) - [0.15_r8,0.3_r8,0.35_r8,0.25_r8])) < 1.0e-14_r8, &
      'wrong accepted volume fractions')

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (.not.condition) error stop 'FAIL: ' // message
  end subroutine

end program test_flow_2d_material_layout
