program test_flow_2d_mesh_factory

  use parallel_communication
  use parameter_list_type
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  implicit none

  type(parameter_list) :: params
  type(unstr_2d_mesh), pointer :: mesh
  character(512) :: path
  character(:), allocatable :: errmsg
  integer :: stat

  call init_parallel_communication
  call get_command_argument(1, path)
  call params%set('mesh-file', trim(path))
  call params%set('exodus-block-modulus', 1000)
  mesh => new_unstr_2d_mesh(params, stat, errmsg)

  if (stat /= 0 .or. .not.associated(mesh)) then
    if (is_IOP) print '(2a)', 'FAIL: ', errmsg
    call halt_parallel_communication
    stop 1
  end if
  if (size(mesh%cell_set_id) /= 1 .or. mesh%cell_set_id(1) /= 1) then
    if (is_IOP) print '(a)', 'FAIL: block IDs were not normalized to cell set 1'
    call halt_parallel_communication
    stop 1
  end if

  if (is_IOP) print '(a)', 'PASS: block IDs normalized to cell set 1'
  call halt_parallel_communication

end program test_flow_2d_mesh_factory
