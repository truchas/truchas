program test_region_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use region_func_type
  use parameter_list_type
  use parameter_list_json
  use unstr_2d_mesh_type
  use string_utilities, only: i_to_c
  implicit none

  type(unstr_2d_mesh), pointer :: mesh
  integer :: status

  !! Create a mesh for testing (only accessed by test_cell_set_region)
  !! The mesh component is parallel so we must initialize MPI, although
  !! we will only run in serial.
  block
    use unstr_2d_mesh_factory
    use parallel_communication
    use truchas_env, only: prefix, overwrite_output
    use truchas_logging_services
    call init_parallel_communication
    prefix = 'run'
    overwrite_output = .true.
    call TLS_initialize
    call TLS_set_verbosity(TLS_VERB_SILENT)
    mesh => new_unstr_2d_mesh([-1.0_r8, -1.0_r8], [1.0_r8, 1.0_r8], [4,4])
    call config_mesh_cell_sets(mesh)
  end block

  call config_mesh_cell_sets(mesh)
  call test1

  if (status /= 0) stop 1

contains

  !! Modify the default cell set data (all cells belong to a single cell set 1)
  !! such that cells to the left of the y-axis belong to cell set 1 and those
  !! to the right to cell set 2. Also add a cell set 3 to which no cell belongs.

  subroutine config_mesh_cell_sets(mesh)
    type(unstr_2d_mesh), intent(inout) :: mesh
    integer :: j
    call mesh%init_cell_centroid
    mesh%cell_set_id = [1,2,3]
    do j = 1, mesh%ncell
      if (mesh%cell_centroid(1,j) > 0) mesh%cell_set_mask(j) = ibset(0,2)
    end do
  end subroutine

  subroutine test1

    type(parameter_list), pointer :: params
    integer :: stat
    character(:), allocatable :: string, errmsg
    type(region_func) :: rfunc
    real(r8), allocatable :: point(:,:)
    integer, allocatable :: cellid(:), regid(:), bitmask(:)
    integer :: j

    string = '{"r1":{"type":"box","lower-corner":[0.25,0.25],"upper-corner":[0.75,0.5]},&
              &"r2":{"type":"disk","center":[0.0,0.0],"radius":0.75},&
              &"r3":{"type":"cell-set","cell-set-ids":[2]},&
              &"r4":{"type":"half-plane","point":[0.0,0.0],"normal":[0.0,1.0]},&
              &"r5":{"type":"background"}}'

    call parameter_list_from_json_string(string, params, errmsg)
    if (.not.associated(params)) then
      call write_fail('test1: ' // errmsg)
      return
    end if

    call rfunc%init(mesh, params, stat, errmsg)
    if (stat /= 0) then
      call write_fail('test1: ' // errmsg)
      return
    end if

    cellid = [1,4,7,10,13,16] ! sampling of cells
    regid  = [4,3,2,2,5,3]     ! the corresponding regions containing their center
    point  = mesh%cell_centroid(:,cellid)
    bitmask = mesh%cell_set_mask(cellid)

    do j = 1, size(cellid)
      if (rfunc%region_index(point(:,j), bitmask(j)) /= regid(j)) &
          call write_fail('wrong region for point ' // i_to_c(j))
    end do

    if (rfunc%region_index([0.5_r8,0.375_r8], 0) /= 1) call write_fail('wrong region for box point')

  end subroutine

  subroutine write_fail(errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
    character(*), intent(in) :: errmsg
    status = 1
    write(error_unit,'(a)') errmsg
  end subroutine

end program test_region_func_type
