program test_2d_regions

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use region_class
  use region_factory
  use parameter_list_type
  use parameter_list_json
  use unstr_2d_mesh_type
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
    mesh => new_unstr_2d_mesh([-1.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [4,2])
    call config_mesh_cell_sets(mesh)
  end block

  status = 0

  call test_box_region
  call test_disk_region
  call test_half_plane_region
  call test_cell_set_region

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

  subroutine test_box_region

    class(region), allocatable :: r
    type(parameter_list) :: params
    real(r8), parameter :: lower(*) = [1, 1], upper(*)=[3,1]
    character(:), allocatable :: errmsg
    integer :: stat, unused

    call params%set('type', 'box')
    call params%set('lower-corner', lower)
    call params%set('upper-corner', upper)
    call alloc_region(r, mesh, params, stat, errmsg)
    if (stat /= 0) then
      call write_fail('test_box_region: '//errmsg)
      return
    end if

    if (.not.r%encloses(lower,unused)) call write_fail('test_box_region: failed test 1')
    if (r%encloses(lower-epsilon(lower),unused)) call write_fail('test_box_region: failed test 2')
    if (.not.r%encloses(upper,unused)) call write_fail('test_box_region: failed test 3')
    if (r%encloses(lower+epsilon(lower),unused)) call write_fail('test_box_region: failed test 4')

    call params%set('complement', .true.)
    call alloc_region(r, mesh, params, stat, errmsg)
    if (stat /= 0) then
      call write_fail('test_box_region: '//errmsg)
      return
    end if

    if (r%encloses(lower,unused)) call write_fail('test_box_region: failed test 5')
    if (.not.r%encloses(lower-epsilon(lower),unused)) call write_fail('test_box_region: failed test 6')
    if (r%encloses(upper,unused)) call write_fail('test_box_region: failed test 7')
    if (.not.r%encloses(lower+epsilon(lower),unused)) call write_fail('test_box_region: failed test 8')

  end subroutine

  subroutine test_disk_region

    class(region), allocatable :: r
    type(parameter_list) :: params
    character(:), allocatable :: errmsg
    integer :: stat, unused

    call params%set('type', 'disk')
    call params%set('center', [1.0_r8, 1.0_r8])
    call params%set('radius', 1.0_r8)
    call alloc_region(r, mesh, params, stat, errmsg)
    if (stat /= 0) then
      call write_fail('test_disk_region: '//errmsg)
      return
    end if

    if (.not.r%encloses([1.0_r8,0.0_r8],unused)) call write_fail('test_disk_region: failed test 1')
    if (r%encloses([0.0_r8,0.0_r8],unused)) call write_fail('test_disk_region: failed test 2')

    call params%set('complement', .true.)
    call alloc_region(r, mesh, params, stat, errmsg)
    if (stat /= 0) then
      call write_fail('test_disk_region: '//errmsg)
      return
    end if

    if (r%encloses([1.0_r8,0.0_r8],unused)) call write_fail('test_disk_region: failed test 3')
    if (.not.r%encloses([0.0_r8,0.0_r8],unused)) call write_fail('test_disk_region: failed test 4')

  end subroutine

  subroutine test_half_plane_region

    class(region), allocatable :: r
    type(parameter_list) :: params
    real(r8), parameter :: point(*) = [2,0], normal(*) = [1,2]
    character(:), allocatable :: errmsg
    integer :: stat, unused

    call params%set('type', 'half-plane')
    call params%set('point', point)
    call params%set('normal', normal)
    call alloc_region(r, mesh, params, stat, errmsg)
    if (stat /= 0) then
      call write_fail('test_half_plane_region: '//errmsg)
      return
    end if

    if (.not.r%encloses(point,unused)) call write_fail('test_half_plane_region: failed test 1')
    if (r%encloses(point+epsilon(point),unused)) call write_fail('test_half_plane_region: failed test 2')

    call params%set('complement', .true.)
    call alloc_region(r, mesh, params, stat, errmsg)
    if (stat /= 0) then
      call write_fail('test_half_plane_region: '//errmsg)
      return
    end if

    if (r%encloses(point,unused)) call write_fail('test_half_plane_region: failed test 3')
    if (r%encloses(point-epsilon(point),unused)) call write_fail('test_half_plane_region: failed test 4')

  end subroutine

  subroutine test_cell_set_region

    class(region), allocatable :: r
    type(parameter_list) :: params
    real(r8), parameter :: unused(*) = [0,0]
    character(:), allocatable :: errmsg
    integer :: stat, bitmask

    call params%set('type', 'cell-set')
    call params%set('cell-set-ids', [2,3])  ! right of y-axis
    call alloc_region(r, mesh, params, stat, errmsg)
    if (stat /= 0) then
      call write_fail('test_cell_set_region: '//errmsg)
      return
    end if

    bitmask = mesh%cell_set_mask(1) ! cell 1 should be left of the y-axis
    if (r%encloses(unused,bitmask)) call write_fail('test_cell_set_region: failed test 1')
    bitmask = mesh%cell_set_mask(mesh%ncell) ! last cell should be right of the y-axis
    if (.not.r%encloses(unused,bitmask)) call write_fail('test_cell_set_region: failed test 2')

    call params%set('complement', .true.)
    call alloc_region(r, mesh, params, stat, errmsg)
    if (stat /= 0) then
      call write_fail('test_cell_set_region: '//errmsg)
      return
    end if

    bitmask = mesh%cell_set_mask(1) ! cell 1 should be left of the y-axis
    if (.not.r%encloses(unused,bitmask)) call write_fail('test_cell_set_region: failed test 3')
    bitmask = mesh%cell_set_mask(mesh%ncell) ! last cell should be right of the y-axis
    if (r%encloses(unused,bitmask)) call write_fail('test_cell_set_region: failed test 4')

  end subroutine

  subroutine write_fail(errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
    character(*), intent(in) :: errmsg
    status = 1
    write(error_unit,'(a)') errmsg
  end subroutine

end program test_2d_regions
