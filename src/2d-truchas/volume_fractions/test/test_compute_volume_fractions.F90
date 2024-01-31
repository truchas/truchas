program test_region_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use region_func_type
  use parameter_list_type
  use parameter_list_json
  use unstr_2d_mesh_type
  use vol_frac_init_procs
  use string_utilities, only: i_to_c
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use parallel_communication
  implicit none

  type(unstr_2d_mesh), pointer :: mesh => null()
  integer :: status

  !! The mesh component is parallel so we must initialize MPI, although
  !! we will only run in serial.
  call init_parallel_communication
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_SILENT)

  block
  end block

  call create_mesh1(mesh)
  call test1

  call create_mesh2(mesh)
  call test2

  if (status /= 0) stop 1

contains

  !! Create a regular mesh with about half of the hex cells split into 2 tri
  !! cells. Modify the default cell set data (all cells belong to a single
  !! cell set 1) such that cells to the left of the y-axis belong to cell
  !! set 1 and those to the right to cell set 2. Also add a cell set 3 to
  !! which no cell belongs.

  subroutine create_mesh1(mesh)
    use unstr_2d_mesh_factory
    type(unstr_2d_mesh), pointer :: mesh
    integer :: j
    if (associated(mesh)) deallocate(mesh)
    mesh => new_unstr_2d_mesh([-1.0_r8, -1.0_r8], [1.0_r8, 1.0_r8], [4,4], ptri=0.5_r8)
    call mesh%init_cell_centroid
    mesh%cell_set_id = [1,2,3]
    do j = 1, mesh%ncell
      if (mesh%cell_centroid(1,j) > 0) mesh%cell_set_mask(j) = ibset(0,2)
    end do
  end subroutine

  !! This is just like create_mesh1 except that we randomly perturb the node
  !! positions. We also don't modify the cell set data (and wont use it in the
  !! test) as it doesn't really work with a perturbed mesh.

  subroutine create_mesh2(mesh)
    use unstr_2d_mesh_factory
    type(unstr_2d_mesh), pointer :: mesh
    integer :: j
    if (associated(mesh)) deallocate(mesh)
    mesh => new_unstr_2d_mesh([-1.0_r8, -1.0_r8], [1.0_r8, 1.0_r8], [4,4], eps=0.1_r8, ptri=0.5_r8)
  end subroutine

  subroutine test1

    type(parameter_list), pointer :: params
    integer :: stat
    character(:), allocatable :: string, errmsg
    type(region_func) :: rfunc
    real(r8), allocatable :: vol_frac(:,:)
    real(r8) :: v(4), err, vol
    integer :: j, rlev

    string = '{"r1":{"type":"disk","center":[0.0,0.0],"radius":0.75},&
              &"r2":{"type":"cell-set","cell-set-ids":[2]},&
              &"r3":{"type":"half-plane","point":[0.0,0.0],"normal":[0.0,1.0]},&
              &"r4":{"type":"background"}}'

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

    rlev = 15
    allocate(vol_frac(rfunc%num_region(),mesh%ncell))
    call compute_volume_fractions(mesh, rfunc, rlev, vol_frac, stat)
    if (stat /= 0) then
      call write_fail('test1: ' // errmsg)
      return
    end if

    !! Expected volumes
    v(1) = 4*atan(1.0_r8) * 0.75_r8 ** 2
    v(2) = 2 - v(1)/2
    v(3) = v(2)/2
    v(4) = v(3)

    write(*,'(a,i0,a)') 'test1: using ', rlev, ' recursion levels'
    write(*,'(a)') 'test1: expecting single precision accuracy in the total region volumes'
    do j = 1, size(vol_frac,1)
      vol = dot_product(mesh%volume, vol_frac(j,:))
      err = abs(vol-v(j))
      write(*,'(a,g0)') 'test1: volume error for region ' // i_to_c(j) // ': ', err
      if (err > epsilon(1.0)) call write_fail('test1: wrong volume for region ' // i_to_c(j))
    end do

  end subroutine

  subroutine test2

    type(parameter_list), pointer :: params
    integer :: stat
    character(:), allocatable :: string, errmsg
    type(region_func) :: rfunc
    real(r8), allocatable :: vol_frac(:,:)
    real(r8) :: v(3), err, vol
    integer :: j, rlev

    string = '{"r1":{"type":"disk","center":[0.0,0.0],"radius":0.75},&
              &"r2":{"type":"half-plane","point":[0.0,0.0],"normal":[1.0,7.0]},&
              &"r3":{"type":"background"}}'

    call parameter_list_from_json_string(string, params, errmsg)
    if (.not.associated(params)) then
      call write_fail('test2: ' // errmsg)
      return
    end if

    call rfunc%init(mesh, params, stat, errmsg)
    if (stat /= 0) then
      call write_fail('test2: ' // errmsg)
      return
    end if

    rlev = 15
    allocate(vol_frac(rfunc%num_region(),mesh%ncell))
    call compute_volume_fractions(mesh, rfunc, rlev, vol_frac, stat)
    if (stat /= 0) then
      call write_fail('test2: ' // errmsg)
      return
    end if

    !! Expected volumes
    v(1) = 4*atan(1.0_r8) * 0.75_r8 ** 2
    v(2) = 2 - v(1)/2
    v(3) = v(2)

    write(*,'(a,i0,a)') 'test2: using ', rlev, ' recursion levels'
    write(*,'(a)') 'test2: expecting single precision accuracy in the total region volumes'
    do j = 1, size(vol_frac,1)
      vol = dot_product(mesh%volume, vol_frac(j,:))
      err = abs(vol-v(j))
      write(*,'(a,g0)') 'test2: volume error for region ' // i_to_c(j) // ': ', err
      if (err > epsilon(1.0)) call write_fail('test2: wrong volume for region ' // i_to_c(j))
    end do

  end subroutine

  subroutine write_fail(errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
    character(*), intent(in) :: errmsg
    status = 1
    write(error_unit,'(a)') errmsg
  end subroutine

end program test_region_func_type
