program test_region_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size
  use region_func_type
  use parameter_list_type
  use parameter_list_json
  use unstr_2d_mesh_type
  use vol_frac_init_procs
  use string_utilities, only: i_to_c
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use parallel_communication
  use simulation_environment_type
  implicit none

  type(unstr_2d_mesh), pointer :: mesh => null()
  integer :: status, stat
  character(:), allocatable :: errmsg
  type(simulation_environment) :: env

  !! The mesh component is parallel so we must initialize MPI, although
  !! we will only run in serial.
  call init_parallel_communication
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_SILENT)
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  call env%simlog%init(env%comm, 'test_compute_volume_fractions.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) call TLS_fatal('initializing simulation log: ' // errmsg)

  status = 0

  call create_mesh1(mesh)
  call test_cell_set
  call test1

  call create_mesh2(mesh)
  call test2

  if (global_any(status /= 0)) stop 1

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
    mesh => new_unstr_2d_mesh(env, [-1.0_r8, -1.0_r8], [1.0_r8, 1.0_r8], [4,4], ptri=0.5_r8)
    call mesh%init_cell_centroid
    mesh%cell_set_id = [1,2,3]
    do j = 1, mesh%ncell
      if (mesh%cell_centroid(1,j) > 0) mesh%cell_set_mask(j) = ibset(0,2)
    end do
  end subroutine

  subroutine test_cell_set

    type(parameter_list), pointer :: params
    integer :: stat, j
    character(:), allocatable :: string, errmsg
    type(region_func) :: rfunc
    real(r8), allocatable :: vol_frac(:,:)

    string = '{"selected":{"type":"cell-set","cell-set-ids":[2]},&
              &"background":{"type":"background"}}'

    call parameter_list_from_json_string(string, params, errmsg)
    if (.not.associated(params)) then
      call write_fail('test_cell_set: ' // errmsg)
      return
    end if

    call rfunc%init(mesh, params, stat, errmsg)
    if (stat /= 0) then
      call write_fail('test_cell_set: ' // errmsg)
      return
    end if

    allocate(vol_frac(rfunc%num_region(),mesh%ncell_onP))
    call compute_volume_fractions(mesh, rfunc, 0, vol_frac, stat)
    if (stat /= 0) then
      call write_fail('test_cell_set: failed to compute volume fractions')
      return
    end if
    call check_fraction_invariants(vol_frac, 'test_cell_set')

    do j = 1, mesh%ncell_onP
      if (btest(mesh%cell_set_mask(j), 2)) then
        if (vol_frac(1,j) /= 1.0_r8 .or. vol_frac(2,j) /= 0.0_r8) then
          call write_fail('test_cell_set: selected cell is not exact')
        end if
      else
        if (vol_frac(1,j) /= 0.0_r8 .or. vol_frac(2,j) /= 1.0_r8) then
          call write_fail('test_cell_set: background cell is not exact')
        end if
      end if
    end do

  end subroutine test_cell_set

  !! This is just like create_mesh1 except that we randomly perturb the node
  !! positions. We also don't modify the cell set data (and wont use it in the
  !! test) as it doesn't really work with a perturbed mesh.

  subroutine create_mesh2(mesh)
    use unstr_2d_mesh_factory
    type(unstr_2d_mesh), pointer :: mesh
    if (associated(mesh)) deallocate(mesh)
    mesh => new_unstr_2d_mesh(env, [-1.0_r8, -1.0_r8], [1.0_r8, 1.0_r8], [4,4], eps=0.1_r8, ptri=0.5_r8)
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
    allocate(vol_frac(rfunc%num_region(),mesh%ncell_onP))
    call compute_volume_fractions(mesh, rfunc, rlev, vol_frac, stat)
    if (stat /= 0) then
      call write_fail('test1: ' // errmsg)
      return
    end if
    call check_fraction_invariants(vol_frac, 'test1')

    !! Expected volumes
    v(1) = 4*atan(1.0_r8) * 0.75_r8 ** 2
    v(2) = 2 - v(1)/2
    v(3) = v(2)/2
    v(4) = v(3)

    if (is_IOP) write(*,'(a,i0,a)') 'test1: using ', rlev, ' recursion levels'
    if (is_IOP) write(*,'(a)') 'test1: expecting single precision accuracy in the total region volumes'
    do j = 1, size(vol_frac,1)
      vol = global_sum(dot_product(mesh%volume(:mesh%ncell_onP), vol_frac(j,:)))
      err = abs(vol-v(j))
      if (is_IOP) write(*,'(a,g0)') 'test1: volume error for region ' // i_to_c(j) // ': ', err
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
    allocate(vol_frac(rfunc%num_region(),mesh%ncell_onP))
    call compute_volume_fractions(mesh, rfunc, rlev, vol_frac, stat)
    if (stat /= 0) then
      call write_fail('test2: ' // errmsg)
      return
    end if
    call check_fraction_invariants(vol_frac, 'test2')

    !! Expected volumes
    v(1) = 4*atan(1.0_r8) * 0.75_r8 ** 2
    v(2) = 2 - v(1)/2
    v(3) = v(2)

    if (is_IOP) write(*,'(a,i0,a)') 'test2: using ', rlev, ' recursion levels'
    if (is_IOP) write(*,'(a)') 'test2: expecting single precision accuracy in the total region volumes'
    do j = 1, size(vol_frac,1)
      vol = global_sum(dot_product(mesh%volume(:mesh%ncell_onP), vol_frac(j,:)))
      err = abs(vol-v(j))
      if (is_IOP) write(*,'(a,g0)') 'test2: volume error for region ' // i_to_c(j) // ': ', err
      if (err > epsilon(1.0)) call write_fail('test2: wrong volume for region ' // i_to_c(j))
    end do

  end subroutine


  subroutine check_fraction_invariants(vol_frac, label)
    real(r8), intent(in) :: vol_frac(:,:)
    character(*), intent(in) :: label

    if (global_any(any(vol_frac < 0.0_r8) .or. any(vol_frac > 1.0_r8))) then
      call write_fail(label // ': volume fractions are not bounded')
    end if
    if (global_any(any(abs(sum(vol_frac, dim=1) - 1.0_r8) > 16.0_r8 * epsilon(1.0_r8)))) then
      call write_fail(label // ': volume fractions do not sum to one')
    end if
  end subroutine check_fraction_invariants

  subroutine write_fail(errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
    character(*), intent(in) :: errmsg
    status = 1
    write(error_unit,'(a)') errmsg
  end subroutine

end program test_region_func_type
