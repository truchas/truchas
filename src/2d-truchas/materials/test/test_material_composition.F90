program test_material_composition

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size
  use parameter_list_type
  use parameter_list_json
  use unstr_2d_mesh_type
  use material_database_type
  use material_model_type
  use material_factory, only: load_material_database
  use material_composition_type
  use parallel_communication
  use simulation_environment_type
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  implicit none

  type(unstr_2d_mesh), pointer :: mesh => null()
  type(material_database), target :: matl_db
  type(material_model) :: matl_model
  type(material_composition) :: composition
  type(parameter_list), pointer :: matl_params, region_params
  character(:), allocatable :: string, errmsg, names(:)
  integer :: stat
  type(simulation_environment) :: env

  call init_parallel_communication
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_SILENT)
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  call env%simlog%init(env%comm, 'test_material_composition.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) call TLS_fatal('initializing simulation log: ' // errmsg)

  mesh => create_mesh()

  string = '{"mat-a":{"properties":{"conductivity":1.0,"density":1.0,"specific-heat":1.0}},&
            &"mat-b":{"properties":{"conductivity":2.0,"density":1.0,"specific-heat":1.0}},&
            &"mat-c":{"properties":{"conductivity":3.0,"density":1.0,"specific-heat":1.0}}}'
  call parameter_list_from_json_string(string, matl_params, errmsg)
  if (.not.associated(matl_params)) call fail('could not parse materials: ' // errmsg)
  call load_material_database(matl_db, matl_params, stat, errmsg)
  if (stat /= 0) call fail('could not load materials: ' // errmsg)

  string = '{"left-disk":{"material":"mat-a","type":"disk","center":[-0.5,0.0],"radius":0.3},&
            &"right-disk":{"material":"mat-a","type":"disk","center":[0.5,0.0],"radius":0.3},&
            &"lower-half":{"material":"mat-b","type":"half-plane","point":[0.0,0.0],"normal":[0.0,1.0]},&
            &"background":{"material":"mat-c","type":"background"}}'
  call parameter_list_from_json_string(string, region_params, errmsg)
  if (.not.associated(region_params)) call fail('could not parse material regions: ' // errmsg)

  call get_material_region_names(region_params, names, stat, errmsg)
  if (stat /= 0) call fail('could not get material region names: ' // errmsg)
  if (size(names) /= 3 .or. any(names /= ['mat-a', 'mat-b', 'mat-c'])) then
    call fail('wrong material-region name sequence')
  end if

  call matl_model%init(names, matl_db, stat, errmsg)
  if (stat /= 0) call fail('could not initialize material model: ' // errmsg)

  call composition%init(mesh, matl_model, region_params, 15, stat, errmsg)
  if (stat /= 0) call fail('could not initialize material composition: ' // errmsg)

  call check_composition(mesh, matl_model, composition)

contains

  function create_mesh() result(mesh)
    use unstr_2d_mesh_factory
    type(unstr_2d_mesh), pointer :: mesh
    mesh => new_unstr_2d_mesh(env, [-1.0_r8, -1.0_r8], [1.0_r8, 1.0_r8], [8,8], ptri=0.5_r8)
    call mesh%init_cell_centroid
  end function create_mesh


  subroutine check_composition(mesh, matl_model, composition)
    type(unstr_2d_mesh), intent(in) :: mesh
    type(material_model), intent(in) :: matl_model
    type(material_composition), intent(in) :: composition

    real(r8) :: expected(3), volume, radius
    logical :: pure(3)
    integer :: mid(3), i, j

    if (matl_model%nmatl /= 3) call fail('wrong number of model materials')
    if (any(composition%vfrac < 0.0_r8) .or. any(composition%vfrac > 1.0_r8)) then
      call fail('material fractions are not bounded')
    end if
    if (any(abs(sum(composition%vfrac, dim=1) - 1.0_r8) > 16.0_r8 * epsilon(1.0_r8))) then
      call fail('material fractions do not sum to one')
    end if

    mid = matl_model%matl_index(['mat-a', 'mat-b', 'mat-c'])
    if (any(mid == 0)) call fail('material index lookup failed')

    radius = 0.3_r8
    expected(1) = 2.0_r8 * 4.0_r8 * atan(1.0_r8) * radius**2
    expected(2) = 2.0_r8 - 4.0_r8 * atan(1.0_r8) * radius**2
    expected(3) = expected(2)
    do i = 1, size(mid)
      volume = global_sum(dot_product(mesh%volume(:mesh%ncell_onP), composition%vfrac(mid(i),:)))
      if (abs(volume-expected(i)) > epsilon(1.0)) call fail('wrong material volume')
    end do

    pure = .false.
    do j = 1, mesh%ncell_onP
      do i = 1, size(mid)
        pure(i) = pure(i) .or. composition%vfrac(mid(i),j) == 1.0_r8
      end do
    end do
    do i = 1, size(mid)
      if (.not.global_any(pure(i))) call fail('expected a pure cell for every material')
    end do
  end subroutine check_composition


  subroutine fail(message)
    character(*), intent(in) :: message
    if (is_IOP) write(*,'(a)') 'test_material_composition: ' // message
    error stop 1
  end subroutine fail

end program test_material_composition
