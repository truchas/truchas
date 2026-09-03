program test_ns_ht_2d_enthalpy_advector

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size
  use parallel_communication
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use simulation_environment_type
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use material_database_type
  use material_model_type
  use material_factory, only: load_material_database
  use material_utilities, only: add_enthalpy_prop
  use bndry_face_func_type
  use scalar_func_class
  use scalar_func_factories, only: alloc_const_scalar_func
  use parameter_list_type
  use parameter_list_json
  use ns_ht_2d_enthalpy_advector_type
  implicit none

  integer :: status, stat
  character(:), allocatable :: errmsg
  type(simulation_environment) :: env

  call init_parallel_communication
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_SILENT)
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  call env%simlog%init(env%comm, 'test_ns_ht_2d_enthalpy_advector.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) call TLS_fatal('initializing simulation log: ' // errmsg)

  status = 0
  call test_material_flux_enthalpy
  call test_closed_transport_conservation
  call test_multiphase_phase_enthalpy

  call halt_parallel_communication
  stop status

contains

  subroutine test_material_flux_enthalpy
    type(unstr_2d_mesh), pointer :: mesh
    type(material_database) :: database
    type(material_model) :: matl_model
    type(parameter_list), pointer :: matl_params
    type(ns_ht_2d_enthalpy_advector) :: advector
    type(bndry_face_func), target :: inflow_bc
    class(scalar_func), allocatable :: inflow_func
    real(r8), allocatable :: cell_temp(:), flux_volumes(:,:), dQ(:), expected(:)
    integer :: c, f, i, neighbor
    real(r8), parameter :: influx_temp = 7.0_r8
    real(r8), parameter :: q = 1.0e-3_r8

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 4], 0.0_r8, 0.0_r8)
    call mesh%init_face_centroid
    call parameter_list_from_json_string( &
      '{"mat-a":{"properties":{"density":3.0,"specific-heat":4.0}},&
       &"mat-b":{"properties":{"density":5.0,"specific-heat":6.0}}}', matl_params, errmsg)
    if (.not.associated(matl_params)) call fail('could not parse material input: ' // errmsg)
    call load_material_database(database, matl_params, stat, errmsg)
    if (stat /= 0) call fail('could not load material database: ' // errmsg)
    call matl_model%init(['mat-a', 'mat-b'], database, stat, errmsg)
    if (stat /= 0) call fail('could not initialize material model: ' // errmsg)
    call add_enthalpy_prop(matl_model, stat, errmsg)
    if (stat /= 0) call fail('could not add enthalpy property: ' // errmsg)

    call inflow_bc%init(mesh)
    call alloc_const_scalar_func(inflow_func, influx_temp)
    call inflow_bc%add(inflow_func, [1], stat, errmsg)
    if (stat /= 0) call fail('could not add inflow thermal BC: ' // errmsg)
    call inflow_bc%add_complete()
    call advector%init(mesh, matl_model, [1, 2], stat, errmsg, inflow_temperature=inflow_bc)
    if (stat /= 0) call fail('could not initialize enthalpy advector: ' // errmsg)

    allocate(cell_temp(mesh%ncell_onP), flux_volumes(2,size(mesh%cface)), dQ(mesh%ncell_onP), &
        expected(mesh%ncell_onP))
    cell_temp = 1.0_r8 + real(mesh%cell_imap%global_index([(c, c=1,mesh%ncell_onP)]), r8)
    flux_volumes = 0.0_r8
    do c = 1, mesh%ncell_onP
      do i = mesh%cstart(c), mesh%cstart(c+1)-1
        f = mesh%cface(i)
        flux_volumes(:,i) = q*[0.4_r8, 0.6_r8]*mesh%unit_normal(1,f)
        if (btest(mesh%cfpar(c), i-mesh%cstart(c)+1)) flux_volumes(:,i) = -flux_volumes(:,i)
      end do
    end do
    call advector%get_advected_enthalpy(0.0_r8, cell_temp, flux_volumes, dQ)

    expected = 0.0_r8
    do c = 1, mesh%ncell_onP
      do i = mesh%cstart(c), mesh%cstart(c+1)-1
        f = mesh%cface(i)
        if (any(flux_volumes(:,i) > 0.0_r8)) then
          expected(c) = expected(c) - dot_product([12.0_r8, 30.0_r8], flux_volumes(:,i))*cell_temp(c)
        else
          neighbor = mesh%cnhbr(i)
          if (neighbor > 0) then
            expected(c) = expected(c) - dot_product([12.0_r8, 30.0_r8], flux_volumes(:,i))* &
                (1.0_r8 + real(mesh%cell_imap%global_index(neighbor), r8))
          else if (mesh%face_centroid(1,f) == 0.0_r8) then
            expected(c) = expected(c) - dot_product([12.0_r8, 30.0_r8], flux_volumes(:,i))*influx_temp
          else
            expected(c) = expected(c) - dot_product([12.0_r8, 30.0_r8], flux_volumes(:,i))*cell_temp(c)
          end if
        end if
      end do
    end do
    call require(all(abs(dQ - expected) < 1.0e-13_r8), 'incorrect material-resolved enthalpy flux')
  end subroutine


  subroutine test_multiphase_phase_enthalpy
    type(unstr_2d_mesh), pointer :: mesh
    type(material_database) :: database
    type(material_model) :: matl_model
    type(parameter_list), pointer :: matl_params
    type(ns_ht_2d_enthalpy_advector) :: advector
    real(r8), allocatable :: cell_temp(:), flux_volumes(:,:), dQ(:), expected(:)
    real(r8), parameter :: q = 1.0e-3_r8

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 4], 0.0_r8, 0.0_r8)
    call parameter_list_from_json_string( &
      '{"pcm":{"properties":{"density":1.0},' // &
      '"phases":{"solid":{"specific-heat":2.0},' // &
      '"liquid":{"specific-heat":3.0}},' // &
      '"phase-changes":{"solid:liquid":{"solidus-temp":0.0,' // &
      '"liquidus-temp":1.0,"latent-heat":1.0}}}}', matl_params, errmsg)
    if (.not.associated(matl_params)) call fail('could not parse phase material input: ' // errmsg)
    call load_material_database(database, matl_params, stat, errmsg)
    if (stat /= 0) call fail('could not load phase material database: ' // errmsg)
    call matl_model%init(['pcm'], database, stat, errmsg)
    if (stat /= 0) call fail('could not initialize phase material model: ' // errmsg)
    call add_enthalpy_prop(matl_model, stat, errmsg)
    if (stat /= 0) call fail('could not add phase enthalpy property: ' // errmsg)
    call advector%init(mesh, matl_model, [2], stat, errmsg)
    if (stat /= 0) call fail('could not initialize phase enthalpy advector: ' // errmsg)

    allocate(cell_temp(mesh%ncell_onP), flux_volumes(1,size(mesh%cface)), &
        dQ(mesh%ncell_onP), expected(mesh%ncell_onP))
    cell_temp = 2.0_r8
    flux_volumes = q
    call advector%get_advected_enthalpy(0.0_r8, cell_temp, flux_volumes, dQ)
    expected = -4.0_r8*q*6.0_r8
    call require(all(abs(dQ - expected) < 1.0e-13_r8), &
        'phase enthalpy was not used for material-resolved enthalpy flux')
  end subroutine


  subroutine test_closed_transport_conservation
    type(unstr_2d_mesh), pointer :: mesh
    type(material_database) :: database
    type(material_model) :: matl_model
    type(parameter_list), pointer :: matl_params
    type(ns_ht_2d_enthalpy_advector) :: advector
    real(r8), allocatable :: cell_temp(:), flux_volumes(:,:), dQ(:)
    integer :: c, f, i
    real(r8), parameter :: q = 1.0e-3_r8

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 4], 0.0_r8, 0.0_r8)
    call parameter_list_from_json_string( &
      '{"mat-a":{"properties":{"density":3.0,"specific-heat":4.0}},&
       &"mat-b":{"properties":{"density":5.0,"specific-heat":6.0}}}', matl_params, errmsg)
    if (.not.associated(matl_params)) call fail('could not parse material input: ' // errmsg)
    call load_material_database(database, matl_params, stat, errmsg)
    if (stat /= 0) call fail('could not load material database: ' // errmsg)
    call matl_model%init(['mat-a', 'mat-b'], database, stat, errmsg)
    if (stat /= 0) call fail('could not initialize material model: ' // errmsg)
    call add_enthalpy_prop(matl_model, stat, errmsg)
    if (stat /= 0) call fail('could not add enthalpy property: ' // errmsg)

    call advector%init(mesh, matl_model, [1, 2], stat, errmsg)
    if (stat /= 0) call fail('could not initialize enthalpy advector: ' // errmsg)

    allocate(cell_temp(mesh%ncell_onP), flux_volumes(2,size(mesh%cface)), dQ(mesh%ncell_onP))
    cell_temp = 1.0_r8 + real(mesh%cell_imap%global_index([(c, c=1,mesh%ncell_onP)]), r8)
    flux_volumes = 0.0_r8
    do c = 1, mesh%ncell_onP
      do i = mesh%cstart(c), mesh%cstart(c+1)-1
        f = mesh%cface(i)
        if (mesh%fcell(2,f) == 0) cycle
        flux_volumes(:,i) = q*[0.4_r8, 0.6_r8]
        if (btest(mesh%cfpar(c), i-mesh%cstart(c)+1)) flux_volumes(:,i) = -flux_volumes(:,i)
      end do
    end do
    call advector%get_advected_enthalpy(0.0_r8, cell_temp, flux_volumes, dQ)
    call require(abs(global_sum(sum(dQ))) < 1.0e-12_r8, &
        'closed material transport does not conserve total enthalpy')
  end subroutine


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine


  subroutine fail(message)
    character(*), intent(in) :: message

    if (is_IOP) print '("ERROR: ",a)', message
    error stop 1
  end subroutine

end program test_ns_ht_2d_enthalpy_advector
