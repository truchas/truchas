!!
!! HTSD_MODEL_FACTORY
!!
!! This module provides a method for initializing an instance of an HTSD_MODEL
!! variable.  It depends essentially on Truchas-specific input data, and
!! encapsulates all the gory details of how the components of the HTSD_MODEL
!! object are defined.  This is also where the bulk of the HT/SD model
!! verification is done.
!!
!! PROGRAMMING INTERFACE
!!
!!  CREATE_HTSD_MODEL (DISC, MMF, STAT, ERRMSG) returns a pointer to a new,
!!    fully-defined, HTSD_MODEL-type object.  DISC is a pointer to the MFD
!!    discretization and MMF is a pointer to a corresponding material mesh
!     function.  The integerSTAT returns a non-zero value if an error is
!!    encountered, and an explanatory message is returned in the character
!!    string ERRMSG.
!!

#include "f90_assert.fpp"

module HTSD_model_factory

  use HTSD_model_type
  use dist_mesh_type
  use mfd_disc_type
  use material_mesh_function
  implicit none
  private

  public :: create_HTSD_model

contains

  function create_HTSD_model (disc, mmf, stat, errmsg) result (model)

    use diffusion_solver_data, only: heat_eqn, num_species
    use property_data_module, only: void_temperature
    use material_interop, only: void_material_index

    type(mfd_disc), intent(in), target :: disc
    type(mat_mf), intent(in), target :: mmf
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg
    type(HTSD_model), pointer :: model

    type(HT_model), pointer :: htmodel => null()
    type(SD_model), pointer :: sdmodel(:) => null()
#ifdef INTEL_COMPILER_WORKAROUND
    type(dist_mesh), pointer :: fubar

    fubar => mmf_mesh(mmf)
    ASSERT(associated(disc%mesh,fubar))
#else

    ASSERT(associated(disc%mesh, mmf_mesh(mmf)))
#endif

    if (heat_eqn) then
      htmodel => create_HT_model(disc%mesh, mmf, stat, errmsg)
      if (stat /= 0) return
    endif

    if (num_species > 0) then
      sdmodel => create_SD_model(disc%mesh, mmf, stat, errmsg)
      if (stat /= 0) return
    end if

    allocate(model)
    call HTSD_model_init (model, disc, htmodel, sdmodel)
    
    if (heat_eqn .and. void_material_index > 0) model%void_temp = void_temperature(void_material_index)

  end function create_HTSD_model

  function create_HT_model (mesh, mmf, stat, errmsg) result (model)

   use ER_driver

   type(dist_mesh), intent(in), target :: mesh
    type(mat_mf), intent(in), target :: mmf
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg
    type(HT_model), pointer :: model

    allocate(model)

    !! Initializes the VF_RAD_PROB components of MODEL.
    model%vf_rad_prob => create_vf_rad_prob(mesh, stat, errmsg)
    if (stat /= 0) return

    !! Defines the equation parameter components of MODEL.
    call define_system_parameters (mesh, mmf, model, stat, errmsg)
    if (stat /= 0) return

    !! Defines the boundary condition components of MODEL.
    call define_system_bc (mesh, model, stat, errmsg)
    if (stat /= 0) return

  contains

    subroutine define_system_parameters (mesh, mmf, model, stat, errmsg)

      use phase_property_table
      use material_utilities
      use material_mesh_function
      use property_mesh_function
      use ds_source_input, only: define_external_source
      use parallel_communication, only: global_any

      type(dist_mesh), intent(in), target :: mesh
      type(mat_mf), intent(in), target :: mmf
      type(HT_model), intent(inout) :: model
      integer, intent(out) :: stat
      character(len=*), intent(out) :: errmsg

      integer, pointer :: matid(:)

      !! Retrieve a list of all the material IDs that may be encountered.
      call mmf_get_all_matid (mmf, matid, drop_void=.true.)

      !! Enthalpy density.
      call required_property_check (matid, 'enthalpy density', stat, errmsg)
      if (stat /= 0) return
      call pmf_create (model%H_of_T, mmf, ppt_property_id('enthalpy density'), stat, errmsg)
      if (global_any(stat /= 0)) then
        stat = -1
        errmsg = 'unexpected error defining H_of_T: ' // trim(errmsg)
        return
      end if

      !! Thermal conductivity.
      call required_property_check (matid, 'conductivity', stat, errmsg)
      if (stat /= 0) return
      call pmf_create (model%conductivity, mmf, ppt_property_id('conductivity'), stat, errmsg)
      !call pmf_set_harmonic_average (model%conductivity)
      if (global_any(stat /= 0)) then
        stat = -1
        errmsg = 'unexpected error defining conductivity: ' // trim(errmsg)
        return
      end if

      !! External heat source.
      call define_external_source (mesh, 'temperature', model%source)

      deallocate(matid)

      stat = 0
      errmsg = ''

    end subroutine define_system_parameters

    subroutine define_system_bc (mesh, model, stat, errmsg)

      use ds_boundary_condition_input, only: get_boundary_data
      use ds_interface_condition_input, only: get_interface_data
      use index_partitioning, only: gather_boundary
      use bitfield_type, only: btest
      use physical_constants, only: stefan_boltzmann, absolute_zero
      use parallel_communication, only: global_any

      type(dist_mesh), intent(in), target :: mesh
      type(HT_model), intent(inout) :: model
      integer, intent(out) :: stat
      character(len=*), intent(out) :: errmsg

      integer :: j
      logical, allocatable :: mask(:), rmask(:)
      integer, allocatable :: setids(:)

      allocate(mask(mesh%nface))

      mask = .false.

      !! Define the internal HTC interface conditions.
      call get_interface_data (mesh, 'temperature', 'HTC', 1, model%ic_htc)
      mask(model%ic_htc%faces(1,:)) = .true.
      mask(model%ic_htc%faces(2,:)) = .true.
      call gather_boundary (mesh%face_ip, mask)
      
      !! Define the gap radiation interface conditions;
      !! may be superimposed with HTC conditions.
      call get_interface_data (mesh, 'temperature', 'radiation', 1, model%ic_rad)
      mask(model%ic_rad%faces(1,:)) = .true.
      mask(model%ic_rad%faces(2,:)) = .true.
      call gather_boundary (mesh%face_ip, mask)

      !! Define the external HTC boundary conditions.
      call get_boundary_data (mesh, 'temperature', 'HTC', 2, model%bc_htc)
      if (global_any(mask(model%bc_htc%faces))) then
        stat = -1
        errmsg = 'temperature HTC boundary condition overlaps with preceding conditions'
        return
      end if
      mask(model%bc_htc%faces) = .true. ! mark the HTC faces

      !! Tag all the enclosure radation (view factor) faces.  They were already
      !! verified to be boundary faces and non-overlapping with each other.
      allocate(rmask(mesh%nface))
      rmask = .false.
      if (associated(model%vf_rad_prob)) then
        do j = 1, size(model%vf_rad_prob)
          rmask(model%vf_rad_prob(j)%faces) = .true.
        end do
        call gather_boundary (mesh%face_ip, rmask)
      end if

      !! Define the (simple) radiation boundary conditions.
      call get_boundary_data (mesh, 'temperature', 'radiation', 2, model%bc_rad)
      if (global_any(rmask(model%bc_rad%faces))) then
        stat = -1
        errmsg = 'temperature radiation boundary condition overlaps with enclosure radiation'
        return
      end if
      rmask(model%bc_rad%faces) = .true. ! mark the radiation faces
      model%sbconst = stefan_boltzmann
      model%abszero = absolute_zero

      !! Radiation BC may be superimposed with the HTC flux conditions.
      mask = mask .or. rmask
      deallocate(rmask)

      !! Define the Dirichlet boundary conditions.
      call get_boundary_data (mesh, 'temperature', 'dirichlet', 1, model%bc_dir)
      if (global_any(mask(model%bc_dir%faces))) then
        stat = -1
        errmsg = 'temperature dirichlet boundary condition overlaps with preceding conditions'
        return
      end if
      mask(model%bc_dir%faces) = .true. ! mark the dirichlet faces

      !! Define the flux boundary conditions.
      call get_boundary_data (mesh, 'temperature', 'flux', 1, model%bc_flux)
      if (global_any(mask(model%bc_flux%faces))) then
        stat = -1
        errmsg = 'temperature flux boundary condition overlaps with preceding conditions'
        return
      end if
      mask(model%bc_flux%faces) = .true. ! mark the flux faces

      !! Finally verify that a condition has been applied to every boundary face.
      !! TODO: THIS DOESN'T WORK PROPERLY IN PARALLEL
      if (global_any(mask.neqv.btest(mesh%face_set_mask,0))) then
        mask = mask .neqv. btest(mesh%face_set_mask,0)
        call mesh%get_face_set_IDs (pack([(j,j=1,mesh%nface)], mask), setids)
        stat = -1
        write(errmsg,'(a,99(:,1x,i0))') &
          'incomplete temperature boundary/interface condition specification; ' // &
          'remaining boundary faces belong to face sets', setids
        deallocate(setids)
        return
      end if

      deallocate(mask)

      stat = 0
      errmsg = ''

    end subroutine define_system_bc

  end function create_HT_model

  function create_vf_rad_prob (mesh, stat, errmsg) result (vf_rad_prob)

    use ER_driver
    use ER_input
    use bitfield_type, only: btest
    use parallel_communication, only: global_any, global_all

    type(dist_mesh), intent(in) :: mesh
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg
    type(ERD_problem), pointer :: vf_rad_prob(:)

    integer :: n, j
    logical, allocatable :: mask(:)
    character(len=31), allocatable :: encl_name(:)

    stat = 0
    errmsg = ''

    !! Initialize the enclosure radiation (view factor) problems, if any.
    n = ERI_num_enclosures()
    if (n > 0) then
      allocate(vf_rad_prob(n), encl_name(n))
      call ERI_get_names (encl_name)
      do j = 1, n
        call ERD_problem_init (vf_rad_prob(j), mesh, encl_name(j))
        !! Verify that these enclosure faces are boundary faces.
        if (.not.global_all(btest(mesh%face_set_mask(vf_rad_prob(j)%faces),0))) then
          stat = -1
          errmsg = 'some enclosure faces are not boundary faces'
          exit
        else if (n > 1) then  ! multiple enclosures
          if (j == 1) then
            allocate(mask(mesh%nface))
            mask = .false.
          !! Verify that these faces don't overlap with other enclosure faces.
          else if (global_any(mask(vf_rad_prob(j)%faces))) then
            stat = -1
            errmsg = 'some enclosure faces belong to other enclosures'
            exit
          !! Tag the faces belonging to this enclosure.
          else
            mask(vf_rad_prob(j)%faces) = .true.
          end if
        end if
      end do
      if (allocated(mask)) deallocate(mask)
      deallocate(encl_name)
    else
      vf_rad_prob => null()
    end if

  end function create_vf_rad_prob

  function create_SD_model (mesh, mmf, stat, errmsg) result (model)

    use diffusion_solver_data, only: num_species, heat_eqn
    use phase_property_table
    use property_mesh_function
    use material_utilities
    use ds_source_input, only: define_external_source
    use ds_boundary_condition_input, only: get_boundary_data
    use bitfield_type, only: btest
    use index_partitioning
    use parallel_communication, only: global_any

    type(dist_mesh), intent(in), target :: mesh
    type(mat_mf), intent(in), target :: mmf
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg
    type(SD_model), pointer :: model(:)

    integer :: n, j
    logical :: mask(mesh%nface)
    character(len=31) :: property, variable
    integer, allocatable :: setids(:)
    integer, pointer :: matids(:)

    allocate(model(num_species))

    !! Define the equation parameter components of MODEL.
    call mmf_get_all_matid (mmf, matids, drop_void=.true.)
    do n = 1, num_species
      write(variable,'(a,i0)') 'concentration', n
      write(property,'(a,i0)') 'diffusivity', n
      call required_property_check (matids, property, stat, errmsg)
      if (stat /= 0) return
      call pmf_create (model(n)%diffusivity, mmf, ppt_property_id(property), stat, errmsg)
      if (global_any(stat /= 0)) then
        stat = -1
        errmsg = 'unexpected error defining diffusivity: ' // trim(errmsg)
        return
      end if
      call define_external_source (mesh, variable, model(n)%source)
      !! Define the optional Soret effect coefficients.
      if (heat_eqn) then
        write(property,'(a,i0)') 'soret', n
        call optional_property_check (matids, property, stat, errmsg)
        if (stat < 0) return
        if (stat == 0) then ! this property was defined
          allocate(model(n)%soret)
          call pmf_create (model(n)%soret, mmf, ppt_property_id(property), stat, errmsg)
          if (stat /= 0) return
        end if
      end if
    end do
    deallocate(matids)

    !! Define the boundary condition components of MODEL.
    do n = 1, num_species
      write(variable,'(a,i0)') 'concentration', n
      mask = .false.  ! used to tag faces where a BC has been applied
      !! Define the Dirichlet BC object for this concentration component.
      call get_boundary_data (mesh, variable, 'dirichlet', 1, model(n)%bc_dir)
      mask(model(n)%bc_dir%faces) = .true.  ! tag the Dirichlet BC faces
      !! Define the flux BC object for this concentration component.
      call get_boundary_data (mesh, variable, 'flux', 1, model(n)%bc_flux)
      !! Verify the flux BC faces don't overlap with preceding BC faces.
      if (global_any(mask(model(n)%bc_flux%faces))) then
        stat = -1
        errmsg = trim(variable) // ': flux BC overlaps with other BC!'
        return
      end if
      mask(model(n)%bc_flux%faces) = .true. ! tag the flux BC faces
      !! Finally verify that a BC has been applied to every boundary face.
      if (global_any(mask.neqv.btest(mesh%face_set_mask,0))) then
        mask = mask .neqv. btest(mesh%face_set_mask,0)
        call mesh%get_face_set_IDs (pack([(j,j=1,mesh%nface)], mask), setids)
        stat = -1
        !! The assumption here is that no bad faces are internal; any errors of
        !! that kind should have been caught when creating the bd_data objects.
        write(errmsg,'(a,99(:,1x,i0))') &
          trim(variable) // ': incomplete BC specification; ' // &
          'remaining boundary faces belong to face sets', setids
        deallocate(setids)
        return
      end if
    end do
    
    stat = 0
    errmsg = ''

  end function create_SD_model

end module HTSD_model_factory
