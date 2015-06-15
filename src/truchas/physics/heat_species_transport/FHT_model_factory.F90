!!
!! FHT_MODEL_FACTORY
!!
!! This module provides a method for initializing an instance of an FHT_MODEL
!! variable.  It depends essentially on Truchas-specific input data.
!!
!! PROGRAMMING INTERFACE
!!
!!  MODEL => CREATE_FHT_MODEL (MESH, MMF, STAT, ERRMSG) allocates and fully
!!    defines the components of the FHT_MODEL-type variable MODEL.  MESH is
!!    a pointer to the distributed mesh and MMF is a pointer to a corresponding
!!    material mesh function.  The integer STAT returns a non-zero value if an
!!    error is encountered, and an explanatory message is returned in the
!!    character string ERRMSG.
!!
!! IMPLEMENTATION NOTES
!!

#include "f90_assert.fpp"

module FHT_model_factory

  use kinds, only: r8
  use FHT_model_type
  use base_mesh_class
  use mfd_disc_type
  use material_mesh_function
  use rad_problem_type
  implicit none
  private
  
  public :: create_FHT_model

contains

  function create_FHT_model (disc, mmf, stat, errmsg) result (model)
  
    type(mfd_disc), intent(in), target :: disc
    type(mat_mf), intent(in), target :: mmf
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg
    type(FHT_model), pointer :: model
    
    class(base_mesh), pointer :: mesh

    !select type (disc_mesh => disc%mesh)
    !type is (dist_mesh)
    !  mesh => disc_mesh
    !class default
    !  INSIST(.false.)
    !end select
    mesh => disc%mesh
    
    allocate(model)
    
    !! Initializes the VF_RAD_PROB components of MODEL.
    call vf_rad_init (mesh, model, stat, errmsg)
    if (stat /= 0) return
    
    !! Defines the heat equation parameter components of MODEL.
    call define_system_parameters (mesh, mmf, model, stat, errmsg)
    if (stat /= 0) return
    
    !! Defines the boundary condition components of MODEL.
    call define_system_bc (mesh, model, stat, errmsg)
    if (stat /= 0) return
    
    !! Perform the final initialization of MODEL.
    call FHT_model_init (model, disc)
   
  end function create_FHT_model
  

  subroutine vf_rad_init (mesh, model, stat, errmsg)
  
    use ER_input
    use bitfield_type, only: btest
    use parallel_communication, only: global_any, global_all
  
    class(base_mesh), intent(in) :: mesh
    type(FHT_model), intent(inout) :: model
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg
    
    integer :: n, j
    logical, allocatable :: mask(:)
    character(len=31), allocatable :: encl_name(:)
    
    stat = 0
    errmsg = ''
    
    !! Initialize the enclosure radiation (view factor) problems, if any.
    n = ERI_num_enclosures()
    if (n > 0) then
      allocate(model%vf_rad_prob(n), encl_name(n))
      call ERI_get_names (encl_name)
      do j = 1, n
        call model%vf_rad_prob(j)%init (mesh, encl_name(j))
        !! Verify that these enclosure faces are boundary faces.
        if (.not.global_all(btest(mesh%face_set_mask(model%vf_rad_prob(j)%faces),0))) then
          stat = -1
          errmsg = 'some enclosure faces are not boundary faces'
          exit
        else if (n > 1) then  ! multiple enclosures
          if (j == 1) then
            allocate(mask(mesh%nface))
            mask = .false.
          !! Verify that these faces don't overlap with other enclosure faces.
          else if (global_any(mask(model%vf_rad_prob(j)%faces))) then
            stat = -1
            errmsg = 'some enclosure faces belong to other enclosures'
            exit
          !! Tag the faces belonging to this enclosure.
          else
            mask(model%vf_rad_prob(j)%faces) = .true.
          end if
        end if
      end do
      if (allocated(mask)) deallocate(mask)
      deallocate(encl_name)
    end if

  end subroutine vf_rad_init


  subroutine define_system_parameters (mesh, mmf, model, stat, errmsg)

    use phase_property_table
    use material_utilities
    use material_mesh_function
    use property_mesh_function
    use ds_source_input, only: define_external_source
    use parallel_communication, only: global_any

    class(base_mesh), intent(in), target :: mesh
    type(mat_mf), intent(in), target :: mmf
    type(FHT_model), intent(inout) :: model
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer, pointer :: matid(:)

    !! Retrieve a list of all the material IDs that may be encountered.
    call mmf_get_all_matid (mmf, matid, drop_void=.true.)

    !! Enthalpy density.
    allocate(model%H_of_T)
    call required_property_check (matid, 'enthalpy density', stat, errmsg)
    if (stat /= 0) return
    call pmf_create (model%H_of_T, mmf, ppt_property_id('enthalpy density'), stat, errmsg)
    if (global_any(stat /= 0)) then
      stat = -1
      errmsg = 'unexpected error defining H_of_T: ' // trim(errmsg)
      return
    end if

    !! Thermal conductivity.
    allocate(model%conductivity)
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
    allocate(model%q)
    call define_external_source (mesh, 'temperature', model%q)

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

    class(base_mesh), intent(in), target :: mesh
    type(FHT_model), intent(inout) :: model
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: j
    logical, allocatable :: mask(:), rmask(:)
    integer, allocatable :: setids(:)

    allocate(mask(mesh%nface))

    mask = .false.

    !! Define the internal HTC interface conditions.
    allocate(model%ic_htc)
    call get_interface_data (mesh, 'temperature', 'HTC', 1, model%ic_htc)
    mask(model%ic_htc%faces(1,:)) = .true.
    mask(model%ic_htc%faces(2,:)) = .true.
    call gather_boundary (mesh%face_ip, mask)

    !! Define the gap radiation interface conditions;
    !! may be superimposed with HTC conditions.
    allocate(model%ic_rad)
    call get_interface_data (mesh, 'temperature', 'radiation', 1, model%ic_rad)
    mask(model%ic_rad%faces(1,:)) = .true.
    mask(model%ic_rad%faces(2,:)) = .true.
    call gather_boundary (mesh%face_ip, mask)

    !! Define the external HTC boundary conditions.
    allocate(model%bc_htc)
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
    allocate(model%bc_rad)
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
    allocate(model%bc_dir)
    call get_boundary_data (mesh, 'temperature', 'dirichlet', 1, model%bc_dir)
    if (global_any(mask(model%bc_dir%faces))) then
      stat = -1
      errmsg = 'temperature dirichlet boundary condition overlaps with preceding conditions'
      return
    end if
    mask(model%bc_dir%faces) = .true. ! mark the dirichlet faces

    !! Define the flux boundary conditions.
    allocate(model%bc_flux)
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
      call mesh%get_face_set_IDs (pack((/(j,j=1,mesh%nface)/), mask), setids)
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

end module FHT_model_factory
