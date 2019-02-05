!!
!! FHT_MODEL_FACTORY
!!
!! This module provides a method for initializing an instance of an FHT_MODEL
!! variable.  It depends essentially on Truchas-specific input data.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  use unstr_mesh_type
  use mfd_disc_type
  use matl_mesh_func_type
  use rad_problem_type
  implicit none
  private
  
  public :: create_FHT_model

contains

  function create_FHT_model (disc, mmf, stat, errmsg) result (model)
  
    type(mfd_disc), intent(in), target :: disc
    type(matl_mesh_func), intent(in), target :: mmf
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg
    type(FHT_model), pointer :: model
    
    type(unstr_mesh), pointer :: mesh

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
  
    type(unstr_mesh), intent(in) :: mesh
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
    use matl_mesh_func_type
    use property_mesh_function
    use ds_source_input, only: define_external_source
    use parallel_communication, only: global_any

    type(unstr_mesh), intent(in), target :: mesh
    type(matl_mesh_func), intent(in), target :: mmf
    type(FHT_model), intent(inout) :: model
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer, allocatable :: matid(:)

    !! Retrieve a list of all the material IDs that may be encountered.
    call mmf%get_all_matl(matid, drop_void=.true.)

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

    stat = 0
    errmsg = ''

  end subroutine define_system_parameters

  subroutine define_system_bc (mesh, model, stat, errmsg)

    use ds_boundary_condition_input, only: get_boundary_data
    use ds_interface_condition_input, only: get_interface_data
    use index_partitioning, only: gather_boundary
    use bitfield_type
    use physical_constants, only: stefan_boltzmann, absolute_zero
    use parallel_communication, only: global_any, global_count
    use string_utilities, only: i_to_c

    type(unstr_mesh), intent(in), target :: mesh
    type(FHT_model), intent(inout) :: model
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: j, n
    logical, allocatable :: mask(:), rmask(:), fmask(:)
    integer, allocatable :: setids(:)
    character(len(errmsg)) :: string1, string2
    type(bitfield) :: bitmask

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

    !! Flux type boundary conditions can be superimposed.
    allocate(fmask(mesh%nface))
    fmask = .false.

    !! Define the simple flux boundary conditions.
    allocate(model%bc_flux)
    call get_boundary_data (mesh, 'temperature', 'flux', 1, model%bc_flux)
    if (global_any(mask(model%bc_flux%faces))) then
      stat = -1
      errmsg = 'temperature flux boundary condition overlaps with interface conditions'
      return
    end if
    fmask(model%bc_flux%faces) = .true. ! mark the simple flux faces

    !! Define the external HTC boundary conditions.
    allocate(model%bc_htc)
    call get_boundary_data (mesh, 'temperature', 'HTC', 2, model%bc_htc)
    if (global_any(mask(model%bc_htc%faces))) then
      stat = -1
      errmsg = 'temperature HTC boundary condition overlaps with interface conditions'
      return
    end if
    fmask(model%bc_htc%faces) = .true. ! mark the HTC faces

    !! Tag all the enclosure radation (view factor) faces.  They were already
    !! verified to be boundary faces and non-overlapping with each other.
    allocate(rmask(mesh%nface))
    rmask = .false.
    if (associated(model%vf_rad_prob)) then
      do j = 1, size(model%vf_rad_prob)
        rmask(model%vf_rad_prob(j)%faces) = .true.
        fmask(model%vf_rad_prob(j)%faces) = .true.
      end do
      call gather_boundary (mesh%face_ip, rmask)
      call gather_boundary (mesh%face_ip, fmask)
    end if

    !! Define the (simple) radiation boundary conditions.
    allocate(model%bc_rad)
    call get_boundary_data (mesh, 'temperature', 'radiation', 2, model%bc_rad)
    if (global_any(rmask(model%bc_rad%faces))) then
      stat = -1
      errmsg = 'temperature radiation boundary condition overlaps with enclosure radiation'
      return
    end if
    deallocate(rmask)
    if (global_any(mask(model%bc_rad%faces))) then
      stat = -1
      errmsg = 'temperature radiation boundary condition overlaps with interface conditions'
      return
    end if
    fmask(model%bc_rad%faces) = .true. ! mark the radiation faces
    model%sbconst = stefan_boltzmann
    model%abszero = absolute_zero

    !! Merge the flux mask into the main mask.
    mask = mask .or. fmask
    deallocate(fmask)

    !! Define the Dirichlet boundary conditions.
    allocate(model%bc_dir)
    call get_boundary_data (mesh, 'temperature', 'dirichlet', 1, model%bc_dir)
    if (global_any(mask(model%bc_dir%faces))) then
      stat = -1
      errmsg = 'temperature dirichlet boundary condition overlaps with other conditions'
      return
    end if
    mask(model%bc_dir%faces) = .true. ! mark the dirichlet faces

    !! Finally verify that a condition has been applied to every boundary face.
    mask = mask .neqv. btest(mesh%face_set_mask,0)
    if (global_any(mask)) then
      call mesh%get_face_set_ids(pack([(j,j=1,mesh%nface)], mask), setids)
      if (size(setids) == 0) then
        string1 = ' (none)'
      else
        write(string1,'(i0,*(:,", ",i0))') setids
      end if
      call mesh%get_link_set_ids(mask, setids)
      if (size(setids) == 0) then
        string2 = ' (none)'
      else
        write(string2,'(i0,*(:,", ",i0))') setids
      end if
      errmsg = 'incomplete temperature boundary/interface specification;' // &
          ' remaining boundary faces belong to face sets ' // trim(string1) // &
          '; and interface link sets ' // trim(string2)
      bitmask = ibset(ZERO_BITFIELD, 0)
      mask = mask .and. (mesh%face_set_mask == bitmask)
      mask(mesh%lface(1,:)) = .false.
      mask(mesh%lface(2,:)) = .false.
      n = global_count(mask(:mesh%nface_onP))
      if (n > 0) errmsg = trim(errmsg) // '; ' // i_to_c(n) // ' faces belong to neither'
      stat = -1
      return
    end if

    stat = 0
    errmsg = ''

  end subroutine define_system_bc

end module FHT_model_factory
