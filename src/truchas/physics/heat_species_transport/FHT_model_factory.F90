!!
!! FHT_MODEL_FACTORY
!!
!! This module provides a method for initializing an instance of an FHT_MODEL
!! variable.  It depends essentially on Truchas-specific input data.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
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

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use FHT_model_type
  use unstr_mesh_type
  use mfd_disc_type
  use matl_mesh_func_type
  use rad_problem_type
  use thermal_bc_factory_class
  use thermal_source_factory_type
  implicit none
  private

  public :: create_FHT_model

contains

  function create_FHT_model (tinit, disc, mmf, tbc_fac, tsrc_fac, stat, errmsg) result (model)

    use diffusion_solver_data, only: void_temperature

    real(r8), intent(in) :: tinit
    type(mfd_disc), intent(in), target :: disc
    type(matl_mesh_func), intent(in), target :: mmf
    class(thermal_bc_factory), intent(inout) :: tbc_fac
    type(thermal_source_factory), intent(inout) :: tsrc_fac
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg
    character(:), allocatable :: errmsg2
    type(FHT_model), pointer :: model

    type(unstr_mesh), pointer :: mesh

    mesh => disc%mesh

    allocate(model)

    !! Initializes the VF_RAD_PROB components of MODEL.
    call vf_rad_init (tinit, mesh, model, stat, errmsg)
    if (stat /= 0) return

    !! Defines the heat equation parameter components of MODEL.
    call define_system_parameters (mesh, mmf, tsrc_fac, model, stat, errmsg)
    if (stat /= 0) return

    !! Defines the boundary condition components of MODEL.
    call define_system_bc (mesh, tbc_fac, model, stat, errmsg2)
    if (stat /= 0) then
      errmsg = errmsg2
      return
    end if

    !! Perform the final initialization of MODEL.
    call FHT_model_init (model, disc)

    model%void_temp = void_temperature

  end function create_FHT_model


  subroutine vf_rad_init (tinit, mesh, model, stat, errmsg)

    use enclosure_radiation_namelist, only: params
    use bitfield_type, only: btest
    use parallel_communication, only: global_any, global_all
    use parameter_list_type

    real(r8), intent(in) :: tinit
    type(unstr_mesh), intent(in) :: mesh
    type(FHT_model), intent(inout) :: model
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: n, j
    logical, allocatable :: mask(:)
    character(len=31), allocatable :: encl_name(:)
    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist

    stat = 0
    errmsg = ''

    !! Initialize the enclosure radiation (view factor) problems, if any.
    piter = parameter_list_iterator(params, sublists_only=.true.)
    n = piter%count()
    if (n > 0) then
      allocate(model%vf_rad_prob(n), encl_name(n))
      do j = 1, n
        plist => piter%sublist()
        call model%vf_rad_prob(j)%init (mesh, piter%name(), plist, tinit)
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
        call piter%next
      end do
      if (allocated(mask)) deallocate(mask)
    end if

  end subroutine vf_rad_init


  subroutine define_system_parameters (mesh, mmf, src_fac, model, stat, errmsg)

    use matl_mesh_func_type
    use ds_source_input, only: define_external_source
    use parallel_communication, only: global_any
    use material_model_driver, only: matl_model
    use material_utilities

    type(unstr_mesh), intent(in), target :: mesh
    type(matl_mesh_func), intent(in), target :: mmf
    type(thermal_source_factory), intent(inout) :: src_fac
    type(FHT_model), intent(inout) :: model
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    !integer, allocatable :: matid(:)
    character(:), allocatable :: errmsg2

    !! Retrieve a list of all the material IDs that may be encountered.
    !call mmf%get_all_matl(matid, drop_void=.true.)

    !! Enthalpy density.
    allocate(model%H_of_T)
    call required_property_check(matl_model, 'enthalpy', stat, errmsg2)
    if (stat /= 0) then
      errmsg = errmsg2
      return
    end if
    call model%H_of_T%init(mmf, 'enthalpy', stat, errmsg)
    if (global_any(stat /= 0)) then
      stat = -1
      errmsg = 'unexpected error defining H_of_T: ' // trim(errmsg)
      return
    end if

    !! Thermal conductivity.
    allocate(model%conductivity)
    call required_property_check(matl_model, 'conductivity', stat, errmsg2)
    if (stat /= 0) then
      errmsg = errmsg2
      return
    end if
    call model%conductivity%init(mmf, 'conductivity', stat, errmsg)
    !call pmf_set_harmonic_average (model%conductivity)
    if (global_any(stat /= 0)) then
      stat = -1
      errmsg = 'unexpected error defining conductivity: ' // trim(errmsg)
      return
    end if

    !! External heat source.
    allocate(model%q)
    call define_external_source (mesh, 'temperature', model%q)

    !! Additional heat sources
    call src_fac%alloc_source_func1(model%src, stat, errmsg2)
    if (stat /= 0) then
      errmsg = errmsg2
      return
    end if

    stat = 0
    errmsg = ''

  end subroutine define_system_parameters

  subroutine define_system_bc(mesh, bc_fac, model, stat, errmsg)

    use bitfield_type
    use parallel_communication, only: global_any, global_count
    use string_utilities, only: i_to_c
    use f08_intrinsics, only: findloc

    type(unstr_mesh), intent(in), target :: mesh
    class(thermal_bc_factory), intent(inout) :: bc_fac
    type(FHT_model), intent(inout) :: model
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, n
    logical, allocatable :: mask(:), rmask(:), fmask(:)
    integer, allocatable :: setids(:)
    character(160) :: string1, string2
    type(bitfield) :: bitmask

    allocate(mask(mesh%nface))

    mask = .false.

    !! Define the internal HTC interface conditions.
    call bc_fac%alloc_htc_ic(model%ic_htc, stat, errmsg)
    if (stat /= 0) return
    if (allocated(model%ic_htc)) then
      mask(model%ic_htc%index(1,:)) = .true.
      mask(model%ic_htc%index(2,:)) = .true.
      call mesh%face_ip%gather_offp(mask)
    end if

    !! Define the gap radiation interface conditions;
    !! may be superimposed with HTC conditions.
    call bc_fac%alloc_rad_ic(model%ic_rad, stat, errmsg)
    if (stat /= 0) return
    if (allocated(model%ic_rad)) then
      mask(model%ic_rad%index(1,:)) = .true.
      mask(model%ic_rad%index(2,:)) = .true.
      call mesh%face_ip%gather_offp(mask)
    end if

    !! Flux type boundary conditions can be superimposed.
    allocate(fmask(mesh%nface))
    fmask = .false.

    !! Define the simple flux boundary conditions.
    call bc_fac%alloc_flux_bc(model%bc_flux, stat, errmsg)
    if (stat /= 0) return
      if (allocated(model%bc_flux)) then
      if (global_any(mask(model%bc_flux%index))) then
        stat = -1
        errmsg = 'temperature flux boundary condition overlaps with interface conditions'
        return
      end if
      fmask(model%bc_flux%index) = .true. ! mark the simple flux faces
    end if

    !! Define the oriented flux boundary conditions.
    call bc_fac%alloc_vflux_bc(model%bc_vflux, stat, errmsg)
    if (stat /= 0) return
      if (allocated(model%bc_vflux)) then
      if (global_any(mask(model%bc_vflux%index))) then
        stat = -1
        errmsg = 'temperature flux boundary condition overlaps with interface conditions'
        return
      end if
      fmask(model%bc_vflux%index) = .true. ! mark the oriented flux faces
    end if

    !! Define the external HTC boundary conditions.
    call bc_fac%alloc_htc_bc(model%bc_htc, stat, errmsg)
    if (stat /= 0) return
    if (allocated(model%bc_htc)) then
      if (global_any(mask(model%bc_htc%index))) then
        stat = -1
        errmsg = 'temperature HTC boundary condition overlaps with interface conditions'
        return
      end if
      fmask(model%bc_htc%index) = .true. ! mark the HTC faces
    end if

    !! Tag all the enclosure radation (view factor) faces.  They were already
    !! verified to be boundary faces and non-overlapping with each other.
    allocate(rmask(mesh%nface))
    rmask = .false.
    if (associated(model%vf_rad_prob)) then
      do j = 1, size(model%vf_rad_prob)
        rmask(model%vf_rad_prob(j)%faces) = .true.
        !fmask(model%vf_rad_prob(j)%faces) = .true.
      end do
      call mesh%face_ip%gather_offp(rmask)
      !call mesh%face_ip%gather_offp(fmask)
    end if

    !! Define the (simple) radiation boundary conditions.
    call bc_fac%alloc_rad_bc(model%bc_rad, stat, errmsg)
    if (stat /= 0) return
    if (allocated(model%bc_rad)) then
      if (global_any(rmask(model%bc_rad%index))) then
        stat = -1
        errmsg = 'temperature radiation boundary condition overlaps with enclosure radiation'
        return
      end if
    end if
    if (allocated(model%bc_rad)) then
      if (global_any(mask(model%bc_rad%index))) then
        stat = -1
        errmsg = 'temperature radiation boundary condition overlaps with interface conditions'
        return
      end if
      fmask(model%bc_rad%index) = .true. ! mark the radiation faces
    end if

    !! Merge the flux mask into the main mask.
    mask = mask .or. fmask
    deallocate(fmask)

    !! Define the Dirichlet boundary conditions.
    call bc_fac%alloc_dir_bc(model%bc_dir, stat, errmsg)
    if (stat /= 0) return
    if (allocated(model%bc_dir)) then
      if (global_any(mask(model%bc_dir%index))) then
        stat = -1
        errmsg = 'temperature dirichlet boundary condition overlaps with other conditions'
        return
      end if
      mask(model%bc_dir%index) = .true. ! mark the dirichlet faces
    end if

    !! Update enclosure radiation face mask for temperature Dirichlet faces.
    if (associated(model%vf_rad_prob) .and. allocated(model%bc_dir)) then
      do n = 1, size(model%vf_rad_prob)
        associate (faces => model%vf_rad_prob(n)%faces, fmask => model%vf_rad_prob(n)%fmask)
        do j = 1, size(faces)
          if (findloc(model%bc_dir%index, faces(j), dim=1) > 0) fmask(j) = .false.
        end do
        end associate
      end do
    end if
    mask = mask .or. rmask

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
