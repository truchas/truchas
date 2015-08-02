!!
!! MESH_MANAGER
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! July 2015
!!

#include "f90_assert.fpp"

module mesh_manager

  use parameter_list_type
  use truchas_logging_services
  use simpl_mesh_type
  implicit none
  private
  
  public :: enable_mesh, named_mesh_ptr, unstr_mesh_ptr, simpl_mesh_ptr
  public :: peek_truchas_mesh_namelists, init_mesh_manager
  
  public :: simpl_mesh ! re-export; do not want this -- FIXME
  
  interface init_mesh_manager
    module procedure init_mesh_manager, init_mesh_manager_params
  end interface
  
  type(parameter_list), target, save :: meshes
  
  !! The value of the 'mesh' parameter will be of this private type.
  !! A PARAMETER_LIST stores (shallow) copies of parameter values, and thus
  !! cannot properly store a mesh object.  Instead we have it store a pointer
  !! to a mesh object by boxing the pointer in a TYPE(ANY_MESH) variable and
  !! and storing that variable.
  type :: any_mesh
    class(*), pointer :: mesh => null()
  end type any_mesh

contains

  !! NNC, Jun 2014.  A new initialization routine that takes a parameter list
  !! to specify the mesh files and parameters.  This is an alternative to
  !! using READ_MESH_NAMELISTS to read the data from an input file, and was
  !! added to simplify the development of unit tests requiring a mesh.  The
  !! expected syntax of the parameter list is a list of sublists.  The name
  !! of a sublist is taken as the name of the mesh, and the sublist has the
  !! following parameters:
  !!
  !!                 'mesh-file' : string (required)
  !!        'coord-scale-factor' : real-scalar (optional; default 1.0)
  !!    'interface-side-set-ids' : integer-array (optional)

  subroutine init_mesh_manager_params (params)

    use kinds, only: r8
    use string_utilities, only: raise_case

    type(parameter_list) :: params

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist1, plist2
    real(r8) :: coord_scale_factor
    character(:), allocatable :: mesh_file
#ifdef INTEL_COMPILER_WORKAROUND
    character(:), allocatable :: name
#endif
    integer, allocatable :: side_sets(:)

    !! Copy the input parameter list to the modules private parameter list.
    !! We want case-insensitive mesh names in the module, so we convert them
    !! to upper case when copying; later searching will use the upper-cased
    !! version of the specified name.
    piter = parameter_list_iterator(params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist1 => piter%sublist()
#ifdef INTEL_COMPILER_WORKAROUND
      !Intel tracking ID: DPD200357444
      name = piter%name()
      name = raise_case(name)
      plist2 => meshes%sublist(name)
#else
      plist2 => meshes%sublist(raise_case(piter%name()))
#endif
      call plist1%get ('mesh-file', mesh_file)
      call plist2%set ('mesh-file', mesh_file)
      call plist2%set ('enabled', .true.)
      
      if (plist1%is_parameter('coord-scale-factor')) then
        call plist1%get ('coord-scale-factor', coord_scale_factor)
        call plist2%set ('coord-scale-factor', coord_scale_factor)
      end if
      
      if (plist1%is_parameter('interface-side-set-ids')) then
        call plist1%get ('interface-side-set-ids', side_sets)
        call plist2%set ('interface-side-set-ids', side_sets)
      end if
      
      call piter%next
    end do

    !! Now call the original initialization procedure.
    call init_mesh_manager

  end subroutine init_mesh_manager_params

  subroutine init_mesh_manager
  
    use unstr_mesh_factory
    use simpl_mesh_factory
  
    integer :: stat
    logical :: enabled, em_mesh
    character(:), allocatable :: errmsg
    type(parameter_list), pointer :: plist
    type(parameter_list_iterator) :: piter
    type(unstr_mesh), pointer :: umesh
    type(simpl_mesh), pointer :: smesh

    piter = parameter_list_iterator(meshes)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      call plist%get ('enabled', enabled, default=.false.)
      if (enabled) then
        call TLS_info ('')
        call TLS_info ('Initializing mesh "' // piter%name() // '" ...')
        call plist%get ('em-mesh', em_mesh, default=.false.)
        if (em_mesh) then
          smesh => new_simpl_mesh(plist, stat, errmsg)
          if (stat /= 0) call TLS_fatal (errmsg)
          call plist%set ('mesh', any_mesh(smesh))
          call smesh%write_profile
        else
          umesh => new_unstr_mesh(plist, stat, errmsg)
          if (stat /= 0) call TLS_fatal (errmsg)
          call plist%set ('mesh', any_mesh(umesh))
          call umesh%write_profile
        end if
        call TLS_info ('  Mesh "' // trim(piter%name()) // '" initialized')
      end if
      call piter%next
    end do

  end subroutine init_mesh_manager

  subroutine peek_truchas_mesh_namelists

    use mesh_input_module, only: mesh_file, mesh_file_format, coordinate_scale_factor, &
                                 gap_element_blocks, interface_side_sets
    use altmesh_input, only: altmesh_exists, altmesh_file, altmesh_coordinate_scale_factor
    use string_utilities, only: raise_case
  
    integer, allocatable :: iarray(:)
    type(parameter_list), pointer :: plist

    !! The MESH namelist: the "MAIN" mesh
    if (raise_case(mesh_file_format) == 'EXODUSII') then
      plist => meshes%sublist('MAIN')
      call plist%set ('mesh', any_mesh())
      call plist%set ('mesh-file', trim(mesh_file))
      call plist%set ('coord-scale-factor', coordinate_scale_factor)
      iarray = pack(gap_element_blocks, mask=(gap_element_blocks > 0))
      if (size(iarray) > 0) call plist%set ('gap-element-block-ids', iarray)
      iarray = pack(interface_side_sets, mask=(interface_side_sets > 0))
      if (size(iarray) > 0) call plist%set ('interface-side-set-ids', iarray)
    end if
    
    !! ALTMESH namelist. the "ALT" mesh used by EM
    if (altmesh_exists) then
      plist => meshes%sublist('ALT')
      call plist%set ('mesh', any_mesh())
      call plist%set ('mesh-file', trim(altmesh_file))
      call plist%set ('coord-scale-factor', altmesh_coordinate_scale_factor)
      call plist%set ('em-mesh', .true.)
    end if
    
  end subroutine peek_truchas_mesh_namelists
  

  subroutine enable_mesh (name, exists)
    use string_utilities, only: raise_case
    character(*), intent(in) :: name
    logical, intent(out) :: exists
    type(parameter_list), pointer :: plist
    exists = meshes%is_sublist(raise_case(name))
    if (exists) then
      plist => meshes%sublist(raise_case(name))
      call plist%set ('enabled', .true.)
    end if
  end subroutine enable_mesh
    

  !! Returns a pointer to the CLASS(BASE_MESH) mesh object that corresponds to
  !! NAME. A null pointer is returned if NAME is not recognized or if the mesh
  !! that corresponds to NAME is not of class BASE_MESH.
  function named_mesh_ptr (name)
    use base_mesh_class
    character(*), intent(in) :: name
    class(base_mesh), pointer :: named_mesh_ptr
    class(*), pointer :: csptr
    csptr => mesh_ptr(name)
    named_mesh_ptr => null()
    if (.not.associated(csptr)) return
    select type (csptr)
    class is (base_mesh)
      named_mesh_ptr => csptr
    end select
  end function named_mesh_ptr

  !! Returns a pointer to the TYPE(UNSTR_MESH) mesh object that corresponds to
  !! NAME. A null pointer is returned if NAME is not recognized or if the mesh
  !! that corresponds to NAME is not of type UNSTR_MESH.
  function unstr_mesh_ptr (name)
    use unstr_mesh_type
    character(*), intent(in) :: name
    type(unstr_mesh), pointer :: unstr_mesh_ptr
    class(*), pointer :: csptr
    csptr => mesh_ptr(name)
    unstr_mesh_ptr => null()
    if (.not.associated(csptr)) return
    select type (csptr)
    class is (unstr_mesh)
      unstr_mesh_ptr => csptr
    end select
  end function unstr_mesh_ptr

  !! Returns a pointer to the TYPE(SIMPL_MESH) mesh object that corresponds to
  !! NAME. A null pointer is returned if NAME is not recognized or if the mesh
  !! that corresponds to NAME is not of type SIMPL_MESH.
  function simpl_mesh_ptr (name)
    use simpl_mesh_type
    character(*), intent(in) :: name
    type(simpl_mesh), pointer :: simpl_mesh_ptr
    class(*), pointer :: csptr
    csptr => mesh_ptr(name)
    simpl_mesh_ptr => null()
    if (.not.associated(csptr)) return
    select type (csptr)
    class is (simpl_mesh)
      simpl_mesh_ptr => csptr
    end select
  end function simpl_mesh_ptr

  !! Auxiliary function returns a CLASS(*) pointer to the arbitrary mesh object
  !! that corresponds to NAME.  NAME is case-insensitive.  A null pointer is
  !! returned if NAME is not recognized.
  function mesh_ptr (name)
    use string_utilities, only: raise_case
    character(*), intent(in) :: name
    class(*), pointer :: mesh_ptr
    type(parameter_list), pointer :: plist
    class(*), allocatable :: cs
    mesh_ptr => null()
    if (meshes%is_sublist(raise_case(name))) then
      plist => meshes%sublist(raise_case(name))
      call plist%get_any ('mesh', cs)
      select type (cs)  ! unbox the mesh pointer
      type is (any_mesh)
        mesh_ptr => cs%mesh
      class default
        INSIST(.false.)
      end select
    end if
  end function mesh_ptr

end module mesh_manager
