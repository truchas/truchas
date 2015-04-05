!!
!! MESH_BROKER
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Last revised 21 March 2006
!!

#include "f90_assert.fpp"

#define DEBUG_WRITE_MESH

module mesh_broker

  use kinds, only: r8
  use parallel_communication
  use dist_mesh_type

  implicit none
  private

  public :: dist_mesh

  public :: read_mesh_namelists, enable_mesh, init_mesh_broker, named_mesh_ptr
  public :: peek_truchas_mesh_namelists

  !!  This module maintains a collection of meshes using a recursive linked-
  !!  list data structure.

  integer, parameter :: MAX_NAME_LEN = 32, MAX_FILE_LEN = 256

  type :: mesh_desc
    character(len=MAX_NAME_LEN) :: name ! User-specified name for the mesh
    character(len=MAX_FILE_LEN) :: file ! mesh file path
    real(kind=r8) :: coord_scale_factor ! scaling factor for the mesh node coordinates
    integer, pointer :: gap_blocks(:) => null()
    integer, pointer :: iface_sset_id(:) => null()
    logical :: enabled = .false.
    logical :: em_mesh = .false. ! whether this mesh is for EM (HACK)
    type(dist_mesh), pointer :: mesh => null()
  end type mesh_desc

  type :: mesh_list
    type(mesh_list_item), pointer :: first => null()
  end type mesh_list

  type :: mesh_list_item
    type(mesh_desc) :: mesh
    type(mesh_list) :: rest
  end type mesh_list_item

  type(mesh_list), save :: meshes
  
  interface init_mesh_broker
    module procedure init_mesh_broker, init_mesh_broker_plist
  end interface

contains

  function named_mesh_ptr (name)
    use string_utilities, only: raise_case
    character(len=*), intent(in) :: name
    type(dist_mesh), pointer :: named_mesh_ptr
    type(mesh_list) :: l
    named_mesh_ptr => null()
    l = meshes
    do while (associated(l%first))
      if (l%first%mesh%enabled) then
        if (l%first%mesh%name == raise_case(name)) then
          named_mesh_ptr => l%first%mesh%mesh
          exit
        end if
      end if
      l = l%first%rest
    end do
  end function named_mesh_ptr

  subroutine enable_mesh (name, exists)
    use string_utilities, only: raise_case
    character(len=*), intent(in) :: name
    logical, intent(out) :: exists
    type(mesh_list) :: l
    exists = .false.
    l = meshes
    do while (associated(l%first))
      if (l%first%mesh%name == raise_case(name)) then
        l%first%mesh%enabled = .true.
        exists = .true.
        exit
      end if
      l = l%first%rest
    end do
  end subroutine enable_mesh

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! INIT_MESH_BROKER
 !!

  subroutine init_mesh_broker ()

    use mesh_importer
    use mesh_modification, only: convert_cells_to_links, create_internal_interfaces
    use dist_mesh_factory, only: new_dist_mesh
    use distributed_tet_mesh, only: create_dist_tet_mesh
#ifdef DEBUG_WRITE_MESH
    use dist_mesh_gmv
#endif

    integer :: stat
    character(len=127) :: errmsg
    type(mesh_list) :: l
    type(external_mesh) :: xmesh

    l = meshes
    do while (associated(l%first))
      if (l%first%mesh%enabled) then
        call write_info ('')
        call write_info ('Initializing distributed mesh ' // trim(l%first%mesh%name) // ' ...')
        call write_info ('  Reading ExodusII mesh file ' // trim(l%first%mesh%file))
        call import_exodus_mesh (l%first%mesh%file, xmesh)
        if (l%first%mesh%coord_scale_factor /= 1.0_r8) &
            xmesh%x = l%first%mesh%coord_scale_factor * xmesh%x
        if (is_IOP) call convert_cells_to_links (xmesh, l%first%mesh%gap_blocks, stat, errmsg)
        call broadcast (stat)
        if (stat /= 0) call halt ('INIT_MESH_BROKER: ' // trim(errmsg))
        if (is_IOP) call create_internal_interfaces (xmesh, l%first%mesh%iface_sset_id, stat, errmsg)
        call broadcast (stat)
        if (stat /= 0) call halt ('INIT_MESH_BROKER: ' // trim(errmsg))
        if (l%first%mesh%em_mesh) then
          select case (xmesh%mesh_type)
          case ('TET')
            allocate (l%first%mesh%mesh)
            call create_dist_tet_mesh (l%first%mesh%mesh, xmesh)
          case default
            INSIST( .false. )
          end select
        else
          select case (xmesh%mesh_type)
          case ('TET', 'HEX')
            l%first%mesh%mesh => new_dist_mesh(xmesh)
          case default
            INSIST( .false. )
          end select
        end if
        call l%first%mesh%mesh%write_profile
        call destroy (xmesh)
        call write_info ('  Distributed mesh ' // trim(l%first%mesh%name) // ' initialized.')
#ifdef DEBUG_WRITE_MESH
        call gmv_open (trim(l%first%mesh%name) // '.gmv')
        call gmv_write_dist_mesh (l%first%mesh%mesh)
        call gmv_close ()
#endif
      end if
      l = l%first%rest
    end do

  end subroutine init_mesh_broker

  subroutine write_info (message, skip)
    use truchas_logging_services
    character(len=*), intent(in) :: message
    logical, intent(in), optional :: skip
    if (present(skip)) then
      if (skip) call TLS_info ('')
    end if
    call TLS_info (message)
  end subroutine write_info
  
  subroutine halt (errmsg)
    use truchas_logging_services
    character(len=*), intent(in) :: errmsg
    call TLS_fatal (errmsg)
  end subroutine halt

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_MESH_NAMELISTS
 !!

  subroutine read_mesh_namelists (lun)

    use string_utilities, only: i_to_c, raise_case
    use input_utilities, only: seek_to_namelist

    integer, intent(in) :: lun

    logical :: found, file_exists
    integer :: stat, count
    character(len=MAX_NAME_LEN) :: name
    character(len=MAX_FILE_LEN) :: file
    real(kind=r8) :: coord_scale_factor
    type(mesh_list) :: old_meshes

    namelist /mesh/ name, file, coord_scale_factor

    !!! Check that a file is opened for reading on LUN? !!!

    if (is_IOP) rewind(lun)
    count = 0

    do  ! until all mesh namelists have been read or an error is found.

      !! Seek to the next instance of the MESH namelist.
      if (is_IOP) call seek_to_namelist (lun, 'MESH', found, iostat=stat)

      call broadcast (stat)
      if (stat /= 0) call abort ('error reading file connected to unit ' // i_to_c(lun))

      call broadcast (found)
      if (.not.found) return

      count = 1 + count

      !! Read the namelist, assigning default values first.
      if (is_IOP) then
        name = char(0)
        file = char(0)
        coord_scale_factor = 1.0_r8
        read(lun,nml=mesh,iostat=stat)
      end if

      call broadcast (stat)
      if (stat /= 0) call abort ('Error reading MESH namelist ' // i_to_c(count))

      !! Check that the values are valid.
      if (is_IOP) then
        !! Check NAME
        if (name == char(0)) call input_error ('NAME must be assigned a value')
        name = raise_case(name)
        if (mesh_name_exists(meshes, name)) call input_error ('a mesh with this name already exists: ' // trim(name))
        !! Check FILE
        if (file == char(0)) call input_error ('FILE must be assigned a value')
        inquire(file=file,exist=file_exists)
        if (.not.file_exists) call input_error ('no such FILE: ' // trim(file))
        !! Check COORD_SCALE_FACTOR
        if (coord_scale_factor <= 0.0_r8) call input_error ('COORD_SCALE_FACTOR must be positive')
      end if

      call broadcast (stat)
      if (stat /= 0) call abort ('Bad values for MESH namelist ' // i_to_c(count))

      !! Broadcast the mesh namelist variables.
      call broadcast (name)
      call broadcast (file)
      call broadcast (coord_scale_factor)

      !! Add this mesh to the existing collection of meshes (ALL).
      old_meshes = meshes
      allocate(meshes%first)
      meshes%first%rest = old_meshes
      meshes%first%mesh%name = name
      meshes%first%mesh%file = file
      meshes%first%mesh%coord_scale_factor = coord_scale_factor
      !! These are unused from this interface.
      allocate(meshes%first%mesh%gap_blocks(0))
      allocate(meshes%first%mesh%iface_sset_id(0))

    end do

  contains

    subroutine input_error (message)
      character(len=*), intent(in) :: message
      write(*,fmt='(2a)') 'FATAL: ', message
      stat = 1
    end subroutine input_error

    subroutine abort (message)
      character(len=*), intent(in) :: message
      write(*,fmt='(2a)') 'FATAL: READ_MESH_NAMELISTS: ', message
      call halt_parallel_communication ()
      stop
    end subroutine abort

  end subroutine read_mesh_namelists

  logical function mesh_name_exists (list, name)
    type(mesh_list),  intent(in) :: list
    character(len=*), intent(in) :: name
    type(mesh_list) :: l
    l = list
    mesh_name_exists = .true.
    do while (associated(l%first))
      if (l%first%mesh%name == name) return
      l = l%first%rest
    end do
    mesh_name_exists = .false.
  end function mesh_name_exists

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! PEEK_TRUCHAS_MESH_NAMELISTS
 !!

  subroutine peek_truchas_mesh_namelists ()

    use mesh_input_module, only: mesh_file, mesh_file_format, coordinate_scale_factor, &
                                 gap_element_blocks, interface_side_sets
    use altmesh_input, only: altmesh_exists, altmesh_file, altmesh_coordinate_scale_factor
    use string_utilities, only: raise_case

    integer :: n
    type(mesh_list) :: old_meshes

    !! MESH namelist: if it's a true Exodus/Genesis mesh file add it to the list.
    if (raise_case(mesh_file_format) == 'EXODUSII') then
      !! Add this mesh to the existing collection of meshes (ALL).
      old_meshes = meshes
      allocate(meshes%first)
      meshes%first%rest = old_meshes
      meshes%first%mesh%name = 'MAIN'
      meshes%first%mesh%file = mesh_file
      meshes%first%mesh%coord_scale_factor = coordinate_scale_factor
      n = count(gap_element_blocks > 0)
      allocate(meshes%first%mesh%gap_blocks(n))
      meshes%first%mesh%gap_blocks = pack(gap_element_blocks, mask=(gap_element_blocks > 0))
      n = count(interface_side_sets > 0)
      allocate(meshes%first%mesh%iface_sset_id(n))
      meshes%first%mesh%iface_sset_id = pack(interface_side_sets, mask=(interface_side_sets > 0))
    end if

    !! ALTMESH namelist.
    if (altmesh_exists) then
      !! Add this mesh to the existing collection of meshes (ALL).
      old_meshes = meshes
      allocate(meshes%first)
      meshes%first%rest = old_meshes
      meshes%first%mesh%name = 'ALT'
      meshes%first%mesh%file = altmesh_file
      meshes%first%mesh%em_mesh = .true.
      meshes%first%mesh%coord_scale_factor = altmesh_coordinate_scale_factor
      allocate(meshes%first%mesh%gap_blocks(0), meshes%first%mesh%iface_sset_id(0))
    end if

  end subroutine peek_truchas_mesh_namelists

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! PEEK_TRUCHAS_PARTITION
 !!

  subroutine peek_truchas_partition (cell_bsize, cell_perm, node_bsize, node_perm)

    use parameter_module, only: ncells, nnodes, ncells_tot, nnodes_tot
    use mesh_module, only: unpermute_mesh_vector, unpermute_vertex_vector

    integer, pointer :: cell_bsize(:), cell_perm(:), node_bsize(:), node_perm(:)

    call allocate_collated_array (cell_bsize, nPE)
    call allocate_collated_array (cell_perm, ncells_tot)

    call allocate_collated_array (node_bsize, nPE)
    call allocate_collated_array (node_perm, nnodes_tot)

    call collate (cell_bsize, ncells)
    call collate (cell_perm, unpermute_mesh_vector)

    call collate (node_bsize, nnodes)
    call collate (node_perm, unpermute_vertex_vector)

  end subroutine peek_truchas_partition

  !! NNC, Jun 2014.  A new initialization routine that takes a parameter list
  !! to specify the mesh files and parameters.  This is an alternative to
  !! using READ_MESH_NAMELISTS to read the data from an input file, and was
  !! added to simplify the development of unit tests requiring a mesh.  The
  !! expected syntax of the parameter list is a list of sublists.  The name
  !! of a sublist is taken as the name of the mesh, and the sublist has the
  !! following parameters:
  !!
  !!              'mesh-file' : string (required)
  !!     'coord-scale-factor' : real-scalar (optional; default 1.0)
  !!    'interface-side-sets' : integer-array (optional)

  subroutine init_mesh_broker_plist (plist)

    use parameter_list_type
    use string_utilities, only: raise_case

    type(parameter_list) :: plist

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: params
    type(mesh_list) :: old_meshes
    character(:), allocatable :: file
#ifdef INTEL_COMPILER_WORKAROUND
    character(:), allocatable :: name
#endif
    integer, allocatable :: side_sets(:)

    !! Generate the internal mesh specification list ...
    piter = parameter_list_iterator(plist, sublists_only=.true.)
    do while (.not.piter%at_end())
      old_meshes = meshes
      allocate(meshes%first)
      meshes%first%rest = old_meshes
#ifdef INTEL_COMPILER_WORKAROUND
      !Intel tracking ID: DPD200357444
      name = piter%name()
      meshes%first%mesh%name = raise_case(name)
#else
      meshes%first%mesh%name = raise_case(piter%name())
#endif
      params => piter%sublist()
      call params%get ('mesh-file', file)
      meshes%first%mesh%file = file
      call params%get ('coord-scale-factor', meshes%first%mesh%coord_scale_factor, default=1.0_r8)
      call params%get ('interface-side-sets', side_sets, default=[integer::])
      allocate(meshes%first%mesh%iface_sset_id(size(side_sets)))
      meshes%first%mesh%iface_sset_id = side_sets
      allocate(meshes%first%mesh%gap_blocks(0))
      meshes%first%mesh%enabled = .true.
      call piter%next
    end do

    !! and then call the original initialization procedure.
    call init_mesh_broker

  end subroutine init_mesh_broker_plist

end module mesh_broker
