!!
!! MESH_IMPORTER
!!
!! Neil N. Carlson <nnc@newmexico.com> 28 Sep 2004
!! Neil N. Carlson <nnc@lanl.gov> 21 Feb 2006
!!
!! This module provides a procedure for reading a mesh from an Exodus II
!! format mesh file and a data structure for encapsulating the mesh data.
!!
!! Instances of the mesh data structure created by this module are intended
!! to be short-lived, serving only as a conduit from the outside world to
!! some initialization procedure that uses the data to create the actual
!! distributed mesh data structures used by the application code.
!!
!! PROGRAMMING INTERFACE
!!
!! The module provides the derived data type EXTERNAL_MESH, which encapsulates
!! the mesh data, and the following subroutines.
!!
!!  CALL IMPORT_EXODUS_MESH (PATH, MESH, STAT)
!!    CHARACTER(LEN=*),    INTENT(IN)  :: PATH
!!    TYPE(EXTERNAL_MESH), INTENT(OUT) :: MESH
!!    INTEGER, OPTIONAL,   INTENT(IN)  :: STAT
!!
!! The Exodus II format mesh is read from PATH and returned in MESH.
!! Currently, only non-hybrid meshes composed of either 4-node tets or 8-node
!! hexes are handled.  When present, STAT returns a positive value if an errror
!! occurred while reading the file, -1 if the mesh type can't be handled, and
!! otherwise 0 if the import was successful.
!!
!!  CALL DESTROY (MESH)
!!    TYPE(EXTERNAL_MESH), INTENT(INOUT) :: MESH
!!
!! All array components of MESH are deallocated, and MESH is returned to its
!! default initialization state.
!!

#include "f90_assert.fpp"

module mesh_importer

  use kinds
  use exodus_mesh_type
  use exodus_mesh_io, only: read_exodus_mesh
  use parallel_communication
  implicit none
  private

  public :: destroy, import_exodus_mesh
  public :: dump_external_mesh ! for debugging use
  public :: collect_main_mesh
  public :: side_set_node_list

  public :: side_set, node_set  ! export these exodus types

  type, public :: external_mesh
    integer :: ncell=0, nnode=0                     ! number of cells and nodes
    integer, allocatable :: cnode(:,:)              ! cell node array (nodes-of-cell)
    real(kind=r8), allocatable :: x(:,:)            ! node positions

    !! Attributes
    character(len=8) :: mesh_type = 'UNKNOWN'       ! tet, hex, etc.

    !! Exodus element block data
    integer :: nblock=0                             ! number of element blocks
    integer, allocatable :: block_id(:)             ! list of block IDs
    integer, allocatable :: cell_block(:)           ! block index for each cell

    !! Exodus side set data; just use the exodus-provided type.
    type(side_set), allocatable :: sset(:)

    !! Exodus node set data; just use the exodus-provided type.
    type(node_set), allocatable :: nset(:)
    
    !! Link data
    integer :: nlink=0, nlblock=0                   ! number of link blocks
    integer, pointer :: lnode(:,:) => null()        ! link node array (nodes-of-link)
    integer, pointer :: link_block(:) => null()     ! block ID for each link
    integer, pointer :: link_block_id(:) => null() ! list of unique block IDs
  end type external_mesh

!! NOTE: cell_block(:) are block IDs (by cells) and block_ids the list of
!! unique IDs.  cell_block(:) ought to index into the block_ids array instead,
!! with the index serving as the internal block number.  Same considerations
!! with the link blocks.  This needs to be fixed and all the using code (there
!! isn't much right now) modified accordingly.  NNC 6/17/2009

  interface destroy
    module procedure destroy_mesh
  end interface

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! IMPORT_EXODUS_MESH
 !!
 !! This procedure reads the mesh from the Exodus II file PATH, and
 !! initializes the data structure MESH with the mesh it describes.
 !!
 !! Only non-hybrid 3-D meshes composed of 4-node tets or 8-node hexes are
 !! currently handled.  If present, the optional argument STAT returns 0
 !! if the import was successful, a positive value if an error was encountered
 !! while reading the file, and -1 if the type of the mesh read is not
 !! recognized.  A message is written and execution halted if an error
 !! occurs and STAT is not present.
 !!

  subroutine import_exodus_mesh (path, mesh)
#ifdef NAG_COMPILER

    use,intrinsic :: f90_unix, only: exit
#endif

    character(len=*),    intent(in)  :: path
    type(external_mesh), intent(out) :: mesh

    integer :: stat, dimen, nvert, ssize
    character(:), allocatable :: errmsg

    !! Import the mesh on the IO processor.
    if (is_IOP) then
      call import_exodus_mesh_iop (path, mesh, stat, errmsg)
      if (stat > 0) then
        write(0,'(a)') 'IMPORT_EXODUS_MESH: Error reading file: ' // errmsg
      else if (stat < 0) then
        write(0,'(a)') 'IMPORT_EXODUS_MESH: Error: mesh type not recognized; must be hex or tet.'
      end if
    end if

    call broadcast (stat)
    if (stat /= 0) then
      call halt_parallel_communication ()
      call exit (1)
    end if

    !! Create 0-sized companion meshes on all the other processors.
    call broadcast (mesh%mesh_type)
    
    select case (mesh%mesh_type)
    case ('TET')
      dimen = 3
      nvert = 4
      ssize = 3
    case ('HEX')
      dimen = 3
      nvert = 8
      ssize = 4
    case default
      INSIST( .false. ) ! we shouldn't be here
    end select
    
    allocate(mesh%lnode(2*ssize,0), mesh%link_block(0), mesh%link_block_id(0))

    if (.not.is_IOP) then
      allocate(mesh%x(dimen,0))
      allocate(mesh%cnode(nvert,0))
      allocate(mesh%block_id(0), mesh%cell_block(0))
      allocate(mesh%sset(0), mesh%nset(0))
    end if

  end subroutine import_exodus_mesh

  subroutine import_exodus_mesh_iop (path, mesh, stat, errmsg)

    character(len=*),    intent(in)  :: path
    type(external_mesh), intent(out) :: mesh
    integer,             intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, n, offset, nvert
    type(exodus_mesh) :: exo_mesh

    !! Read the Exodus II mesh file.
    call read_exodus_mesh (path, exo_mesh, stat, errmsg)
    if (stat /= 0) then
      stat = 1
      return ! status value is positive
    end if

    !INSIST( defined(exo_mesh) )

    call identify_mesh_type (exo_mesh, mesh%mesh_type)
    select case (mesh%mesh_type)
    case ('TET')
      nvert = 4
    case ('HEX')
      nvert = 8
    case default
      stat = -1
      return
    end select

    mesh%nnode  = exo_mesh%num_node
    mesh%ncell  = exo_mesh%num_elem
    mesh%nblock = exo_mesh%num_eblk

    !! Take the node positions.
    call move_alloc (exo_mesh%coord, mesh%x)

    !! Allocate the remaining array components.
    allocate(mesh%cnode(nvert,mesh%ncell), mesh%block_id(mesh%nblock), mesh%cell_block(mesh%ncell))

    mesh%block_id = exo_mesh%eblk%ID

    !! Copy the cell information, block by block.
    offset = 0
    do n = 1, exo_mesh%num_eblk
      mesh%cnode(:,offset+1:offset+exo_mesh%eblk(n)%num_elem) = exo_mesh%eblk(n)%connect(:,:)
      mesh%cell_block(offset+1:offset+exo_mesh%eblk(n)%num_elem) = exo_mesh%eblk(n)%ID
      offset = offset + exo_mesh%eblk(n)%num_elem
    end do

    !! Take the side set data.
    call move_alloc (exo_mesh%sset, mesh%sset)
    if (.not.allocated(mesh%sset)) allocate(mesh%sset(0))

    !! Take the node set data.
    call move_alloc (exo_mesh%nset, mesh%nset)
    if (.not.allocated(mesh%nset)) allocate(mesh%nset(0))

  end subroutine import_exodus_mesh_iop


  subroutine identify_mesh_type (exo_mesh, mesh_type)

    use string_utilities

    type(exodus_mesh), intent(in) :: exo_mesh
    character(len=*), intent(out) :: mesh_type

    integer :: n
    character(len(exo_mesh%eblk(1)%elem_type)) :: elem_type

    mesh_type = 'UNKNOWN'

    select case (exo_mesh%num_dim)
    case (3)
      elem_type = raise_case(exo_mesh%eblk(1)%elem_type)
      select case (elem_type(1:3))
      case ('TET')
        do n = 1, exo_mesh%num_eblk
          if (raise_case(exo_mesh%eblk(n)%elem_type(1:3)) /= 'TET') return
          if (size(exo_mesh%eblk(n)%connect,dim=1) /= 4) return
        end do
        mesh_type = 'TET'
      case ('HEX')
        do n = 1, exo_mesh%num_eblk
          if (raise_case(exo_mesh%eblk(n)%elem_type(1:3)) /= 'HEX') return
          if (size(exo_mesh%eblk(n)%connect,dim=1) /= 8) return
        end do
        mesh_type = 'HEX'
      end select
    end select

  end subroutine identify_mesh_type


  subroutine destroy_mesh (mesh)
    type(external_mesh), intent(inout) :: mesh
    type(external_mesh) :: default
    if (associated(mesh%lnode)) deallocate(mesh%lnode)
    if (associated(mesh%link_block)) deallocate(mesh%link_block)
    if (associated(mesh%link_block_id)) deallocate(mesh%link_block_id)
    mesh = default  ! set default initialization values
  end subroutine destroy_mesh

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! COLLECT_MAIN_MESH
 !!
 !! This subroutine creates an EXTERNAL_MESH object corresponding to the main
 !! mesh, rummaging about in various Truchas modules to collect the data.
 !! This is a temporary bridge to the old mesh data structures needed by
 !! EM_DATA_PROXY to create the grid-to-grid mappings.  Note that the side set
 !! info is not defined as it is not required by EM_DATA_PROXY.
 !!

  subroutine collect_main_mesh (mesh)

    use parameter_module, only: ncells_tot, nnodes_tot
    use mesh_module, only: MMesh => Mesh, Vertex, CELL_TET, mesh_has_cblockid_data

    type(external_mesh), intent(out) :: mesh

    integer :: k

    if (is_IOP) then
      mesh%ncell = ncells_tot
      mesh%nnode = nnodes_tot
    end if

    if (global_all(MMesh(:)%Cell_Shape == CELL_TET)) then

      !! Tet-cell node array. Note the funky way Truchas describes a tet as a degenerate hex.
      mesh%mesh_type = 'TET'
      allocate(mesh%cnode(4,mesh%ncell))
      call collate (mesh%cnode(1,:), MMesh(:)%Ngbr_Vrtx_Orig(1))
      call collate (mesh%cnode(2,:), MMesh(:)%Ngbr_Vrtx_Orig(3))
      call collate (mesh%cnode(3,:), MMesh(:)%Ngbr_Vrtx_Orig(4))
      call collate (mesh%cnode(4,:), MMesh(:)%Ngbr_Vrtx_Orig(5))

    else

      !! Hex-cell node array.
      mesh%mesh_type = 'HEX'
      allocate(mesh%cnode(8,mesh%ncell))
      do k = 1, 8
        call collate (mesh%cnode(k,:), MMesh(:)%Ngbr_Vrtx_Orig(k))
      end do

    end if

    !! Node positions.
    allocate(mesh%x(3,mesh%nnode))
    do k = 1, 3
      call collate (mesh%x(k,:), vertex(:)%Coord(k))
    end do

    !! Cell block IDs.
    if (mesh_has_cblockid_data) then
      allocate(mesh%cell_block(mesh%ncell))
      call collate (mesh%cell_block, MMesh(:)%CBlockID)
    end if

    !! The old Truchas mesh structure doesn't save the actual list of block
    !! IDs when it reads the main exodus mesh initially.  This information
    !! can be recovered from the cell block ID array (it's just the list of
    !! distinct values) but since we have no need for it presently, we'll
    !! just leave that component (and the number of blocks) undefined for now.

    !! When the main mesh is read into the old Truchas mesh structure the
    !! compact Exodus side set data is unpacked into a large rank-3 array.
    !! Since the only user of this routine doesn't require this info we don't
    !! go to the effort of recreating the side set structures.

  end subroutine collect_main_mesh

  subroutine dump_external_mesh (mesh, unit)
    type(external_mesh), intent(in) :: mesh
    integer, intent(in) :: unit
    integer :: j
    write(unit,fmt='(a)') 'EXTERNAL_MESH('
    write(unit,fmt='(t4,a,i7)') 'NCELL=', mesh%ncell
    write(unit,fmt='(t4,a,i7)') 'NNODE=', mesh%nnode
    if (allocated(mesh%cnode)) then
      write(unit,fmt='(t4,a,t15,8i8)') 'CNODE=', mesh%cnode(:,1)
      write(unit,fmt='((t15,8i8))') mesh%cnode(:,2:)
    else
      write(unit,fmt='(t4,a)') 'CNODE => NULL()'
    end if
    if (allocated(mesh%x)) then
      write(unit,fmt='(t4,a,(t15,3es22.14))') 'X=', mesh%x
    else
      write(unit,fmt='(t4,a)') 'X => NULL()'
    end if
    write(unit,fmt='(t4,a)') 'MESH_TYPE= "' // trim(mesh%mesh_type) // '"'
    write(unit,fmt='(t4,a,i4)') 'NBLOCK=', mesh%nblock
    if (allocated(mesh%block_id)) then
      write(unit,fmt='(t4,a,(t15,10i5))') 'BLOCK_ID=', mesh%block_id
    else
      write(unit,fmt='(t4,a)') 'BLOCK_ID => NULL()'
    end if
    if (allocated(mesh%cell_block)) then
      write(unit,fmt='(t4,a,(t15,10i5))') 'CELL_BLOCK=', mesh%cell_block
    else
      write(unit,fmt='(t4,a)') 'CELL_BLOCK => NULL()'
    end if
    if (allocated(mesh%sset)) then
      write(unit,fmt='(t4,a)') 'SSET= ...'
      do j = 1, size(mesh%sset)
        call mesh%sset(j)%dump(unit)
      end do
    else
      write(unit,fmt='(t4,a)') 'SSET => NULL()'
    end if
    if (allocated(mesh%nset)) then
      write(unit,fmt='(t4,a)') 'NSET= ...'
      do j = 1, size(mesh%nset)
        call mesh%nset(j)%dump(unit)
      end do
    else
      write(unit,fmt='(t4,a)') 'NSET => NULL()'
    end if
    write(unit,fmt='(a)') ')'
  end subroutine dump_external_mesh

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SIDE_SET_NODE_LIST
 !!
 !! This function returns a pointer to a rank-1 integer array holding the side
 !! set node list for side set N in MESH.  Currently, only 4-node tets and
 !! 8-node hexes are handled.  MESH must be well-defined; if it is not, the
 !! behavior of the function will be unpredictable. A null pointer is returned
 !! if anything not understood is encountered.  The caller is responsible for
 !! deallocating the returned pointer.
 !!
 !! NB: The only proper use of the function is as the target of a pointer
 !! assignment; if used otherwise (in an expression, e.g.) a memory leak will
 !! result.
 !!

  function side_set_node_list (mesh, n) result (list)

    use cell_topology, only: TETRA4_FACE_VERT, HEX8_FACE_VERT

    type(external_mesh), intent(in) :: mesh
    integer, intent(in) :: n
    integer, pointer :: list(:)

    integer :: i, j, k, fsize, offset
    integer, pointer :: fvert(:,:)

    list => null()

    if (n < 1 .or. n > size(mesh%sset)) return

    select case (mesh%mesh_type)
    case ('TET')
      fvert => TETRA4_FACE_VERT
    case ('HEX')
      fvert => HEX8_FACE_VERT
    case default
      return
    end select

    fsize = size(fvert,1)

    allocate(list(fsize*mesh%sset(n)%num_side))

    offset = 0
    do i = 1, mesh%sset(n)%num_side
      j = mesh%sset(n)%elem(i)
      k = mesh%sset(n)%face(i)
      list(offset+1:offset+fsize) = mesh%cnode(fvert(:,k),j)
      offset = offset + fsize
    end do

  end function side_set_node_list

end module mesh_importer
