!!
!! EXODUS_MESH_READER
!!
!! Neil N. Carlson <nnc@lanl.gov> 5 Jul 2004
!! Last revised 10 Nov 2004
!!
!! This module provides a procedure to read mesh data from an Exodus II file.
!! The Fortran 90 interface to the netCDF library is used to read data directly
!! from the file, rather than using Sandia's official Exodus II library, which
!! is proprietary.  The Exodus II manual [1] (out-of-date), direct inspection
!! of Exodus II files produced by Cubit, and private communication with
!! maintainers have been used to determine what netCDF entities may be present
!! in the file.
!!
!! This is a private module; application code should only use the top-level
!! module EXODUS.
!!
!! Not all the possible mesh data present in an Exodus II file is read; only
!! that which is needed to fill out the mesh data structure.  In particular,
!! the following optional data objects are currently ignored, if present:
!!  o The node and element number maps which specify a mapping from the
!!    internal, 1-based, contiguous numbering to a user-space numbering.
!!  o The element order map which specifies a 'good' order to process the
!!    elements.
!!  o Distribution factors for side and node sets.
!!
!! [1] L.A.Schoof and V.R.Yarberry, "Exodus II: A Finite Element Data Model",
!!     Sandia report SAND92-2137.  This can be obtained at
!!     http://endo.sandia.gov/SEACAS/Documentation/exodusII.pdf
!!

module exodus_mesh_reader

  use netcdf
  use string_utilities
  use exodus_mesh_type
  use exodus_errors
  implicit none
  private

  public :: read_exodus_mesh

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_EXODUS_MESH
 !!
 !! This routine reads mesh data from the Exodus II file PATH, returning it
 !! in the mesh data structure MESH.  If an error occurs, a nonzero value is
 !! returned in the optional argurment STAT, if present, otherwise an error
 !! message is written to stdout and the program is halted.
 !!

  subroutine read_exodus_mesh (path, mesh, stat)

    character(len=*),  intent(in)  :: path
    type(exodus_mesh), intent(out) :: mesh
    integer, optional, intent(out) :: stat

    integer :: ncid, status, n

    !! Open the exodus mesh file (a netCDF dataset).
    status = nf90_open(path, NF90_NOWRITE, ncid)
    if (status /= NF90_NOERR) then
      call set_stat (EX_EOPENR, stat)
      return
    end if

    !! Verify that the dataset isn't a 'large model'.
    status = nf90_get_att(ncid, NF90_GLOBAL, 'file_size', n)
    if (status == NF90_NOERR) then
      if (n /= 0) then
        call set_stat (EX_ELFILE, stat)
        return
      end if
    end if

    !! Read the dataset title.
    status = nf90_get_att(ncid, NF90_GLOBAL, 'title', mesh%title)
    if (status /= NF90_NOERR) then
      call set_stat (EX_ERGLOB, stat)
      return
    end if

    !! Read the node position data.
    call read_nodes (ncid, mesh, status)
    if (status /= 0) then
      call set_stat (EX_ERNODE, stat)
      return
    end if

    !! Read the cell and cell block ID data.
    call read_element_blocks (ncid, mesh, status)
    if (status /= 0) then
      call set_stat (EX_EREBLK, stat)
      return
    end if

    !! Read the node set data.
    call read_node_sets (ncid, mesh, status)
    if (status /= 0) then
      call set_stat (EX_ERNSET, stat)
      return
    end if

    !! Read the side set data.
    call read_side_sets (ncid, mesh, status)
    if (status /= 0) then
      call set_stat (EX_ERSSET, stat)
      return
    end if

    !! Close the exodus mesh file.
    status = nf90_close(ncid)
    if (status /= NF90_NOERR) then
      call set_stat (EX_ECLOSE, stat)
      return
    end if
    
    if (present(stat)) stat = 0

  end subroutine read_exodus_mesh

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_NODES
 !!
 !! This auxillary procedure reads the node data from the Exodus II file
 !! associated with the netCDF file handle NCID, and initializes the
 !! corresponding components of the derived type MESH.  If an error occurs,
 !! the netCDF error code is returned in STATUS; otherwise STATUS is returned
 !! with the value NF90_NOERR (=0).
 !!
 !! NOTES
 !! o The coordinate name array is ignored.
 !!

  subroutine read_nodes (ncid, mesh, status)

    integer,           intent(in)    :: ncid
    type(exodus_mesh), intent(inout) :: mesh
    integer,           intent(out)   :: status

    integer :: dimid, varid, n

    !! Get the number of nodes.
    status = nf90_inq_dimid(ncid, 'num_nodes', dimid)
    if (status /= NF90_NOERR) return
    status = nf90_inquire_dimension(ncid, dimid, len=mesh%num_node)
    if (status /= NF90_NOERR) return

    !! Get the dimension of the node coordinate space.
    status = nf90_inq_dimid(ncid, 'num_dim', dimid)
    if (status /= NF90_NOERR) return
    status = nf90_inquire_dimension(ncid, dimid, len=mesh%num_dim)
    if (status /= NF90_NOERR) return

    !! Allocate the node coordinate array.
    allocate(mesh%coord(mesh%num_dim,mesh%num_node))
#ifdef NAG_COMPILER
    !! When compiled with -nan, something in nf90_get_var below copies this array,
    !! triggering an arithmetic exception on some x86 platforms only.  Not really
    !! supposed to use -nan on just part of the program for this very reason.  NNC 6/24/2005.
    mesh%coord = 0.0
#endif

    !! Read the node positions, coordinate by coordinate.
    status = nf90_inq_varid(ncid, 'coord', varid)
    if (status /= NF90_NOERR) return
    do n = 1, mesh%num_dim
      status = nf90_get_var(ncid, varid, mesh%coord(n,:), start=(/1,n/))!, count=(/mesh%num_node,1/))
      if (status /= NF90_NOERR) return
    end do

  end subroutine read_nodes

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_ELEMENT_BLOCKS
 !!
 !! This auxillary procedure reads the element block data from the Exodus II
 !! file associated with the netCDF file handle NCID, and initializes the
 !! corresponding components of the derived type MESH.  If an error occurs,
 !! the netCDF error code is returned in STATUS; otherwise STATUS is returned
 !! with the value NF90_NOERR (=0).
 !!
 !! NOTES
 !! o The status array is ignored.
 !! o Element attributes, if any, are ignored.
 !!

  subroutine read_element_blocks (ncid, mesh, status)

    integer,           intent(in)    :: ncid
    type(exodus_mesh), intent(inout) :: mesh
    integer,           intent(out)   :: status

    integer :: dimid, varid, n, nodes_per_elem

    !! Read the total number of elements.
    status = nf90_inq_dimid(ncid, 'num_elem', dimid)
    if (status /= NF90_NOERR) return
    status = nf90_inquire_dimension(ncid, dimid, len=mesh%num_elem)
    if (status /= NF90_NOERR) return

    !! Read the number of element blocks.
    status = nf90_inq_dimid(ncid, 'num_el_blk', dimid)
    if (status /= NF90_NOERR) return
    status = nf90_inquire_dimension(ncid, dimid, len=mesh%num_eblk)
    if (status /= NF90_NOERR) return

    allocate(mesh%eblk(mesh%num_eblk))

    !! Read the list of element block IDs.
    status = nf90_inq_varid(ncid, 'eb_prop1', varid)
    if (status /= NF90_NOERR) return
    status = nf90_get_var(ncid, varid, mesh%eblk%ID)
    if (status /= NF90_NOERR) return

    do n = 1, mesh%num_eblk
      !! Read the number of elements in this block.
      status = nf90_inq_dimid(ncid, 'num_el_in_blk'//i_to_c(n), dimid)
      if (status /= NF90_NOERR) return
      status = nf90_inquire_dimension(ncid, dimid, len=mesh%eblk(n)%num_elem)
      if (status /= NF90_NOERR) return

      !! Read the number of nodes per element for this block.
      status = nf90_inq_dimid(ncid, 'num_nod_per_el'//i_to_c(n), dimid)
      if (status /= NF90_NOERR) return
      status = nf90_inquire_dimension(ncid, dimid, len=nodes_per_elem)
      if (status /= NF90_NOERR) return

      allocate(mesh%eblk(n)%connect(nodes_per_elem,mesh%eblk(n)%num_elem))

      !! Read the element type for this block.
      status = nf90_inq_varid(ncid, 'connect'//i_to_c(n), varid)
      if (status /= NF90_NOERR) return
      status = nf90_get_att(ncid, varid, 'elem_type', mesh%eblk(n)%elem_type)
      if (status /= NF90_NOERR) return
      
      !! Remove any null-termination that may be present.
      mesh%eblk(n)%elem_type = ctrim(mesh%eblk(n)%elem_type)

      !! Read the connections array for this block.
      status = nf90_get_var(ncid, varid, mesh%eblk(n)%connect)
      if (status /= NF90_NOERR) return
    end do
    
  end subroutine read_element_blocks

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_NODE_SETS
 !!
 !! This auxillary procedure reads the node set data from the Exodus II file
 !! associated with the netCDF file handle NCID, and initializes the node set
 !! components of the derived type MESH.  If an error occurs, the netCDF error
 !! code is returned in STATUS; otherwise STATUS is returned with the value
 !! NF90_NOERR (=0).  If the Exodus II file contains no node set data, NUM_NSET
 !! is assigned the value 0, and the remaining node set components should be
 !! ignored.
 !!
 !! NOTES
 !! o Distribution factors are ignored.
 !! o The status array is ignored.
 !!

  subroutine read_node_sets (ncid, mesh, status)

    integer,           intent(in)    :: ncid
    type(exodus_mesh), intent(inout) :: mesh
    integer,           intent(out)   :: status

    integer :: dimid, varid, n

    !! Read the number of node sets, if any.
    status = nf90_inq_dimid(ncid, 'num_node_sets', dimid)
    if (status /= NF90_NOERR) then  ! we conclude there are no node sets.
      mesh%num_nset = 0
      status = NF90_NOERR
      return
    end if
    status = nf90_inquire_dimension(ncid, dimid, len=mesh%num_nset)
    if (status /= NF90_NOERR) return

    allocate(mesh%nset(mesh%num_nset))

    !! Read the node set ID array.
    status = nf90_inq_varid(ncid, 'ns_prop1', varid)
    if (status /= NF90_NOERR) return
    status = nf90_get_var(ncid, varid, mesh%nset%ID)
    if (status /= NF90_NOERR) return

    do n = 1, mesh%num_nset
      !! Read the number of nodes in this node set.
      status = nf90_inq_dimid(ncid, 'num_nod_ns'//i_to_c(n), dimid)
      if (status /= NF90_NOERR) return
      status = nf90_inquire_dimension(ncid, dimid, len=mesh%nset(n)%num_node)
      if (status /= NF90_NOERR) return

      allocate(mesh%nset(n)%node(mesh%nset(n)%num_node))

      !! Read the node array for this node set.
      status = nf90_inq_varid(ncid, 'node_ns'//i_to_c(n), varid)
      if (status /= NF90_NOERR) return
      status = nf90_get_var(ncid, varid, mesh%nset(n)%node)
      if (status /= NF90_NOERR) return
    end do

  end subroutine read_node_sets

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_SIDE_SETS
 !!
 !! This auxillary procedure reads the side set data from the Exodus II file
 !! associated with the netCDF file handle NCID, and initializes the side set
 !! components of the derived type MESH.  If an error occurs, the netCDF error
 !! code is returned in STATUS; otherwise STATUS is returned with the value
 !! NF90_NOERR (=0).  If the Exodus II file contains no side set data, NUM_SSET
 !! is assigned the value 0, and the remaining side set components should be
 !! ignored.
 !!
 !! NOTES
 !! o Distribution factors are ignored.
 !! o The status array is ignored.
 !!

  subroutine read_side_sets (ncid, mesh, status)

    integer,           intent(in)    :: ncid
    type(exodus_mesh), intent(inout) :: mesh
    integer,           intent(out)   :: status

    integer :: dimid, varid, n

    !! Read the number of side sets, if any.
    status = nf90_inq_dimid(ncid, 'num_side_sets', dimid)
    if (status /= NF90_NOERR) then  ! we conclude there are no side sets.
      mesh%num_sset = 0
      status = NF90_NOERR
      return
    end if
    status = nf90_inquire_dimension(ncid, dimid, len=mesh%num_sset)
    if (status /= NF90_NOERR) return

    allocate(mesh%sset(mesh%num_sset))

    !! Read the side set ID array.
    status = nf90_inq_varid(ncid, 'ss_prop1', varid)
    if (status /= NF90_NOERR) return
    status = nf90_get_var(ncid, varid, mesh%sset%ID)
    if (status /= NF90_NOERR) return

    do n = 1, mesh%num_sset
      !! Read the number of sides in this side set.
      status = nf90_inq_dimid(ncid, 'num_side_ss'//i_to_c(n), dimid)
      if (status /= NF90_NOERR) return
      status = nf90_inquire_dimension(ncid, dimid, len=mesh%sset(n)%num_side)
      if (status /= NF90_NOERR) return

      allocate(mesh%sset(n)%elem(mesh%sset(n)%num_side))
      allocate(mesh%sset(n)%face(mesh%sset(n)%num_side))

      !! Read the element/face arrays for this side set.
      status = nf90_inq_varid(ncid, 'elem_ss'//i_to_c(n), varid)
      if (status /= NF90_NOERR) return
      status = nf90_get_var(ncid, varid, mesh%sset(n)%elem)
      if (status /= NF90_NOERR) return
      status = nf90_inq_varid(ncid, 'side_ss'//i_to_c(n), varid)
      if (status /= NF90_NOERR) return
      status = nf90_get_var(ncid, varid, mesh%sset(n)%face)
      if (status /= NF90_NOERR) return
    end do

  end subroutine read_side_sets

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SET_STAT
 !!
 !! This is the error handler.  If the optional output argument STAT is present
 !! then this routine simply sets its value to the value of the input argument
 !! ERRNO, which specifies the error number, and then returns.  If STAT is not
 !! present, then the error message corresponding to the value of ERRNO is
 !! written to stdout, and the program halted.
 !!

  subroutine set_stat (errno, stat)
    integer, intent(in) :: errno
    integer, optional, intent(out) :: stat
    if (present(stat)) then
      stat = errno
    else
      write(unit=0,fmt='(2a)') 'READ_EXODUS_MESH: fatal error: ', exo_err_str(errno)
      stop
    end if
  end subroutine set_stat

end module exodus_mesh_reader
