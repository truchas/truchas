!!
!! EXODUS_MESH_WRITER
!!
!! Neil N. Carlson <nnc@lanl.gov> 20 Sep 2004
!! Last revised 10 Nov 2004
!!
!! This module provides a procedure to write a mesh to an Exodus II format
!! file.  The Fortran 90 interface to the netCDF library is used to write
!! the data directly to the file, rather than Sandia's official Exodus II
!! library, which is proprietary.  The Exodus II manual [1] (out-of-date),
!! direct inspection of Exodus II files produced by Cubit, and private
!! communication with maintainers have been used to determine what netCDF
!! entities are required to form a valid Exodus II file.
!!
!! This is a private module; application code should only use the top-level
!! module EXODUS.
!!
!! [1] L.A.Schoof and V.R.Yarberry, "Exodus II: A Finite Element Data Model",
!!     Sandia report SAND92-2137.  This can be obtained at
!!     http://endo.sandia.gov/SEACAS/Documentation/exodusII.pdf
!!

module exodus_mesh_writer

  use netcdf
  use string_utilities
  use exodus_mesh_type
  use exodus_errors
  implicit none
  private

  public :: write_exodus_mesh

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! WRITE_EXODUS_MESH
 !!
 !! This routine writes the mesh data structure MESH to the Exodus II file PATH.
 !! The input character arguments CREATOR and VERSION should contain the name
 !! and version of the program responsible for creating the mesh.  If an error
 !! occurs, a nonzero value is returned in the optional argurment STAT, if
 !! present, otherwise an error message is written to stdout and the program
 !! is halted.
 !!

  subroutine write_exodus_mesh (path, mesh, creator, version, stat)

    character(len=*),  intent(in)  :: path
    type(exodus_mesh), intent(in)  :: mesh
    character(len=*),  intent(in)  :: creator
    character(len=*),  intent(in)  :: version
    integer, optional, intent(out) :: stat

    integer :: ncid, status

    !! Verify that the mesh is well-enough defined to write.
    if (.not.defined(mesh)) then
      call set_stat (EX_BADMESH, stat)
      return
    end if
print *, '$'

    !! Create the exodus mesh file (a netCDF data set).
    status = nf90_create(path, NF90_CLOBBER, ncid)
    if (status /= NF90_NOERR) then
      call set_stat (EX_EOPENW, stat)
      return
    end if
print *, '$$'

    call write_global_data (ncid, mesh, status)
    if (status /= 0) then
      call set_stat (EX_EWGLOB, stat)
      return
    end if
print *, '$$$'

    call write_nodes (ncid, mesh, status)
    if (status /= 0) then
      call set_stat (EX_EWNODE, stat)
      return
    end if
print *, '$$$$'

    call write_element_blocks (ncid, mesh, status)
    if (status /= 0) then
      call set_stat (EX_EWEBLK, stat)
      return
    end if
print *, '$$$$$'

    call write_node_sets (ncid, mesh, status)
    if (status /= 0) then
      call set_stat (EX_EWNSET, stat)
      return
    end if
print *, '$$$$$$'

    call write_side_sets (ncid, mesh, status)
    if (status /= 0) then
      call set_stat (EX_EWSSET, stat)
      return
    end if
print *, '$$$$$$$'

    call write_qa_record (ncid, creator, version, status)
    if (status /= 0) then
      call set_stat (EX_EWQA, stat)
      return
    end if
print *, '$$$$$$$$'

    !! Close the exodus mesh file.
    status = nf90_close(ncid)
    if (status /= NF90_NOERR) then
      call set_stat (EX_ECLOSE, stat)
      return
    end if
    
    if (present(stat)) stat = 0

  end subroutine write_exodus_mesh

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! WRITE_NODES
 !!
 !! This auxillary procedure writes the node data from the mesh data structure
 !! MESH to the Exodus II file associated with the NetCDF file handle NCID.  If
 !! an error occurs, the netCDF error code is returned in STATUS; otherwise
 !! STATUS is returned with the value CF90_NOERR (=0).
 !!
 !! NOTES
 !! o The names of the coordinates seem to be a required element of the format.
 !!   Here they've been hardwired to 'x', 'y', and 'z'; we don't really care.
 !!

  subroutine write_nodes (ncid, mesh, status)

    integer,           intent(in)  :: ncid
    type(exodus_mesh), intent(in)  :: mesh
    integer,           intent(out) :: status

    integer :: n, dim1, dim2, namesID, coordID
#ifdef PGI_COMPILER_WORKAROUND
    character(len=2) :: NAMES(3)
    NAMES = (/ 'x', 'y', 'z'/) // char(0)
#else
    character(len=2), parameter :: NAMES(3) = (/ 'x', 'y', 'z'/) // char(0)
#endif

    !! Define the dimensions.
    status = nf90_def_dim(ncid, 'num_dim',   mesh%num_dim,  dim2)
    if (status /= NF90_NOERR) return
    status = nf90_def_dim(ncid, 'num_nodes', mesh%num_node, dim1)
    if (status /= NF90_NOERR) return

    !! Define the coordinate array.
    status = nf90_def_var(ncid, 'coord', NF90_DOUBLE, (/dim1,dim2/), coordID)
    if (status /= NF90_NOERR) return

    !! Define the coordinate names array.
    status = nf90_inq_dimid(ncid, 'len_string', dim1)
    if (status /= NF90_NOERR) return
    status = nf90_def_var(ncid, 'coor_names', NF90_CHAR, (/dim1,dim2/), namesID)
    if (status /= NF90_NOERR) return

    !! Put the dataset into data mode.
    status = nf90_enddef(ncid)
    if (status /= NF90_NOERR) return

    !! Write the coordinate data.
    status = nf90_put_var(ncid, coordID, transpose(mesh%coord))
    if (status /= NF90_NOERR) return

    !! Write the coordinate names array.
    do n = 1, mesh%num_dim
      status = nf90_put_var(ncid, namesID, NAMES(n), start=(/1,n/))
      if (status /= NF90_NOERR) return
    end do

    !! Return the dataset to define mode.
    status = nf90_redef(ncid)

  end subroutine write_nodes

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! WRITE_ELEMENT_BLOCKS
 !!
 !! This auxillary procedure writes the element block data from the mesh data
 !! structure MESH to the Exodus II file associated with the NetCDF file handle
 !! NCID.  If an error occurs, the netCDF error code is returned in STATUS;
 !! otherwise STATUS is returned with the value NF90_NOERR (=0).
 !!
 !! NOTES
 !! o Element attributes are not supported.
 !! o A mask array 'eb_status' is written that indicates all blocks are
 !!   active (value 1).  This pertains to a parallel database, and may or may
 !!   not be a required element of the format, though Cubit does produce this
 !!   array in its Exodus output.
 !!

  subroutine write_element_blocks (ncid, mesh, status)

    integer,           intent(in)  :: ncid
    type(exodus_mesh), intent(in)  :: mesh
    integer,           intent(out) :: status

    integer :: n, dim1, dim2
    integer :: idID, statusID, connectID(mesh%num_eblk)

    !! Define dimensions for the element blocks.
    status = nf90_def_dim(ncid, 'num_elem',   mesh%num_elem, dim1)
    if (status /= NF90_NOERR) return
    status = nf90_def_dim(ncid, 'num_el_blk', mesh%num_eblk, dim1)
    if (status /= NF90_NOERR) return

    !! Define the element block status array.
    status = nf90_def_var(ncid, 'eb_status', NF90_INT, dim1, statusID)
    if (status /= NF90_NOERR) return

    !! Define the element block ID array.
    status = nf90_def_var(ncid, 'eb_prop1', NF90_INT, dim1, idID)
    if (status /= NF90_NOERR) return
    status = nf90_put_att(ncid, idID, 'name', 'ID')
    if (status /= NF90_NOERR) return

    do n = 1, mesh%num_eblk
      !! Define the dimensions for this element block.
      status = nf90_def_dim(ncid, 'num_el_in_blk'//i_to_c(n),  size(mesh%eblk(n)%connect,dim=2), dim2)
      if (status /= NF90_NOERR) return
      status = nf90_def_dim(ncid, 'num_nod_per_el'//i_to_c(n), size(mesh%eblk(n)%connect,dim=1), dim1)
      if (status /= NF90_NOERR) return

      !! Define the connections array for this element block.
      status = nf90_def_var(ncid, 'connect'//i_to_c(n), NF90_INT, (/dim1, dim2/), connectID(n))
      if (status /= NF90_NOERR) return
      status = nf90_put_att(ncid, connectID(n), 'elem_type', trim(mesh%eblk(n)%elem_type))
      if (status /= NF90_NOERR) return
    end do

    !! Put the dataset into data mode.
    status = nf90_enddef(ncid)
    if (status /= NF90_NOERR) return

    !! Write the element block status array (all 1's)
    status = nf90_put_var(ncid, statusID, spread(1,dim=1,ncopies=mesh%num_eblk))
    if (status /= NF90_NOERR) return

    !! Write the element block ID array data
    status = nf90_put_var(ncid, idID, mesh%eblk%ID)
    if (status /= NF90_NOERR) return

print *, '>'
    !! Write the connections array data.
    do n = 1, mesh%num_eblk
      status = nf90_put_var(ncid, connectID(n), mesh%eblk(n)%connect)
      if (status /= NF90_NOERR) return
    end do
print *, '>>'

    !! Return the dataset to define mode.
    status = nf90_redef(ncid)

  end subroutine write_element_blocks

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! WRITE_NODE_SETS
 !!
 !! This auxillary procedure writes the node set data from the mesh data
 !! structure MESH to the Exodus II file associated with the netCDF file handle
 !! NCID.  If an error occurs, the netCDF error code is returned in STATUS;
 !! otherwise STATUS is returned with the value NF90_NOERR (=0).
 !!
 !! NOTES
 !! o Distribution factors are not supported.
 !! o A mask array 'ns_status' is written that indicates all node sets are
 !!   active (value 1).  This pertains to a parallel database, and may or may
 !!   not be a required element of the format, though Cubit does produce this
 !!   array in its Exodus output.
 !!

  subroutine write_node_sets (ncid, mesh, status)

    integer,           intent(in)  :: ncid
    type(exodus_mesh), intent(in)  :: mesh
    integer,           intent(out) :: status

    integer :: n, dimid
    integer :: idID, statusID, nodeID(mesh%num_nset)

    if (mesh%num_nset <= 0) then
      status = NF90_NOERR
      return
    end if

    !! Define the node set dimension.
    status = nf90_def_dim(ncid, 'num_node_sets', mesh%num_nset, dimid)
    if (status /= NF90_NOERR) return

    !! Define the node set status array.
    status = nf90_def_var(ncid, 'ns_status', NF90_INT, dimid, statusID)
    if (status /= NF90_NOERR) return

    !! Define the node set ID array.
    status = nf90_def_var(ncid, 'ns_prop1', NF90_INT, dimid, idID)
    if (status /= NF90_NOERR) return
    status = nf90_put_att(ncid, idID, 'name', 'ID')
    if (status /= NF90_NOERR) return

    do n = 1, mesh%num_nset
      !! Define the dimension for this node set.
      status = nf90_def_dim(ncid, 'num_nod_ns'//i_to_c(n), mesh%nset(n)%num_node, dimid)
      if (status /= NF90_NOERR) return

      !! Define the node array for this node set.
      status = nf90_def_var(ncid, 'node_ns'//i_to_c(n), NF90_INT, dimid, nodeID(n))
      if (status /= NF90_NOERR) return
    end do

    !! Put the dataset into data mode.
    status = nf90_enddef(ncid)
    if (status /= NF90_NOERR) return

    !! Write the node set status array (all 1's)
    status = nf90_put_var(ncid, statusID, spread(1,dim=1,ncopies=mesh%num_nset))
    if (status /= NF90_NOERR) return

    !! Write the node set ID array data.
    status = nf90_put_var(ncid, idID, mesh%nset%ID)
    if (status /= NF90_NOERR) return

    !! Write the node set node array data.
    do n = 1, mesh%num_nset
      status = nf90_put_var(ncid, nodeID(n), mesh%nset(n)%node)
      if (status /= NF90_NOERR) return
    end do

    !! Return the dataset to define mode.
    status = nf90_redef(ncid)

  end subroutine write_node_sets

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! WRITE_SIDE_SETS
 !!
 !! This auxillary procedure writes the side set data from the mesh data
 !! structure MESH to the Exodus II file associated with the netCDF file handle
 !! NCID.  If an error occurs, the netCDF error code is returned in STATUS;
 !! otherwise STATUS is returned with the value NF90_NOERR (=0).
 !!
 !! NOTES
 !! o Distribution factors are not supported.
 !! o A mask array 'ss_status' is written that indicates all side sets are
 !!   active (value 1).  This pertains to a parallel database, and may or may
 !!   not be a required element of the format, though Cubit does produce this
 !!   array in its Exodus output.
 !!

  subroutine write_side_sets (ncid, mesh, status)

    integer,           intent(in)  :: ncid
    type(exodus_mesh), intent(in)  :: mesh
    integer,           intent(out) :: status

    integer :: n, dimid
    integer :: idID, statusID, elemID(mesh%num_sset), faceID(mesh%num_sset)

    if (mesh%num_sset <= 0) then
      status = NF90_NOERR
      return
    end if

    !! Define the side set dimension.
    status = nf90_def_dim(ncid, 'num_side_sets', mesh%num_sset, dimid)
    if (status /= NF90_NOERR) return

    !! Define the side set status array.
    status = nf90_def_var(ncid, 'ss_status', NF90_INT, dimid, statusID)
    if (status /= NF90_NOERR) return

    !! Define the side set ID array.
    status = nf90_def_var(ncid, 'ss_prop1', NF90_INT, dimid, idID)
    if (status /= NF90_NOERR) return
    status = nf90_put_att(ncid, idID, 'name', 'ID')
    if (status /= NF90_NOERR) return

    do n = 1, mesh%num_sset
      !! Define the dimension for this side set.
      status = nf90_def_dim(ncid, 'num_side_ss'//i_to_c(n), mesh%sset(n)%num_side, dimid)
      if (status /= NF90_NOERR) return

      !! Define the element and face arrays for this side set.
      status = nf90_def_var(ncid, 'elem_ss'//i_to_c(n), NF90_INT, dimid, elemID(n))
      if (status /= NF90_NOERR) return
      status = nf90_def_var(ncid, 'side_ss'//i_to_c(n), NF90_INT, dimid, faceID(n))
      if (status /= NF90_NOERR) return
    end do

    !! Put the dataset into data mode.
    status = nf90_enddef(ncid)
    if (status /= NF90_NOERR) return

    !! Write the side set status array (all 1's)
    status = nf90_put_var(ncid, statusID, spread(1,dim=1,ncopies=mesh%num_sset))
    if (status /= NF90_NOERR) return

    !! Write the side set ID array data.
    status = nf90_put_var(ncid, idID, mesh%sset%ID)
    if (status /= NF90_NOERR) return

    !! Write the side set element/face array data.
    do n = 1, mesh%num_sset
      status = nf90_put_var(ncid, elemID(n), mesh%sset(n)%elem)
      if (status /= NF90_NOERR) return
      status = nf90_put_var(ncid, faceID(n), mesh%sset(n)%face)
      if (status /= NF90_NOERR) return
    end do

    !! Return the dataset to define mode.
    status = nf90_redef(ncid)

  end subroutine write_side_sets

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! WRITE_GLOBAL_DATA
 !!
 !! This auxillary procedure writes some required global data to the Exodus II
 !! file associated with the netCDF file handle NCID.  If an error occurs, the
 !! netCDF error code is returned in STATUS; otherwise STATUS returns the value
 !! NF90_NOERR (=0).
 !!
 !! NOTES
 !! o The Exodus format requires that the database version number (DB_VERSION) and
 !!   the API version number (API_VERSION) be stored in the file.  These are
 !!   numbers associated with the official Exodus II library, which we obviously
 !!   aren't using.  The numbers written here are those found in the output
 !!   from Cubit 8.0.  We don't care, but perhaps some exodus utility does.
 !! o We use 64-bit reals for the floating point data stored in the mesh data
 !!   structure, so here we hardwire the floating point word size stored in
 !!   the file to 8.
 !! o Based on an inspection of the exodus output produced by Cubit the
 !!   character length dimensions defined here have be extended by 1, presumably
 !!   to account for an additional null character to terminate the string.
 !!

  subroutine write_global_data (ncid, mesh, status)

    integer,           intent(in)  :: ncid
    type(exodus_mesh), intent(in)  :: mesh
    integer,           intent(out) :: status

    integer :: dimid

    real, parameter :: API_VERSION = 3.04
    real, parameter :: DB_VERSION  = 2.03
    integer,parameter :: FP_WORD_SIZE = 8   ! we'll be writing 64-bit floating point data

    !! Define some common dimensions for the dataset.
    status = nf90_def_dim(ncid, 'len_string', MAX_STR_LENGTH+1,  dimid)
    if (status /= NF90_NOERR) return
    status = nf90_def_dim(ncid, 'len_line',   MAX_LINE_Length+1, dimid)
    if (status /= NF90_NOERR) return
    
    status = nf90_def_dim(ncid, 'time_step', NF90_UNLIMITED, dimid)
    if (status /= NF90_NOERR) return

    !! Write the global attributes for the dataset.
    status = nf90_put_att(ncid, NF90_GLOBAL, 'api_version', API_VERSION)
    if (status /= NF90_NOERR) return
    status = nf90_put_att(ncid, NF90_GLOBAL, 'version', DB_VERSION)
    if (status /= NF90_NOERR) return
    status = nf90_put_att(ncid, NF90_GLOBAL, 'floating_point_word_size', FP_WORD_SIZE)
    if (status /= NF90_NOERR) return
    status = nf90_put_att(ncid, NF90_GLOBAL, 'title', mesh%title)

  end subroutine write_global_data

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! WRITE_QA_RECORD
 !!
 !! This auxillary routine writes a 'quality assurance' record to the Exodus II
 !! file associated to the netCDF file handle NCID.  This probably isn't a
 !! required element of the format, but it seems useful to include it anyway.
 !! The character arguments CREATOR and VERSION should describe the name and
 !! version of the program that created the mesh being written.  The date and
 !! time, which are also written, are automatically generated.
 !!

  subroutine write_qa_record (ncid, creator, version, status)

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: creator, version
    integer, intent(out) :: status

    integer :: j, dim1, dim2, dim3, varid
    character(len=MAX_STR_LENGTH+1) :: qa(4)
    character(len=8)  :: date
    character(len=10) :: time

    qa(1) = trim(adjustl(creator)) // char(0)
    qa(2) = trim(adjustl(version)) // char(0)

    call date_and_time (date, time)

    qa(3) = date(5:6) // '/' // date(7:8) // '/' // date(1:4) // char(0)
    qa(4) = time(1:2) // ':' // time(3:4) // ':' // time(5:6) // char(0)

    status = nf90_inq_dimid (ncid, 'len_string', dim1)
    if (status /= NF90_NOERR) return
    status = nf90_def_dim (ncid, 'four', 4, dim2)
    if (status /= NF90_NOERR) return
    status = nf90_def_dim (ncid, 'num_qa_rec', 1, dim3)
    if (status /= NF90_NOERR) return
    status = nf90_def_var (ncid, 'qa_records', NF90_CHAR, (/ dim1, dim2, dim3 /), varid)
    if (status /= NF90_NOERR) return

    status = nf90_enddef(ncid)
    if (status /= NF90_NOERR) return
    do j = 1, 4
      status = nf90_put_var (ncid, varid, qa(j), start=(/1,j,1/), count=(/len_trim(qa(j)),1,1/))
      if (status /= NF90_NOERR) return
    end do
    status = nf90_redef(ncid)

  end subroutine write_qa_record

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
      write(unit=0,fmt='(2a)') 'WRITE_EXODUS_MESH: fatal error: ', exo_err_str(errno)
      stop
    end if
  end subroutine set_stat

end module exodus_mesh_writer
