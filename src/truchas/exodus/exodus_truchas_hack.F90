!!
!! EXODUS_TRUCHAS_HACK
!!
!! Neil N. Carlson <nnc@lanl.gov> 25 Sep 2004
!! Revised 25 Mar 2010
!!
!! This module provides procedures to read certain data from an ExodusII file.
!! These are intended solely to accomodate specific needs in Truchas as it is
!! currently designed.  Please DO NOT use these for any other purpose.
!!

module exodus_truchas_hack

  use netcdf
  use string_utilities
  implicit none
  private

  public :: read_exodus_mesh_size, read_exodus_side_sets

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_EXODUS_MESH_SIZE
 !!
 !! This routine reads the number of nodes and elements from the ExodusII
 !! file PATH.  If an error occurs, a nonzero value is returned in STATUS.
 !!

  subroutine read_exodus_mesh_size (path, num_node, num_elem, stat)

    character(len=*), intent(in) :: path
    integer, intent(out) :: num_node, num_elem
    integer, intent(out), optional :: stat

    integer :: status, ncid, dimid

    if (present(stat)) stat = 0

    !! Open the exodus mesh file (a netCDF dataset).
    status = nf90_open(path, NF90_NOWRITE, ncid)
    if (status /= NF90_NOERR) goto 666

    !! Get the number of nodes.
    status = nf90_inq_dimid(ncid, 'num_nodes', dimid)
    if (status /= NF90_NOERR) goto 666
    status = nf90_inquire_dimension(ncid, dimid, len=num_node)
    if (status /= NF90_NOERR) goto 666

    !! Read the number of elements.
    status = nf90_inq_dimid(ncid, 'num_elem', dimid)
    if (status /= NF90_NOERR) goto 666
    status = nf90_inquire_dimension(ncid, dimid, len=num_elem)
    if (status /= NF90_NOERR) goto 666

    !! Close the exodus mesh file.
    status = nf90_close(ncid)
    if (status /= NF90_NOERR) goto 666

    return

    !! Handle NetCDF errors.
666 if (present(stat)) then ! caller handles the error
      stat = status
      return
    else  ! we deal with it
      write(0,fmt='(2a)') 'READ_EXODUS_MESH_SIZE: FATAL: ', nf90_strerror(status)
      stop
    end if

  end subroutine read_exodus_mesh_size

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_EXODUS_SIDE_SETS
 !!
 !! This routine reads the side set information from the ExodusII file PATH.
 !! The data is returned in the supplied rank-3 MESH_FACE_SET_TOT array; its
 !! shape should be number of side sets, by number of faces per cell (6), by
 !! number of elements.  The old Truchas mesh data structure stores tet and
 !! wedge elements as specific forms of degenerate hexes, and it also uses a
 !! different numbering of the sides of a hex.  (See MESH_READ from the
 !! MESH_INPUT_MODULE for details.)  This routine takes care of translating
 !! the ExodusII side info appropriately.
 !!

  subroutine read_exodus_side_sets (path, mesh_face_set_tot, status)

    character(len=*), intent(in) :: path
    integer, pointer :: mesh_face_set_tot(:,:,:)
    integer, intent(out) :: status

    integer :: ncid, dimid, varid, num_elem, num_sset, num_eblk, num_side, n, j, tside
    integer, allocatable :: eblk_size(:), eblk_type(:), sset_id(:), elem(:), side(:)

    !! Exodus side numbering to (old) Truchas side numbering.
    integer, parameter :: TET_SIDE_MAP(4)=(/4,1,3,5/)
    integer, parameter :: WEDGE_SIDE_MAP(5)=(/5,1,2,3,4/)
    integer, parameter :: HEX_SIDE_MAP(6)=(/2,4,1,3,5,6/)

    !! Element types identified by their number of nodes.
    integer, parameter :: TET=4, WEDGE=6, HEX=8

    if (associated(mesh_face_set_tot)) then
      status = 1
      return
    end if

    !! Open the exodus mesh file (a NetCDF dataset).
    status = nf90_open(path, NF90_NOWRITE, ncid)
    if (status /= NF90_NOERR) return

    !! Get the total number of elements.
    status = nf90_inq_dimid(ncid, 'num_elem', dimid)
    if (status /= NF90_NOERR) return
    status = nf90_inquire_dimension(ncid, dimid, len=num_elem)
    if (status /= NF90_NOERR) return

    !! Read the number of side sets, if any.
    status = nf90_inq_dimid(ncid, 'num_side_sets', dimid)
    if (status /= NF90_NOERR) then  ! we conclude there are no side sets.
      status = NF90_NOERR
      num_sset = 0
    else
      status = nf90_inquire_dimension(ncid, dimid, len=num_sset)
      if (status /= NF90_NOERR) return
    end if

    if (num_sset == 0) return ! leaving array unallocated
    allocate(mesh_face_set_tot(num_sset,6,num_elem))

    !! Get the number of element blocks.
    status = nf90_inq_dimid(ncid, 'num_el_blk', dimid)
    if (status /= NF90_NOERR) return
    status = nf90_inquire_dimension(ncid, dimid, len=num_eblk)
    if (status /= NF90_NOERR) return

    !! For each block get the number of elements and the element "type".
    allocate(eblk_size(num_eblk), eblk_type(num_eblk))
    do n = 1, num_eblk
      status = nf90_inq_dimid(ncid, 'num_el_in_blk'//i_to_c(n), dimid)
      if (status /= NF90_NOERR) return
      status = nf90_inquire_dimension(ncid, dimid, len=eblk_size(n))
      if (status /= NF90_NOERR) return
      status = nf90_inq_dimid(ncid, 'num_nod_per_el'//i_to_c(n), dimid)
      if (status /= NF90_NOERR) return
      status = nf90_inquire_dimension(ncid, dimid, len=eblk_type(n))
      if (status /= NF90_NOERR) return
    end do

    !! Read the side set ID array.
    allocate(sset_id(num_sset))
    status = nf90_inq_varid(ncid, 'ss_prop1', varid)
    if (status /= NF90_NOERR) return
    status = nf90_get_var(ncid, varid, sset_id)
    if (status /= NF90_NOERR) return

    mesh_face_set_tot = 0

    do n = 1, num_sset
      !! Read the number of sides in this side set.
      status = nf90_inq_dimid(ncid, 'num_side_ss'//i_to_c(n), dimid)
      if (status /= NF90_NOERR) return
      status = nf90_inquire_dimension(ncid, dimid, len=num_side)
      if (status /= NF90_NOERR) return

      !! Read the element/side arrays for this side set.
      allocate(elem(num_side), side(num_side))
      status = nf90_inq_varid(ncid, 'elem_ss'//i_to_c(n), varid)
      if (status /= NF90_NOERR) return
      status = nf90_get_var(ncid, varid, elem)
      if (status /= NF90_NOERR) return
      status = nf90_inq_varid(ncid, 'side_ss'//i_to_c(n), varid)
      if (status /= NF90_NOERR) return
      status = nf90_get_var(ncid, varid, side)
      if (status /= NF90_NOERR) return

      !! Unpack the info into the Truchas array, translating the side indices.
      do j = 1, num_side
        tside = truchas_side(elem(j),side(j))
        if (tside == 0) then
          status = 2
          return
        end if
        mesh_face_set_tot(n,tside,elem(j)) = sset_id(n)
      end do
      deallocate(elem, side)
    end do

    deallocate(sset_id, eblk_type, eblk_size)

    !! Close the exodus mesh file.
    status = nf90_close(ncid)

  contains

    integer function truchas_side (j, k)

      integer, intent(in) :: j, k

      integer :: b, l

      !! Locate element J as element L of block B.
      b = 1
      l = j
      do while (l > eblk_size(b))
        l = l - eblk_size(b)
        b = b + 1
      end do

      !! Map the ExodusII side index to the corresponding Truchas index.
      select case (eblk_type(b))
      case (TET)
        truchas_side = TET_SIDE_MAP(k)
      case (WEDGE)
        truchas_side = WEDGE_SIDE_MAP(k)
      case (HEX)
        truchas_side = HEX_SIDE_MAP(k)
      case default
        truchas_side = 0
      end select

    end function truchas_side

  end subroutine read_exodus_side_sets

end module exodus_truchas_hack
