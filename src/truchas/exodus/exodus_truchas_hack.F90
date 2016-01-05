!!
!! EXODUS_TRUCHAS_HACK
!!
!! This module provides procedures to read certain data from an ExodusII file.
!! These are intended solely to accomodate specific needs in Truchas as it is
!! currently designed.  Please DO NOT use these for any other purpose.
!!
!! Neil N. Carlson <nnc@lanl.gov> 25 Sep 2004
!! Revised February 2015
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
!!  CALL READ_EXODUS_MESH_SIZE (PATH, NUM_NODE, NUM_ELEM [,STAT])
!!
!!    Reads the number of nodes and elements from the ExodusII file PATH.
!!    If an error occurs, a nonzero value is returned in STAT, if present;
!!    otherwise an error message is written to the preconnected error unit
!!    and the program is stopped.
!!
!!  CALL READ_EXODUS_SIDE_SETS (PATH, MESH_FACE_SET_TOT, STATUS)
!!
!!    Reads the side set information from the ExodusII file PATH.  The data
!!    is returned in the supplied rank-3 MESH_FACE_SET_TOT array; its shape
!!    should be number of side sets, by number of faces per cell (6), by number
!!    of elements.  The old Truchas mesh data structure stores tet and wedge
!!    elements as specific forms of degenerate hexes, and it also uses a
!!    different numbering of the sides of a hex.  (See MESH_READ from the
!!    MESH_INPUT_MODULE for details.)  This routine takes care of translating
!!    the ExodusII side info appropriately.  If an error occurs, a nonzero
!!    value is returned in STAT.
!!
!! NOTES
!!
!! The internal implementation of the procedures was completely rewritten
!! to use the official ExodusII library to read the file instead of making
!! direct NetCDF calls.  The procedure interfaces are unchanged.
!!
!! Trailing white space in the path string is significant to the ExodusII
!! library function ex_open that underlies the EXODUS_FILE:OPEN method, which
!! does not modify the path string other than to append a null character to
!! form a proper C string.  The EXODUS_FILE:OPEN caller is thus expected to
!! pass the precise path.  Instead of propagating this requirement up to the
!! caller of these procedures, these procedures trim the received path to
!! remove any trailing spaces, which are thus deemed to be insignificant.
!!

module exodus_truchas_hack

  implicit none
  private

  public :: read_exodus_mesh_size, read_exodus_side_sets

contains

  subroutine read_exodus_mesh_size (path, num_node, num_elem, stat)

    use exodus_file_type
    use,intrinsic :: iso_fortran_env, only: error_unit
#ifdef NAG_COMPILER
    use,intrinsic :: f90_unix, only: exit
#endif

    character(*), intent(in) :: path
    integer, intent(out) :: num_node, num_elem
    integer, intent(out), optional :: stat

    integer :: status
    type(exodus_file) :: file
    character(:), allocatable :: errmsg

    call file%open (trim(path), status, errmsg)
    if (status /= 0) then
      if (present(stat)) then
        stat = status
        return
      else
        write(error_unit,'(2a)') 'READ_EXODUS_MESH_SIZE: error opening file: ', errmsg
        call exit (1)
      end if
    end if

    num_node = file%num_node
    num_elem = file%num_elem

    if (present(stat)) stat = 0

  end subroutine read_exodus_mesh_size


  subroutine read_exodus_side_sets (path, mesh_face_set_tot, status)

    use exodus_file_type

    character(len=*), intent(in) :: path
    integer, pointer :: mesh_face_set_tot(:,:,:)
    integer, intent(out) :: status

    integer :: num_side, n, j, tside
    integer, allocatable :: elem(:), side(:)
    type(exodus_file) :: file

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

    !! Open the exodus mesh file.
    call file%open (trim(path), status)
    if (status /= 0) return

    if (file%num_side_sets == 0) return ! leaving array unallocated
    allocate(mesh_face_set_tot(file%num_side_sets,6,file%num_elem))

    mesh_face_set_tot = 0

    do n = 1, file%num_side_sets
      !! Read the element/side arrays for this side set.
      num_side = file%side_set(n)%num_side_in_set
      allocate(elem(num_side), side(num_side))
      call file%get_side_set (file%side_set(n)%ID, elem, side, status)
      if (status /= 0) return
      !! Unpack the info into the Truchas array, translating the side indices.
      do j = 1, num_side
        tside = truchas_side(elem(j),side(j))
        if (tside == 0) then
          status = 2
          return
        end if
        mesh_face_set_tot(n,tside,elem(j)) = file%side_set(n)%ID
      end do
      deallocate(elem, side)
    end do

  contains

    integer function truchas_side (j, k)

      integer, intent(in) :: j, k

      integer :: b, l

      associate (eblk_size => file%elem_blk%num_elem, &
                 eblk_type => file%elem_blk%num_nodes_per_elem)
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
      end associate

    end function truchas_side

  end subroutine read_exodus_side_sets

end module exodus_truchas_hack
