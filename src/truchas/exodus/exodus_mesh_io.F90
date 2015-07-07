!!
!! EXODUS_MESH_IO
!!
!! This module provides procedures to read/write an EXODUS_MESH object
!! from/to an ExodusII mesh file.
!!
!! Neil N. Carlson <nnc@lanl.gov> 5 Jul 2004
!! Last revised February 2015
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL READ_EXODUS_MESH (PATH, MESH [,STAT [,ERRMSG]])
!!
!!    This subroutine initializes the EXODUS_MESH object MESH with data read
!!    from an ExodusII file at the given PATH.  Trailing blanks in PATH are
!!    significant, so the caller should trim them from the actual argument if
!!    needed. MESH is intent(out) and is finalized on entry to the subroutine.
!!    If the optional integer argument STAT is present, it is assigned a
!!    nonzero value if an error occurs, and the deferred-length allocatable
!!    character string ERRMSG, if present, is assigned an explanatory message.
!!    If STAT is not present and an error occurs, the error message is written
!!    to the preconnected error unit and the program is stopped.
!!
!!  CALL WRITE_EXODUS_MESH (PATH, MESH, CREATOR, VERSION [,STAT [,ERRMSG]])
!!
!!    This subroutine writes the mesh from the EXODUS_MESH object MESH to a
!!    new ExodusII file at the given PATH.  Trailing blanks in PATH are
!!    significant, so the caller should trim them from the actual argument if
!!    needed.  MESH is intent(in) and must be well-defined.  The file PATH
!!    must not exist.  CREATOR and VERSION are intent(in) character strings
!!    that identify the program creating the file and its version.  These
!!    are written to the QA records of the ExodusII file.  If the optional
!!    integer argument STAT is present, it is assigned a nonzero value if an
!!    error occurs, and the deferred-length allocatable character string
!!    ERRMSG, if present, is assigned an explanatory message.  If STAT is not
!!    present and an error occurs, the error message is written to the
!!    preconnected error unit and the program is stopped.
!!
!! NOTES
!!
!! This complete rewrite of the original version, which interfaced with the
!! file via direct NetCDF calls, now uses the official ExodusII library to
!! read and write from the file.
!!
!! Not all the possible mesh data present in an ExodusII file is read; only
!! that which is needed to fill out the EXODUS_MESH data structure.
!!
!! The read procedure is built on an EXODUS_FILE object and its methods are
!! used to get data from the file; the ExodusII library is used indirectly.
!!
!! The write procedure, in contrast, uses direct ExodusII library calls to
!! create the file and put the mesh data into it, as the EXODUS_FILE object
!! provides only methods for reading data from a file.
!!

module exodus_mesh_io

  use exodus_mesh_type
  use exodus_c_binding
  use string_utilities, only: i_to_c
  use,intrinsic :: iso_c_binding, only: c_int, c_ptr, c_null_ptr, c_loc, c_char, c_null_char
  implicit none
  private

  public :: read_exodus_mesh, write_exodus_mesh

contains

  subroutine read_exodus_mesh (path, mesh, stat, errmsg)

    use exodus_file_type

    character(*), intent(in)  :: path
    class(exodus_mesh), intent(out) :: mesh
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg

    integer :: ierr, n
    type(exodus_file) :: file
    character(:), allocatable :: errstr

    !! Open the exodus mesh file.
    call file%open (trim(path), ierr, errstr)
    if (ierr /= 0) then
      call error_handler ('READ_EXODUS_MESH', 'error opening file: '//errstr, stat, errmsg)
      return
    end if

    mesh%title    = file%title
    mesh%num_dim  = file%num_dim
    mesh%num_node = file%num_node
    mesh%num_elem = file%num_elem
    mesh%num_eblk = file%num_elem_blk
    mesh%num_nset = file%num_node_sets
    mesh%num_sset = file%num_side_sets

    !! Read the node coordinate data.
    allocate(mesh%coord(mesh%num_dim,mesh%num_node))
#ifdef NAG_COMPILER
    !! When compiled with -nan, some underlying NetCDF call below copies this array,
    !! triggering an arithmetic exception on some x86 platforms only.  Not really
    !! supposed to use -nan on just part of the program for this very reason.  NNC 6/24/2005.
    mesh%coord = 0.0
#endif
    call file%get_coord (mesh%coord, ierr, errstr)
    if (ierr /= 0) then
      call error_handler ('READ_EXODUS_MESH', &
          'error reading coordinate data: '//errstr, stat, errmsg)
      return
    end if

    !! Read the element block data.
    allocate(mesh%eblk(mesh%num_eblk))
    mesh%eblk%ID = file%elem_blk%ID
    mesh%eblk%num_elem = file%elem_blk%num_elem
    mesh%eblk%num_nodes_per_elem = file%elem_blk%num_nodes_per_elem
    do n = 1, mesh%num_eblk
      mesh%eblk(n)%elem_type = file%elem_blk(n)%elem_type
      allocate(mesh%eblk(n)%connect(mesh%eblk(n)%num_nodes_per_elem,mesh%eblk(n)%num_elem))
      call file%get_connect (mesh%eblk(n)%ID, mesh%eblk(n)%connect, ierr, errstr)
      if (ierr /= 0) then
        call error_handler ('READ_EXODUS_MESH', &
            'error reading connectivity data: '//errstr, stat, errmsg)
        return
      end if
    end do

    !! Read the node set data.
    allocate(mesh%nset(mesh%num_nset))
    if (mesh%num_nset > 0) then
      mesh%nset%ID = file%node_set%ID
      mesh%nset%num_node = file%node_set%num_node_in_set
      do n = 1, mesh%num_nset
        allocate(mesh%nset(n)%node(mesh%nset(n)%num_node))
        call file%get_node_set (mesh%nset(n)%ID, mesh%nset(n)%node, ierr, errstr)
        if (ierr /= 0) then
          call error_handler ('READ_EXODUS_MESH', &
              'error reading node set data: '//errstr, stat, errmsg)
          return
        end if
      end do
    end if

    !! Read the side set data.
    allocate(mesh%sset(mesh%num_sset))
    if (mesh%num_sset > 0) then
      mesh%sset%ID = file%side_set%ID
      mesh%sset%num_side = file%side_set%num_side_in_set
      do n = 1, mesh%num_sset
        allocate(mesh%sset(n)%elem(mesh%sset(n)%num_side))
        allocate(mesh%sset(n)%face(mesh%sset(n)%num_side))
        call file%get_side_set (mesh%sset(n)%ID, mesh%sset(n)%elem, mesh%sset(n)%face, ierr, errstr)
        if (ierr /= 0) then
          call error_handler ('READ_EXODUS_MESH', &
              'error reading side set data: '//errstr, stat, errmsg)
          return
        end if
      end do
    end if

    call file%close

    if (present(stat)) stat = 0

  end subroutine read_exodus_mesh


  subroutine write_exodus_mesh (path, mesh, creator, version, stat, errmsg)

    character(*), intent(in) :: path
    class(exodus_mesh), intent(in) :: mesh
    character(*), intent(in)  :: creator, version
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg

    integer :: ierr
    integer(c_int) :: exoid, comp_ws, io_ws
    character(:), allocatable :: string

    comp_ws = 8 ! will pass 8-byte floating point variables to exodus
    io_ws = 8   ! and have exodus store it as 8-byte floating point data
    exoid = ex_create (trim(path)//c_null_char, EX_CLOBBER, comp_ws, io_ws)
    if (exoid < 0) then
      call error_handler ('WRITE_EXODUS_MESH', &
          'error creating exodus file; ex_create returned '//i_to_c(exoid), stat, errmsg)
      return
    end if

    ierr = ex_put_init(exoid, mesh%title//c_null_char, mesh%num_dim, mesh%num_node, &
                       mesh%num_elem, mesh%num_eblk, mesh%num_nset, mesh%num_sset)
    if (ierr /= 0) then
      call error_handler ('WRITE_EXODUS_MESH', &
          'error writing initialization params; ex_put_init returned '//i_to_c(ierr), stat, errmsg)
      return
    end if

    call put_coord (exoid, mesh, ierr, string)
    if (ierr /= 0) then
      call error_handler ('WRITE_EXODUS_MESH', string, stat, errmsg)
      return
    end if

    call put_elem_blocks (exoid, mesh, ierr, string)
    if (ierr /= 0) then
      call error_handler ('WRITE_EXODUS_MESH', string, stat, errmsg)
      return
    end if

    call put_node_sets (exoid, mesh, ierr, string)
    if (ierr /= 0) then
      call error_handler ('WRITE_EXODUS_MESH', string, stat, errmsg)
      return
    end if

    call put_side_sets (exoid, mesh, ierr, string)
    if (ierr /= 0) then
      call error_handler ('WRITE_EXODUS_MESH', string, stat, errmsg)
      return
    end if

    call put_qa_record (exoid, creator, version, ierr, string)
    if (ierr /= 0) then
      call error_handler ('WRITE_EXODUS_MESH', string, stat, errmsg)
      return
    end if

    ierr = ex_close(exoid)
    if (ierr /= 0) then
      call error_handler ('WRITE_EXODUS_MESH', &
          'error closing file; ex_close returned ' // i_to_c(ierr), stat, errmsg)
      return
    end if

    if (present(stat)) stat = 0

  end subroutine write_exodus_mesh

  !! Auxiliary error handler.  Defines the optional STAT and ERRMSG arguments
  !! if present; otherwise it writes the provided error message to the error
  !! unit and calls the system exit.

  subroutine error_handler (proc, string, stat, errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
#ifdef NAG_COMPILER
    use,intrinsic :: f90_unix, only: exit
#endif
    character(*), intent(in) :: proc, string
    integer, intent(out), optional :: stat
    character(:), allocatable, optional :: errmsg
    if (present(stat)) then
      stat = 1
      if (present(errmsg)) errmsg = proc // ': ' // string
    else
      write(error_unit,'(3a)') proc, ': ', string
      call exit (1)
    end if
  end subroutine error_handler

  !! Auxiliary procedures used by WRITE_EXODUS_MESH follow !!!!!!!!!!!!!!!!!!!!!

  subroutine put_coord (exoid, mesh, stat, errmsg)

    integer(c_int), intent(in)  :: exoid
    type(exodus_mesh), intent(in) :: mesh
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(kind(mesh%coord)), target :: array(mesh%num_node)

    !! Write the x coordinates.
    array = mesh%coord(1,:)
    stat = ex_put_coord (exoid, c_loc(array), c_null_ptr, c_null_ptr)
    if (stat /= 0) then
      errmsg = 'error writing x coord data; ex_put_coord returned ' // i_to_c(stat)
      return
    end if
    if (mesh%num_dim == 1) return

    !! Write the y coordinates.
    array = mesh%coord(2,:)
    stat = ex_put_coord (exoid, c_null_ptr, c_loc(array), c_null_ptr)
    if (stat /= 0) then
      errmsg = 'error writing y coord data; ex_put_coord returned ' // i_to_c(stat)
      return
    end if
    if (mesh%num_dim == 2) return

    !! Write the z coordinates.
    array = mesh%coord(3,:)
    stat = ex_put_coord (exoid, c_null_ptr, c_null_ptr, c_loc(array))
    if (stat /= 0) then
      errmsg = 'error writing z coord data; ex_put_coord returned ' // i_to_c(stat)
      return
    end if

  end subroutine put_coord


  subroutine put_elem_blocks (exoid, mesh, stat, errmsg)

    integer(c_int), intent(in) :: exoid
    type(exodus_mesh), intent(in) :: mesh
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    target :: mesh  ! for the C_LOC call

    integer :: n

    !! Write the element block parameters.
    do n = 1, mesh%num_eblk
      stat = ex_put_elem_block(exoid, mesh%eblk(n)%ID, &
                               trim(mesh%eblk(n)%elem_type)//c_null_char, &
                               mesh%eblk(n)%num_elem, size(mesh%eblk(n)%connect,1), 0)
      if (stat /= 0) then
        errmsg = 'error writing element block params for ID=' // i_to_c(mesh%eblk(n)%ID) // &
                 '; ex_put_elem_block returned ' // i_to_c(stat)
        return
      end if
    end do

    !! Write the element block connectivity.
    do n = 1, mesh%num_eblk
      stat = ex_put_elem_conn(exoid, mesh%eblk(n)%ID, c_loc(mesh%eblk(n)%connect))
      if (stat /= 0) then
        errmsg = 'error writing element block connectivity for ID=' // i_to_c(mesh%eblk(n)%ID) // &
                 '; ex_put_elem_block returned ' // i_to_c(stat)
        return
      end if
    end do

  end subroutine put_elem_blocks


  subroutine put_node_sets (exoid, mesh, stat, errmsg)

    integer(c_int), intent(in) :: exoid
    type(exodus_mesh), intent(in) :: mesh
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    target :: mesh  ! for the C_LOC call

    integer :: n

    !! Write the node set parameters.
    do n = 1, mesh%num_nset
      stat = ex_put_node_set_param(exoid, mesh%nset(n)%ID, mesh%nset(n)%num_node, 0)
      if (stat /= 0) then
        errmsg = 'error writing node set params for ID=' // i_to_c(mesh%nset(n)%ID) // &
                 ' ex_put_node_set_param returned ' // i_to_c(stat)
        return
      end if
    end do

    !! Write the nodes sets.
    do n = 1, mesh%num_nset
      stat = ex_put_node_set(exoid, mesh%nset(n)%ID, c_loc(mesh%nset(n)%node))
      if (stat /= 0) then
        errmsg = 'error writing node set data for ID=' // i_to_c(mesh%nset(n)%ID) // &
                 ' ex_put_node_set returned ' // i_to_c(stat)
        return
      end if
    end do

    stat = 0  ! must define; there need not be any node sets

  end subroutine put_node_sets


  subroutine put_side_sets (exoid, mesh, stat, errmsg)

    integer(c_int), intent(in) :: exoid
    type(exodus_mesh), intent(in) :: mesh
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    target :: mesh  ! for the C_LOC call

    integer :: n

    !! Write the side set parameters.
    do n = 1, mesh%num_sset
      stat = ex_put_side_set_param(exoid, mesh%sset(n)%ID, mesh%sset(n)%num_side, 0)
      if (stat /= 0) then
        errmsg = 'error writing side set params for ID=' // i_to_c(mesh%sset(n)%ID) // &
                 ' ex_put_side_set_param returned ' // i_to_c(stat)
        return
      end if
    end do

    !! Write the side sets.
    do n = 1, mesh%num_sset
      stat = ex_put_side_set(exoid, mesh%sset(n)%ID, c_loc(mesh%sset(n)%elem), c_loc(mesh%sset(n)%face))
      if (stat /= 0) then
        errmsg = 'error writing side set data for ID=' // i_to_c(mesh%sset(n)%ID) // &
                 ' ex_put_side_set returned ' // i_to_c(stat)
        return
      end if
    end do

    stat = 0  ! must define; there need not be any side sets

  end subroutine put_side_sets


  subroutine put_qa_record (exoid, creator, version, stat, errmsg)

    integer, intent(in) :: exoid
    character(*), intent(in) :: creator, version
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j
    character(MAX_STR_LENGTH+1,kind=c_char), target :: qa(4)
    character(8)  :: date
    character(10) :: time
    type(c_ptr) :: qa_records(4,1)

    qa(1) = creator(:min(len_trim(creator),MAX_STR_LENGTH)) // c_null_char
    qa(2) = version(:min(len_trim(creator),MAX_STR_LENGTH)) // c_null_char

    call date_and_time (date, time)

    qa(3) = date(5:6) // '/' // date(7:8) // '/' // date(1:4) // c_null_char
    qa(4) = time(1:2) // ':' // time(3:4) // ':' // time(5:6) // c_null_char

    do j = 1, 4
      call aux (qa(j), qa_records(j,1))
    end do

    stat = ex_put_qa(exoid, 1, qa_records)
    if (stat /= 0) then
      errmsg = 'error writing QA records; ex_put_qa returned ' // i_to_c(stat)
      return
    end if

  contains

    subroutine aux (string, c_loc_string)
      character(kind=c_char), intent(in), target :: string(*)
      type(c_ptr) :: c_loc_string
      c_loc_string = c_loc(string)
    end subroutine

  end subroutine put_qa_record

end module exodus_mesh_io
