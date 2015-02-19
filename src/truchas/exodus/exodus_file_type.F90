!!
!! EXODUS_FILE_TYPE
!!
!! A high-level interface to reading data from an Exodus II mesh file.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2014; updated February 2015.
!!
!! PROGRAMMING INTERFACE
!!
!!  This module defines a minimal high-level interface for reading mesh data
!!  from an Exodus II file.  The derived type EXODUS_FILE holds basic mesh
!!  metadata as public (read-only) components that are populated when the file
!!  is opened. The data components names correspond to variable names described
!!  in the Exodus II document.  In addition there are type bound procedures for
!!  opening and closing the file, and for reading the actual mesh data.
!!
!!  In the following functions, if the optional integer argument STAT is
!!  present, it is assigned the value 0 if no error was encountered; otherwise
!!  it is assigned a non-zero value.  In the latter case, the allocatable
!!  deferred-length character string ERRMSG, if present, is assigned an
!!  explanatory message.  If STAT is not present and an error is encountered,
!!  the error message is written to the preconnected error unit and the
!!  program is terminated.
!!
!!  OPEN (PATH [,STAT [,ERRMSG]]) opens the Exodus II mesh file located at the
!!    specified PATH, and assigns values to all the data components of the
!!    object.
!!
!!  CLOSE () closes the Exodus II mesh file associated with the object and
!!    deallocates any storage associated with the object, returning it to its
!!    default initialization state.  It is not necessary to call this; this
!!    happens automatically when the object goes out of scope.
!!
!!  GET_COORD (COORD [,STAT [,ERRMSG]]) returns the node coordinates in the
!!    rank-2 array COORD.  The first dimension corresponds to the coordinate,
!!    and the second to node index.
!!
!!  GET_CONNECT (ELEM_BLK_ID, CONNECT [,STAT [,ERRMSG]]) returns the element
!!    connectivity data for the specified element block ID in the array CONNECT.
!!
!!  GET_NODE_SET (NODE_SET_ID, NODE_LIST [,STAT [,ERRMSG]]) returns the list of
!!    node numbers for the nodes belonging to the specified node set.
!!
!!  GET_SIDE_SET (SIDE_SET_ID, ELEM_LIST, SIDE_LIST [,STAT [,ERRMSG]]) returns
!!    the list of element numbers and corresponding local element side indices
!!    for the sides belonging to the specified side set.
!!

module exodus_file_type

  use exodus_c_binding
  use,intrinsic :: iso_c_binding
  implicit none
  private

  public :: ex_get_side_set

  type :: exo_elem_blk_info
    integer(c_int) :: ID
    integer(c_int) :: num_elem
    integer(c_int) :: num_nodes_per_elem
    character(:), allocatable :: elem_type
  end type exo_elem_blk_info

  type :: exo_node_set_info
    integer(c_int) :: ID
    integer(c_int) :: num_node_in_set
  end type exo_node_set_info

  type :: exo_side_set_info
    integer(c_int) :: ID
    integer(c_int) :: num_side_in_set
  end type exo_side_set_info

  type, public :: exodus_file
    integer(c_int) :: exoid = -1
    integer(c_int) :: num_dim = 0
    integer(c_int) :: num_node = 0
    integer(c_int) :: num_elem = 0
    integer(c_int) :: num_elem_blk = 0
    integer(c_int) :: num_node_sets = 0
    integer(c_int) :: num_side_sets = 0
    character(:), allocatable :: title
    type(exo_elem_blk_info), allocatable :: elem_blk(:)
    type(exo_node_set_info), allocatable :: node_set(:)
    type(exo_side_set_info), allocatable :: side_set(:)
  contains
    procedure :: open  => exodus_file_open
    procedure :: close => exodus_file_close
    procedure :: get_coord => exodus_file_get_coord
    procedure :: get_connect => exodus_file_get_connect
    procedure :: get_node_set => exodus_file_get_node_set
    procedure :: get_side_set => exodus_file_get_side_set
    procedure :: dump
    final :: exodus_file_final
  end type exodus_file

contains

 !! Final subroutine for EXODUS_FILE objects. Closes the Exodus file if it is
 !! open.  An open file is indicated by a nonnegative EXOID.  A valid EXOID
 !! is probably positive, though a strict reading of the spec allows for 0.
 !! Exodus library errors are ignored.

  subroutine exodus_file_final (this)
    type(exodus_file), intent(inout) :: this
    integer :: error
    if (this%exoid >= 0) error = ex_close(this%exoid)
  end subroutine exodus_file_final

 !! Closes the Exodus file and default initializes the components of the
 !! object.  The subroutine appears to do nothing, however the passed object
 !! is INTENT(OUT) and this causes the object to be finalized and its
 !! components default initialized on entry to the subroutine.

  subroutine exodus_file_close (this)
    class(exodus_file), intent(out) :: this
  end subroutine exodus_file_close

 !! Opens the specified Exodus file and reads the relevant mesh meta data
 !! (sizes, element block IDs, side set IDs, etc.), storing it in public
 !! data components of the object.

  subroutine exodus_file_open (this, path, stat, errmsg)

    class(exodus_file), intent(out) :: this
    character(*), intent(in) :: path
    integer, intent(out), optional :: stat
    character(:), allocatable, optional :: errmsg

    integer(c_int) :: comp_ws, io_ws, error, idummy
    real(c_float) :: version
    character(MAX_LINE_LENGTH,kind=c_char) :: title
    character(MAX_STR_LENGTH,kind=c_char) :: elem_type
    integer :: n

    if (present(stat)) stat = 0

    !! Open the Exodus file.
    comp_ws = 8 ! code will use 8-byte reals
    io_ws = 0   ! let ex_open tell us the floating point data size in the file
    this%exoid = ex_open(path//c_null_char, EX_READ, comp_ws, io_ws, version)
    if (this%exoid < 0) then
      call exo_error_handler ('EXODUS_FILE:OPEN', &
          'ex_open returned an error', this%exoid, stat, errmsg)
      return
    end if
    !TODO: should probably emit a warning if comp_ws /= io_ws

    !! Get the mesh sizes.
    error = ex_get_init(this%exoid, title, this%num_dim, this%num_node, this%num_elem, &
                        this%num_elem_blk, this%num_node_sets, this%num_side_sets)
    if (error /= 0) then
      call exo_error_handler ('EXODUS_FILE:OPEN', &
          'ex_get_init returned an error', error, stat, errmsg)
      return
    end if
    this%title = ctrim(title,MAX_LINE_LENGTH)

    !! Read the element block info.
    allocate(this%elem_blk(this%num_elem_blk))
    error = ex_get_elem_blk_ids(this%exoid, this%elem_blk%ID)
    if (error /= 0) then
      call exo_error_handler ('EXODUS_FILE:OPEN', &
          'ex_get_elem_blk_ids returned an error', error, stat, errmsg)
      return
    end if
    do n = 1, this%num_elem_blk
      error = ex_get_elem_block (this%exoid, this%elem_blk(n)%ID, elem_type, &
          this%elem_blk(n)%num_elem, this%elem_blk(n)%num_nodes_per_elem, idummy)
      if (error /= 0) then
        call exo_error_handler ('EXODUS_FILE:OPEN', &
            'ex_get_elem_block returned an error', error, stat, errmsg)
        return
      end if
      this%elem_blk(n)%elem_type = ctrim(elem_type,MAX_STR_LENGTH)
    end do

    !! Read the node set info.
    allocate(this%node_set(this%num_node_sets))
    if (this%num_node_sets > 0) then
      error = ex_get_node_set_ids(this%exoid, this%node_set%ID)
      if (error /= 0) then
        call exo_error_handler ('EXODUS_FILE:OPEN', &
            'ex_get_node_set_ids returned an error', error, stat, errmsg)
        return
      end if
      do n = 1, this%num_node_sets
        error = ex_get_node_set_param(this%exoid, this%node_set(n)%ID, &
            this%node_set(n)%num_node_in_set, idummy)
        if (error /= 0) then
          call exo_error_handler ('EXODUS_FILE:OPEN', &
              'ex_get_node_set_param returned an error', error, stat, errmsg)
          return
        end if
      end do
    end if

    !! Read the side set info.
    allocate(this%side_set(this%num_side_sets))
    if (this%num_side_sets > 0) then
      error = ex_get_side_set_ids(this%exoid, this%side_set%ID)
      if (error /= 0) then
        call exo_error_handler ('EXODUS_FILE:OPEN', &
            'ex_get_side_set_ids returned an error', error, stat, errmsg)
        return
      end if
      do n = 1, this%num_side_sets
        error = ex_get_side_set_param(this%exoid, this%side_set(n)%ID, &
            this%side_set(n)%num_side_in_set, idummy)
        if (error /= 0) then
          call exo_error_handler ('EXODUS_FILE:OPEN', &
              'ex_get_side_set_param returned an error', error, stat, errmsg)
          return
        end if
      end do
    end if

  end subroutine exodus_file_open

 !! Reads the connectivity array for the specified element block.  In the
 !! C interface, the array as passed as a void pointer, which will accept
 !! either a 4-byte or 8-byte integer array (but only one will be valid for
 !! a particular file).  Here we assume a 4-byte (INTEGER(C_INT)) -- the
 !! 8-byte interface is not documented.  This is a likely source of breakage.

  subroutine exodus_file_get_connect (this, elem_blk_id, connect, stat, errmsg)

    use string_utilities, only: i_to_c

    class(exodus_file), intent(in) :: this
    integer(c_int), intent(in)  :: elem_blk_id
    integer(c_int), intent(out) :: connect(:,:)
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg

    integer :: n, ierr

    if (present(stat)) stat = 0

    !! Check the shape of the CONNECT array.
    do n = 1, size(this%elem_blk)
      if (this%elem_blk(n)%ID == elem_blk_id) exit
    end do
    if (n > size(this%elem_blk)) then
      call error_handler ('EXODUS_FILE:GET_CONNECTIVITY', &
          'no such element block ID: ' // i_to_c(elem_blk_id), stat, errmsg)
      return
    end if
    if (size(connect,dim=1) /= this%elem_blk(n)%num_nodes_per_elem .or. &
        size(connect,dim=2) /= this%elem_blk(n)%num_elem) then
      call error_handler ('EXODUS_FILE:GET_CONNECTIVITY', &
          'invalid connect array shape', stat, errmsg)
    end if

    !! Read the connectivity data for this element block.
    call get_connect_aux (this%exoid, elem_blk_id, connect, ierr)
    if (ierr /= 0) then
      call exo_error_handler ('EXODUS_FILE:GET_CONNECT', &
          'ex_get_connect returned an error', ierr, stat, errmsg)
      return
    end if

  contains

    !! The library function assumes a contiguous CONNECT array.  The point of
    !! this auxillary function is to force, if necessary, a copy-out to CONNECT
    !! from a contiguous temporary passed to the library function.
    subroutine get_connect_aux (exoid, elem_blk_id, connect, ierr)
      integer(c_int), intent(in) :: exoid, elem_blk_id
      integer(c_int), intent(out), target :: connect(*)
      integer, intent(out) :: ierr
      ierr = ex_get_elem_conn (exoid, elem_blk_id, c_loc(connect))
    end subroutine

  end subroutine exodus_file_get_connect

 !! Reads the node coordinates.  Handles the juggling from the C interface
 !! that works with each coordinate as a separate contiguous array, and the
 !! desired interface where the coordinates are merged into a single array.

  subroutine exodus_file_get_coord (this, coord, stat, errmsg)

    use,intrinsic :: iso_fortran_env, only: real64  ! 8-byte reals (matches open params)

    class(exodus_file), intent(in) :: this
    real(real64), intent(out) :: coord(:,:)
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg

    real(real64), target :: array(this%num_node)
    integer(c_int) :: error

    if (present(stat)) stat = 0

    !! Check the shape of the COORD array.
    if (size(coord,1) /= this%num_dim .or. size(coord,2) /= this%num_node) then
      call error_handler ('EXODUS_FILE:GET_COORD', &
          'invalid coord array shape', stat, errmsg)
      return
    end if

    !! Read the x coordinates.
    error = ex_get_coord(this%exoid, c_loc(array), c_null_ptr, c_null_ptr)
    if (error /= 0) then
      call exo_error_handler ('EXODUS_FILE:GET_COORD', &
          ' ex_get_coord returned an error for X', error, stat, errmsg)
      return
    end if
    coord(1,:) = array
    if (this%num_dim == 1) return

    !! Read the y coordinates.
    error = ex_get_coord(this%exoid, c_null_ptr, c_loc(array), c_null_ptr)
    if (error /= 0) then
      call exo_error_handler ('EXODUS_FILE:GET_COORD', &
          ' ex_get_coord returned an error for Y', error, stat, errmsg)
      return
    end if
    coord(2,:) = array
    if (this%num_dim == 2) return

    !! Read the z coordinates.
    error = ex_get_coord(this%exoid, c_null_ptr, c_null_ptr, c_loc(array))
    if (error /= 0) then
      call exo_error_handler ('EXODUS_FILE:GET_COORD', &
          ' ex_get_coord returned an error for Z', error, stat, errmsg)
      return
    end if
    coord(3,:) = array

  end subroutine exodus_file_get_coord

 !! Reads the node set array for the specified node set.  In the
 !! C interface, the array is passed as a void pointer, which will accept
 !! either a 4-byte or 8-byte integer array (but only one will be valid for
 !! a particular file).  Here we assume a 4-byte (INTEGER(C_INT)) -- the
 !! 8-byte interface is not documented.  This is a likely source of breakage.

  subroutine exodus_file_get_node_set (this, node_set_id, node_list, stat, errmsg)

    use string_utilities, only: i_to_c

    class(exodus_file), intent(in) :: this
    integer(c_int), intent(in)  :: node_set_id
    integer(c_int), intent(out) :: node_list(:)
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg

    integer :: n, ierr

    if (present(stat)) stat = 0

    !! Check the shape of the NODE_LIST array.
    do n = 1, size(this%node_set)
      if (this%node_set(n)%ID == node_set_id) exit
    end do
    if (n > size(this%node_set)) then
      call error_handler ('EXODUS_FILE:GET_NODE_SET', &
          'no such node set ID: ' // i_to_c(node_set_id), stat, errmsg)
      return
    end if
    if (size(node_list) /= this%node_set(n)%num_node_in_set) then
      call error_handler ('EXODUS_FILE:GET_NODE_SET', &
          'invalid node_list array size', stat, errmsg)
    end if

    !! Read the data for this node set.
    call get_node_set_aux (this%exoid, node_set_id, node_list, ierr)
    if (ierr /= 0) then
      call exo_error_handler ('EXODUS_FILE:GET_NODE_SET', &
          'ex_get_node_set returned an error', ierr, stat, errmsg)
      return
    end if

  contains

    !! The library function assumes contiguous arrays.  The point of this
    !! auxillary function is to force, if necessary, a copy-out to the output
    !! arrays from contiguous temporaries passed to the library function.
    subroutine get_node_set_aux (exoid, node_set_id, node_list, ierr)
      integer(c_int), intent(in) :: exoid, node_set_id
      integer(c_int), intent(out), target :: node_list(*)
      integer, intent(out) :: ierr
      ierr = ex_get_node_set (exoid, node_set_id, c_loc(node_list))
    end subroutine

  end subroutine exodus_file_get_node_set

 !! Reads the side set arrays for the specified side set.  In the
 !! C interface, the arrays are passed as void pointers, which will accept
 !! either a 4-byte or 8-byte integer array (but only one will be valid for
 !! a particular file).  Here we assume a 4-byte (INTEGER(C_INT)) -- the
 !! 8-byte interface is not documented.  This is a likely source of breakage.

  subroutine exodus_file_get_side_set (this, side_set_id, elem_list, side_list, stat, errmsg)

    use string_utilities, only: i_to_c

    class(exodus_file), intent(in) :: this
    integer(c_int), intent(in)  :: side_set_id
    integer(c_int), intent(out) :: elem_list(:), side_list(:)
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg

    integer :: n, ierr

    if (present(stat)) stat = 0

    !! Check the shape of the ELEM_LIST and SIDE_LIST arrays.
    do n = 1, size(this%side_set)
      if (this%side_set(n)%ID == side_set_id) exit
    end do
    if (n > size(this%side_set)) then
      call error_handler ('EXODUS_FILE:GET_SIDE_SET', &
          'no such side set ID: ' // i_to_c(side_set_id), stat, errmsg)
      return
    end if
    if (size(elem_list) /= this%side_set(n)%num_side_in_set) then
      call error_handler ('EXODUS_FILE:GET_SIDE_SET', &
          'invalid elem_list array size', stat, errmsg)
    end if
    if (size(side_list) /= size(elem_list)) then
      call error_handler ('EXODUS_FILE:GET_SIDE_SET', &
          'invalid side_list array size', stat, errmsg)
    end if

    !! Read the data for this side set.
    call get_side_set_aux (this%exoid, side_set_id, elem_list, side_list, ierr)
    if (ierr /= 0) then
      call exo_error_handler ('EXODUS_FILE:GET_SIDE_SET', &
          'ex_get_side_set returned an error', ierr, stat, errmsg)
      return
    end if

  contains

    !! The library function assumes contiguous arrays.  The point of this
    !! auxillary function is to force, if necessary, a copy-out to the output
    !! arrays from contiguous temporaries passed to the library function.
    subroutine get_side_set_aux (exoid, side_set_id, elem_list, side_list, ierr)
      integer(c_int), intent(in) :: exoid, side_set_id
      integer(c_int), intent(out), target :: elem_list(*), side_list(*)
      integer, intent(out) :: ierr
      ierr = ex_get_side_set (exoid, side_set_id, c_loc(elem_list), c_loc(side_list))
    end subroutine

  end subroutine exodus_file_get_side_set

!!!! AUXILLARY PROCEDURES FOLLOW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Exodus returns strings in user-allocated storage.  Although the API
  !! specifies a maximum length for these strings, they appear to be null
  !! terminated.  On the Fortran side we need to strip off this null
  !! termination.  This auxillary function does this.

  function ctrim (string, maxlen) result (s)
    character(*,kind=c_char), intent(in) :: string
    integer, intent(in) :: maxlen
    character(:), allocatable :: s
    integer :: n
    do n = 1, maxlen
      if (string(n:n) == c_null_char) exit
    end do
    n = n - 1
    s = string(:n)
  end function ctrim

  !! Auxillary exodus error handler.  Defines the optional STAT and ERRMSG
  !! arguments if present; otherwise it writes an error message to the error
  !! unit and calls the system exit.  The error message is constructed from
  !! PROC (calling procedure name), STRING (error description), and ERROR
  !! (the exodus error code).  It's unfortunate that the exodus library
  !! provides no means for acquiring its message for the error code, only a
  !! method that writes it to stderr.

  subroutine exo_error_handler (proc, string, error, stat, errmsg)
    use string_utilities, only: i_to_c
    use,intrinsic :: iso_fortran_env, only: error_unit
#ifdef NAG_COMPILER
    use,intrinsic :: f90_unix, only: exit
#endif
    character(*), intent(in) :: proc, string
    integer(c_int), intent(in) :: error
    integer, intent(out), optional :: stat
    character(:), allocatable, optional :: errmsg
    if (present(stat)) then
      stat = 1
      if (present(errmsg)) errmsg = string // ' (EXOERR=' // i_to_c(error) // ')'
    else
      write(error_unit,'(4a,i0,a)') proc, ': ', string, ' (EXOERR=', error, ')'
      call exit (1)
    end if
  end subroutine exo_error_handler

  !! Auxillary error handler.  Defines the optional STAT and ERRMSG arguments
  !! if present; otherwise it writes an error message to the error unit and
  !! calls the system exit.  The error message is constructed from PROC
  !! (calling procedure name) and STRING (error description).

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
      if (present(errmsg)) errmsg = string
    else
      write(error_unit,'(3a)') proc, ': ', string
      call exit (1)
    end if
  end subroutine error_handler

  subroutine dump (this)
    class(exodus_file), intent(in) :: this
    integer :: j
    write(*,'(a,i0)') 'exoid=', this%exoid
    write(*,'(a,i0)') 'num_dim=', this%num_dim
    write(*,'(a,i0)') 'num_node=', this%num_node
    write(*,'(a,i0)') 'num_elem=', this%num_elem
    write(*,'(a,i0)') 'num_elem_blk=', this%num_elem_blk
    write(*,'(a,i0)') 'num_node_sets=', this%num_node_sets
    write(*,'(a,i0)') 'num_side_sets=', this%num_side_sets
    write(*,'(2a)')   'title=', this%title
    if (allocated(this%elem_blk)) then
      do j = 1, size(this%elem_blk)
        write(*,'(3(a,i0),2a)') 'elem_blk: ID=', this%elem_blk(j)%ID, &
                   ', num_elem=', this%elem_blk(j)%num_elem, &
         ', num_nodes_per_elem=', this%elem_blk(j)%num_nodes_per_elem, &
                  ', elem_type=', this%elem_blk(j)%elem_type
      end do
    end if
    if (allocated(this%node_set)) then
      do j = 1, size(this%node_set)
        write(*,'(2(a,i0))') 'node_set: ID=', this%node_set(j)%ID, &
                   ', num_node_in_set=', this%node_set(j)%num_node_in_set
      end do
    end if
    if (allocated(this%side_set)) then
      do j = 1, size(this%side_set)
        write(*,'(2(a,i0))') 'side_set: ID=', this%side_set(j)%ID, &
                   ', num_side_in_set=', this%side_set(j)%num_side_in_set
      end do
    end if
  end subroutine dump

end module exodus_file_type
