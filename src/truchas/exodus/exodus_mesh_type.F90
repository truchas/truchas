!!
!! EXODUS_MESH_TYPE
!!
!! This module provides derived data types for encapsulating the mesh data
!! described by an Exodus II file, and several basic procedures that operate
!! on those types.
!!
!! Neil N. Carlson <nnc@lanl.gov> 16 Sep 2004
!! Major revisions February 2015
!!
!! PROGRAMMING INTERFACE
!!
!! The module defines the derived type EXODUS_MESH and subsidiary derived
!! types ELEM_BLK, NODE_SET, and SIDE_SET for encapsulating the mesh data.
!! These were designed to closely emulate the structure of the data as it
!! is described by the ExodusII mesh format [1]. This facilitates the reading
!! and writing of mesh data by ensuring an equivalence between the in-memory
!! and on-disk information content.  The types have the following public
!! data components.
!!
!!  EXODUS_MESH
!!      num_dim - The number of spatial coordinates per node (1, 2, or 3).
!!     num_node - The number of nodes
!!     num_elem - The number of elements
!!     num_eblk - The number of element blocks
!!     num_nset - The number of node sets (possibly 0)
!!     num_sset - The number of side sets (possibly 0)
!!         eblk - An allocatable rank-1 array of type ELEM_BLK holding the
!!                list of elements.  The size of the array equals num_eblk.
!!         nset - An allocatable rank-1 array of type NODE_SET holding the
!!                list of node sets.  The array size equals num_nset.
!!         sset - An allocatable rank-1 array of type SIDE_SET holding the
!!                list of side sets.  The array size equals num_sset.
!!        coord - An allocatable rank-2 real array containing the coordinates
!!                of the nodes.  The real kind is REAL64.  The first dimension
!!                is the coordinate index and has extent num_dim.  The second
!!                dimension is the node index and has extent num_node.
!!        title - A deferred-length allocatable character holding the title
!!                for the data set.
!!
!!  ELEM_BLK
!!           id - An arbitrary, unique, positive integer that identifies the
!!                particular element block.
!!     num_elem - The number of elements in the element block.
!!    elem_type - A deferred-length allocatable character string that
!!                identifies the type of element contained in the element
!!                block.  The interpretation of this value is left to the
!!                client codes that use or produce the mesh data, however
!!                it is recommended to adhere to the conventions described
!!                in the ExodusII manual insofar as possible.
!!      connect - An allocatable rank-2 integer array holding the list of nodes
!!                which define each element in the element block.  The first
!!                dimension is the node index and its extent is the number of
!!                nodes per element, which will depend on the element type of
!!                the block.  The second dimension is the element index and its
!!                extent is num_elem.
!!  NODE_SET
!!           id - An arbitrary, unique, positive integer that identifies the
!!                particular node set.
!!     num_node - The number of nodes in the node set.
!!         node - An allocatable rank-1 integer array holding the list of nodes
!!                in the node set.  Its size equals num_node.
!!
!!  SIDE_SET
!!           id - An arbitrary, unique, positive integer that identifies the
!!                particular side set.
!!     num_side - The number of sides in the node set.
!!         elem - An allocatable rank-1 integer array holding the list of
!!                elements in the side set.  Its size equals num_side.
!!         face - An allocatable rank-1 integer array holding the list of local
!!                face numbers for the elements in the side set.  Its size
!!                equals num_side.
!!         node - An allocatable rank-1 integer array holding the list of nodes
!!                in the node set.  Its size equals num_node.
!!
!!  The labeling of nodes and faces of each element type should adhere to the
!!  conventions described in the ExodusII manual.
!!
!!  Storage associated with objects of these types is automatically freed when
!!  the object is deallocated or otherwise goes out of scope.
!!
!!  The EXODUS_MESH derived type has the following type bound procedures:
!!
!!    SIDE_SET_NODE_LIST(N) returns a pointer to a rank-1 default integer array
!!      containing the side set node list for side set index N.
!!      NB: N is the sequential index of the side set in the %SSET component
!!      of the mesh object, *not* the side set ID! (FIXME)
!!      This function is similar to the ExodusII library procedure
!!      ex_get_side_set_node_list (exgssn).  A null pointer is returned if
!!      anything not understood is encountered.
!!      NB: The only proper use of the function is as the target of a pointer
!!      assignment; if used otherwise (in an expression, e.g.) a memory leak
!!      will occur.  The caller is responsible for deallocating the array.
!!
!!    SIDE_NODE_LIST(ELEM, SIDE [,REVERSE]) returns a pointer to a rank-1
!!      default integer array containing the indices of the nodes on one side
!!      of an element.  The element and local side indices are specified by
!!      ELEM and SIDE.  The returned node list begins with the lowest numbered
!!      node, and by default continues counter-clockwise about the side with
!!      respect to the outward orientation of that side.  If the optional
!!      argument REVERSE is present with value true, the nodes are instead
!!      ordered with respect to the inward orientation of the side.  A null
!!      pointer is returned if the ELEM or SIDE values are out of range, or
!!      if an element of unknown type is encountered.
!!      NB: The only proper use of the function is as the target of a pointer
!!      assignment; if used otherwise (in an expression, e.g.) a memory leak
!!      will occur.  The caller is responsible for deallocating the array.
!!
!!    SIDE_SIZE_LIST(ELEM_TYPE) returns a pointer to a rank-1 default integer
!!      array containing the sizes of (i.e., the number of nodes on) the sides
!!      of an element of the specified type.  Recognized values for ELEM_TYPE
!!      are a subset of the element type strings described by the ExodusII
!!      manual: 'TETRA', 'TETRA4', 'WEDGE', 'WEDGE6', 'HEX', and 'HEX8'.  The
!!      sizes are listed according to the ExodusII face ordering convention.
!!      A null pointer is returned if the given ELEM_TYPE is not recognized.
!!      NB: The target of the pointer result is static storage and consequently
!!      should never be altered, nor should the pointer be deallocated.
!!      Note that this function makes no use whatsoever of the object data.
!!
!!  Each type also has the following type bound procedures that are useful
!!  only for debugging purposes:
!!
!!    DEFINED() returns true if the object is defined; otherwise it returns
!!      false.  Defined means that the components are properly defined and
!!      superficially reasonable.
!!
!!    DUMP(LUN) writes the object data to the logical unit LUN.
!!
!! NOTES
!!
!! The ExodusII format may include additional data that we are not using,
!! including: optional node and element number maps that specify a mapping
!! from the internal numbering to a user-space numbering; an optional element
!! order map that specifies a 'good' order to process the elements, and others.
!!
!! [1] L.A.Schoof and V.R.Yarberry, "Exodus II: A Finite Element Data Model",
!!     Sandia report SAND92-2137.  This can be obtained at
!!     http://endo.sandia.gov/SEACAS/Documentation/exodusII.pdf
!!

module exodus_mesh_type

  use,intrinsic :: iso_fortran_env, only: real64
  implicit none
  private

  type, public :: elem_blk
    integer :: id = 0, num_elem = 0, num_nodes_per_elem = 0
    character(:), allocatable :: elem_type
    integer, allocatable :: connect(:,:)
  contains
    !! Debugging methods
    procedure :: defined => elem_blk_defined
    procedure :: dump => elem_blk_dump
  end type elem_blk

  type, public :: node_set
    integer :: id = 0, num_node = 0
    integer, allocatable :: node(:)
  contains
    !! Debugging methods
    procedure :: defined => node_set_defined
    procedure :: dump => node_set_dump
  end type node_set

  type, public :: side_set
    integer :: id = 0, num_side = 0
    integer, allocatable :: elem(:), face(:)
  contains
    !! Debugging methods
    procedure :: defined => side_set_defined
    procedure :: dump => side_set_dump
  end type side_set

  type, public :: exodus_mesh
    integer :: num_dim  = 0   ! spatial dimension of the mesh
    integer :: num_node = 0   ! number of nodes in the mesh
    integer :: num_elem = 0   ! number of elements in the mesh
    integer :: num_eblk = 0   ! number of element blocks
    integer :: num_nset = 0   ! number of node sets
    integer :: num_sset = 0   ! number of side sets
    type(elem_blk), allocatable :: eblk(:)    ! list of element blocks
    type(node_set), allocatable :: nset(:)    ! list of node sets
    type(side_set), allocatable :: sset(:)    ! list of side sets
    real(real64),   allocatable :: coord(:,:) ! node coordinates
    character(:),   allocatable :: title
  contains
    procedure :: side_node_list
    procedure :: side_set_node_list
    procedure :: get_concat_elem_conn
    procedure, nopass :: side_size_list
    !! Debugging methods
    procedure :: defined => exodus_mesh_defined
    procedure :: dump => exodus_mesh_dump
  end type exodus_mesh

  !! Element specification data (static) as described in the ExodusII manual.
  !! SIDES(XSIDE(K):XSIDE(K+1)-1) is the list of local nodes of side K for a
  !! particular type of element.  SSIZE(K) = XSIDE(K+1)-XSIDE(K) is the size
  !! (number of nodes) of side K of the element.

  integer, target, private :: TETRA4_XSIDE(5), TETRA4_SIDES(12), TETRA4_SSIZE(4)
  data TETRA4_XSIDE/1,4,7,10,13/
  data TETRA4_SIDES/1,2,4, 2,3,4, 1,4,3, 1,3,2/
  data TETRA4_SSIZE/3,3,3,3/

  integer, target, private :: PYRAMID5_XSIDE(6), PYRAMID5_SIDES(16), PYRAMID5_SSIZE(5)
  data PYRAMID5_XSIDE/1,4,7,10,13,17/
  data PYRAMID5_SIDES/1,2,5, 2,3,5, 3,4,5, 1,5,4, 1,4,3,2/
  data PYRAMID5_SSIZE/3,3,3,3,4/

  integer, target, private :: WEDGE6_XSIDE(6), WEDGE6_SIDES(18), WEDGE6_SSIZE(5)
  data WEDGE6_XSIDE/1,5,9,13,16,19/
  data WEDGE6_SIDES/1,2,5,4, 2,3,6,5, 1,4,6,3, 1,3,2, 4,5,6/
  data WEDGE6_SSIZE/4,4,4,3,3/

  integer, target, private :: HEX8_XSIDE(7), HEX8_SIDES(24), HEX8_SSIZE(6)
  data HEX8_XSIDE/1,5,9,13,17,21,25/
  data HEX8_SIDES/1,2,6,5, 2,3,7,6, 3,4,8,7, 1,5,8,4, 1,4,3,2, 5,6,7,8/
  data HEX8_SSIZE/4,4,4,4,4,4/

contains

  function side_node_list (this, elem, side, reverse) result (list)

    class(exodus_mesh), intent(in) :: this
    integer, intent(in) :: elem, side
    logical, intent(in), optional :: reverse
    integer, pointer :: list(:)

    integer :: b, l, n
    integer, pointer :: xside(:), sides(:)

    list => null()

    if (elem < 1 .or. elem > this%num_elem) return

    !! Locate the element as element L of block B.
    b = 1
    l = elem
    do while (l > this%eblk(b)%num_elem)
      l = l - this%eblk(b)%num_elem
      b = b + 1
    end do

    select case (this%eblk(b)%elem_type(1:3))
    case ('TET')
      xside => TETRA4_XSIDE
      sides => TETRA4_SIDES
    case ('PYR')
      xside => PYRAMID5_XSIDE
      sides => PYRAMID5_SIDES
    case ('WED')
      xside => WEDGE6_XSIDE
      sides => WEDGE6_SIDES
    case ('HEX')
      xside => HEX8_XSIDE
      sides => HEX8_SIDES
    case default
      return  ! unknown element type
    end select

    if (side < 1 .or. side > size(xside)-1) return  ! bad side index

    allocate(list(xside(side+1)-xside(side)))
    list = this%eblk(b)%connect(sides(xside(side):xside(side+1)-1),l)

    !! Normalize the list by rotating the smallest value to the initial position.
    list = cshift(list, shift=minloc(list,dim=1)-1)

    !! If specified, reverse the orientation of the side.
    if (present(reverse)) then
      if (reverse) then
        select case (size(list))
        case (3)
          n = list(2)
          list(2) = list(3)
          list(3) = n
        case (4)
          n = list(2)
          list(2) = list(4)
          list(4) = n
        end select
      end if
    end if

  end function side_node_list


  function side_set_node_list (this, n) result (list)

    class(exodus_mesh), intent(in) :: this
    integer, intent(in) :: n
    integer, pointer :: list(:)

    integer :: i, j, k, b, l, len
    integer, pointer :: xside(:), sides(:)

    list => null()

    if (n < 1 .or. n > this%num_sset) return

    !! Pass 1: determine the length of the node list.
    len = 0
    do j = 1, this%sset(n)%num_side

      !! Locate the element in the element block structure.
      b = 1
      l = this%sset(n)%elem(j)
      do while (l > this%eblk(b)%num_elem)
        l = l - this%eblk(b)%num_elem
        b = b + 1
      end do

      select case (trim(this%eblk(b)%elem_type(1:3)))
      case ('TET')
        xside => TETRA4_XSIDE
      case ('PYR')
        xside => PYRAMID5_XSIDE
      case ('WED')
        xside => WEDGE6_XSIDE
      case ('HEX')
        xside => HEX8_XSIDE
      case default
        return  ! unknown element type
      end select

      k = this%sset(n)%face(j)
      if (k < 1 .or. k > size(xside)-1) return  ! bad side index
      len = len + xside(k+1) - xside(k)
    end do

    allocate(list(len))

    !! Pass 2: fill the node list.
    len = 0
    do j = 1, this%sset(n)%num_side

      !! Locate the element in the element block structure.
      b = 1
      l = this%sset(n)%elem(j)
      do while (l > this%eblk(b)%num_elem)
        l = l - this%eblk(b)%num_elem
        b = b + 1
      end do

      select case (trim(this%eblk(b)%elem_type(1:3)))
      case ('TET')
        xside => TETRA4_XSIDE
        sides => TETRA4_SIDES
      case ('PYR')
        xside => PYRAMID5_XSIDE
        sides => PYRAMID5_SIDES
      case ('WED')
        xside => WEDGE6_XSIDE
        sides => WEDGE6_SIDES
      case ('HEX')
        xside => HEX8_XSIDE
        sides => HEX8_SIDES
      end select

      k = this%sset(n)%face(j)
      do i = xside(k), xside(k+1)-1
        len = len + 1
        list(len) = this%eblk(b)%connect(sides(i),l)
      end do

    end do

  end function side_set_node_list


  function side_size_list (elem_type) result (list)

    character(*), intent(in) :: elem_type
    integer, pointer :: list(:)

    select case (elem_type(1:3))
    case ('TET')
      list => TETRA4_SSIZE
    case ('PYR')
      list => PYRAMID5_SSIZE
    case ('WED')
      list => WEDGE6_SSIZE
    case ('HEX')
      list => HEX8_SSIZE
    case default
      list => null()
    end select

  end function side_size_list

  !! This ExodusII-like procedure concatenates the mesh connectivity data that
  !! is stored within an array of element blocks, into a packed mixed-element
  !! ragged-array data structure: CONNECT(XCONNECT(j):XCONNECT(j+1)-1) is the
  !! connectivity of cell j. As a convenience to the caller, the procedure
  !! allocates the allocatable output arrays XCONNECT and CONNECT, primarily
  !! because the required size of CONNECT is not immediately known.
  !! NB: The procedure handles a default-initialized mesh object gracefully,
  !! allocating and defining XCONNECT and CONNECT for a mesh with 0 elements.
  !! This is useful in a parallel context where the mesh object is only
  !! default initialized on all but the IO process. 
    
  subroutine get_concat_elem_conn (this, xconnect, connect)
  
    use,intrinsic :: iso_c_binding, only: c_loc, c_f_pointer
    
    class(exodus_mesh), intent(in), target :: this
    integer, allocatable, intent(out) :: xconnect(:), connect(:)
    
    integer :: j, k, n, offset
    integer, pointer :: flat_conn(:)

    if (this%num_elem == 0) then
      allocate(xconnect(1), connect(0))
      xconnect(1) = 1
      return
    end if
    
    n = sum(this%eblk%num_nodes_per_elem * this%eblk%num_elem)
    allocate(xconnect(this%num_elem+1), connect(n))
    n = 0
    offset = 0
    xconnect(1) = 1
    do j = 1, this%num_eblk
      do k = 1, this%eblk(j)%num_elem
        n = n + 1
        xconnect(n+1) = xconnect(n) + this%eblk(j)%num_nodes_per_elem
      end do
      associate (conn => this%eblk(j)%connect)
        call c_f_pointer (c_loc(conn), flat_conn, shape=[size(conn)])
        connect(offset+1:offset+size(conn)) = flat_conn
        offset = offset + size(conn)
      end associate
    end do
    
  end subroutine get_concat_elem_conn
    
!!!! TYPE BOUND DEBUGGING PROCEDURES FOLLOW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  elemental logical function exodus_mesh_defined (this)
    class(exodus_mesh), intent(in) :: this
    integer :: n
    exodus_mesh_defined = .false.
    !! Check the node coordinate info.
    if (.not.allocated(this%coord)) return
    if (this%num_dim /= size(this%coord,dim=1)) return
    if (this%num_dim < 1 .or. this%num_dim > 3) return
    if (this%num_node /= size(this%coord,dim=2)) return
    !! Check the element block info.
    if (.not.allocated(this%eblk)) return
    if (.not.all(this%eblk%defined())) return
    if (.not.unique(this%eblk%id)) return
    if (this%num_elem /= sum(this%eblk%num_elem)) return
    do n = 1, size(this%eblk) ! Check that the CONNECT values are in-range.
      if (minval(this%eblk(n)%connect) < 1) return
      if (maxval(this%eblk(n)%connect) > this%num_node) return
    end do
    !! If present, check the node set info.
    if (allocated(this%nset)) then
      if (this%num_nset /= size(this%nset)) return
      if (.not.all(this%nset%defined())) return
      if (.not.unique(this%nset%id)) return
      do n = 1, size(this%nset) ! Check that the NODE values are in-range.
        if (minval(this%nset(n)%node) < 1) return
        if (maxval(this%nset(n)%node) > this%num_node) return
      end do
    else
      if (this%num_nset /= 0) return
    end if
    !! If present, check the side set info.
    if (allocated(this%sset)) then
      if (this%num_sset /= size(this%sset)) return
      if (.not.all(this%sset%defined())) return
      if (.not.unique(this%sset%id)) return
      do n = 1, size(this%sset) ! Check that the ELEM and FACE values are in-range.
        if (minval(this%sset(n)%elem) < 1) return
        if (maxval(this%sset(n)%elem) > this%num_elem) return
        if (minval(this%sset(n)%face) < 1) return
        !if (maxval(this%sset(n)%face) > ???) return
      end do
    else
      if (this%num_sset /= 0) return
    end if
    exodus_mesh_defined = .true.
  contains
    pure logical function unique (list)
      integer, intent(in) :: list(:)
      integer :: i, j
      unique = .false.
      do i = 1, size(list)-1
        do j = i+1, size(list)
          if (list(i) == list(j)) return
        end do
      end do
      unique = .true.
    end function unique
  end function exodus_mesh_defined

  elemental logical function elem_blk_defined (this)
    class(elem_blk), intent(in) :: this
    elem_blk_defined = .false.
    if (this%id <= 0) return
    if (.not.allocated(this%connect)) return
    if (this%num_nodes_per_elem /= size(this%connect,dim=1)) return
    if (this%num_elem /= size(this%connect,dim=2)) return
    elem_blk_defined = .true.
  end function elem_blk_defined

  elemental logical function node_set_defined (this)
    class(node_set), intent(in) :: this
    node_set_defined = .false.
    if (this%id <= 0) return
    if (.not.allocated(this%node)) return
    if (this%num_node /= size(this%node)) return
    node_set_defined = .true.
  end function node_set_defined

  elemental logical function side_set_defined (this)
    class(side_set), intent(in) :: this
    side_set_defined = .false.
    if (this%id <= 0) return
    if (.not.allocated(this%elem)) return
    if (.not.allocated(this%face)) return
    if (size(this%elem) /= size(this%face)) return
    if (this%num_side /= size(this%elem)) return
    side_set_defined = .true.
  end function side_set_defined

  subroutine exodus_mesh_dump (this, lun)
    class(exodus_mesh), intent(in) :: this
    integer, intent(in) :: lun
    integer :: j
    write(lun,'(a)') 'EXODUS_MESH('
    write(lun,'(t4,a,i6)') 'NUM_DIM= ', this%num_dim
    write(lun,'(t4,a,i6)') 'NUM_NODE=', this%num_node
    write(lun,'(t4,a,i6)') 'NUM_ELEM=', this%num_elem
    write(lun,'(t4,a,i6)') 'NUM_EBLK=', this%num_eblk
    write(lun,'(t4,a,i6)') 'NUM_NSET=', this%num_nset
    write(lun,'(t4,a,i6)') 'NUM_SSET=', this%num_sset
    if (allocated(this%eblk)) then
      write(lun,'(t4,a)') 'EBLK= ...'
      do j = 1, size(this%eblk)
        call this%eblk(j)%dump(lun)
      end do
    else
      write(lun,'(t4,a)') 'EBLK unallocated'
    end if
    if (allocated(this%nset)) then
      write(lun,'(t4,a)') 'NSET= ...'
      do j = 1, size(this%nset)
        call this%nset(j)%dump(lun)
      end do
    else
      write(lun,'(t4,a)') 'NSET unallocated'
    end if
    if (allocated(this%sset)) then
      write(lun,'(t4,a)') 'SSET= ...'
      do j = 1, size(this%sset)
        call this%sset(j)%dump(lun)
      end do
    else
      write(lun,'(t4,a)') 'SSET unallocated'
    end if
    if (allocated(this%coord)) then
      write(lun,'(t4,a,(t12,3(1x,es22.14)))') 'COORD=', this%coord(:,1)
      do j = 2, size(this%coord,dim=2)
        write(lun,'((t12,3(1x,es22.14)))') this%coord(:,j)
      end do
    else
      write(lun,'(t4,a)') 'COORD unallocated'
    end if
  end subroutine exodus_mesh_dump

  subroutine node_set_dump (this, lun)
    class(node_set), intent(in) :: this
    integer, intent(in) :: lun
    write(lun,'(a)') 'NODE_SET('
    write(lun,'(t4,a,i6)') 'ID=', this%id
    write(lun,'(t4,a,i6)') 'NUM_NODE=', this%num_node
    if (allocated(this%node)) then
      write(lun,'(t4,a,(t10,10(1x,i6)))') 'NODE=', this%node
    else
      write(lun,'(t4,a)') 'NODE unallocated'
    end if
    write(lun,'(t4,a)') ')'
  end subroutine node_set_dump

  subroutine side_set_dump (this, lun)
    class(side_set), intent(in) :: this
    integer, intent(in) :: lun
    write(lun,'(a)') 'SIDE_SET('
    write(lun,'(t4,a,i6)') 'ID=', this%id
    write(lun,'(t4,a,i6)') 'NUM_SIDE=', this%num_side
    if (allocated(this%elem)) then
      write(lun,'(t4,a,(t10,10(1x,i6)))') 'ELEM=', this%elem
    else
      write(lun,'(t4,a)') 'ELEM unallocated'
    end if
    if (allocated(this%face)) then
      write(lun,'(t4,a,(t10,10(1x,i6)))') 'FACE=', this%face
    else
      write(lun,'(t4,a)') 'FACE unallocated'
    end if
    write(lun,'(t4,a)') ')'
  end subroutine side_set_dump

  subroutine elem_blk_dump (this, lun)
    class(elem_blk), intent(in) :: this
    integer, intent(in) :: lun
    integer :: j
    write(lun,'(a)') 'ELEMENT_BLOCK('
    write(lun,'(t4,a,i6)') 'ID=', this%id
    write(lun,'(t4,a,i6)') 'NUM_ELEM=', this%num_elem
    write(lun,'(t4,a,i6)') 'NUM_NODES_PER_ELEM=', this%num_nodes_per_elem
    write(lun,'(t4,a)')  'ELEM_TYPE= "' // this%elem_type // '"'
    if (allocated(this%connect)) then
      write(lun,'(t4,a,(t12,8(1x,i6)))') 'CONNECT=', this%connect(:,1)
      do j = 2, size(this%connect,dim=2)
        write(lun,'((t12,8(1x,i6)))') this%connect(:,j)
      end do
    else
      write(lun,'(t4,a)') 'CONNECT unallocated'
    end if
    write(lun,'(t4,a)') ')'
  end subroutine elem_blk_dump

end module exodus_mesh_type
