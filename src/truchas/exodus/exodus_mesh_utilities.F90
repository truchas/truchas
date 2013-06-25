!!
!! EXODUS_MESH_UTILITIES
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 28 Oct 2005
!!
!! This module provides a collection of useful mesh procedures.  These
!! are provisional/experimental.
!!

module exodus_mesh_utilities

  use exodus_mesh_type
#ifdef PGI_COMPILER_WORKAROUND
  ! bizarrely, this works around a bug in mesh_importer when use'ing exodus
#else
  implicit none
#endif
  private

  public :: side_node_list, side_size_list, side_set_node_list

  !! Element specification data (static) as described in the Exodus II manual.
  !! SIDES(XSIDE(K):XSIDE(K+1)-1) is the list of local nodes of side K for a
  !! particular type of element.  SSIZE(K) = XSIDE(K+1)-XSIDE(K) is the size
  !! (number of nodes) of side K of the element.

  integer, target, private :: TETRA4_XSIDE(5), TETRA4_SIDES(12), TETRA4_SSIZE(4)
  data TETRA4_XSIDE/1,4,7,10,13/
  data TETRA4_SIDES/1,2,4, 2,3,4, 1,4,3, 1,3,2/
  data TETRA4_SSIZE/3,3,3,3/

  integer, target, private :: WEDGE6_XSIDE(6), WEDGE6_SIDES(18), WEDGE6_SSIZE(5)
  data WEDGE6_XSIDE/1,5,9,13,16,19/
  data WEDGE6_SIDES/1,2,5,4, 2,3,6,5, 1,4,6,3, 1,3,2, 4,5,6/
  data WEDGE6_SSIZE/4,4,4,3,3/

  integer, target, private :: HEX8_XSIDE(7), HEX8_SIDES(24), HEX8_SSIZE(6)
  data HEX8_XSIDE/1,5,9,13,17,21,25/
  data HEX8_SIDES/1,2,6,5, 2,3,7,6, 3,4,8,7, 1,5,8,4, 1,4,3,2, 5,6,7,8/
  data HEX8_SSIZE/4,4,4,4,4,4/

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SIDE_NODE_LIST
 !!
 !! This function returns a pointer to a rank-1 integer array containing the
 !! nodes on one side of an element in the mesh data structure MESH.  The side
 !! and element indices are specified by SIDE and ELEM.  By default the nodes
 !! are listed CCW about the side with respect to the outward orientation of
 !! the side.  If the optional argument REVERSE is present with value true,
 !! the nodes are listed with respect to the inward orientation of the side.
 !! The node list begins with the lowest numered node.  MESH must be defined;
 !! if it is not, the behavior of the function will be unpredictable. A null
 !! pointer is returned If ELEM or SIDE are out-of-range, or if the element is
 !! of unknown type.  Currently, only 4-node tets, 6-node wedges, and 8-node
 !! hexes are handled.
 !!
 !! NB: The only proper use of the function is as the target of a pointer
 !! assignment; if used otherwise (in an expression, e.g.) a memory leak will
 !! result.
 !!

  function side_node_list (mesh, elem, side, reverse) result (list)

    type(exodus_mesh), intent(in) :: mesh
    integer, intent(in) :: elem, side
    logical, intent(in), optional :: reverse
    integer, pointer :: list(:)

    integer :: b, l, n
    integer, pointer :: xside(:), sides(:)

    list => null()

    if (elem < 1 .or. elem > mesh%num_elem) return

    !! Locate the element as element L of block B.
    b = 1
    l = elem
    do while (l > mesh%eblk(b)%num_elem)
      l = l - mesh%eblk(b)%num_elem
      b = b + 1
    end do

    select case (mesh%eblk(b)%elem_type)
    case ('TETRA', 'TETRA4')
      xside => TETRA4_XSIDE
      sides => TETRA4_SIDES
    case ('WEDGE', 'WEDGE6')
      xside => WEDGE6_XSIDE
      sides => WEDGE6_SIDES
    case ('HEX', 'HEX8')
      xside => HEX8_XSIDE
      sides => HEX8_SIDES
    case default
      return  ! unknown element type
    end select

    if (side < 1 .or. side > size(xside)-1) return  ! bad side index

    allocate(list(xside(side+1)-xside(side)))
    list = mesh%eblk(b)%connect(sides(xside(side):xside(side+1)-1),l)

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

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SIDE_SIZE_LIST
 !!
 !! This function returns a pointer to a rank-1 integer array containing the
 !! sizes of the sides of an element whose type is specified by the argument.
 !! The target of this pointer should NEVER be altered.
 !!

  function side_size_list (elem_type) result (list)

    character(len=*), intent(in) :: elem_type
    integer, pointer :: list(:)

    select case (elem_type)
    case ('TETRA', 'TETRA4')
      list => TETRA4_SSIZE
    case ('WEDGE', 'WEDGE6')
      list => WEDGE6_SSIZE
    case ('HEX', 'HEX8')
      list => HEX8_SSIZE
    case default
      list => null()
    end select

  end function side_size_list

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SIDE_SET_NODE_LIST
 !!
 !! This function returns a pointer to a rank-1 integer array holding the side
 !! set node list for side set N in the Exodus MESH; c.f. EXGSSN described in
 !! section 3.11.4 in the Exodus II reference manual.  This data is written by
 !! the utility program exotxt, but it does not actually exist in the binary
 !! Exodus file.  Currently, only 4-node tets, 6-node wedges, and 8-node hexes
 !! are handled.  MESH must be well-defined; if it is not, the behavior of the
 !! function will be unpredictable.  A null pointer is returned if anything
 !! not understood is encountered.
 !!
 !! NB: The only proper use of the function is as the target of a pointer
 !! assignment; if used otherwise (in an expression, e.g.) a memory leak will
 !! result.
 !!

  function side_set_node_list (mesh, n) result (list)

    type(exodus_mesh), intent(in) :: mesh
    integer, intent(in) :: n
    integer, pointer :: list(:)

    integer :: i, j, k, b, l, len
    integer, pointer :: xside(:), sides(:)

    list => null()

    if (n < 1 .or. n > mesh%num_sset) return

    !! Pass 1: determine the length of the node list.
    len = 0
    do j = 1, mesh%sset(n)%num_side

      !! Locate the element in the element block structure.
      b = 1
      l = mesh%sset(n)%elem(j)
      do while (l > mesh%eblk(b)%num_elem)
        l = l - mesh%eblk(b)%num_elem
        b = b + 1
      end do

      select case (trim(mesh%eblk(b)%elem_type))
      case ('TETRA', 'TETRA4')
        xside => TETRA4_XSIDE
      case ('WEDGE', 'WEDGE6')
        xside => WEDGE6_XSIDE
      case ('HEX', 'HEX8')
        xside => HEX8_XSIDE
      case default
        return  ! unknown element type
      end select

      k = mesh%sset(n)%face(j)
      if (k < 1 .or. k > size(xside)-1) return  ! bad side index
      len = len + xside(k+1) - xside(k)
    end do

    allocate(list(len))

    !! Pass 2: fill the node list.
    len = 0
    do j = 1, mesh%sset(n)%num_side

      !! Locate the element in the element block structure.
      b = 1
      l = mesh%sset(n)%elem(j)
      do while (l > mesh%eblk(b)%num_elem)
        l = l - mesh%eblk(b)%num_elem
        b = b + 1
      end do

      select case (trim(mesh%eblk(b)%elem_type))
      case ('TETRA', 'TETRA4')
        xside => TETRA4_XSIDE
        sides => TETRA4_SIDES
      case ('WEDGE', 'WEDGE6')
        xside => WEDGE6_XSIDE
        sides => WEDGE6_SIDES
      case ('HEX', 'HEX8')
        xside => HEX8_XSIDE
        sides => HEX8_SIDES
      end select

      k = mesh%sset(n)%face(j)
      do i = xside(k), xside(k+1)-1
        len = len + 1
        list(len) = mesh%eblk(b)%connect(sides(i),l)
      end do

    end do

  end function side_set_node_list

end module exodus_mesh_utilities
