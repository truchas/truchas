!!
!! CELL_TOPOLOGY
!!
!! Neil Carlson <nnc@lanl.gov>
!! 16 Feb 2006
!!
!! PROGRAMMING INTERFACE
!!
!!  LIST => TET_FACE_NODES(CNODES, K[, NORMALIZE][, REVERSE])
!!  LIST => HEX_FACE_NODES(CNODES, K[, NORMALIZE][, REVERSE])
!!
!!    These functions return a pointer to a rank-1 integer array containing
!!    the nodes of one side of a cell.  CNODES is a rank-1 integer array
!!    containing the cell node indices, properly oriented, and K is the side
!!    index of the cell.  CNODES must have size 4 for the tet cell function
!!    and 8 for the hex cell function.  By default the side nodes are listed
!!    CCW about the side with respect to the outward orientation of the side.
!!    If the optional argument REVERSE is present with value true, the nodes
!!    are listed with respect to the inward orientation of the side.  If the
!!    optional argument NORMALIZE is present with value true, the returned
!!    node list begins with the lowest numbered node.  A null pointer is
!!    returned if the side index is out of range.
!!
!!    The only proper use of these functions is as the target of a pointer
!!    assignment; if used otherwise (in an expression, for example) a memory
!!    leak will result.  The caller is responsible for deallocating the
!!    pointer when it is no longer needed.
!!
!!    Note that the returned node list for a side is uniquely determined for
!!    an orientation if normalized, but only when the cell is not degenerate.
!!
!!  CALL REVERSE_FACET (LIST) permutes the node list LIST defining an oriented
!!    facet so as to reverse the orientation of the facet.  The facet is
!!    interpreted as an edge if the size of LIST is 2, in which case the pair
!!    of nodes is swapped.  Otherwise the facet is interpreted as a polygon
!!    with the nodes listed in CCW order about the facet with respect to one
!!    of the two possible orientations.  In this latter case, the reversal of
!!    the orientation leaves the initial node of the list unchanged (and so
!!    this operation commutes with NORMALIZE_FACET).
!!
!!  CALL NORMALIZE_FACET (LIST) permutes the node list LIST defining an
!!    oriented polygonal facet so that the initial node is the smallest of
!!    the nodes while preserving the orientation of the facet.  As long as
!!    the facet is not degenerate this ordering is unique.  For uniformity,
!!    this routine may be passed an edge facet (SIZE(LIST) == 2); in this
!!    case LIST is returned unchanged.
!!

#include "f90_assert.fpp"

module cell_topology

  implicit none
  private

  public :: get_tet_face_nodes, tet_face_nodes
  public :: get_hex_face_nodes, hex_face_nodes, hex_edge_nodes
  public :: parity
  public :: normalize_facet, reverse_facet
  public :: link_face_nodes

  integer, target, public :: TETRA4_XFACE(5), TETRA4_FACES(12), TETRA4_FSIZE(4), TETRA4_FACE_SIG(4)
  data TETRA4_XFACE/1,4,7,10,13/
  data TETRA4_FACES/2,3,4, 1,4,3, 1,2,4, 1,3,2/
  data TETRA4_FSIZE/3,3,3,3/
  data TETRA4_FACE_SIG/b'1110', b'1101', b'1011', b'0111'/

  integer, target, public :: HEX8_XFACE(7), HEX8_FACES(24), HEX8_FSIZE(6), HEX8_FACE_SIG(6)
  data HEX8_XFACE/1,5,9,13,17,21,25/
  data HEX8_FACES/3,4,8,7, 1,2,6,5, 1,5,8,4, 2,3,7,6, 1,4,3,2, 5,6,7,8/
  data HEX8_FSIZE/4,4,4,4,4,4/
  data HEX8_FACE_SIG/b'11001100', b'00110011', b'10011001', b'01100110', b'00001111', b'11110000'/

  integer, private :: HEX8_EDGES(2,12)
  data HEX8_EDGES/1,2, 1,4, 1,5, 2,3, 2,6, 3,4, 3,7, 4,8, 5,6, 5,8, 6,7, 7,8/
  
  integer, target, public :: HEX8_FACE_VERT(4,6)
  data HEX8_FACE_VERT/3,4,8,7, 1,2,6,5, 1,5,8,4, 2,3,7,6, 1,4,3,2, 5,6,7,8/
  
  integer, target, public :: HEX8_VERT_FACE(3,8)
  data HEX8_VERT_FACE/2,3,5, 2,4,5, 1,4,5, 1,3,5, 2,3,6, 2,4,6, 1,4,6, 1,3,6/
  
  integer, target, public :: TETRA4_FACE_VERT(3,4)
  data TETRA4_FACE_VERT/2,3,4, 1,4,3, 1,2,4, 1,3,2/
  
  integer, target, public :: TETRA4_VERT_FACE(3,4)
  data TETRA4_VERT_FACE/2,3,4, 1,3,4, 1,2,4, 1,2,3/
  
  integer :: TRI_LINK_FACE_VERT(3,2)
  data TRI_LINK_FACE_VERT/1,2,3, 4,6,5/
  
  integer :: QUAD_LINK_FACE_VERT(4,2)
  data QUAD_LINK_FACE_VERT/1,2,3,4, 5,8,7,6/

contains

  subroutine get_tet_face_nodes (cnodes, k, fnodes, normalize, reverse)

    integer, intent(in) :: cnodes(:)
    integer, intent(in) :: k
    integer, intent(out) :: fnodes(:)
    logical, intent(in), optional :: normalize
    logical, intent(in), optional :: reverse

    integer :: n

    ASSERT(size(fnodes) == 3)
    ASSERT(size(cnodes) == 4)
    ASSERT(k >= 1 .and. k <= 4)

    fnodes = cnodes(TETRA4_FACES(TETRA4_XFACE(k):TETRA4_XFACE(k+1)-1))

    !! If specified, rotate the smallest value to the initial position.
    if (present(normalize)) then
      if (normalize) fnodes = cshift(fnodes, shift=minloc(fnodes,dim=1)-1)
    end if

    !! If specified, reverse the orientation of the face.
    if (present(reverse)) then
      if (reverse) then
        n = fnodes(2)
        fnodes(2) = fnodes(3)
        fnodes(3) = n
      end if
    end if

  end subroutine get_tet_face_nodes

  function tet_face_nodes (cnodes, k, normalize, reverse) result (list)

    integer, intent(in) :: cnodes(:)
    integer, intent(in) :: k
    logical, intent(in), optional :: normalize
    logical, intent(in), optional :: reverse
    integer, pointer :: list(:)

    integer :: n

    list => null()
    if (k < 1 .or. k > 4 .or. size(cnodes) /= 4) return
    allocate(list(3))
    list = cnodes(TETRA4_FACES(TETRA4_XFACE(k):TETRA4_XFACE(k+1)-1))

    !! If specified, rotate the smallest value to the initial position.
    if (present(normalize)) then
      if (normalize) then
        list = cshift(list, shift=minloc(list,dim=1)-1)
      end if
    end if

    !! If specified, reverse the orientation of the face.
    if (present(reverse)) then
      if (reverse) then
        n = list(2)
        list(2) = list(3)
        list(3) = n
      end if
    end if

  end function tet_face_nodes


  subroutine get_hex_face_nodes (cnodes, k, fnodes, normalize, reverse)

    integer, intent(in) :: cnodes(:)
    integer, intent(in) :: k
    integer, intent(out) :: fnodes(:)
    logical, intent(in), optional :: normalize
    logical, intent(in), optional :: reverse

    integer :: n

    ASSERT(size(fnodes) == 4)
    ASSERT(size(cnodes) == 8)
    ASSERT(k >= 1 .and. k <= 6)

    fnodes = cnodes(HEX8_FACES(HEX8_XFACE(k):HEX8_XFACE(k+1)-1))

    !! If specified, rotate the smallest value to the initial position.
    if (present(normalize)) then
      if (normalize) fnodes = cshift(fnodes, shift=minloc(fnodes,dim=1)-1)
    end if

    !! If specified, reverse the orientation of the face.
    if (present(reverse)) then
      if (reverse) then
        n = fnodes(2)
        fnodes(2) = fnodes(4)
        fnodes(4) = n
      end if
    end if

  end subroutine get_hex_face_nodes


  function hex_face_nodes (cnodes, k, normalize, reverse) result (list)

    integer, intent(in) :: cnodes(:)
    integer, intent(in) :: k
    logical, intent(in), optional :: normalize
    logical, intent(in), optional :: reverse
    integer, pointer :: list(:)

    integer :: n

    list => null()
    if (k < 1 .or. k > 6 .or. size(cnodes) /= 8) return
    allocate(list(4))
    list = cnodes(HEX8_FACES(HEX8_XFACE(k):HEX8_XFACE(k+1)-1))

    !! If specified, rotate the smallest value to the initial position.
    if (present(normalize)) then
      if (normalize) then
        list = cshift(list, shift=minloc(list,dim=1)-1)
      end if
    end if

    !! If specified, reverse the orientation of the face.
    if (present(reverse)) then
      if (reverse) then
        n = list(2)
        list(2) = list(4)
        list(4) = n
      end if
    end if

  end function hex_face_nodes


  function hex_edge_nodes (cnodes, k, reverse) result (list)

    integer, intent(in) :: cnodes(:)
    integer, intent(in) :: k
    logical, intent(in), optional :: reverse
    integer, pointer :: list(:)

    integer :: n

    list => null()
    if (k < 1 .or. k > 12 .or. size(cnodes) /= 8) return
    allocate(list(2))
    list = cnodes(HEX8_EDGES(:,k))

    !! If specified, reverse the orientation of the edge.
    if (present(reverse)) then
      if (reverse) then
        n = list(1)
        list(1) = list(2)
        list(2) = n
      end if
    end if

  end function hex_edge_nodes

  !! This function return a pointer to the list of nodes defining the specified
  !! oriented link face K in {1,2}. The faces are oriented outward with respect
  !! to the cells they belongs to; if we regard the link nodes as defining a
  !! wedge or hex element, the faces are inward oriented with respect to it.
  !! The pointer target is allocated by the function and the caller is
  !! responsible for deallocating it when it is no longer needed.

  pure function link_face_nodes (lnodes, k, normalize, reverse) result (fnodes)

    integer, intent(in) :: lnodes(:), k
    logical, intent(in), optional :: normalize, reverse
    integer, pointer :: fnodes(:)

    fnodes => null()
    if (k < 1 .or. k > 2) return
    select case (size(lnodes))
    case (6) ! triangular link faces
      allocate(fnodes(3))
      fnodes = lnodes(TRI_LINK_FACE_VERT(:,k))
    case (8) ! quadrilateral link faces
      allocate(fnodes(4))
      fnodes = lnodes(QUAD_LINK_FACE_VERT(:,k))
    end select
    if (.not.associated(fnodes)) return
    
    !! If specified, rotate the smallest value to the initial position.
    if (present(normalize)) then
      if (normalize) call normalize_facet (fnodes)
    end if
    
    !! If specified, reverse the orientation of the face.
    if (present(reverse)) then
      if (reverse) call reverse_facet (fnodes)
    end if

  end function link_face_nodes


  pure integer function parity (f, g)
  
    integer, intent(in) :: f(:), g(:)
    
    integer :: n, j
    
    parity = 0
    
    if (size(f) /= size(g)) return
    
    select case (size(f))
    
    case (2)  ! edge facets
    
      if (f(1) == g(1)) then
        if (f(2) == g(2)) parity = 1
      else if (f(1) == g(2)) then
        if (f(2) == g(1)) parity = -1
      end if
      
    case (3)  ! tri facets
      
      if (f(1) == g(1)) then
        if (f(2) == g(2)) then
          if (f(3) == g(3)) parity = 1
        else if (f(2) == g(3)) then
          if (f(3) == g(2)) parity = -1
        end if
      end if
    
    case (4)  ! quad facets
    
      if (f(1) == g(1)) then
        if (f(2) == g(2)) then
          if (f(3) /= g(3)) return
          if (f(4) /= g(4)) return
          parity = 1
        else if (f(2) == g(4)) then
          if (f(3) /= g(3)) return
          if (f(4) /= g(2)) return
          parity = -1
        end if
      end if
      
    case (5:) ! more general facets
    
      n = size(f)
      if (f(1) == g(1)) then
        if (f(2) == g(2)) then
          do j = 3, n
            if (f(j) /= g(j)) return
          end do
          parity = 1
        else if (f(2) == g(n)) then
          do j = 3, n
            if (f(j) /= g(n+2-j)) return
          end do
          parity = -1
        end if
      end if
      
    end select
    
  end function parity

  pure subroutine normalize_facet (list)
    integer, intent(inout) :: list(:)
    if (size(list) > 2) list = cshift(list, shift=minloc(list,dim=1)-1)
  end subroutine normalize_facet
    
  pure subroutine reverse_facet (list)
    integer, intent(inout) :: list(:)
    integer :: j1, j2, n
    select case (size(list))
    case (2)
      n = list(1)
      list(1) = list(2)
      list(2) = n
    case (3)
      n = list(2)
      list(2) = list(3)
      list(3) = n
    case (4)
      n = list(2)
      list(2) = list(4)
      list(4) = n
    case (5:)
      j1 = 2
      j2 = size(list)
      do while (j1 < j2)
        n = list(j1)
        list(j1) = list(j2)
        list(j2) = n
        j1 = j1 + 1
        j2 = j2 - 1
      end do
    end select
  end subroutine reverse_facet

end module cell_topology
