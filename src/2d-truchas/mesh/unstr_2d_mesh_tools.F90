! TODO:
! This contains versions of procedures from UNSTR_MESH_TOOLS specialized for
! 2D meshes. They are very nearly identical -- can this be merged with it?

#include "f90_assert.fpp"

module unstr_2d_mesh_tools

  implicit none
  private

  public :: get_cell_neighbor_array, label_mesh_faces

  interface get_cell_neighbor_array
    procedure get_cell_neighbor_array_1
  end interface

contains

  subroutine get_cell_neighbor_array_1(cstart, cnode, cnhbr, stat)

    use face_neighbor_table_type

    integer, intent(in) :: cstart(:), cnode(:)
    integer, allocatable, intent(out) :: cnhbr(:)
    integer, intent(out) :: stat

    type(face_neighbor_table) :: nhbr_table

    call nhbr_table%init(cstart, cnode, dim=2)
    call get_cnhbr_aux(nhbr_table, cstart, cnode, cnhbr, stat)

  end subroutine get_cell_neighbor_array_1

  subroutine get_cnhbr_aux(nhbr_table, cstart, cnode, cnhbr, stat)

    use face_neighbor_table_type
    use cell_topology_2d, only: get_face_nodes

    type(face_neighbor_table), intent(in) :: nhbr_table
    integer, intent(in) :: cstart(:), cnode(:)
    integer, allocatable, intent(out) :: cnhbr(:)
    integer, intent(out) :: stat

    integer :: i, j, k, jj, kk, nmatch, bad_faces, ncell, offset
    integer, allocatable :: face(:)
    type(face_neighbor), allocatable :: nhbrs(:)

    ncell = size(cstart) - 1

    allocate(cnhbr(cstart(ncell+1)-1))

    cnhbr = 0
    bad_faces = 0
    stat = 0

    do j = 1, ncell
      associate (cell => cnode(cstart(j):cstart(j+1)-1))
        offset = cstart(j)-1  ! for getting local face index
        do k = cstart(j), cstart(j+1)-1
          if (cnhbr(k) /= 0) cycle  ! info already assigned
          !! Get a face and the list of its neighbor cells.
          call get_face_nodes(cell, k-offset, face)
          call nhbr_table%get_neighbors(face, nhbrs)
          !! Locate the cell neighbor, but scan all for valid topology.
          jj = 0
          kk = 0
          nmatch = 0
          do i = 1, size(nhbrs)
            if (nhbrs(i)%j == j) cycle  ! found myself
            nmatch = nmatch + 1
            if (nhbrs(i)%p == -1) then  ! a good match (if only one)
              jj = nhbrs(i)%j
              kk = nhbrs(i)%k
            end if
          end do
          !! Store the neighbor information.
          if (nmatch == 1 .and. jj /= 0) then
            !! Found a unique neighbor; assign the neighbor data.
            cnhbr(k) = jj ! my neighbor, and
            cnhbr(kk) = j ! my neighbor's neighbor (me!)
          else if (nmatch /= 0) then
            !! Bad mesh topology; tag the faces involved.
            bad_faces = 1 + bad_faces
            do i = 1, size(nhbrs)
              cnhbr(nhbrs(i)%k) = -bad_faces
            end do
            stat = -1
          end if
        end do
      end associate
    end do

  end subroutine get_cnhbr_aux


  subroutine label_mesh_faces(cstart, cnode, nface, cface)

    use cell_topology_2d
    use facet_table_type

    integer, intent(in)  :: cstart(:), cnode(:)
    integer, intent(out) :: nface
    integer, allocatable, intent(out) :: cface(:)

    integer :: j, k, n, offset, ncell, nlink, max_face, node_max
    integer, allocatable :: face(:)
    type(facet_table) :: table

    ncell = size(cstart) - 1

    allocate(cface(cstart(ncell+1)-1))

    ASSERT(minval(cnode) > 0)
    ASSERT(size(cface) == size(cnode))

    max_face = size(cface)  ! worst case; realistically, closer to half this
    node_max = maxval(cnode)
    call table%init(max_face, node_max)

    do j = 1, ncell
      associate (cell => cnode(cstart(j):cstart(j+1)-1))
        offset = cstart(j)-1  ! for getting local face index
        do k = cstart(j), cstart(j+1)-1
          call get_face_nodes(cell, k-offset, face)
          call table%get_facet_label(face, cface(k), insert=.true.)
        end do
      end associate
    end do
    nface = table%number_of_facets()
    INSIST(all(cface /= 0))

  end subroutine label_mesh_faces

end module unstr_2d_mesh_tools
