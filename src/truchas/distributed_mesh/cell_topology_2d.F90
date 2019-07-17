! TODO:
! This contains versions of procedures from CELL_TOPOLOGY specialized
! for 2D cell types. Can this be reasonably merged into that module?

module cell_topology_2d

  implicit none
  private

  public :: get_face_nodes

contains

  !! Returns the list of nodes defining the specified face of the given cell.
  !! CNODES is the list of nodes defining the cell and INDEX the face of the
  !! cell.  The result is returned in the allocatable array FNODES and is
  !! oriented outward with respect to the cell. If REVERSE is present with
  !! value true, the returned FNODES list is oriented inward with respect to
  !! the cell. This works for any polygonal cell, not just tri and quad.
  !! The optional NORMALIZE argument is unused but retained in order to keep
  !! the same interface as the 3D version

  pure subroutine get_face_nodes(cnodes, index, fnodes, normalize, reverse)

    integer, intent(in) :: cnodes(:), index
    integer, allocatable, intent(inout) :: fnodes(:)
    logical, intent(in), optional :: normalize, reverse

    integer :: n

    if (index == size(cnodes)) then
      fnodes = [cnodes(index), cnodes(1)]
    else
      fnodes = cnodes(index:index+1)
    end if

    if (present(reverse)) then
      if (reverse) then
        n = fnodes(1)
        fnodes(1) = fnodes(2)
        fnodes(2) = n
      end if
    end if

  end subroutine get_face_nodes

end module cell_topology_2d
