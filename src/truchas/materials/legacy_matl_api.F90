#include "f90_assert.fpp"

module legacy_matl_api

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use matl_module
  use matl_utilities
  use graph_type
  implicit none
  private

  !! Export procedures/data from matl_module and matl_utilities
  public :: nmat, mat_slot, matl, gather_vof
  public :: matl_get_vof, matl_set_vof, read_matl_data, matl_get_cell_vof, update_matl

  !! New initialization procedures
  public :: matl_redef, matl_enddef, matl_add_cell_mat, matl_set_cell_vof
  public :: matl_normalize_vof

  type(graph), allocatable :: g

contains

  subroutine matl_redef
    use legacy_mesh_api, only: ncells
    integer :: n
    if (allocated(g)) deallocate(g)
    allocate(g)
    n = max(ncells, nmat)
    call g%init(n, directed=.true., self_edge=.true.)
  end subroutine matl_redef

  subroutine matl_add_cell_mat(cellid, matid)
    use legacy_mesh_api, only: ncells
    integer, intent(in) :: cellid, matid
    ASSERT(allocated(g))
    ASSERT(cellid > 0 .and. cellid <= ncells)
    ASSERT(matid > 0 .and. matid <= nmat)
    call g%add_edge(cellid, matid)
  end subroutine matl_add_cell_mat

  subroutine matl_enddef

    use legacy_mesh_api, only: ncells

    integer :: i, j
    integer, allocatable :: xadj(:), adjncy(:)

    ASSERT(allocated(g))

    call g%get_adjacency(xadj, adjncy)
    deallocate(g)

    !! Allocate the required number of slots and zero it all out
    mat_slot_new = maxval(xadj(2:ncells+1) - xadj(1:ncells))
    if (mat_slot_new > mat_slot) call slot_increase(matl, mat_slot, mat_slot_new)
    do i = 1, mat_slot
      call slot_set(matl, i)
    end do

    !! Definte the material ids
    do j = 1, ncells
      associate (list => adjncy(xadj(j):xadj(j+1)-1))
        do i = 1, size(list)
          matl(i)%cell(j)%id = list(i)
        end do
      end associate
    end do

  end subroutine matl_enddef

  subroutine matl_set_cell_vof(n, m, vfrac)
    integer, intent(in) :: n, m
    real(r8), intent(in) :: vfrac
    integer :: s
    do s = 1, mat_slot
      if (matl(s)%cell(n)%id == m) then
        matl(s)%cell(n)%vof = vfrac
        return
      end if
    end do
    INSIST(.false.)
  end subroutine matl_set_cell_vof

  !! Possibly called as a final fixup from init_module::matl_init.
  !! Should not be needed.
  subroutine matl_normalize_vof
    use legacy_mesh_api, only: ncells
    integer :: s
    real(r8) :: total(ncells)
    total = 0.0_r8
    do s = 1, mat_slot
      total = total + matl(s)%cell%vof
    end do
    do s = 1, mat_slot
      matl(s)%cell%vof = matl(s)%cell%vof / total
    end do
  end subroutine matl_normalize_vof

end module
