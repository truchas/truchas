!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  unsplit_geometry_volume_tracker_type
!! 
!!  Author: Robert Chiodi (robertchiodi@lanl.gov)
!!  June 2019
!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module traversal_tracker_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64, i8 => int64

  integer, parameter, public :: indiv_cell_capacity = 200

  type, public :: traversal_tracker
     integer, private :: cells_encountered(indiv_cell_capacity)
     integer, private :: size     
     integer, private :: number_of_cells_visited
   contains
     procedure, public :: init
     procedure, public :: still_cells_to_visit
     procedure, public :: add_cell
     procedure, public :: get_next_cell
     procedure, public :: cell_not_encountered
     procedure, public :: cells_left_to_visit
     procedure, public :: get_number_of_visited_cells
     procedure, public :: get_visited_cell_index
  end type traversal_tracker
  
 
contains

  subroutine init(this, a_starting_index)

    class(traversal_tracker), intent(out) :: this
    integer, intent(in) :: a_starting_index(:)

    integer :: n
    
    this%cells_encountered = -1
    this%number_of_cells_visited = 0
    this%size = 0
    do n = 1, size(a_starting_index)
      this%size = this%size + 1
      this%cells_encountered(this%size) = a_starting_index(n)
    end do
    
  end subroutine init

  function cells_left_to_visit(this) result(a_cells_to_visit)

    class(traversal_tracker), intent(in) :: this
    integer :: a_cells_to_visit

    a_cells_to_visit = this%size - this%number_of_cells_visited
    return
  end function cells_left_to_visit
  
  function still_cells_to_visit(this) result(keep_going)

    class(traversal_tracker), intent(in) :: this
    logical :: keep_going

    keep_going = (this%cells_left_to_visit()) /= 0
    return
  end function still_cells_to_visit

  subroutine add_cell(this, a_index)

    class(traversal_tracker), intent(inout) :: this
    integer, intent(in) :: a_index

    ASSERT(this%size < indiv_cell_capacity)
    this%size = this%size + 1
    this%cells_encountered(this%size) = a_index

  end subroutine add_cell

  ! When you get this next cell ID, it will be
  ! marked as visited, so make sure you actually
  ! use it!
  function get_next_cell(this) result(a_cell_index)

    class(traversal_tracker), intent(inout) :: this
    integer :: a_cell_index

    ASSERT(this%cells_left_to_visit() > 0)
    this%number_of_cells_visited = this%number_of_cells_visited + 1    
    a_cell_index = this%cells_encountered(this%number_of_cells_visited)
    return
  end function get_next_cell

  function cell_not_encountered(this, a_index) result(a_not_encountered)

    class(traversal_tracker), intent(in) :: this
    integer, intent(in) :: a_index
    logical :: a_not_encountered

    integer :: n

    a_not_encountered = .true.
    do n = 1, this%size
       if(this%cells_encountered(n) == a_index) then
          a_not_encountered = .false.
          return
       end if
    end do
    
    return
  end function cell_not_encountered

  function get_number_of_visited_cells(this) result(a_size)
    
    class(traversal_tracker), intent(in) :: this
    integer :: a_size

    a_size = this%number_of_cells_visited
    return
    
  end function get_number_of_visited_cells

  function get_visited_cell_index(this, a_index) result(a_cell_index)
    class(traversal_tracker), intent(in) :: this
    integer, intent(in) :: a_index
    integer :: a_cell_index

    ASSERT(a_index <= this%number_of_cells_visited)
    a_cell_index = this%cells_encountered(a_index)

  end function get_visited_cell_index

end module traversal_tracker_type
