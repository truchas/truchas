!!
!! Neil N. Carlson <nnc@lanl.gov> March 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NB: This is incomplete!
!!

program test_graph_module

#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  use graph_type
  implicit none
  
  integer :: stat = 0
  
  call test_0_nodes
  call test_basic
  call test_basic_self
  call test_comp1
  call test_comp2
  
  call exit (stat)
  
contains

  subroutine test_0_nodes
    type(graph) :: g
    integer, allocatable :: xadj(:), adjncy(:)
    call g%init (0)
    call g%get_adjacency (xadj, adjncy)
    if (allocated(xadj)) then
      if (size(xadj) /= 1 .or.  xadj(1) /= 1) then
        call write_fail ('test_0_nodes: wrong xadj size/value')
      end if
    else
      call write_fail ('test_0_nodes: xadj not allocated')
    end if
    if (allocated(adjncy)) then
      if (size(adjncy) /= 0) then
        call write_fail ('test_0_nodes: wrong adjncy size')
      end if
    else
      call write_fail ('test_0_nodes: adjncy not allocated')
    end if
  end subroutine test_0_nodes
  
  subroutine test_basic
    type(graph), allocatable :: g
    integer, allocatable :: xadj(:), adjncy(:)
    allocate(g)
    call g%init (3)
    call g%get_adjacency (xadj, adjncy)
    if (size(xadj) /= 4) call write_fail ('test_basic: xadj wrong size')
    if (any(xadj /= 1)) call write_fail ('test_basic: xadj wrong values')
    if (size(adjncy) /= 0) call write_fail ('test_basic: wrong adjncy size')
    call g%add_clique ([3,2,1])
    call g%add_edge (1, [1,2,3]) ! duplicates with self edge
    call g%add_edge (2, 3)  ! duplicate
    call g%get_adjacency (xadj, adjncy)
    if (size(xadj) /= 4) call write_fail ('test_basic: xadj wrong size')
    if (any(xadj /= [1,3,5,7])) call write_fail ('test_basic: xadj wrong values')
    if (any(adjncy /= [2,3,1,3,1,2])) call write_fail ('test_basic: xadj wrong values')
  end subroutine test_basic
  
  subroutine test_basic_self
    type(graph), allocatable :: g
    integer, allocatable :: xadj(:), adjncy(:)
    allocate(g)
    call g%init (3, self_edge=.true.)
    call g%get_adjacency (xadj, adjncy)
    if (size(xadj) /= 4) call write_fail ('test_basic_self: xadj wrong size')
    if (any(xadj /= 1)) call write_fail ('test_basic_self: xadj wrong values')
    if (size(adjncy) /= 0) call write_fail ('test_basic_self: wrong adjncy size')
    call g%add_clique ([3,2])
    call g%add_edge (1, [2,3])
    call g%add_edge (2, 3)  ! duplicate
    call g%get_adjacency (xadj, adjncy)
    if (size(xadj) /= 4) call write_fail ('test_basic_self: xadj wrong size')
    if (any(xadj /= [1,3,6,9])) call write_fail ('test_basic_self: xadj wrong values')
    if (any(adjncy /= [2,3,1,2,3,1,2,3])) call write_fail ('test_basic_self: xadj wrong values')
  end subroutine test_basic_self
    
  subroutine test_comp1
    type(graph), allocatable :: g
    integer :: ncomp
    integer, allocatable :: xcomp(:), comp(:)
    allocate(g)
    call g%init (4)
    call g%get_components (ncomp, xcomp, comp)
    if (ncomp /= 4) call write_fail ('test_comp1: wrong ncomp value')
    if (size(xcomp) /= 5) call write_fail ('test_comp1: wrong xcomp size')
    if (any(xcomp /= [1,2,3,4,5])) call write_fail ('test_comp1: wrong xcomp values')
    if (size(comp) /= 4) call write_fail ('test_comp1: wrong comp size')
    if (any(comp /= [1,2,3,4])) call write_fail ('test_comp1: wring comp values')
  end subroutine test_comp1
  
  subroutine test_comp2
    type(graph), allocatable :: g
    integer :: ncomp
    integer, allocatable :: xcomp(:), comp(:)
    allocate(g)
    call g%init (9)
    call g%add_clique ([1,3,5])
    call g%add_clique ([4,6,8,9])
    call g%get_components (ncomp, xcomp, comp)
    if (ncomp /= 4) call write_fail ('test_comp2: wrong ncomp value')
    if (size(xcomp) /= 5) call write_fail ('test_comp2: wrong xcomp size')
    if (any(xcomp /= [1,4,5,9,10])) call write_fail ('test_comp2: wrong xcomp values')
    if (size(comp) /= 9) call write_fail ('test_comp2: wrong comp size')
    if (any(comp /= [1,3,5,2,4,6,8,9,7])) call write_fail ('test_comp2: wrong comp values')
    call g%add_edge (2, [5,7])
    call g%add_edge (4, 7)
    call g%get_components (ncomp, xcomp, comp)
    if (ncomp /= 1) call write_fail ('test_comp2: wrong ncomp value')
    if (size(xcomp) /= 2) call write_fail ('test_comp2: wrong xcomp size')
    if (any(xcomp /= [1,10])) call write_fail ('test_comp2: wrong xcomp values')
    if (size(comp) /= 9) call write_fail ('test_comp2: wrong comp size')
    if (any(comp /= [1,2,3,4,5,6,7,8,9])) call write_fail ('test_comp2: wrong comp values')
  end subroutine test_comp2

  subroutine write_fail (errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
    character(*), intent(in) :: errmsg
    stat = 1
    write(error_unit,'(a)') errmsg
  end subroutine

end program test_graph_module
