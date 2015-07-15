!!
!! Neil N. Carlson <nnc@lanl.gov> March 2015
!!
!! NB: This is incomplete!
!!

program test_graph_module

#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  use graph_type
  !use truchas_logging_services  ! because asserts in graph_type use it
  implicit none
  
  integer :: stat = 0
  
  !call 
  
  call test_0_nodes
  call test_basic
  call test_basic_self
  call test1
  
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
    
  subroutine test1
    type(graph) :: g
    call g%init (3)
    call g%add_clique ([1, 2, 3])
  end subroutine test1


  subroutine write_fail (errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
    character(*), intent(in) :: errmsg
    stat = 1
    write(error_unit,'(a)') errmsg
  end subroutine

end program test_graph_module
