program test_block_pcsr_bridge

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use index_map_type
  use pbsr_matrix_type
  use pcsr_matrix_type
  use block_pcsr_bridge_type
  implicit none

  integer :: status

  call init_parallel_communication
  status = 0
  call test_bridge
  call halt_parallel_communication
  stop status

contains

  subroutine test_bridge
    type(index_map), target :: imap
    type(pcsr_graph), pointer :: graph
    type(pbsr_matrix) :: block_matrix
    type(block_pcsr_bridge), target :: bridge
    type(pcsr_matrix), pointer :: scalar_matrix
    real(r8) :: diagonal(2,2), upper(2,2), lower(2,2)
    real(r8), allocatable :: x_block(:,:), y_block(:,:), x_scalar(:), y_scalar(:)
    integer :: i, offp(1)

    allocate(graph)
    if (nPE == 1) then
      call imap%init(3)
    else
      offp = modulo(this_PE, nPE)*3 + 1
      call imap%init(3, offp)
    end if
    call graph%init(imap)
    call graph%add_edge(1, 1)
    call graph%add_edge(2, 2)
    call graph%add_edge(3, 3)
    call graph%add_edge(1, 2)
    call graph%add_edge(2, 1)
    call graph%add_edge(2, 3)
    call graph%add_edge(3, 2)
    if (nPE > 1) call graph%add_edge(3, 4)
    call graph%add_complete()
    call block_matrix%init(2, graph, take_graph=.true.)

    diagonal = reshape([3.0_r8, 1.0_r8, 1.0_r8, 4.0_r8], [2,2])
    upper = reshape([-1.0_r8, 2.0_r8, 0.5_r8, -2.0_r8], [2,2])
    lower = transpose(upper)
    do i = 1, 3
      call block_matrix%set(i, i, diagonal)
    end do
    call block_matrix%set(1, 2, upper)
    call block_matrix%set(2, 1, lower)
    call block_matrix%set(2, 3, upper)
    call block_matrix%set(3, 2, lower)

    call bridge%update(block_matrix)
    scalar_matrix => bridge%matrix()
    allocate(x_block(2,imap%local_size), y_block(2,imap%onp_size))
    x_block = reshape([(real(i, r8), i=1,size(x_block))], shape(x_block))
    call imap%gather_offp(x_block)
    call block_matrix%matvec(x_block, y_block)
    allocate(x_scalar(size(x_block)), y_scalar(size(y_block)))
    x_scalar = reshape(x_block, shape(x_scalar))
    call scalar_matrix%matvec(x_scalar, y_scalar)
    call require(maxval(abs(y_scalar - reshape(y_block, shape(y_scalar)))) < 1.0e-12_r8, &
        'scalar matrix does not reproduce block matrix product')

    diagonal = 2.0_r8*diagonal
    do i = 1, 3
      call block_matrix%set(i, i, diagonal)
    end do
    call bridge%update(block_matrix)
    call block_matrix%matvec(x_block, y_block)
    call scalar_matrix%matvec(x_scalar, y_scalar)
    call require(maxval(abs(y_scalar - reshape(y_block, shape(y_scalar)))) < 1.0e-12_r8, &
        'bridge did not refresh scalar matrix values')
  end subroutine


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine

end program test_block_pcsr_bridge
