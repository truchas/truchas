program test_polygon_sides

  use exodus_mesh_type
  use exodus_mesh_io
  implicit none

  type(exodus_mesh) :: mesh
  integer, pointer :: sizes(:), list(:)
  integer :: status

  status = 0
  call test_quad(mesh, status)
  call test_tri(mesh, status)
  if (command_argument_count() > 0) call test_file(mesh, status)
  if (status /= 0) stop 1
  print '(a)', 'PASS'

contains

  subroutine test_quad(mesh, status)
    type(exodus_mesh), intent(out) :: mesh
    integer, intent(inout) :: status

    integer :: side

    mesh%num_dim = 2
    mesh%num_node = 4
    mesh%num_elem = 1
    mesh%num_eblk = 1
    mesh%num_sset = 1
    allocate(mesh%eblk(1), mesh%sset(1))
    mesh%eblk(1)%num_elem = 1
    mesh%eblk(1)%num_nodes_per_elem = 4
    mesh%eblk(1)%elem_type = 'QUAD4'
    allocate(mesh%eblk(1)%connect(4,1)); mesh%eblk(1)%connect(:,1) = [1,2,3,4]
    mesh%sset(1)%num_side = 4
    allocate(mesh%sset(1)%elem(4), mesh%sset(1)%face(4))
    mesh%sset(1)%elem = 1
    mesh%sset(1)%face = [1,2,3,4]

    sizes => mesh%side_size_list('QUAD4')
    call require(all(sizes == [2,2,2,2]), 'QUAD4 side sizes', status)
    do side = 1, 4
      list => mesh%side_node_list(1, side)
      call require(size(list) == 2, 'QUAD4 side node count', status)
      deallocate(list)
    end do
    list => mesh%side_set_node_list(1)
    call require(all(list == [1,2,2,3,3,4,4,1]), 'QUAD4 side-set nodes', status)
    deallocate(list)
  end subroutine


  subroutine test_file(mesh, status)
    type(exodus_mesh), intent(out) :: mesh
    integer, intent(inout) :: status

    character(256) :: filename
    integer :: io_status
    integer, pointer :: list(:)

    call get_command_argument(1, filename)
    call read_exodus_mesh(trim(filename), mesh, io_status)
    call require(io_status == 0, 'read Exodus mesh', status)
    if (io_status /= 0) return

    call require(mesh%num_dim == 2, '2D Exodus mesh dimension', status)
    call require(mesh%num_eblk == 1, '2D Exodus element block count', status)
    call require(mesh%eblk(1)%num_nodes_per_elem == 4, &
                 '2D Exodus QUAD4 node count', status)
    call require(mesh%eblk(1)%elem_type(1:4) == 'QUAD', &
                 '2D Exodus QUAD element type', status)

    list => mesh%side_node_list(1, 1)
    call require(size(list) == 2, 'read QUAD4 side node count', status)
    deallocate(list)
    list => mesh%side_set_node_list(1)
    call require(size(list) == 2 * mesh%sset(1)%num_side, &
                 'read QUAD4 side-set node count', status)
    deallocate(list)
  end subroutine


  subroutine test_tri(mesh, status)
    type(exodus_mesh), intent(out) :: mesh
    integer, intent(inout) :: status

    integer :: side

    mesh%num_dim = 2
    mesh%num_node = 3
    mesh%num_elem = 1
    mesh%num_eblk = 1
    mesh%num_sset = 1
    allocate(mesh%eblk(1), mesh%sset(1))
    mesh%eblk(1)%num_elem = 1
    mesh%eblk(1)%num_nodes_per_elem = 3
    mesh%eblk(1)%elem_type = 'TRI3'
    allocate(mesh%eblk(1)%connect(3,1)); mesh%eblk(1)%connect(:,1) = [1,2,3]
    mesh%sset(1)%num_side = 3
    allocate(mesh%sset(1)%elem(3), mesh%sset(1)%face(3))
    mesh%sset(1)%elem = 1
    mesh%sset(1)%face = [1,2,3]

    sizes => mesh%side_size_list('TRI3')
    call require(all(sizes == [2,2,2]), 'TRI3 side sizes', status)
    do side = 1, 3
      list => mesh%side_node_list(1, side)
      call require(size(list) == 2, 'TRI3 side node count', status)
      deallocate(list)
    end do
    list => mesh%side_set_node_list(1)
    call require(all(list == [1,2,2,3,3,1]), 'TRI3 side-set nodes', status)
    deallocate(list)
  end subroutine


  subroutine require(condition, message, status)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    integer, intent(inout) :: status

    if (.not.condition) then
      print '(a)', 'FAIL: ' // message
      status = 1
    end if
  end subroutine

end program test_polygon_sides
