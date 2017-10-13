module flow_operators

  use unstr_mesh_type
  implicit none
  private

contains

  ! gradient of cell-centered scalar `x` evaluated at cell centers
  subroutine gradient_cc(m, gx, gy, gz, x, w_node0, w_node1, w_face)
    type(unstr_mesh), intent(in) :: m
    real(r8), intent(out) :: gx(:), gy(:), gz(:), w_node(:), w_face(:)
    real(r8), intent(in) :: x(:)

    integer :: i, j

    ! node average data stored in w_node0 workspace array
    call node_avg(m, x, w_node0, w_node1)
    call gather_boundary(m%node_ip, w_node0)

    do i = 1, m%nface
      associate (fn => m%fnode(m%xfnode(i):m%xfnode(i+1)-1))
        select case (size(fn))
        case (3,4)
          ! linear interpolation of vertex values to face centroid is arithmetic mean
          ! for triangle and quadrilateral faces
          w_face(i) = sum(w_node0(fn))/real(size(fn),r8)
        case default
          call TLS_fatal("wrong number of faces in gradient_cc")
        end select
      end associate
    end do


  end subroutine gradient_cc


  subroutine node_avg(m, x_cell, x_node, w_node)
    type(unstr_mesh), intent(in) :: m
    real(r8), intent(in) :: x_cell(:)
    real(r8), intent(out) :: x_node(:), w_node(:)

    integer :: i, j

    x_node = 0.0_r8
    w_node = 0.0_r8

    do i = 1, m%ncell
      associate (cn => m%cnode(m%xcnode(i):m%xcnode(i+1)-1))
        do j = 1, size(cn)
          x_node(cn(j)) = x_node(cn(j)) + m%volume(i)*x_cell(i)
          w_node(cn(j)) = w_node(cn(j)) + m%volume(i)
        end do
      end associate
    end do

    do i = 1, m%nnode_onP
      x_node(i) = x_node(i)/w_node(i)
    end do
  end subroutine node_avg

end module flow_operators
