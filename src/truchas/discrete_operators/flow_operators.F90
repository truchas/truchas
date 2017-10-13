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
        case (3)
          w_face(i) = lerp3(m%x(:,fn), w_node0(fn), sum(m%x(:,fn),dim=2)/3.0_r8)
        case (4)

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


  function lerp3(x, v, p)
    real(r8), intent(in) :: x(:,:), v(:), p(:)
    real(r8) :: lerp3

    real(r8) :: t0(3), t1(3), s(3), c0, c1

    t0 = x(:,2) - x(:,1)
    t1 = x(:,3) - x(:,1)
    s = p - x(:,1)

    ! write s in terms of t0, t1
    c0 = (s(1)*(t0(1)*t1(2)**2 + t0(1)*t1(3)**2 - t0(2)*t1(1)*t1(2) - t0(3)*t1(1)*t1(3)) &
        + s(2)*(-t0(1)*t1(1)*t1(2) + t0(2)*t1(1)**2 + t0(2)*t1(3)**2 - t0(3)*t1(2)*t1(3)) &
        + s(3)*(-t0(1)*t1(1)*t1(3) - t0(2)*t1(2)*t1(3) + t0(3)*t1(1)**2 + t0(3)*t1(2)**2))/ &
        (t0(1)**2*t1(2)**2 + t0(1)**2*t1(3)**2 - 2*t0(1)*t0(2)*t1(1)*t1(2) - &
        2*t0(1)*t0(3)*t1(1)*t1(3) + t0(2)**2*t1(1)**2 + t0(2)**2*t1(3)**2 - &
        2*t0(2)*t0(3)*t1(2)*t1(3) + t0(3)**2*t1(1)**2 + t0(3)**2*t1(2)**2)
    c1 = -(s(1)*(t0(1)*t0(2)*t1(2) + t0(1)*t0(3)*t1(3) - t0(2)**2*t1(1) - t0(3)**2*t1(1)) &
        + s(2)*(-t0(1)**2*t1(2) + t0(1)*t0(2)*t1(1) + t0(2)*t0(3)*t1(3) - t0(3)**2*t1(2)) &
        + s(3)*(-t0(1)**2*t1(3) + t0(1)*t0(3)*t1(1) - t0(2)**2*t1(3) + t0(2)*t0(3)*t1(2))) / &
        (t0(1)**2*t1(2)**2 + t0(1)**2*t1(3)**2 - 2*t0(1)*t0(2)*t1(1)*t1(2) - &
        2*t0(1)*t0(3)*t1(1)*t1(3) + t0(2)**2*t1(1)**2 + t0(2)**2*t1(3)**2 - &
        2*t0(2)*t0(3)*t1(2)*t1(3) + t0(3)**2*t1(1)**2 + t0(3)**2*t1(2)**2)

    lerp3 = c0*v(0) + c0*v(1) - c1*v(0) + c1*v(2) + v(0)

  end function lerp3
end module flow_operators
