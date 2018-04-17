module flow_operators

  use kinds, only: r8
  use truchas_logging_services
  use unstr_mesh
  use flow_mesh_type
  use bndry_func_class
  use index_partitioning
  implicit none
  private

  public :: gradient_cc, gradient_fc, interpolate_fc, interpolate_cf

  type :: flow_operator
    type(flow_mesh), pointer :: mesh => null()
    real(r8), allocatable, target :: ds(:,:) ! dxyz/||dxyz||^2
  end type flow_operator

  type(flow_operator), allocatable :: this

contains

  subroutine flow_operators_init(mesh)
    type(flow_mesh), pointer, intent(in) :: mesh
    !
    type(unstr_mesh), pointer :: m
    integer :: j, c1, c2

    if (allocated(this)) return

    allocate(this)
    this%mesh => mesh
    m => mesh%mesh
    allocate(ds(3, size(m%nface)))

    this%ds = 0.0_r8

    do j = 1, m%nface_onP
      c1 = mesh%fcell(1,j)
      c2 = mesh%fcell(2,j)
      if (c1 == 0) then
        this%ds(:,j) = mesh%face_centroid(:,j) - mesh%cell_centroid(:,c2)
        this%ds(:,j) = this%ds(:,j)/sum(this%ds(:,j)**2)
      else
        this%ds(:,j) = mesh%cell_centroid(:,c1) - mesh%cell_centroid(:,c2)
        this%ds(:,j) = this%ds(:,j)/sum(this%ds(:,j)**2)
      end if
    end do

    do j = 1, 3
      call gather_boundary(m%face_ip, this%ds(j,:))
    end do
  end subroutine flow_operators_init

  function flow_gradient_coefficients() result(p)
    real(r8), pointer :: p(:,:)
    p => this%ds
  end function flow_gradient_coefficients

  ! gradient of cell-centered scalar `x` evaluated at cell centers
  subroutine gradient_cc(gx, gy, gz, x, w_node0, w_node1, w_face)
    real(r8), intent(out) :: gx(:), gy(:), gz(:), w_node0(:), w_node1(:), w_face(:)
    real(r8), intent(in) :: x(:)

    type(unstr_mesh), pointer :: m
    integer :: i, j, k

    m => this%mesh%mesh
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

    do i = 1, m%ncell_onP
      gx(i) = 0.0_r8
      gy(i) = 0.0_r8
      gz(i) = 0.0_r8
      associate (fi => m%cface(m%xcface(i):m%xcface(i+1)-1))
        do j = 1, size(fi)
          k = fi(j)
          if (btest(m%cfpar(i),pos=j)) then ! true if normal points inward
            ! note the thae normal associate with the mesh object has already been scaled
            ! by the face area, so we do not do it again
            gx(i) = gx(i) - m%normal(1,k)*w_face(k)
            gy(i) = gy(i) - m%normal(2,k)*w_face(k)
            gz(i) = gz(i) - m%normal(3,k)*w_face(k)
          else
            gx(i) = gx(i) + m%normal(1,k)*w_face(k)
            gy(i) = gy(i) + m%normal(2,k)*w_face(k)
            gz(i) = gz(i) + m%normal(3,k)*w_face(k)
          end if
        end do
      end associate
      gx(i) = gx(i)/m%volume(i)
      gy(i) = gy(i)/m%volume(i)
      gz(i) = gz(i)/m%volume(i)
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


  ! gradient of cell-centered scalar `x` evaluated at face centers
  ! result only valid on nface_onP
  subroutine gradient_cf(g, x, t, normal_flux_bc, dirichlet_bc, &
      inactive_faces, inactive_default, gravity)
    real(r8), intent(out) :: g(:,:)
    real(r8), intent(in) :: x(:), t
    class(bndry_func), optional, intent(inout) :: normal_flux_bc, dirichlet_bc
    integer, intent(in), optional :: inactive_faces(:)
    real(r8), intent(in), optional :: inactive_default
    real(r8), intent(in), optional :: gravity(:,:) ! for dynamic pressure grad

    integer :: j,i

    g = 0.0_r8

    associate (m => this%mesh%mesh, f=> this%mesh)
      if (present(gravity)) then
        do j = 1, m%nface_onP
          associate (n => f%fcell(:,j), y => gravity(:,j))
            if (n(1) > 0) g(:,j) = this%ds(:,j) * (x(n(1))+y(1)-x(n(2))-y(2))
          end associate
        end do
      else
        do j = 1, m%nface_onP
          associate (n => f%fcell(:,j))
            if (n(1) > 0) g(:,j) = this%ds(:,j) * (x(n(1))-x(n(2)))
          end associate
        end do
      end if

      if (present(normal_flux_bc)) then
        call normal_flux_bc%compute(t)
        associate (index => normal_flux_bc%index, value => normal_flux_bc%value)
          g(:,index) = value*m%normal(:,index)/m%area(index)
        end associate
      end if

      if (present(dirichlet_bc)) then
        call dirichlet_bc%compute(t)
        associate (index => dirichlet_bc%index, value => dirichlet_bc%value)
          if (present(gravity)) then
            do i = 1, size(index)
              j = index(i)
              associate(n => f%fcell(2,j), y => gravity(2,j))
                g(:,j) = this%ds(:,j)*(value(i)-x(n)-y)
              end associate
            end do
          else
            do i = 1, size(index)
              j = index(i)
              associate(n => f%fcell(2,j))
                g(:,j) = this%ds(:,j)*(value(i)-x(n))
              end associate
            end do
          end if
        end associate
      end if


      if (present(inactive_faces) .and. present(inactive_default)) then
        do j = 1, m%nface_onP
          if (inactive_faces(j) > 0) g(:,j) = inactive_default
        end do
      end if

    end associate

  end subroutine gradient_cf

  ! interpolation of vector quantitiy xf to faces
  ! result only valid on ncell_onP
  subroutine interpolate_fc(ic, xf)
    real(r8), intent(out) :: ic(:,:)
    real(r8), intent(in) :: xf(:,:)

    integer :: i, dim


    associate (m => this%mesh%mesh, f=> this%mesh)
      do i = 1, m%ncell_onP
        associate (fn => m%cface(m%xcface(i):m%xcface(i+1)-1))
          do dim = 1, size(ic,dim=1)
            g(dim,i) = dot_product(abs(m%normal(dim,fn)), xf(dim,fn))/sum(abs(m%normal(dim,fn)))
          end do
        end associate
      end do
    end associate
  end subroutine interpolate_fc


  ! normal component of interpolation of vector cell centered vector
  ! quantity 'x' to face centers, using weights 'w'.  Dirichlet
  ! boundary conditions on faces may be supplied as well a list of
  ! inactive faces and a default value for xf at inactive faces.
  ! result only valid on nface_onP
  subroutine interpolate_cf(xf, x, w, bc, inactive_faces, inactive_default)
    real(r8), intent(out) :: xf(:)
    real(r8), intent(in) :: x(:,:), w(:,:)
    class(bndry_func), optional, intent(inout) :: bc
    integer, intent(in), optional :: inactive_faces(:)
    real(r8), intent(in), optional :: inactive_default

    integer :: j,i

    associate (m => this%mesh%mesh, f=> this%mesh)
      do j = 1, m%nface_onP
        associate (n => f%fcell(:,j))
          if (n(1) > 0) then
            xf(j) = (w(1,j)*dot_product(m%normal(:,j),x(i,n(1))) + &
                w(2,j)*dot_product(m%normal(:,j),x(i,n(2)))) / m%area(j)
          else
            xf(j) = dot_product(m%normal(:,j),xf(:,n(2)))/m%area(j)
          end if
        end associate
      end do

      if (present(bc)) then
        associate (index => bc%index, value => bc%value)
          do i = 1, size(index)
            j = index(i)
            xf(j) = dot_product(m%normal(:,j),value(:,j))/m%area(j)
          end do
        end associate
      end if

      if (present(inactive_faces) .and. present(inactive_default)) then
        do j = 1, m%nface_onP
          if (inactive_faces(j) > 0) xf(j) = inactive_default
        end do
      end if
    end associate
  end subroutine interpolate_cf


end module flow_operators
