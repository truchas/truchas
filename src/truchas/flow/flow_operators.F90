!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module flow_operators

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  use unstr_mesh_type
  use flow_domain_types
  use bndry_func_class
  use bndry_vfunc_class
  use index_partitioning
  implicit none
  private

  public :: gradient_cc, gradient_cf, interpolate_fc, interpolate_cf, &
      flow_gradient_coefficients, flow_operators_init

  interface gradient_cf
    module procedure gradient_cf_scalar, gradient_cf_vector
  end interface gradient_cf

  type :: flow_operator
    type(unstr_mesh), pointer :: mesh => null()
    ! ds = dxyz/||dxyz||^2 -- coefficients for computing face gradients
    ! for a face `f` adjacent to cells `i`,`j` oriented such that the face normal
    ! points into cell `i`, the face gradient is given by ds(:,f)*(P(i)-P(j))
    real(r8), pointer :: ds(:,:)
    ! face workspace
    real(r8), allocatable :: work(:)
  end type flow_operator

  type(flow_operator), allocatable :: this

contains

  subroutine flow_operators_init(mesh)
    type(unstr_mesh), intent(in), target :: mesh
    !
    integer :: j, c1, c2

    if (allocated(this)) return

    allocate(this)
    this%mesh => mesh
    allocate(this%ds(3,mesh%nface))
    allocate(this%work(mesh%nface))
    this%ds = 0.0_r8

    do j = 1, mesh%nface_onP
      c1 = mesh%fcell(2,j) ! in
      c2 = mesh%fcell(1,j) ! out
      if (c1 == 0) then
        this%ds(:,j) = mesh%face_centroid(:,j) - mesh%cell_centroid(:,c2)
        this%ds(:,j) = this%ds(:,j)/sum(this%ds(:,j)**2)
      else
        this%ds(:,j) = mesh%cell_centroid(:,c1) - mesh%cell_centroid(:,c2)
        this%ds(:,j) = this%ds(:,j)/sum(this%ds(:,j)**2)
      end if
    end do

    do j = 1, 3
      call gather_boundary(mesh%face_ip, this%ds(j,:))
    end do
  end subroutine flow_operators_init

  function flow_gradient_coefficients() result(p)
    real(r8), pointer :: p(:,:)
    p => this%ds
  end function flow_gradient_coefficients

  ! gradient of cell-centered scalar `x` evaluated at cell centers
  subroutine gradient_cc(gx, gy, gz, x, w_node0, w_node1)
    real(r8), intent(out) :: gx(:), gy(:), gz(:), w_node0(:), w_node1(:)
    real(r8), intent(in) :: x(:)

    integer :: i, j, k

    ! initialization routines are currently handled by flow.  This routine is called by
    ! volume tracker so there is a chance of calling things out of sync.
    INSIST(allocated(this))

    associate (w_face => this%work)
      ! node average data stored in w_node0 workspace array
      call node_avg(x, w_node0, w_node1)
      call gather_boundary(this%mesh%node_ip, w_node0)

      do i = 1, this%mesh%nface
        associate (fn => this%mesh%fnode(this%mesh%xfnode(i):this%mesh%xfnode(i+1)-1))
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

      do i = 1, this%mesh%ncell_onP
        gx(i) = 0.0_r8
        gy(i) = 0.0_r8
        gz(i) = 0.0_r8
        associate (fi => this%mesh%cface(this%mesh%xcface(i):this%mesh%xcface(i+1)-1))
          do j = 1, size(fi)
            k = fi(j)
            if (btest(this%mesh%cfpar(i),pos=j)) then ! true if normal points inward
              ! note that the normal associated with the mesh object has already been scaled
              ! by the face area, so we do not do it again
              gx(i) = gx(i) - this%mesh%normal(1,k)*w_face(k)
              gy(i) = gy(i) - this%mesh%normal(2,k)*w_face(k)
              gz(i) = gz(i) - this%mesh%normal(3,k)*w_face(k)
            else
              gx(i) = gx(i) + this%mesh%normal(1,k)*w_face(k)
              gy(i) = gy(i) + this%mesh%normal(2,k)*w_face(k)
              gz(i) = gz(i) + this%mesh%normal(3,k)*w_face(k)
            end if
          end do
        end associate
        gx(i) = gx(i)/this%mesh%volume(i)
        gy(i) = gy(i)/this%mesh%volume(i)
        gz(i) = gz(i)/this%mesh%volume(i)
      end do
    end associate
  end subroutine gradient_cc


  subroutine node_avg(x_cell, x_node, w_node)

    real(r8), intent(in) :: x_cell(:)
    real(r8), intent(out) :: x_node(:), w_node(:)

    integer :: i, j

    x_node = 0.0_r8
    w_node = 0.0_r8

    do i = 1, this%mesh%ncell
      associate (cn => this%mesh%cnode(this%mesh%xcnode(i):this%mesh%xcnode(i+1)-1))
        do j = 1, size(cn)
          x_node(cn(j)) = x_node(cn(j)) + this%mesh%volume(i)*x_cell(i)
          w_node(cn(j)) = w_node(cn(j)) + this%mesh%volume(i)
        end do
      end associate
    end do

    do i = 1, this%mesh%nnode_onP
      x_node(i) = x_node(i)/w_node(i)
    end do

  end subroutine node_avg

  ! gradient of cell-centered scalar `x` evaluated at face centers
  ! result only valid on nface_onP

  subroutine gradient_cf_scalar(g, x, normal_flux_bc, dirichlet_bc, &
      face_t, non_regular_default, gravity)

    real(r8), intent(out) :: g(:,:)
    real(r8), intent(in) :: x(:)
    class(bndry_func), optional, intent(inout) :: normal_flux_bc, dirichlet_bc
    integer, intent(in), optional :: face_t(:)
    real(r8), intent(in), optional :: non_regular_default
    real(r8), intent(in), optional :: gravity(:,:) ! for dynamic pressure grad

    integer :: j,i

    g = 0.0_r8

    if (present(gravity)) then
      do j = 1, this%mesh%nface_onP
        associate (n => this%mesh%fcell(:,j), y => gravity(:,j))
          if (n(2) > 0) g(:,j) = this%ds(:,j) * (x(n(2))+y(2)-x(n(1))-y(1))
        end associate
      end do
    else
      do j = 1, this%mesh%nface_onP
        associate (n => this%mesh%fcell(:,j))
          if (n(2) > 0) g(:,j) = this%ds(:,j) * (x(n(2))-x(n(1)))
        end associate
      end do
    end if

    if (present(normal_flux_bc)) then
      associate (faces => normal_flux_bc%index, value => normal_flux_bc%value)
        do i = 1, size(faces)
          j = faces(i)
          g(:,j) = value(i)*this%mesh%normal(:,j)/this%mesh%area(j)
        end do
      end associate
    end if

    if (present(dirichlet_bc)) then
      associate (faces => dirichlet_bc%index, value => dirichlet_bc%value)
        if (present(gravity)) then
          do i = 1, size(faces)
            j = faces(i)
            associate(n => this%mesh%fcell(1,j), y => gravity(1,j))
              g(:,j) = this%ds(:,j)*(value(i)-x(n)-y)
            end associate
          end do
        else
          do i = 1, size(faces)
            j = faces(i)
            associate(n => this%mesh%fcell(1,j))
              g(:,j) = this%ds(:,j)*(value(i)-x(n))
            end associate
          end do
        end if
      end associate
    end if


    if (present(face_t) .and. present(non_regular_default)) then
      do j = 1, this%mesh%nface_onP
        if (face_t(j) > regular_t) g(:,j) = non_regular_default
      end do
    end if

  end subroutine gradient_cf_scalar

  ! gradient of cell-centered vector `x` evaluated at face centers
  ! result only valid on nface_onP.  This routine assumes that all
  ! components of the gradient are zero on neumann walls

  subroutine gradient_cf_vector(g, x, zero_normal_bc, dirichlet_bc, &
      face_t, non_regular_default)

    real(r8), intent(out) :: g(:,:,:)
    real(r8), intent(in) :: x(:,:)
    class(bndry_func), optional, intent(in) :: zero_normal_bc
    class(bndry_vfunc), optional, intent(in) :: dirichlet_bc
    integer, intent(in), optional :: face_t(:)
    real(r8), intent(in), optional :: non_regular_default

    real(r8) :: v_normal, v(3)
    integer :: j,i,d

    g = 0.0_r8

    do j = 1, this%mesh%nface_onP
      associate (n => this%mesh%fcell(:,j))
        if (n(2) > 0) then
          do d = 1, size(x,dim=1)
            g(:,d,j) = this%ds(:,j) * (x(d,n(2))-x(d,n(1)))
          end do
        end if
      end associate
    end do

    if (present(zero_normal_bc)) then
      associate (faces => zero_normal_bc%index)
        do i = 1, size(faces)
          j = faces(i)
          associate(n => this%mesh%fcell(1,j))
            ! double check the signs of variables here
            v_normal = dot_product(x(:,n), this%mesh%normal(:,j))/this%mesh%area(j)
            v = x(:,n) - v_normal*this%mesh%normal(:,j)/this%mesh%area(j)
            do d = 1, size(x,dim=1)
              g(:,d,j) = this%ds(:,j)*(v(d)-x(d,n))
            end do
          end associate
        end do
      end associate
    end if

    if (present(dirichlet_bc)) then
      associate (faces => dirichlet_bc%index, value => dirichlet_bc%value)
        do i = 1, size(faces)
          j = faces(i)
          associate(n => this%mesh%fcell(1,j))
            do d = 1, size(x,dim=1)
              g(:,d,j) = this%ds(:,j)*(value(d,i)-x(d,n))
            end do
          end associate
        end do
      end associate
    end if


    if (present(face_t) .and. present(non_regular_default)) then
      do j = 1, this%mesh%nface_onP
        if (face_t(j) > regular_t) g(:,:,j) = non_regular_default
      end do
    end if

  end subroutine gradient_cf_vector

  ! interpolation of vector quantitiy xf to faces
  ! result only valid on ncell_onP
  ! inactive_faces and extra_ignore_faces do not participate in averaging

  subroutine interpolate_fc(ic, xf, face_t, extra_ignore_faces)

    real(r8), intent(out) :: ic(:,:)
    real(r8), intent(in) :: xf(:,:)
    integer, intent(in) , optional :: face_t(:)
    integer, intent(in) , optional :: extra_ignore_faces(:)

    integer :: i, dim
    real(r8) :: tmp

    this%work = 1.0_r8
    if (present(face_t)) then
      do i = 1, this%mesh%nface
        if (face_t(i) > regular_t) this%work(i) = 0.0_r8
      end do
    end if
    if (present(extra_ignore_faces)) then
      this%work(extra_ignore_faces) = 0.0_r8
    end if

    do i = 1, this%mesh%ncell_onP
      associate (fn => this%mesh%cface(this%mesh%xcface(i):this%mesh%xcface(i+1)-1))
        if (sum(this%work(fn)) >= 1.0_r8) then
          do dim = 1, size(ic,dim=1)
            tmp = sum(this%work(fn)*abs(this%mesh%normal(dim,fn)))
            if (tmp /= 0.0_r8) then
              ic(dim,i) = dot_product(this%work(fn)*abs(this%mesh%normal(dim,fn)), xf(dim,fn)) / tmp
            else
              ic(dim,i) = 0.0_r8
            end if
          end do
        else
          ic(:,i) = 0.0_r8
        end if
      end associate
    end do

  end subroutine interpolate_fc

  ! normal component of interpolation of vector cell centered vector
  ! quantity 'x' to face centers, using weights 'w'.  Dirichlet
  ! boundary conditions on faces may be supplied as well a list of
  ! inactive faces and a default value for xf at inactive faces.
  ! result only valid on nface_onP

  subroutine interpolate_cf(xf, x, w, bc, bc_zn, face_t, non_regular_default)

    real(r8), intent(out) :: xf(:)
    real(r8), intent(in) :: x(:,:), w(:,:)
    class(bndry_vfunc), optional, intent(in) :: bc
    class(bndry_func), optional, intent(in) :: bc_zn
    integer, intent(in), optional :: face_t(:)
    real(r8), intent(in), optional :: non_regular_default

    integer :: j,i

    do j = 1, this%mesh%nface_onP
      associate (n => this%mesh%fcell(:,j))
        if (n(2) > 0) then
          xf(j) = (w(1,j)*dot_product(this%mesh%normal(:,j),x(:,n(1))) + &
              w(2,j)*dot_product(this%mesh%normal(:,j),x(:,n(2)))) / this%mesh%area(j)
        else
          xf(j) = dot_product(this%mesh%normal(:,j),x(:,n(1)))/this%mesh%area(j)
        end if
      end associate
    end do

    if (present(bc)) then
      associate (faces => bc%index, value => bc%value)
        do i = 1, size(faces)
          j = faces(i)
          xf(j) = dot_product(this%mesh%normal(:,j),value(:,i))/this%mesh%area(j)
        end do
      end associate
    end if

    if (present(bc_zn)) then
      associate (faces => bc_zn%index)
        do i = 1, size(faces)
          j = faces(i)
          xf(j) = 0.0_r8
        end do
      end associate
    end if

    if (present(face_t) .and. present(non_regular_default)) then
      do j = 1, this%mesh%nface_onP
        if (face_t(j) > regular_t) xf(j) = non_regular_default
      end do
    end if

  end subroutine interpolate_cf

end module flow_operators
