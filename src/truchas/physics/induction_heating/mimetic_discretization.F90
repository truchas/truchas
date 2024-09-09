!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module mimetic_discretization

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simpl_mesh_type
  implicit none
  private

  !! Discrete natural operators and their transposes
  public :: grad, grad_t, curl, curl_t, div, div_t

  !! Inner product matrices
  public :: W1_matrix_HS, W2_matrix_HS
  public :: W1_matrix_WE, W2_matrix_WE
  public :: w1_face_matrix

  !! Interpolation procedures
  public :: w1_vector_on_cells, w2_vector_on_cells, w3_scalar_on_cells
  public :: eval_w0_interp_coef, eval_w1_interp_coef, eval_w2_interp_coef, eval_w3_interp_coef

  !! Gradient operator on a tetrahedral cell
  real(r8), parameter, public :: cell_grad(6,4) = &
      reshape([-1, -1, -1,  0,  0,  0, &
                1,  0,  0, -1, -1,  0, &
                0,  1,  0,  1,  0, -1, &
                0,  0,  1,  0,  1,  1], shape(cell_grad))

  !! Curl operator on a tetrahedral cell
  real(r8), parameter, public :: cell_curl(4,6) = &
      reshape([0,  0,  1,  1, &
               0,  1,  0, -1, &
               0, -1, -1,  0, &
               1,  0,  0,  1, &
              -1,  0,  1,  0, &
               1,  1,  0,  0], shape=shape(cell_curl))

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! THE DISCRETE NATURAL OPERATORS:
 !!
 !!      Grad    Curl    Div
 !!   W0 ---> W1 ---> W2 --> W3
 !!

  subroutine grad (mesh, u, v, increment, mask)

    type(simpl_mesh), intent(in) :: mesh
    real(kind=r8),   intent(in) :: u(:)
    real(kind=r8), intent(inout) :: v(:)
    logical, intent(in), optional :: increment, mask(:)

    integer :: j
    logical :: increment_

    ASSERT( size(u) == mesh%nnode )
    ASSERT( size(v) == mesh%nedge )
    ASSERT( allocated(mesh%enode) )

    increment_ = .false.
    if (present(increment)) increment_ = increment

    if (present(mask)) then
      do j = 1, mesh%nedge
        if (.not.increment_) v(j) = 0
        if (mask(j)) v(j) = v(j) + u(mesh%enode(2,j)) - u(mesh%enode(1,j))
      end do
    else
      do j = 1, mesh%nedge
        if (.not.increment_) v(j) = 0
        v(j) = v(j) + u(mesh%enode(2,j)) - u(mesh%enode(1,j))
      end do
    end if

  end subroutine grad

  subroutine grad_t (mesh, u, v, local, mask)

    type(simpl_mesh), intent(in) :: mesh
    real(kind=r8),   intent(in) :: u(:)
    real(kind=r8), intent(out) :: v(:)
    logical, intent(in), optional :: local
    logical, intent(in), optional :: mask(:)

    integer :: j

    ASSERT( size(u) == mesh%nedge )
    ASSERT( size(v) == mesh%nnode )
    ASSERT( allocated(mesh%enode) )

    if (present(mask)) then
      v = 0.0_r8
      do j = 1, mesh%nedge
        if (.not.mask(j)) cycle
        v(mesh%enode(1,j)) = v(mesh%enode(1,j)) - u(j)
        v(mesh%enode(2,j)) = v(mesh%enode(2,j)) + u(j)
      end do
    else
      v = 0.0_r8
      do j = 1, mesh%nedge
        v(mesh%enode(1,j)) = v(mesh%enode(1,j)) - u(j)
        v(mesh%enode(2,j)) = v(mesh%enode(2,j)) + u(j)
      end do
    end if

    if (present(local)) then
      if (local) return
    end if

    call mesh%node_imap%gather_offp(v)

  end subroutine grad_t

  function curl (mesh, u) result (v)

    type(simpl_mesh), intent(in) :: mesh
    real(kind=r8),   intent(in) :: u(:)
    real(kind=r8) :: v(mesh%nface)

    integer :: j

    ASSERT( size(u) == mesh%nedge )
    ASSERT( allocated(mesh%fedge) )

    do j = 1, mesh%nface
      v(j) = u(mesh%fedge(1,j)) - u(mesh%fedge(2,j)) + u(mesh%fedge(3,j))
    end do

  end function curl

  function curl_t (mesh, u) result (v)

    type(simpl_mesh), intent(in) :: mesh
    real(kind=r8),   intent(in) :: u(:)
    real(kind=r8) :: v(mesh%nedge)

    integer :: j

    ASSERT( size(u) == mesh%nface )
    ASSERT( allocated(mesh%fedge) )

    v = 0.0_r8
    do j = 1, mesh%nface
      v(mesh%fedge(1,j)) = v(mesh%fedge(1,j)) + u(j)
      v(mesh%fedge(2,j)) = v(mesh%fedge(2,j)) - u(j)
      v(mesh%fedge(3,j)) = v(mesh%fedge(3,j)) + u(j)
    end do

    call mesh%edge_imap%gather_offp(v)

  end function curl_t

  function div (mesh, u) result (v)

    type(simpl_mesh), intent(in) :: mesh
    real(kind=r8),   intent(in) :: u(:)
    real(kind=r8) :: v(mesh%ncell)

    integer :: j

    ASSERT( size(u) == mesh%nface )
    ASSERT( allocated(mesh%cface) )

    do j = 1, mesh%ncell
      v(j) = u(mesh%cface(1,j)) - u(mesh%cface(2,j)) + u(mesh%cface(3,j)) - u(mesh%cface(4,j))
    end do

  end function div

  function div_t (mesh, u) result (v)

    type(simpl_mesh), intent(in) :: mesh
    real(kind=r8),   intent(in) :: u(:)
    real(kind=r8) :: v(mesh%nface)

    integer :: j

    ASSERT( size(u) == mesh%ncell )
    ASSERT( allocated(mesh%cface) )

    v = 0.0_r8
    do j = 1, mesh%ncell
      v(mesh%cface(1,j)) = v(mesh%cface(1,j)) + u(j)
      v(mesh%cface(2,j)) = v(mesh%cface(2,j)) - u(j)
      v(mesh%cface(3,j)) = v(mesh%cface(3,j)) + u(j)
      v(mesh%cface(4,j)) = v(mesh%cface(4,j)) - u(j)
    end do

    call mesh%face_imap%gather_offp(v)

  end function div_t

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! W1_MATRIX_HS
 !!
 !! Returns the Hyman-Shashkov form of the W2 inner product matrix for a cell.
 !! Result is a symmetric 6x6 matrix is stored in upper-packed storage mode:
 !! the upper triangle packed by column into a one-dimensional array.
 !!

  function W1_matrix_HS (mesh, cell) result (matrix)

    use simplex_geometry, only: tet_face_normal

    type(simpl_mesh), intent(in) :: mesh
    integer,         intent(in) :: cell
    real(kind=r8) :: matrix(21)

    real(kind=r8) :: c, p(3,4), pp(4,4)

    !ASSERT( defined(mesh) )
    ASSERT( (cell >= 1) .and. (cell <= mesh%ncell) )

    p = tet_face_normal(mesh%x(:,mesh%cnode(:,cell)))
    pp = matmul(transpose(p),p)  ! inner products (WASTEFUL-FIXME)
    c = 1.0_r8 / (36.0_r8 * abs(mesh%volume(cell)))

    matrix( 1) = c*(pp(1,1) + pp(2,2))
    matrix( 3) = c*(pp(1,1) + pp(3,3))
    matrix( 6) = c*(pp(1,1) + pp(4,4))
    matrix(10) = c*(pp(2,2) + pp(3,3))
    matrix(15) = c*(pp(2,2) + pp(4,4))
    matrix(21) = c*(pp(3,3) + pp(4,4))
    matrix( 2) = c*pp(2,3)
    matrix( 4) = c*pp(2,4)
    matrix( 7) = c*(-pp(1,3))
    matrix(11) = c*(-pp(1,4))
    matrix(16) = 0.0_r8
    matrix( 5) = c*pp(3,4)
    matrix( 8) = c*pp(1,2)
    matrix(12) = 0.0_r8
    matrix(17) = c*(-pp(1,4))
    matrix( 9) = 0.0_r8
    matrix(13) = c*pp(1,2)
    matrix(18) = c*pp(1,3)
    matrix(14) = c*pp(3,4)
    matrix(19) = c*(-pp(2,4))
    matrix(20) = c*pp(2,3)

  end function W1_matrix_HS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! W1_MATRIX_WE
 !!
 !! Returns the Whitney element form of the W2 inner product matrix for a cell.
 !! Result is a symmetric 6x6 matrix is stored in upper-packed storage mode:
 !! the upper triangle packed by column into a one-dimensional array.
 !!

  function W1_matrix_WE (mesh, cell) result (matrix)

    use simplex_geometry, only: tet_face_normal

    type(simpl_mesh), intent(in) :: mesh
    integer,         intent(in) :: cell
    real(kind=r8) :: matrix(21)

    real(kind=r8) :: c, p(3,4), pp(4,4)

    !ASSERT( defined(mesh) )
    ASSERT( (cell >= 1) .and. (cell <= mesh%ncell) )

    p = tet_face_normal(mesh%x(:,mesh%cnode(:,cell)))
    pp = matmul(transpose(p),p)  ! inner products (WASTEFUL-FIXME)
    c = 1.0_r8 / (90.0_r8 * abs(mesh%volume(cell)))

    matrix( 1) = c*(pp(1,1) - pp(1,2) + pp(2,2))
    matrix( 3) = c*(pp(1,1) - pp(1,3) + pp(3,3))
    matrix( 6) = c*(pp(1,1) - pp(1,4) + pp(4,4))
    matrix(10) = c*(pp(2,2) - pp(2,3) + pp(3,3))
    matrix(15) = c*(pp(2,2) - pp(2,4) + pp(4,4))
    matrix(21) = c*(pp(3,3) - pp(3,4) + pp(4,4))
    matrix( 2) = c*(0.5_r8*(pp(1,1) - pp(1,2) - pp(1,3)) + pp(2,3))
    matrix( 4) = c*(0.5_r8*(pp(1,1) - pp(1,2) - pp(1,4)) + pp(2,4))
    matrix( 7) = c*(0.5_r8*(pp(1,2) - pp(2,2) + pp(2,3)) - pp(1,3))
    matrix(11) = c*(0.5_r8*(pp(1,2) - pp(2,2) + pp(2,4)) - pp(1,4))
    matrix(16) = c*(0.5_r8*(pp(1,3) - pp(2,3) - pp(1,4) + pp(2,4)))
    matrix( 5) = c*(0.5_r8*(pp(1,1) - pp(1,3) - pp(1,4)) + pp(3,4))
    matrix( 8) = c*(0.5_r8*(pp(3,3) - pp(1,3) - pp(2,3)) + pp(1,2))
    matrix(12) = c*(0.5_r8*(pp(1,2) - pp(2,3) - pp(1,4) + pp(3,4)))
    matrix(17) = c*(0.5_r8*(pp(1,3) - pp(3,3) + pp(3,4)) - pp(1,4))
    matrix( 9) = c*(0.5_r8*(pp(1,2) - pp(1,3) - pp(2,4) + pp(3,4)))
    matrix(13) = c*(0.5_r8*(pp(4,4) - pp(1,4) - pp(2,4)) + pp(1,2))
    matrix(18) = c*(0.5_r8*(pp(4,4) - pp(1,4) - pp(3,4)) + pp(1,3))
    matrix(14) = c*(0.5_r8*(pp(2,2) - pp(2,3) - pp(2,4)) + pp(3,4))
    matrix(19) = c*(0.5_r8*(pp(2,3) - pp(3,3) + pp(3,4)) - pp(2,4))
    matrix(20) = c*(0.5_r8*(pp(4,4) - pp(2,4) - pp(3,4)) + pp(2,3))

  end function W1_matrix_WE

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! W2_MATRIX_HS
 !!
 !! Returns the Hyman-Shashkov form of the W2 inner product matrix for a cell.
 !! Result is a symmetric 4x4 matrix is stored in upper-packed storage mode:
 !! the upper triangle packed by column into a one-dimensional array.
 !!

  function W2_matrix_HS (mesh, cell) result (matrix)

    type(simpl_mesh),  intent(in)  :: mesh
    integer,          intent(in)  :: cell
    real(kind=r8) :: matrix(10)

    real(kind=r8) :: c, lsq(6)

    !ASSERT( defined(mesh) )
    ASSERT( (cell >= 1) .and. (cell <= mesh%ncell) )

    lsq = mesh%length(mesh%cedge(:,cell))**2
    c = 1.0_r8 / (36.0_r8 * abs(mesh%volume(cell)))

    matrix( 1) = c*(lsq(1) + lsq(2) + lsq(3))
    matrix( 3) = c*(lsq(1) + lsq(4) + lsq(5))
    matrix( 6) = c*(lsq(2) + lsq(4) + lsq(6))
    matrix(10) = c*(lsq(3) + lsq(5) + lsq(6))
    matrix( 2) = -c*(0.5_r8*(lsq(2) + lsq(3) + lsq(4) + lsq(5)) - lsq(1))
    matrix( 4) =  c*(0.5_r8*(lsq(1) + lsq(3) + lsq(4) + lsq(6)) - lsq(2))
    matrix( 7) = -c*(0.5_r8*(lsq(1) + lsq(2) + lsq(5) + lsq(6)) - lsq(3))
    matrix( 5) = -c*(0.5_r8*(lsq(1) + lsq(2) + lsq(5) + lsq(6)) - lsq(4))
    matrix( 8) =  c*(0.5_r8*(lsq(1) + lsq(3) + lsq(4) + lsq(6)) - lsq(5))
    matrix( 9) = -c*(0.5_r8*(lsq(2) + lsq(3) + lsq(4) + lsq(5)) - lsq(6))

  end function W2_matrix_HS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! W2_MATRIX_WE
 !!
 !! Returns the Whitney element form of the W2 inner product matrix for a cell.
 !! Result is a symmetric 4x4 matrix is stored in upper-packed storage mode:
 !! the upper triangle packed by column into a one-dimensional array.
 !!

  function W2_matrix_WE (mesh, cell) result (matrix)

    type(simpl_mesh),  intent(in)  :: mesh
    integer,          intent(in)  :: cell
    real(kind=r8) :: matrix(10)

    real(kind=r8) :: c, lsq(6)

    !ASSERT( defined(mesh) )
    ASSERT( (cell >= 1) .and. (cell <= mesh%ncell) )

    lsq = mesh%length(mesh%cedge(:,cell))**2
    c = 1.0_r8 / (90.0_r8 * abs(mesh%volume(cell)))

    matrix( 1) = c*(2.0_r8*(lsq(1)+lsq(2)+lsq(3)) - 0.5_r8*(lsq(4)+lsq(5)+lsq(6)))
    matrix( 3) = c*(2.0_r8*(lsq(1)+lsq(4)+lsq(5)) - 0.5_r8*(lsq(2)+lsq(3)+lsq(6)))
    matrix( 6) = c*(2.0_r8*(lsq(2)+lsq(4)+lsq(6)) - 0.5_r8*(lsq(1)+lsq(3)+lsq(5)))
    matrix(10) = c*(2.0_r8*(lsq(3)+lsq(5)+lsq(6)) - 0.5_r8*(lsq(1)+lsq(2)+lsq(4)))
    matrix( 2) = -c*(0.75_r8*(lsq(2)+lsq(3)+lsq(4)+lsq(5)) - 3.0_r8*lsq(1) - 0.5_r8*lsq(6))
    matrix( 4) =  c*(0.75_r8*(lsq(1)+lsq(3)+lsq(4)+lsq(6)) - 3.0_r8*lsq(2) - 0.5_r8*lsq(5))
    matrix( 7) = -c*(0.75_r8*(lsq(1)+lsq(2)+lsq(5)+lsq(6)) - 3.0_r8*lsq(3) - 0.5_r8*lsq(4))
    matrix( 5) = -c*(0.75_r8*(lsq(1)+lsq(2)+lsq(5)+lsq(6)) - 3.0_r8*lsq(4) - 0.5_r8*lsq(3))
    matrix( 8) =  c*(0.75_r8*(lsq(1)+lsq(3)+lsq(4)+lsq(6)) - 3.0_r8*lsq(5) - 0.5_r8*lsq(2))
    matrix( 9) = -c*(0.75_r8*(lsq(2)+lsq(3)+lsq(4)+lsq(5)) - 3.0_r8*lsq(6) - 0.5_r8*lsq(1))

  end function W2_matrix_WE

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! VECTOR FIELD RECOVERY PROCEDURES
 !!

  function w1_vector_on_cells (mesh, u) result (v)

    use simplex_geometry, only: tet_face_normal

    type(simpl_mesh), intent(in) :: mesh
    real(kind=r8),   intent(in) :: u(:)
    real :: v(3,mesh%ncell)

    integer :: j, k, n, m
    real(kind=r8) :: p(3,4), t(3)

    ASSERT( size(u) == mesh%nedge )

    do j = 1, mesh%ncell

      !! Area-weighted outer face normals
      p = tet_face_normal(mesh%x(:,mesh%cnode(:,j)))

      t = 0.0_r8
      k = 0
      do m = 1, 4       !! NB: run through edges mn in order (index k)
        do n = m+1, 4
          k = k + 1
          t = t + (p(:,m) - p(:,n)) * u(mesh%cedge(k,j))
        end do
      end do
      v(:,j) = t / (12.0_r8 * mesh%volume(j))

    end do

  end function w1_vector_on_cells


  function w2_vector_on_cells (mesh, u) result (v)

    type(simpl_mesh), intent(in) :: mesh
    real(kind=r8),   intent(in) :: u(:)
    real :: v(3,mesh%ncell)

    integer :: j, m, n
    real(kind=r8) :: f(4), t(3), x(3,4)

    ASSERT( size(u) == mesh%nface )

    do j = 1, mesh%ncell

      !! Face fluxes (all outward or all inward depending on sign of volume)
      f = u(mesh%cface(:,j))
      f(2) = -f(2)
      f(4) = -f(4)

      x = mesh%x(:,mesh%cnode(:,j))

      t = 0.0_r8
      do m = 1, 4
        do n = m+1, 4
          t = t + (f(m) - f(n)) * (x(:,n) - x(:,m))
        end do
      end do
      v(:,j) = t / (12.0_r8 * mesh%volume(j))

    end do

  end function w2_vector_on_cells


  function w3_scalar_on_cells (mesh, u) result (v)

    type(simpl_mesh), intent(in) :: mesh
    real(kind=r8),   intent(in) :: u(:)
    real :: v(mesh%ncell)

    ASSERT( size(u) == mesh%ncell )

    v = u / mesh%volume

  end function w3_scalar_on_cells

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! LINEAR INTERPOLATION PROCEDURES
 !!
 !! These procedures compute the multipliers for linear interpolation. It is
 !! unclear where these procedures belong.  These are somewhat generic
 !! procedures that could be grouped with general mesh support  routines,
 !! except for the need to account for the global orientation of the face/edge
 !! values which is specific to the discretization.  On the other hand the
 !! determination of the cell to which the interpolation point belongs,
 !! definitely does not belong here.  This is why the cell is passed as an
 !! argument in addition to the point.
 !!

  subroutine eval_w0_interp_coef (mesh, point, cell, coef, index)

    use simplex_geometry, only: bary_coord

    type(simpl_mesh), intent(in)  :: mesh
    real(kind=r8),   intent(in)  :: point(:)
    integer,         intent(in)  :: cell
    real(kind=r8),   intent(out) :: coef(:)
    integer,         intent(out) :: index(:)

    ASSERT( size(point) == 3 )
    ASSERT( cell > 0 .and. cell <= mesh%ncell )
    ASSERT( size(coef) == 4 .and. size(index) == 4 )

    coef = bary_coord(mesh%x(:,mesh%cnode(:,cell)), point)
    index = mesh%cnode(:,cell)

  end subroutine eval_w0_interp_coef


  subroutine eval_w1_interp_coef (mesh, point, cell, coef, index)

    use simplex_geometry, only: tet_face_normal, bary_coord

    type(simpl_mesh), intent(in)  :: mesh
    real(kind=r8),   intent(in)  :: point(:)    ! interpolation point
    integer,         intent(in)  :: cell        ! index of the containing cell
    real(kind=r8),   intent(out) :: coef(:,:)   ! interpolation coefficients
    integer,         intent(out) :: index(:)    ! interpolation indices

    real(kind=r8) :: p(3,4), lambda(4)

    ASSERT( size(point) == 3 )
    ASSERT( cell >= 1 .and. cell <= mesh%ncell )
    ASSERT( all(shape(coef) == (/3,6/)) )
    ASSERT( size(index) == 6 )

    !! Area-weighted face normals (outer if volume > 0, inner if volume < 0)
    p = tet_face_normal(mesh%x(:,mesh%cnode(:,cell)))

    !! Barycentric coordinates (scaled) of the interpolation point
    lambda = bary_coord(mesh%x(:,mesh%cnode(:,cell)), point) / (3.0_r8 * mesh%volume(cell))

    !! Local interpolation coefficients
    coef(:,1) = lambda(2) * p(:,1) - lambda(1) * p(:,2)
    coef(:,2) = lambda(3) * p(:,1) - lambda(1) * p(:,3)
    coef(:,3) = lambda(4) * p(:,1) - lambda(1) * p(:,4)
    coef(:,4) = lambda(3) * p(:,2) - lambda(2) * p(:,3)
    coef(:,5) = lambda(4) * p(:,2) - lambda(2) * p(:,4)
    coef(:,6) = lambda(4) * p(:,3) - lambda(3) * p(:,4)

    index = mesh%cedge(:,cell)

  end subroutine eval_w1_interp_coef


  subroutine eval_w2_interp_coef (mesh, point, cell, coef, index)

    use simplex_geometry, only: bary_coord

    type(simpl_mesh), intent(in)  :: mesh
    real(kind=r8),   intent(in)  :: point(:)    ! interpolation point
    integer,         intent(in)  :: cell        ! index of the containing cell
    real(kind=r8),   intent(out) :: coef(:,:)   ! interpolation coefficients
    integer,         intent(out) :: index(:)    ! interpolation indices

    real(kind=r8) :: q(3,6), lambda(4)

    ASSERT( size(point) == 3)
    ASSERT( cell >= 1 .and. cell <= mesh%ncell )
    ASSERT( all(shape(coef) == (/3,4/)) )
    ASSERT( size(index) == 4 )

    !! Length-weighted edge vectors
    q(:,1) = mesh%x(:,mesh%cnode(2,cell)) - mesh%x(:,mesh%cnode(1,cell))
    q(:,2) = mesh%x(:,mesh%cnode(3,cell)) - mesh%x(:,mesh%cnode(1,cell))
    q(:,3) = mesh%x(:,mesh%cnode(4,cell)) - mesh%x(:,mesh%cnode(1,cell))
    q(:,4) = mesh%x(:,mesh%cnode(3,cell)) - mesh%x(:,mesh%cnode(2,cell))
    q(:,5) = mesh%x(:,mesh%cnode(4,cell)) - mesh%x(:,mesh%cnode(2,cell))
    q(:,6) = mesh%x(:,mesh%cnode(4,cell)) - mesh%x(:,mesh%cnode(3,cell))

    !! Barycentric coordinates (scaled) of the interpolation point
    lambda = bary_coord(mesh%x(:,mesh%cnode(:,cell)), point) / (3.0_r8 * mesh%volume(cell))

    !! Local interpolation coefficients
    coef(:,1) =   lambda(2)*q(:,1) + lambda(3)*q(:,2) + lambda(4)*q(:,3)
    coef(:,2) = - lambda(1)*q(:,1) + lambda(3)*q(:,4) + lambda(4)*q(:,5)
    coef(:,3) = - lambda(1)*q(:,2) - lambda(2)*q(:,4) + lambda(4)*q(:,6)
    coef(:,4) = - lambda(1)*q(:,3) - lambda(2)*q(:,5) - lambda(3)*q(:,6)

    coef(:,2) = -coef(:,2)
    coef(:,4) = -coef(:,4)

    index = mesh%cface(:,cell)

  end subroutine eval_w2_interp_coef


  subroutine eval_w3_interp_coef (mesh, point, cell, coef, index)

    type(simpl_mesh), intent(in)  :: mesh
    real(kind=r8),   intent(in)  :: point(:)
    integer,         intent(in)  :: cell
    real(kind=r8),   intent(out) :: coef(:)
    integer,         intent(out) :: index(:)

    ASSERT( size(point) == 3 )
    ASSERT( cell > 0 .and. cell <= mesh%ncell )
    ASSERT( size(coef) == 1 .and. size(index) == 1 )

    coef(1) = 1.0_r8 / mesh%volume(cell)
    index(1) = cell

  end subroutine eval_w3_interp_coef


  function w1_face_matrix(mesh, face) result(matrix)
    type(simpl_mesh), intent(in) :: mesh
    integer, intent(in) :: face
    real(r8) :: matrix(6)
    real(r8) :: c, lsq(3)
    c = 1.0_r8 / (48.0_r8 * mesh%area(face))
    lsq = mesh%length(mesh%fedge(:,face))**2
    matrix(1) = c*(lsq(2) + lsq(3) - lsq(1))
    matrix(3) = c*(lsq(1) + lsq(3) - lsq(2))
    matrix(6) = c*(lsq(1) + lsq(2) - lsq(3))
    matrix(2) = -c*(lsq(1) + lsq(2) - 3*lsq(3))
    matrix(4) =  c*(lsq(1) + lsq(3) - 3*lsq(2))
    matrix(5) = -c*(lsq(2) + lsq(3) - 3*lsq(1))
  end function

end module mimetic_discretization
