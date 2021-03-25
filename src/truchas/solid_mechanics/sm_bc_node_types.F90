!!
!! Zach Jibben <zjibben@lanl.gov>
!! February 2021
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_bc_node_types

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use scalar_func_class
  !use sm_bc_face_type
  !use scalar_func_containers, only: scalar_func_box
  implicit none
  private

  !! Private types
  type :: scalar_func_ptr
    class(scalar_func), pointer :: f => null() ! unowned reference
  contains
    procedure :: eval => scalar_func_ptr_eval
  end type scalar_func_ptr


  !! Public types
  type, public :: sm_bc_d1
    private
    integer, allocatable, public :: index(:)
    real(r8), allocatable, public :: value(:,:)
    real(r8), allocatable, public :: normal(:,:)
    type(scalar_func_ptr), allocatable :: displf(:)
  contains
    procedure :: init => sm_bc_d1_init
    procedure :: compute => sm_bc_d1_compute
  end type sm_bc_d1

  type, public :: sm_bc_d2
    private
    integer, allocatable, public :: index(:)
    real(r8), allocatable, public :: value(:,:)
    real(r8), allocatable, public :: tangent(:,:)
    real(r8), allocatable :: normal(:,:,:)
    type(scalar_func_ptr), allocatable :: displf(:,:)
  contains
    procedure :: init => sm_bc_d1_init
    procedure :: compute => sm_bc_d1_compute
  end type sm_bc_d2

  type, public :: sm_bc_d3
    private
    integer, allocatable, public :: index(:)
    real(r8), allocatable, public :: value(:,:)
    real(r8), allocatable :: normal(:,:,:)
    type(scalar_func_ptr), allocatable :: displf(:,:)
  contains
    procedure :: init => sm_bc_d1_init
    procedure :: compute => sm_bc_d1_compute
  end type sm_bc_d3

contains

  !! Initialize the single-normal displacement BC nodes.
  !!
  !! The main difficulty here, and in the other below routines for initializing
  !! node-BCs is identifying the appropriate nodes. Here we take a
  !! face_displacement object which holds BCs on faces, along with some helpful
  !! mapping arrays for looking at nodes connected to those faces, and quickly
  !! searching BC faces neighboring each given node. Thus one private procedure
  !! uses faces surrounding a given node to decide whether that node meets the
  !! conditions for this type.
  !!
  !! Here we're looking for nodes touching only faces with a displacement-BC in
  !! a single normal direction. That normal direction may change smoothly over
  !! a curved surface, so there is a small tolerance for deciding what qualifies
  !! as 'distinct' normal vectors between faces. Currently the tolerance is set
  !! to about 25 degrees.
  subroutine sm_bc_d1_init(this, face_displacement)

    use sm_bc_node_list_type

    class(sm_bc_d1), intent(out) :: this
    type(sm_bc_node_list), intent(in) :: nodebc
    type(sm_bc_list), intent(in), target :: bc

    integer :: nnode, ni

    nnode = 0
    do ni = 1, size(face_displacement%ni_index)
      if (is_single_normal_node(ni)) nnode = nnode + 1
    end do
    allocate(this%index(nnode), this%value(3,nnode), this%normal(3,nnode), this%displf(nnode))

    nnode = 0
    do ni = 1, nodebc%nnode
      if (.not.is_single_normal_node(ni)) cycle
      nnode = nnode + 1
      this%index(nnode) = nodebc%node(ni)
      bcid = nodebc%bcid(nodebc%offset(ni))
      this%displf(nnode)%f => bc%displacement(bcid)%f
      this%normal(:,nnode) = nodebc%normal(:,nodebc%offset(ni))
    end do

  contains

    !! NB: Here we are looking for nodes with only one displacement BC. The
    !! obvious scenario is where all faces share just one BC along a flat
    !! surface. But we also want to include scenarios where surrounding faces
    !! share a single BC along a curved surface.
    !!
    !!   - INCLUDE: surrounding faces have 1 normal regardless of face sets
    !!   - INCLUDE: surrounding faces with multiple normals (curved surface)
    !!              when part of the same face set.
    !!   - EXCLUDE: surrounding faces have multiple BCs OR face sets w/
    !!              different normal vectors
    !!   - EXCLUDE: surrounding faces have non-displacement BCs
    !!
    !! In all of the EXCLUDE cases, we instead apply either a multi-normal
    !! displacment BC, or a displacement+contact combination BC.
    !!
    !! TODO: There are a couple scenarios which may result in unexpected
    !! behavior. Both of these should be uncommon in practice.
    !!
    !!   - Touching displacement BCs in the same direction, but with different
    !!     displacement values. As written this will apply a single-normal BC
    !!     with the first displacement value seen.
    !!   - Touching displacement BCs of the same displacement value, but along a
    !!     curved surface rather than an edge or corner. If the curvature of
    !!     this surface is high enough to surpase the tol parameter, this will
    !!     be treated as a multi-normal BC.
    pure logical function is_single_normal_node(ni)

      integer, intent(in) :: ni

      real(r8), parameter :: tol = 0.9_r8
      integer :: ni, fi, xfi, n, face_set
      real(r8) :: n1(3), n2(3)

      is_single_normal_node = .false.

      xfi = nodebc%offset(ni)
      bc1 = nodebc%bcid(xfi)
      n1 = nodebc%normal(:,xfi)
      do xfi = nodebc%offset(ni)+1, nodebc%offset(ni+1)-1
        bc2 = nodebc%bcid(xfi)
        n2 = nodebc%normal(:,xfi)

        ! See above NB
        if (bc1 /= bc2 .and. dot_product(n1, n2) < tol) return
      end do

      is_single_normal_node = .true.

    end function is_single_normal_node

  end subroutine sm_bc_d1_init


  subroutine sm_bc_d1_compute(this, t)

    class(sm_bc_d1), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: i
    real(r8) :: args(0:size(this%mesh%x,dim=1))

    args(0) = t
    !this%value = 0
    do i = 1, size(this%index)
      args(1:) = this%mesh%x(:,this%index(i))
      this%value(:,i) = this%normal(:,i) * this%displf(i)%eval(args)

      ! ni = this%ni_index(i)
      ! weight = 0
      ! do xfi = xnifi(ni), xnifi(ni+1)-1
      !   fi = nifi(xfi)
      !   weight = weight + this%area_ip(xfi)
      !   this%value(:,i) = this%value(:,i) + this%area_ip(xfi) * this%face_displacement%value(:,fi)
      ! end do
      ! this%value(:,i)  = this%value(:,i) / weight
    end do

  end subroutine sm_bc_d1_compute


  subroutine sm_bc_d2_compute(this, t)

    external dgesv ! LAPACK

    class(sm_bc_d2), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: i, stat
    real(r8) :: args(0:size(this%mesh%x,dim=1)), displ(2), matrix(3,3), ipiv(3)

    args(0) = t
    do i = 1, size(this%index)
      args(1:) = this%mesh%x(:,this%index(i))

      displ(1) = this%displf(1,i)%eval(args) ! associated with this%normal(:,1,i)
      displ(2) = this%displf(2,i)%eval(args) ! associated with this%normal(:,2,i)
      ! !displ(3) = 0 ! associated with this%normal(:,3,i) or this%tangent(:,i)

      ! call dgesv(2, 1, this%normal(:,:,i), 2, ipiv, displ, 1, stat)
      ! INSIST(stat == 0)

      this%value(:,i) = displacement_vector(this%normal(:,:,i), displ)
    end do

  contains

    !! in (d1,d2) displacements in the given normal directions,
    !! out a displacement vector a = b1*n1 + b2*n2 such that
    !!   dot(n1, a) = d1
    !!   dot(n2, a) = d2
    function displacement_vector(normal, displacements) result(a)
      real(r8), intent(in) :: normal(:,:), displ(:)
      real(r8) :: a(3)
      real(r8) :: matrix(2,2)
      matrix(1,1) = 1
      matrix(2,1) = dot_product(normal(:,1),normal(:,2))
      matrix(2,2) = 1
      matrix(1,2) = matrix(2,1)
      call solve2x2(matrix, displacements)
      a = matmul(normal, displacements)
    end function displacement_vector

    !! Solve 2x2 a*x=b linear system.
    subroutine solve2x2(a, x)

      real(r8), intent(inout) :: a(:,:) ! in a out inv(a)
      real(r8), intent(inout) :: x(:) ! in b out x

      real(r8) :: det, s

      ASSERT(size(x) == 2)
      ASSERT(size(a,dim=1) == 2 .and. size(a,dim=2) == 2)

      det = 1 - a(1,2)*a(2,1)
      ASSERT(det /= 0)

      s = a(1,1)
      a(1,1) = a(2,2)
      a(2,2) = s

      a(1,2) = -a(1,2)
      a(2,1) = -a(2,1)

      a = a/det

      x = matmul(a, x)

    end subroutine solve2x2

  end subroutine sm_bc_d2_compute


  !! When given three displacements (d1,d2,d3) associated with linearly
  !! independent normal directions, we compute a displacement vector a
  !! such that:
  !!   dot(n1, a) = d1
  !!   dot(n2, a) = d2
  !!   dot(n3, a) = d3
  subroutine sm_bc_d3_compute(this, t)

    external dgesv ! LAPACK

    class(sm_bc_d3), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: i, stat
    real(r8) :: args(0:size(this%mesh%x,dim=1)), displ(3), matrix(3,3), ipiv(3)

    args(0) = t
    do i = 1, size(this%index)
      args(1:) = this%mesh%x(:,this%index(i))

      displ(1) = this%displf(1,i)%eval(args) ! associated with this%normal(:,1,i)
      displ(2) = this%displf(2,i)%eval(args) ! associated with this%normal(:,2,i)
      displ(3) = this%displf(3,i)%eval(args) ! associated with this%normal(:,3,i)

      call dgesv(3, 1, this%normal(:,:,i), 3, ipiv, displ, 1, stat)
      INSIST(stat == 0)

      this%value(:,i) = displ
    end do

  end subroutine sm_bc_d3_compute


  !! Private methods
  pure real(r8) function scalar_func_ptr_eval(this, x)
    class(scalar_func_ptr), intent(in) :: this
    real(r8), intent(in) :: x(:)
    scalar_func_ptr_eval = this%f%eval(x)
  end function scalar_func_ptr_eval

end module sm_bc_node_types
