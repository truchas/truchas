! TODO: finish documentation
!!
!! MFD_2D_DISC_TYPE
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! January 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module mfd_2d_disc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  implicit none
  private

  type, public :: mfd_2d_disc
    type(unstr_2d_mesh), pointer :: mesh => null()  ! reference only - do not own
    integer, allocatable :: xminv(:)
    real(r8), allocatable :: minv(:)
  contains
    procedure :: init => mfd_2d_disc_init
    procedure :: apply_diff => mfd_2d_disc_apply_diff
  end type mfd_2d_disc

  type, public :: mfd_2d_cell
    integer :: nface
    real(r8), allocatable :: cell_center(:)
    real(r8), allocatable :: face_center(:,:)
    real(r8), allocatable :: face_normal(:,:)
    real(r8), allocatable :: face_area(:)
    real(r8) :: volume
  contains
    procedure :: init => mfd_2d_cell_init
    procedure :: compute_flux_matrix_inv => mfd_2d_cell_compute_flux_matrix_inv
  end type mfd_2d_cell

contains

  !TODO: add optional MINV initialization? See 3D mfd_disc.
  subroutine mfd_2d_disc_init(this, mesh)
    class(mfd_2d_disc), intent(out) :: this
    type(unstr_2d_mesh), intent(in), target :: mesh
    this%mesh => mesh
    call init_minv(this)
  end subroutine mfd_2d_disc_init


  !! Allocates and initializes the THIS%MINV.
  subroutine init_minv(this)

    type(mfd_2d_disc), intent(inout) :: this

    integer :: j, n
    type(mfd_2d_cell), allocatable :: cell

    !! Initialize XMINV indexing array.
    allocate(this%xminv(this%mesh%ncell+1))
    this%xminv(1) = 1
    do j = 1, this%mesh%ncell
      n = this%mesh%cstart(j+1)-this%mesh%cstart(j)  ! number of sides
      this%xminv(j+1) = this%xminv(j) + (n*(n+1))/2
    end do

    !! Populate MINV, where MINV(XMINV(j):XMINV(j+1)-1) is the inverse M matrix
    !! of cell j, stored in upper packed matrix format.
    allocate(this%minv(this%xminv(this%mesh%ncell+1)-1))
    allocate(cell)
    do j = 1, this%mesh%ncell
      associate (minv => this%minv(this%xminv(j):this%xminv(j+1)-1))
          call cell%init(j, this%mesh)
          call cell%compute_flux_matrix_inv(1.0_r8, minv)
      end associate
    end do

  end subroutine init_minv


  !! Applies local MFD diffusion operator
  subroutine mfd_2d_disc_apply_diff (this, coef, ucell, uface, rcell, rface)

    use upper_packed_matrix_procs, only: sym_matmul

    class(mfd_2d_disc), intent(in) :: this
    real(r8), intent(in)  :: coef(:)
    real(r8), intent(in)  :: ucell(:), uface(:)
    real(r8), intent(out) :: rcell(:), rface(:)

    integer :: j
    real(r8), allocatable :: flux(:)

    ASSERT(size(coef) == this%mesh%ncell)
    ASSERT(size(ucell) == size(coef))
    ASSERT(size(rcell) == size(ucell))
    ASSERT(size(uface) == this%mesh%nface)
    ASSERT(size(rface) == size(uface))

    rface = 0.0_r8
    do j = 1, this%mesh%ncell
      associate (cface => this%mesh%cface(this%mesh%cstart(j):this%mesh%cstart(j+1)-1), &
                 minv  => this%minv(this%xminv(j):this%xminv(j+1)-1))
        flux = coef(j) * sym_matmul(minv, ucell(j) - uface(cface))
        rface(cface) = rface(cface) - flux
        rcell(j) = sum(flux)
      end associate
    end do

  end subroutine mfd_2d_disc_apply_diff


  !! MFD_2D_CELL type bound procedures !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mfd_2d_cell_init(this, cellid, mesh)

    use cell_geometry, only: cell_centroid_2d
    use bitfield_type, only: btest

    class(mfd_2d_cell), intent(out) :: this
    integer, intent(in) :: cellid
    type(unstr_2d_mesh), intent(in) :: mesh

    real(r8) :: parity
    integer :: j

    associate (cnode => mesh%cnode(mesh%cstart(cellid):mesh%cstart(cellid+1)-1), &
               cface => mesh%cface(mesh%cstart(cellid):mesh%cstart(cellid+1)-1))
      this%cell_center = cell_centroid_2d(mesh%x(:,cnode))
      this%volume      = mesh%volume(cellid)
      this%nface       = size(cface)
      this%face_area   = mesh%area(cface)
      allocate(this%face_normal(2,this%nface), this%face_center(2,this%nface))
      do j = 1, this%nface
        this%face_center(:,j) = 0.5_r8 * sum(mesh%x(:,mesh%fnode(:,cface(j))), dim=2)
        !! Outward facing normals
        parity = merge(-1.0, 1.0, btest(mesh%cfpar(cellid),j))
        this%face_normal(:,j) = mesh%normal(:,cface(j)) * parity
      end do
    end associate

  end subroutine mfd_2d_cell_init

  !! This is procedure computes inverse of the MFD mass matrix for a polygonal
  !! cell using the standard choice of MFD parameters:
  !! W - is an inverse of mass matrix
  !!
  !! W = (coef/volume)*N*N' + (coef/volume)*(I - Q*Q')
  !!
  !! N - matrix of normals
  !! Q - orthonormal basis of range(R)
  !! R - rows of (face center)-(cell center) for each cell face

  subroutine mfd_2d_cell_compute_flux_matrix_inv(this, coef, matrix)

    class(mfd_2d_cell), intent(in) :: this
    real(r8), intent(in) :: coef
    real(r8), intent(out) :: matrix(:)

    integer :: i, j, n, loc
    real(r8) :: stab_val
    real(r8) :: Q(2,this%nface), R(2,this%nface)

    matrix = 0.0_r8
    n = this%nface

    do i = 1, n
      R(:,i) = (this%face_center(:,i) - this%cell_center(:))*this%face_area(i)
    end do

    !! Compute Q by modified Grammm-Schimdt orthogonalization
    Q(1,:) = R(1,:) / norm2(R(1,:))
    Q(2,:) = R(2,:) - Q(1,:) * dot_product(Q(1,:), R(2,:))
    Q(2,:) = Q(2,:) / norm2(Q(2,:))

    !! Assemble W matrix
    do i = 1, n
       do j = i, n
          loc = i + j*(j-1)/2

          matrix(loc) = dot_product(this%face_normal(:,i), this%face_normal(:,j))

          stab_val = -dot_product(Q(:,i), Q(:,j))
          if (i.eq.j) stab_val = stab_val + 1

          stab_val = stab_val * this%face_area(i) * this%face_area(j)

          matrix(loc) = (matrix(loc) + stab_val)*(coef/this%volume)
       end do
    end do

  end subroutine mfd_2d_cell_compute_flux_matrix_inv

end module mfd_2d_disc_type
