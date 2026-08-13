!!
!! MFD_2D_DISC_TYPE
!!
!! This module defines the MFD_2D_DISC type used by thermal diffusion models to
!! represent a mimetic finite difference (MFD) discretization on an
!! unstructured mesh. Its primary state is the cell-local inverse flux mass
!! matrix MINV, stored in upper packed format. The standard construction
!! follows the generic polyhedral-cell MFD method described in:
!!
!!   Konstantin Lipnikov, Gianmarco Manzini, and Mikhail Shashkov. Mimetic
!!   finite difference method. Journal of Computational Physics, 257:1163-1227,
!!   2014.
!!
!!   V. Gyrya, K. Lipnikov, G. Manzini, and D. Svyatskiy. M-adaptation in the
!!   mimetic finite difference method. Mathematical Models & Methods in Applied
!!   Sciences, 24(8):1621-1663, 2014.
!!
!! The type is used by trusted code; its public data components should be
!! treated as read-only after initialization. It treats its mesh as a local
!! serial mesh, with no ownership or on/off-process semantics. Callers are
!! responsible for ensuring input values are valid on the full local mesh
!! extent they pass in, and for interpreting which entries of any result are
!! meaningful for their distributed algorithm.
!!
!! David Neill-Asanza <davidhneill@gmail.com>, January 2020
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module mfd_2d_disc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  implicit none
  private

  integer, parameter :: MFD_CELL_NFACE_MAX = 4

  type, public :: mfd_2d_disc
    type(unstr_2d_mesh), pointer :: mesh => null()  ! unowned reference
    integer, allocatable :: xminv(:)
    real(r8), allocatable :: minv(:)
    integer, private :: nface_max
  contains
    procedure :: init
    procedure :: apply_diff
    procedure :: add_diff
  end type

  !! Private type for internal use.
  type :: mfd_cell
    integer :: nfaces = 0
    real(r8) :: face_centers(2,MFD_CELL_NFACE_MAX)
    real(r8) :: cell_center(2)
    real(r8) :: face_area(MFD_CELL_NFACE_MAX)
    real(r8) :: face_normals(2,MFD_CELL_NFACE_MAX)
    real(r8) :: volume
  contains
    procedure :: init => mfd_cell_init
    procedure :: compute_flux_matrix_inv
  end type

contains

  subroutine init(this, mesh)

    class(mfd_2d_disc), intent(out) :: this
    type(unstr_2d_mesh), intent(in), target :: mesh

    integer :: j, n
    type(mfd_cell) :: cell

    this%mesh => mesh
    this%nface_max = max_cell_faces(mesh)
    ASSERT(this%nface_max <= MFD_CELL_NFACE_MAX)

    !! Initialize XMINV indexing array.
    allocate(this%xminv(this%mesh%ncell+1))
    this%xminv(1) = 1
    do j = 1, this%mesh%ncell
      n = this%mesh%cstart(j+1)-this%mesh%cstart(j)  ! number of sides
      this%xminv(j+1) = this%xminv(j) + (n*(n+1))/2
    end do

    !! Populate MINV, where MINV(XMINV(j):XMINV(j+1)-1) is the inverse of
    !! the mass matrix M of cell j, stored in upper packed matrix format.
    allocate(this%minv(this%xminv(this%mesh%ncell+1)-1))
    do j = 1, this%mesh%ncell
      associate (minv => this%minv(this%xminv(j):this%xminv(j+1)-1))
        call cell%init(j, this%mesh)
        call cell%compute_flux_matrix_inv(1.0_r8, minv)
      end associate
    end do

  contains

    integer function max_cell_faces(mesh)
      type(unstr_2d_mesh), intent(in) :: mesh
      if (mesh%ncell == 0) then
        max_cell_faces = 0
      else
        max_cell_faces = maxval(mesh%cstart(2:) - mesh%cstart(:mesh%ncell))
      end if
    end function

  end subroutine init

  !! Apply the local MFD diffusion operator as a serial operator over the full
  !! mesh held by this object. No communication is performed; COEF, UCELL, and
  !! UFACE must be valid over their full local extents. RCELL and RFACE are
  !! filled over those same extents, and callers decide which entries to use.

  subroutine apply_diff(this, coef, ucell, uface, rcell, rface)

    !use upper_packed_matrix_procs, only: upm_sym_matvec

    class(mfd_2d_disc), intent(in) :: this
    real(r8), intent(in)  :: coef(:)
    real(r8), intent(in)  :: ucell(:), uface(:)
    real(r8), intent(out) :: rcell(:), rface(:)

    integer :: j, n
    real(r8), allocatable :: du(:), flux(:)

    ASSERT(size(coef) == this%mesh%ncell)
    ASSERT(size(ucell) == size(coef))
    ASSERT(size(rcell) == size(ucell))
    ASSERT(size(uface) == this%mesh%nface)
    ASSERT(size(rface) == size(uface))

    allocate(du(this%nface_max), flux(this%nface_max))
    rface = 0.0_r8
    do j = 1, this%mesh%ncell
      associate (cface => this%mesh%cface(this%mesh%cstart(j):this%mesh%cstart(j+1)-1), &
                 minv  => this%minv(this%xminv(j):this%xminv(j+1)-1))
        n = size(cface)
        du(:n) = ucell(j) - uface(cface)
        call upm_sym_matvec(minv, du(:n), flux(:n))
        flux(:n) = coef(j) * flux(:n)
        rface(cface) = rface(cface) - flux(:n)
        rcell(j) = sum(flux(:n))
      end associate
    end do

  end subroutine apply_diff

  !! Add the local MFD diffusion operator result into RCELL and RFACE. This
  !! has the same input contract as APPLY_DIFF, but accumulates into the output
  !! arrays instead of overwriting them.

  subroutine add_diff(this, coef, ucell, uface, rcell, rface)

    !use upper_packed_matrix_procs, only: upm_sym_matvec

    class(mfd_2d_disc), intent(in) :: this
    real(r8), intent(in)    :: coef(:)
    real(r8), intent(in)    :: ucell(:), uface(:)
    real(r8), intent(inout) :: rcell(:), rface(:)

    integer :: j, n
    real(r8), allocatable :: du(:), flux(:)

    ASSERT(size(coef) == this%mesh%ncell)
    ASSERT(size(ucell) == size(coef))
    ASSERT(size(rcell) == size(ucell))
    ASSERT(size(uface) == this%mesh%nface)
    ASSERT(size(rface) == size(uface))

    allocate(du(this%nface_max), flux(this%nface_max))
    do j = 1, this%mesh%ncell
      associate (cface => this%mesh%cface(this%mesh%cstart(j):this%mesh%cstart(j+1)-1), &
                 minv  => this%minv(this%xminv(j):this%xminv(j+1)-1))
        n = size(cface)
        du(:n) = ucell(j) - uface(cface)
        call upm_sym_matvec(minv, du(:n), flux(:n))
        flux(:n) = coef(j) * flux(:n)
        rface(cface) = rface(cface) - flux(:n)
        rcell(j) = rcell(j) + sum(flux(:n))
      end associate
    end do

  end subroutine add_diff

  !! Multiply an upper-packed symmetric matrix by a vector.
  !! NB: Provided by a newer version of the UPPER_PACKED_MATRIX_PROCS module.

  pure subroutine upm_sym_matvec(a, b, c)
    real(r8), intent(in) :: a(:), b(:)
    real(r8), intent(out) :: c(:)
    integer :: i, j, l
    real(r8) :: bi, ci
    l = 1
    do i = 1, size(b)
      bi = b(i)
      ci = 0.0_r8
      do j = 1, i - 1
        ci = ci + a(l) * b(j)
        c(j) = c(j) + a(l) * bi
        l = l + 1
      end do
      c(i) = ci + a(l) * bi
      l = l + 1
    end do
  end subroutine

  !! MFD_CELL type bound procedures !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mfd_cell_init(this, cellid, mesh)

    use cell_geometry, only: cell_centroid_2d
    use bitfield_type, only: btest

    class(mfd_cell), intent(out) :: this
    integer, intent(in) :: cellid
    type(unstr_2d_mesh), intent(in) :: mesh

    real(r8) :: parity
    integer :: j

    associate (cnode => mesh%cnode(mesh%cstart(cellid):mesh%cstart(cellid+1)-1), &
               cface => mesh%cface(mesh%cstart(cellid):mesh%cstart(cellid+1)-1))
      this%cell_center = cell_centroid_2d(mesh%x(:,cnode))
      this%volume      = mesh%volume(cellid)
      this%nfaces      = size(cface)
      this%face_area(:this%nfaces) = mesh%area(cface)
      do j = 1, this%nfaces
        this%face_centers(:,j) = 0.5_r8 * sum(mesh%x(:,mesh%fnode(:,cface(j))), dim=2)
        !! Outward facing normals
        parity = merge(-1.0, 1.0, btest(mesh%cfpar(cellid),j))
        this%face_normals(:,j) = mesh%normal(:,cface(j)) * parity
      end do
    end associate

  end subroutine mfd_cell_init

  !! Compute the inverse MFD mass matrix for a polyhedral cell using the
  !! standard parameter choice; see references listed above. With N the matrix
  !! of face normals and Q an orthonormal basis for range(R), the inverse is
  !! W = (coef/volume) * (N*N' + I - Q*Q').

  subroutine compute_flux_matrix_inv(this, coef, matrix)

    class(mfd_cell), intent(in) :: this
    real(r8), intent(in) :: coef
    real(r8), intent(out) :: matrix(:)

    integer :: i, j, n, loc
    real(r8) :: stab_val
    real(r8) :: Q(2,this%nfaces), R(2,this%nfaces)

    matrix = 0.0_r8
    n = this%nfaces

    do i = 1, n
      R(:,i) = (this%face_centers(:,i) - this%cell_center(:))*this%face_area(i)
    end do

    !! Compute Q by modified Grammm-Schimdt orthogonalization
    Q(1,:) = R(1,:) / norm2(R(1,:))
    Q(2,:) = R(2,:) - Q(1,:) * dot_product(Q(1,:), R(2,:))
    Q(2,:) = Q(2,:) / norm2(Q(2,:))

    !! Assemble W matrix
    do i = 1, n
       do j = i, n
          loc = i + j*(j-1)/2

          matrix(loc) = dot_product(this%face_normals(:,i), this%face_normals(:,j))

          stab_val = -dot_product(Q(:,i), Q(:,j))
          if (i.eq.j) stab_val = stab_val + 1

          stab_val = stab_val * this%face_area(i) * this%face_area(j)

          matrix(loc) = (matrix(loc) + stab_val)*(coef/this%volume)
       end do
    end do

  end subroutine compute_flux_matrix_inv

end module mfd_2d_disc_type
