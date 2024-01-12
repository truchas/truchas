!!
!! VOL_FRAC_INIT_PROCS
!!
!! This module provides procedures for computing the cell-based volume fractions
!! occupied by the regions of a user-specified disjoint decomposition of the
!! domain of the mesh.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! January 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module vol_frac_init_procs

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use region_func_type
  implicit none
  private

  public :: compute_volume_fractions

  !! Private type implementing the divide-and-conquer method for
  !! computing the volume fractions for a triangular cell.
  type :: tri_cell
    type(region_func), pointer :: rfunc => null()
    real(r8) :: x(2,3)
    integer  :: regid(3)
    integer  :: nreg
    integer  :: bitmask   ! cell set membership of parent cell
  contains
    procedure :: init => tri_cell_init
    procedure :: region_id
    procedure :: centroid_region_id
    procedure :: vol_frac
    procedure :: volume
    procedure :: subdivide
  end type

contains

  subroutine compute_volume_fractions(mesh, rfunc, rlev, vol_frac, stat)

    use unstr_2d_mesh_type

    type(unstr_2d_mesh), intent(in) :: mesh
    type(region_func), intent(in), target :: rfunc
    integer,  intent(in)  :: rlev
    real(r8), intent(out) :: vol_frac(:,:)
    integer,  intent(out) :: stat

    integer :: i, j
    type(tri_cell) :: tri
    type(tri_cell), allocatable :: subtri(:)

    ASSERT(size(vol_frac,dim=1) == rfunc%num_region())
    ASSERT(size(vol_frac,dim=2) == mesh%ncell_onP)

    stat = 0

    do j = 1, mesh%ncell_onP
      associate (cnode => mesh%cnode(mesh%cstart(j):mesh%cstart(j+1)-1))
        if (size(cnode) == 3) then  ! triangular cell
          call tri%init(mesh%x(:,cnode), mesh%cell_set_mask(j), rfunc, stat)
          if (stat /= 0) return
          vol_frac(:,j) = tri%vol_frac(rlev, stat)
          if (stat /= 0) return
        else ! more general polygonal cell
          call triangulate(mesh%x(:,cnode), mesh%cell_set_mask(j), rfunc, subtri, stat)
          if (stat /= 0) return
          vol_frac(:,j) = 0.0_r8
          do i = 1, size(subtri)
            vol_frac(:,j) = vol_frac(:,j) + subtri(i)%volume() * subtri(i)%vol_frac(rlev, stat)
            if (stat /= 0) return
          end do
          vol_frac(:,j) = vol_frac(:,j) / sum(vol_frac(:,j))
        end if
      end associate
    end do

  end subroutine compute_volume_fractions

  !! This auxiliary subroutine subdivides an oriented convex polygon with nodes
  !! x(:,j), j = 1, n, into an array of triangular cells SUBTRI. STAT returns 1
  !! if any of the nodes do not belong to a region; otherwise STAT returns 0.
  !! In any case, SUBTRI is fully initialized.

  subroutine triangulate(x, bitmask, rfunc, subtri, stat)

    real(r8), intent(in) :: x(:,:)
    integer, intent(in) :: bitmask
    type(region_func), intent(in), target :: rfunc
    type(tri_cell), allocatable, intent(out) :: subtri(:)
    integer, intent(out) :: stat

    integer :: n, i
    integer, allocatable :: regid(:)

    n = size(x,dim=2)

    allocate(regid(n))
    do i = 1, size(regid)
      regid(i) = rfunc%region_index(x(:,i), bitmask)
    end do
    stat = merge(1, 0, any(regid == 0))

    allocate(subtri(n-2))
    do i = 1, size(subtri)
      subtri(i)%bitmask = bitmask
      subtri(i)%nreg = rfunc%num_region()
      subtri(i)%rfunc => rfunc
    end do

    call triangulate_aux(x, regid, [(i, i=1,n)], subtri)

  end subroutine

  !! Recursively subdivide a convex polygon into an array of triangular cells.
  !! The nodes of the oriented polygon are x(:,cnode(i)), i = 1, size(cnode).

  recursive subroutine triangulate_aux(x, regid, cnode, subtri)

    real(r8), intent(in) :: x(:,:)
    integer, intent(in) :: regid(:)
    integer, intent(in) :: cnode(:)
    type(tri_cell), intent(inout) :: subtri(:)

    integer :: i, j, n, m, p, q
    integer, allocatable :: subcnode(:)
    real(r8) :: tmp, minlen

    n = size(cnode)
    ASSERT(n >= 3)
    ASSERT(size(subtri) == n-2)

    if (n == 3) then ! a triangle; store its node coord and region IDs and return
      subtri(1)%x = x(:,cnode)
      subtri(1)%regid = regid(cnode)
    else ! split polygon along shortest diagonal and recurse
      minlen = huge(minlen)
      do i = 1, n-2
        do j = i+2, merge(n-1, n, i==1)
          tmp = norm2(x(:,cnode(i))-x(:,cnode(j)))
          if (tmp < minlen) then
            minlen = tmp
            p = i
            q = j
          end if
        end do
      end do
      subcnode = [(cnode(i), i = p, q)]
      m = size(subcnode) - 2
      call triangulate_aux(x, regid, subcnode, subtri(:m))
      subcnode = [(cnode(i), i = 1, p), (cnode(i), i = q, n)]
      call triangulate_aux(x, regid, subcnode, subtri(m+1:))
    end if

  end subroutine triangulate_aux

!!!! TRI_CELL TYPE BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tri_cell_init(this, x, bitmask, rfunc, stat)
    class(tri_cell), intent(out) :: this
    real(r8), intent(in) :: x(:,:)
    integer, intent(in) :: bitmask
    type(region_func), intent(in), target :: rfunc
    integer, intent(out) :: stat
    integer :: i
    this%x = x
    this%bitmask = bitmask
    this%rfunc => rfunc
    this%nreg = this%rfunc%num_region()
    do i = 1, 3
      this%regid(i) = this%rfunc%region_index(this%x(:,i), this%bitmask)
    end do
    stat = merge(1, 0, any(this%regid == 0))
  end subroutine

  !! If the cell vertices all belong to the same region, return its index;
  !! otherwise return 0 to signal a multi-region cell.
  pure integer function region_id(this) result(n)
    class(tri_cell), intent(in) :: this
    integer :: i
    n = 0
    do i = 2, 3
      if (this%regid(i) /= this%regid(1)) return
    end do
    n = this%regid(1)
  end function

  !! Return the region id that contains the cell center.
  pure integer function centroid_region_id(this) result(n)
    class(tri_cell), intent(in) :: this
    n = this%rfunc%region_index(sum(this%x,dim=2)/3.0_r8, this%bitmask)
  end function

  !! Return the volume of the triangular cell.
  pure real(r8) function volume(this)
    use cell_geometry, only: tri_area
    class(tri_cell), intent(in) :: this
    volume = tri_area(this%x)
  end function

  !! Returns a vector of region volume fractions
  recursive function vol_frac(this, rlev, stat)

    class(tri_cell), intent(in) :: this
    integer, intent(in)  :: rlev
    integer, intent(out) :: stat
    real(r8) :: vol_frac(this%nreg)

    integer :: n, i

    vol_frac = 0.0_r8
    n = this%region_id()
    if (n == 0) then ! a multi-region cell
      if (rlev > 0) then ! subdivide and recurse
        block
          type(tri_cell) :: subtri(4)
          call this%subdivide(subtri, stat)
          if (stat /= 0) return
          do i = 1, 4
            vol_frac = vol_frac + subtri(i)%vol_frac(rlev-1, stat)
            if (stat /= 0) return
          end do
          vol_frac = 0.25_r8 * vol_frac
        end block
      else ! apportion it all to the region containing the cell center
        n = this%centroid_region_id()
        if (n == 0) then
          stat = 2
          return
        end if
        vol_frac(n) = 1.0_r8
        !TODO: more accurate apportionment via interface reconstruction as in 3D
      end if
    else
      vol_frac(n) = 1.0_r8
    end if

    stat = 0 ! success; all interpolated points belong to a region

  end function

  !! Subdivide the triangular cell THIS into 4 congruent triangular cells and
  !! return the result in the array SUBTRI. The subcells inherit the values of
  !! the parent components other than the vertex coordinates and their
  !! corresponding region IDs. STAT returns 2 if an interpolated point does not
  !! belong to a region; otherwise STAT returns 0.

  subroutine subdivide(this, subtri, stat)

    class(tri_cell), intent(in) :: this
    type(tri_cell), intent(out) :: subtri(4)
    integer, intent(out) :: stat

    integer :: i

    associate (xmid => subtri(4)%x, rmid => subtri(4)%regid)
      xmid(:,1) = 0.5_r8*(this%x(:,2) + this%x(:,3))
      xmid(:,2) = 0.5_r8*(this%x(:,3) + this%x(:,1))
      xmid(:,3) = 0.5_r8*(this%x(:,1) + this%x(:,2))
      subtri(1)%x(:,1) = this%x(:,1)
      subtri(1)%x(:,2) = xmid(:,3)
      subtri(1)%x(:,3) = xmid(:,2)
      subtri(2)%x(:,1) = xmid(:,3)
      subtri(2)%x(:,2) = this%x(:,2)
      subtri(2)%x(:,3) = xmid(:,1)
      subtri(3)%x(:,1) = xmid(:,2)
      subtri(3)%x(:,2) = xmid(:,1)
      subtri(3)%x(:,3) = this%x(:,3)

      do i = 1, 3
        rmid(i) = this%rfunc%region_index(xmid(:,i), this%bitmask)
      end do
      stat = merge(2, 0, any(rmid == 0))

      subtri(1)%regid(1) = this%regid(1)
      subtri(1)%regid(2) = rmid(3)
      subtri(1)%regid(3) = rmid(2)
      subtri(2)%regid(1) = rmid(3)
      subtri(2)%regid(2) = this%regid(2)
      subtri(2)%regid(3) = rmid(1)
      subtri(3)%regid(1) = rmid(2)
      subtri(3)%regid(2) = rmid(1)
      subtri(3)%regid(3) = this%regid(3)
    end associate

    do i = 1, 4
      subtri(i)%nreg = this%nreg
      subtri(i)%rfunc => this%rfunc
      subtri(i)%bitmask = this%bitmask
    end do

  end subroutine subdivide

end module vol_frac_init_procs
