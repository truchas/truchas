!!
!! CELL_GEOMETRY
!!
!! This module provides some basic geometric primitives for 3D cell types and
!! and their facets.  The module CELL_TOPOLOGY is the companion to this module.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 9 Feb 2006; updated August 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! Some procedures come with important caveats; please read and understand the
!! comments below before using the procedures.  This applies particularly to
!! cells with non-planar faces.
!!

#include "f90_assert.fpp"

module cell_geometry

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  !! Cell volumes
  public :: tet_volume, pyramid_volume, wedge_volume, hex_volume
  public :: cell_volume ! wraps the preceding functions
  public :: eval_hex_volumes  ! special version that returns corner volumes also

  !! Cell centers
  public :: cell_center, cell_centroid
  public :: cell_centroid_2d

  !! Face normals
  public :: tri_face_normal, quad_face_normal
  public :: face_normal ! wraps the preceding functions

  !! Cell face normals
  public :: tet_face_normals, pyramid_face_normals, wedge_face_normals, hex_face_normals
  public :: cell_face_normals ! wraps the preceding functions

  !! Cell face centers
  public :: tet_face_centers, pyramid_face_centers, wedge_face_centers, hex_face_centers
  public :: cell_face_centers
  public :: polygon_center

  !! Algebraic primitives
  public :: cross_product, triple_product, vector_length, tri_area, normalized

  interface tri_area
    procedure tri_area_length, tri_area_coord
  end interface

contains

  !! Returns the volume of the given 3D cell.  X(:,k) are the coordinates of
  !! the kth vertex of the cell.  The cell type is inferred from the number of
  !! vertices.  The function handles tet, pyramid, wedge, and hex cell types;
  !! 0 is returned for anything else.  The volume of a tet is computed using
  !! the common formula involving the triple product of three edges sharing a
  !! vertex.  The problem for the other cell types is much more complex.  A
  !! general hexahedron is the image of a unit cube under the tri-linear mapping
  !! defined by the vertices of the hexahedron, and will not generally have
  !! planar faces.  Its volume is computed using a formula involving the volume
  !! of 10 sub-tetrahedra: O.V. Ushakova, "Conditions of nondegeneracy of
  !! three-dimensional cells.  A formula of a volume of cells.", SIAM J. Sci.
  !! Comput., 23, 2001.  Pyramid and wedge cells are regarded as degenerate
  !! hexahedra, and their volumes computed using the formula for a hex.

  pure function cell_volume (x) result (vol)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: vol
    select case (size(x,dim=1))
    case (2)
      select case (size(x,dim=2))
      case (3)  ! tri
        vol = tri_area(x)
      case (4)  ! quad
        vol = quad_area(x)
      end select
    case (3)
      select case (size(x,dim=2))
      case (4)  ! tet
        vol = tet_volume(x)
      case (5)  ! pyramid
        vol = pyramid_volume(x)
      case (6)  ! wedge
        vol = wedge_volume(x)
      case (8)  ! hex
        vol = hex_volume(x)
      case default
        vol = 0.0_r8
      end select
    end select
  end function cell_volume

  pure function tet_volume (x) result (vol)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: vol
    vol = triple_product(x(:,2)-x(:,1), x(:,3)-x(:,1), x(:,4)-x(:,1)) / 6.0_r8
  end function tet_volume

  pure function pyramid_volume (x) result (vol)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: vol
    vol = 0.5_r8 * (tet_volume(x(:,[1,2,4,5])) &
                  + tet_volume(x(:,[2,3,1,5])) &
                  + tet_volume(x(:,[3,4,2,5])) &
                  + tet_volume(x(:,[4,1,3,5])) )
  end function pyramid_volume

  pure function wedge_volume (x) result (vol)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: vol
    vol = 0.5_r8 * (tet_volume(x(:,[1,2,3,4])) &
                  + tet_volume(x(:,[5,4,6,2])) &
                  + tet_volume(x(:,[2,3,4,6])) &
                  + tet_volume(x(:,[2,3,1,5])) &
                  + tet_volume(x(:,[4,6,5,1])) &
                  + tet_volume(x(:,[1,3,6,5])) )
  end function wedge_volume

  pure function hex_volume (x) result (hvol)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: hvol, cvol(8)
    call eval_hex_volumes (x, hvol, cvol)
  end function hex_volume

  !! Computes the corner volumes CVOL(:) of a hex cell and its total volume HVOL.
  !! This uses the formula of Ushakova cited above.

  pure subroutine eval_hex_volumes (x, hvol, cvol)
    real(r8), intent(in)  :: x(:,:)
    real(r8), intent(out) :: hvol, cvol(:)
    cvol(1) = tet_volume(x(:,[1,2,4,5]))
    cvol(2) = tet_volume(x(:,[2,3,1,6]))
    cvol(3) = tet_volume(x(:,[3,4,2,7]))
    cvol(4) = tet_volume(x(:,[4,1,3,8]))
    cvol(5) = tet_volume(x(:,[5,8,6,1]))
    cvol(6) = tet_volume(x(:,[6,5,7,2]))
    cvol(7) = tet_volume(x(:,[7,6,8,3]))
    cvol(8) = tet_volume(x(:,[8,7,5,4]))
    hvol = 0.5_r8 * (sum(cvol) + tet_volume(x(:,[1,3,8,6])) + tet_volume(x(:,[2,4,5,7])))
  end subroutine eval_hex_volumes

  !! Returns the outward area-weighted normal vectors to the faces of the
  !! given 3D cell. X(:,k) are the coordinates of the kth vertex of the cell.
  !! The cell type is inferred from the number of vertices.  The function
  !! handles tet, pyramid, wedge, and hex cell types; a 0-sized array result
  !! is returned for anything else.  For a planar face (triangles and special
  !! quadrilateral faces) the normal is unambiguous.  For a warped quad face
  !! the normal direction is defined to be that orthogonal to the diagonals of
  !! the face, and the area that of the projection of the face onto a plane
  !! orthogonal to the normal.  This definition coincides with the integral of
  !! the unit normal over the face.

  pure function cell_face_normals (x) result (a)
    real(r8), intent(in) :: x(:,:)
    real(r8), allocatable :: a(:,:)
    select case (size(x,2))
    case (4)
      a = tet_face_normals(x)
    case (5)
      a = pyramid_face_normals(x)
    case (6)
      a = wedge_face_normals(x)
    case (8)
      a = hex_face_normals(x)
    case default
      allocate(a(3,0))
    end select
  end function cell_face_normals

  pure function tet_face_normals (x) result (a)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: a(3,4)
    !! NB: These must be consistent with the TET4 vertex and face labelings
    !! defined in the CELL_TOPOLOGY module.  To avoid a layer of indirection,
    !! its TET4_FACES array was not used here.
    a(:,1) = 0.5_r8 * cross_product(x(:,1)-x(:,4), x(:,2)-x(:,4))
    a(:,2) = 0.5_r8 * cross_product(x(:,2)-x(:,4), x(:,3)-x(:,4))
    a(:,3) = 0.5_r8 * cross_product(x(:,3)-x(:,4), x(:,1)-x(:,4))
    a(:,4) = 0.5_r8 * cross_product(x(:,3)-x(:,1), x(:,2)-x(:,1))
  end function tet_face_normals

  pure function pyramid_face_normals (x) result (a)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: a(3,5)
    !! NB: These must be consistent with the PYR5 vertex and face labelings
    !! defined in the CELL_TOPOLOGY module.  To avoid a layer of indirection,
    !! its PYR5_FACES array was not used here.
    a(:,1) = 0.5_r8 * cross_product(x(:,2)-x(:,1), x(:,5)-x(:,1))
    a(:,2) = 0.5_r8 * cross_product(x(:,3)-x(:,2), x(:,5)-x(:,2))
    a(:,3) = 0.5_r8 * cross_product(x(:,4)-x(:,3), x(:,5)-x(:,3))
    a(:,4) = 0.5_r8 * cross_product(x(:,1)-x(:,4), x(:,5)-x(:,4))
    a(:,5) = 0.5_r8 * cross_product(x(:,3)-x(:,1), x(:,2)-x(:,4))

  end function pyramid_face_normals

  pure function wedge_face_normals (x) result (a)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: a(3,5)
    !! NB: These must be consistent with the WED6 vertex and face labelings
    !! defined in the CELL_TOPOLOGY module.  To avoid a layer of indirection,
    !! its WED6_FACES array was not used here.
    a(:,1) = 0.5_r8 * cross_product(x(:,5)-x(:,1), x(:,4)-x(:,2))
    a(:,2) = 0.5_r8 * cross_product(x(:,6)-x(:,2), x(:,5)-x(:,3))
    a(:,3) = 0.5_r8 * cross_product(x(:,4)-x(:,3), x(:,6)-x(:,1))
    a(:,4) = 0.5_r8 * cross_product(x(:,3)-x(:,1), x(:,2)-x(:,1))
    a(:,5) = 0.5_r8 * cross_product(x(:,5)-x(:,4), x(:,6)-x(:,4))
  end function wedge_face_normals

  pure function hex_face_normals (x) result (a)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: a(3,6)
    !! NB: These must be consistent with the HEX8 vertex and face labelings
    !! defined in the CELL_TOPOLOGY module.  To avoid a layer of indirection,
    !! its HEX8_FACES array was not used here.
    a(:,1) = 0.5_r8 * cross_product(x(:,6)-x(:,1), x(:,5)-x(:,2))
    a(:,2) = 0.5_r8 * cross_product(x(:,7)-x(:,2), x(:,6)-x(:,3))
    a(:,3) = 0.5_r8 * cross_product(x(:,8)-x(:,3), x(:,7)-x(:,4))
    a(:,4) = 0.5_r8 * cross_product(x(:,5)-x(:,4), x(:,8)-x(:,1))
    a(:,5) = 0.5_r8 * cross_product(x(:,3)-x(:,1), x(:,2)-x(:,4))
    a(:,6) = 0.5_r8 * cross_product(x(:,7)-x(:,5), x(:,8)-x(:,6))
  end function hex_face_normals

  !! These funcctions returns the area-weighted normal to the given oriented
  !! triangle or quadrilateral face in 3D.  X(:,k) are the coordinates of the
  !! kth vertex of the face.  In the case of FACE_NORMAL, the face type is
  !! inferred from the number of vertices.  For warped quadrilateral faces
  !! the definition of the normal is as described for CELL_FACE_NORMALS.

  pure function tri_face_normal (x) result (a)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: a(3)
    a = 0.5_r8 * cross_product(x(:,2)-x(:,1), x(:,3)-x(:,2))
  end function tri_face_normal

  pure function quad_face_normal (x) result (a)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: a(3)
    a = 0.5_r8 * cross_product(x(:,3)-x(:,1), x(:,4)-x(:,2))
  end function quad_face_normal

  pure function face_normal (x) result (a)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: a(3)
    select case (size(x,dim=2))
    case (3)  ! triangular face
      a = 0.5_r8 * cross_product(x(:,2)-x(:,1),x(:,3)-x(:,2))
    case (4)  ! quadrilateral face
      a = 0.5_r8 * cross_product(x(:,3)-x(:,1), x(:,4)-x(:,2))
    case default
      a = 0.0_r8
    end select
  end function face_normal

  !! Returns the coordinates of the face centers of the given 3D cell.
  !! X(:,k) are the coordinates of the kth vertex of the cell.  The cell
  !! type is inferred from the number of vertices.  The function handles tet,
  !! pyramid, wedge, and hex cell types; a 0-sized array result is returned
  !! for anything else.  For triangular faces and planar quadrilateral faces
  !! the definition of the center is unambiguous and corresponds to the center
  !! of mass of the face.  Warped quadrilateral faces are approximated by a
  !! triangulated surface formed by adding the algebraic mean of the vertices
  !! as a vertex, and taking the center of mass of that surface as the center
  !! of the warped quadrilateral face.  NB: These must be consistent with the
  !! TET5, PYR5, WED6, HEX8 vertex and face labelings defined in the
  !! CELL_TOPOLOGY module.  They were not used here to avoid indirection.

  pure function cell_face_centers (x) result (xc)
    real(r8), intent(in) :: x(:,:)
    real(r8), allocatable :: xc(:,:)
    select case (size(x,dim=2))
    case (4)  ! tet
      xc = tet_face_centers(x)
    case (5)  ! pyramid
      xc = pyramid_face_centers(x)
    case (6)  ! wedge
      xc = wedge_face_centers(x)
    case (8)  ! hex
      xc = hex_face_centers(x)
    case default
      allocate(xc(3,0))
    end select
  end function cell_face_centers

  pure function tet_face_centers (x) result (xc)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: xc(3,4)
    xc(:,1) = (x(:,1) + x(:,2) + x(:,4)) / 3.0_r8
    xc(:,2) = (x(:,2) + x(:,3) + x(:,4)) / 3.0_r8
    xc(:,3) = (x(:,1) + x(:,3) + x(:,4)) / 3.0_r8
    xc(:,4) = (x(:,1) + x(:,2) + x(:,3)) / 3.0_r8
  end function tet_face_centers

  pure function pyramid_face_centers (x) result (xc)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: xc(3,5)
    xc(:,1) = (x(:,1) + x(:,2) + x(:,5)) / 3.0_r8
    xc(:,2) = (x(:,2) + x(:,3) + x(:,5)) / 3.0_r8
    xc(:,3) = (x(:,3) + x(:,4) + x(:,5)) / 3.0_r8
    xc(:,4) = (x(:,1) + x(:,4) + x(:,5)) / 3.0_r8
    xc(:,5) = polygon_center(x(:,[1,4,3,2]))
  end function pyramid_face_centers

  pure function wedge_face_centers (x) result (xc)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: xc(3,5)
    xc(:,1) = polygon_center(x(:,[1,2,5,4]))
    xc(:,2) = polygon_center(x(:,[2,3,6,5]))
    xc(:,3) = polygon_center(x(:,[1,4,6,3]))
    xc(:,4) = (x(:,1) + x(:,2) + x(:,3)) / 3.0_r8
    xc(:,5) = (x(:,4) + x(:,5) + x(:,6)) / 3.0_r8
  end function wedge_face_centers

  pure function hex_face_centers (x) result (xc)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: xc(3,6)
    xc(:,1) = polygon_center(x(:,[1,2,6,5]))
    xc(:,2) = polygon_center(x(:,[2,3,7,6]))
    xc(:,3) = polygon_center(x(:,[3,4,8,7]))
    xc(:,4) = polygon_center(x(:,[1,5,8,4]))
    xc(:,5) = polygon_center(x(:,[1,4,3,2]))
    xc(:,6) = polygon_center(x(:,[5,6,7,8]))
  end function hex_face_centers

  !! Returns the center of a polygonal face in 3D.  Except for triangles,
  !! a polygonal face is not generally planar and a practical definition
  !! of the center is arguably ambiguous.  This function computes the center
  !! of mass of a triangulated surface that approximates the polygonal face.
  !! The triangulation is formed by connecting the vertices of the polygon
  !! with the algebraic mean of the vertices.

  pure function polygon_center (x) result (xc)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: xc(3), xtri(3,3), asum, atri
    integer :: n, j
    n = size(x,dim=2)
    xtri(:,3) = sum(x,dim=2) / n
    asum = 0.0_r8; xc = 0.0_r8
    do j = 1, n
      xtri(:,1) = x(:,j)
      xtri(:,2) = x(:,modulo(j,n)+1)
      atri = tri_area(xtri)
      asum = asum + atri
      xc = xc + atri*sum(xtri,dim=2)/3.0_r8
    end do
    xc = xc / asum
  end function polygon_center

  !! Returns the center of mass of the given 3D cell.  X(:,k) are the
  !! coordinates of the kth vertiex of the cell.  The cell type is inferred
  !! from the number of vertices.  The function handles tet, pyramid, wedge,
  !! and hex cell types; a 0 result is returned for anything else. The center
  !! of mass of a tetrahedron is simply the algebraic mean of its vertices.
  !! The problem for the other cell types is handled by subdividing them into
  !! tetrahedra.  This is straightforward if they have planar faces, but this
  !! will generally not be the case.  This is where we punt and do something
  !! heuristic based on Ushakova's exact formula for the volume of a general
  !! hexahedron.  Noting that it is the average volume of two different
  !! subdivisions of a hexahedron into tetrahedrons (four corner tets plus the
  !! remaining internal tet -- done two ways), we use the same subdivisions to
  !! compute the center of mass.  In this case the formula is not exact, but
  !! it has the virtue that the errors at a warped quadrilateral face from the
  !! adjacent cells cancel when computing the center of mass for an assembly
  !! of cells. As for volumes, pyramid and wedge cells are regarded as
  !! degenerate hexahedra.

  pure function cell_center (x) result (c)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: c(3), svol, smom(3)
    select case (size(x,dim=2))
    case (4)  ! tet
      c = sum(x,dim=2)/4
    case (5)  ! pyramid
      svol = 0.0_r8; smom = 0.0_r8
      call aux (x(:,[1,2,4,5]), svol, smom)
      call aux (x(:,[2,3,1,5]), svol, smom)
      call aux (x(:,[3,4,2,5]), svol, smom)
      call aux (x(:,[4,1,3,5]), svol, smom)
      c = smom/svol
    case (6)  ! wedge
      svol = 0.0_r8; smom = 0.0_r8
      call aux (x(:,[1,2,3,4]), svol, smom)
      call aux (x(:,[5,4,6,2]), svol, smom)
      call aux (x(:,[2,3,4,6]), svol, smom)
      call aux (x(:,[2,3,1,5]), svol, smom)
      call aux (x(:,[4,6,5,1]), svol, smom)
      call aux (x(:,[1,3,6,5]), svol, smom)
      c = smom/svol
   case (8)  ! hex
      svol = 0.0_r8; smom = 0.0_r8
      call aux (x(:,[1,2,4,5]), svol, smom)
      call aux (x(:,[2,3,1,6]), svol, smom)
      call aux (x(:,[3,4,2,7]), svol, smom)
      call aux (x(:,[4,1,3,8]), svol, smom)
      call aux (x(:,[5,8,6,1]), svol, smom)
      call aux (x(:,[6,5,7,2]), svol, smom)
      call aux (x(:,[7,6,8,3]), svol, smom)
      call aux (x(:,[8,7,5,4]), svol, smom)
      call aux (x(:,[1,3,8,6]), svol, smom)
      call aux (x(:,[2,4,5,7]), svol, smom)
      c = smom/svol
    case default
      c = 0.0_r8
    end select
  contains
    !! Accumulates the moment and volume for the given cell.
    pure subroutine aux (x, svol, smom)
      real(r8), intent(in) :: x(:,:)
      real(r8), intent(inout) :: svol, smom(:)
      svol = svol + tet_volume(x)
      smom = smom + tet_volume(x) * sum(x,dim=2)/4
    end subroutine aux
  end function cell_center

  !! Returns the center of mass of the given 3D cell.  X(:,k) are the
  !! coordinates of the kth vertiex of the cell.  The cell type is inferred
  !! from the number of vertices.  The function handles tet, pyramid, wedge,
  !! and hex cell types; a 0 result is returned for anything else. The center
  !! of mass of a tetrahedron is simply the algebraic mean of its vertices.
  !! The other cell types are treated as degenerate hexahedra and we use the
  !! original Truchas algorithm (described below) which is supposedly exact
  !! for hexahedra with warped faces.

  pure function cell_centroid (x) result (c)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: c(3)
    select case (size(x,dim=2))
    case (4)  ! tet
      c = sum(x,dim=2)/4
    case (5)  ! pyramid
      call hex_centroid(x(:,[1,2,3,4,5,5,5,5]), pyramid_volume(x), c)
    case (6)  ! wedge
      call hex_centroid(x(:,[1,2,3,1,4,5,6,4]), wedge_volume(x), c)
   case (8)  ! hex
      call hex_centroid(x, hex_volume(x), c)
    case default
      c = 0.0_r8
    end select
  end function cell_centroid

  !! This is derived from the original cell_geometry_module::cell_centroid
  !! procedure.  It has been stripped down to operate on a single cell, whose
  !! vertex coordinates are passed, instead of the entire mesh, and code for
  !! handling 2D cells removed.  The algorithm, whose source is unknown, is
  !! itself unaltered and is supposedly exact for hexahedra with warped faces
  !! -- hexahedra that are the image of a reference cube under the canonical
  !! trilinear map. It was applied to both degenerate and non-degenerate hexes.
  !! NB: Recent investigation with mathematica (May 2020) have cast serious
  !! doubt on the exactness claim.

  pure subroutine hex_centroid(xv, vol, xc)

    real(r8), intent(in)  :: xv(:,:)  ! cell vertex coordinates
    real(r8), intent(in)  :: vol      ! cell volume
    real(r8), intent(out) :: xc(:)    ! cell centroid coordinates

    integer :: i, ip1, ip2
    real(r8) :: L(3), M(3), N(3), LxD3(3), MxD2(3), NxD1(3)
    real(r8) :: tmp(3), D1(3), D2(3), D3(3), Dv(3), D1xDv(3), D2xDv(3), D3xDv(3)

    L  = 0.25_r8*(xv(:,1) + xv(:,2) + xv(:,6) + xv(:,5) - xv(:,3) - xv(:,4) - xv(:,8) - xv(:,7))
    M  = 0.25_r8*(xv(:,2) + xv(:,3) + xv(:,7) + xv(:,6) - xv(:,4) - xv(:,1) - xv(:,5) - xv(:,8))
    N  = 0.25_r8*(xv(:,8) + xv(:,5) + xv(:,6) + xv(:,7) - xv(:,3) - xv(:,2) - xv(:,1) - xv(:,4))
    D1 = (xv(:,2) + xv(:,6) + xv(:,4) + xv(:,8) - xv(:,1) - xv(:,5) - xv(:,3) - xv(:,7))
    D2 = (xv(:,3) + xv(:,4) + xv(:,5) + xv(:,6) - xv(:,1) - xv(:,2) - xv(:,7) - xv(:,8))
    D3 = (xv(:,1) + xv(:,4) + xv(:,6) + xv(:,7) - xv(:,2) - xv(:,3) - xv(:,5) - xv(:,8))
    Dv = (xv(:,1) + xv(:,3) + xv(:,6) + xv(:,8) - xv(:,2) - xv(:,4) - xv(:,5) - xv(:,7))

    do i = 1, 3
      ip1 = modulo(i,3) + 1
      ip2 = modulo(ip1,3) + 1

      LxD3(i) = L(ip1)*D3(ip2) - L(ip2)*D3(ip1)
      MxD2(i) = M(ip1)*D2(ip2) - M(ip2)*D2(ip1)
      NxD1(i) = N(ip1)*D1(ip2) - N(ip2)*D1(ip1)

      D1xDv(i) = D1(ip1)*Dv(ip2) - D1(ip2)*Dv(ip1)
      D2xDv(i) = D2(ip1)*Dv(ip2) - D2(ip2)*Dv(ip1)
      D3xDv(i) = D3(ip1)*Dv(ip2) - D3(ip2)*Dv(ip1)
    end do

    tmp = 0.0_r8
    do i = 1, 3
      tmp(1) = tmp(1) + L(i)*(MxD2(i) - NxD1(i)) + (N(i)*D2xDv(i) - M(i)*D1xDv(i))/12
      tmp(2) = tmp(2) + M(i)*(NxD1(i) - LxD3(i)) + (L(i)*D1xDv(i) - N(i)*D3xDv(i))/12
      tmp(3) = tmp(3) + N(i)*(LxD3(i) - MxD2(i)) + (M(i)*D3xDv(i) - L(i)*D2xDv(i))/12
    end do
    tmp = 0.5_r8 + tmp/(24*vol)

    do i = 1, 3
      xc(i) = xv(i,4) + tmp(1)*(xv(i,1) - xv(i,4)) &
                      + tmp(2)*(xv(i,3) - xv(i,4)) &
                      + tmp(1)*tmp(2)*(xv(i,2) + xv(i,4) - xv(i,1) - xv(i,3)) &
                      + tmp(3)*(xv(i,8) - xv(i,4) &
                      + tmp(1)*(xv(i,4) + xv(i,5) - xv(i,1) - xv(i,8)) &
                      + tmp(2)*(xv(i,4) + xv(i,7) - xv(i,3) - xv(i,8)) &
                      + tmp(1)*tmp(2)*Dv(i))
    end do

  end subroutine hex_centroid

  !! Centroid of a convex polygon in the plane.
  pure function cell_centroid_2d(x) result(xc)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: xc(2)
    select case (size(x,dim=2))
    case (3)
      xc = (x(:,1) + x(:,2) + x(:,3)) / 3.0_r8
    case (4)
      block
        real(r8) :: a2, a4, f
        a2 = tri_area(x(:,1:3))
        a4 = tri_area(x(:,[1,3,4]))
        f = a2 / (a2 + a4)
        xc = (x(:,1) + x(:,3) + f*x(:,2) + (1-f)*x(:,4)) / 3.0_r8
      end block
    case (5:)
      block
        integer :: j
        real :: a, aj
        a = 0; xc = 0
        do j = 2, size(x,dim=2)-1
          aj = tri_area(x(:,[1,j,j+1]))
          a  = a + aj
          xc = xc + aj * (x(:,1) + x(:,j) + x(:,j+1))
        end do
        xc = xc / (3*a)
      end block
    case default
      xc = 0.0_r8
    end select
  end function cell_centroid_2d

!!!! ALGEBRAIC PRIMITIVES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Computes A x B
  pure function cross_product (a, b) result (axb)
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: axb(3)
    axb(1) = a(2)*b(3) - a(3)*b(2)
    axb(2) = a(3)*b(1) - a(1)*b(3)
    axb(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_product

  !! Computes A . B X C
  pure function triple_product (a, b, c) result (abc)
    real(r8), intent(in) :: a(:), b(:), c(:)
    real(r8) :: abc
    abc = a(1)*(b(2)*c(3) - b(3)*c(2)) + &
          a(2)*(b(3)*c(1) - b(1)*c(3)) + &
          a(3)*(b(1)*c(2) - b(2)*c(1))
  end function triple_product

  !! Computes the Euclidean length of a vector X(:).  This is really meant
  !! for 2 and 3-vectors where a numerically robust method is used, but it
  !! also works with 1 and n-vectors, n>2, delegating to the NORM2 intrinsic
  !! function in the latter case.

  pure function vector_length (x) result (l)

    real(r8), intent(in) :: x(:)
    real(r8) :: l, a, b, c, t

    select case (size(x))
    case (1)  ! 1-vector
      l = abs(x(1))
    case (2)  ! 2-vector
      a = abs(x(1)); b = abs(x(2))
      !! Swap largest value to A
      if (b > a) then
        t = a; a = b; b = t   ! swap A and B
      end if
      if (a == 0.0_r8) then
        l = 0.0_r8
      else
        l = a * sqrt(1.0_r8 + (b/a)**2)
      end if
    case (3)  ! 3-vector
      a = abs(x(1)); b = abs(x(2)); c = abs(x(3))
      !! Swap largest value to A
      if (b > a) then
        if (c > b) then
          t = a; a = c; c = t ! swap A and C
        else
          t = a; a = b; b = t ! swap A and B
        end if
      else if (c > a) then
        t = a; a = c; c = t   ! swap A and C
      end if
      if (a == 0.0_r8) then
        l = 0.0_r8
      else
        l = a * sqrt(1.0_r8 + ((b/a)**2 + (c/a)**2))
      end if
    case default  ! do nothing special for the general case
      l = norm2(x)
    end select

  end function vector_length

  !! Computes the area of a triangle given the coordinates of its three
  !! vertices: X(:,j) are the coordinates of vertex j. The vertices may
  !! be given in any order and the coordinates n-dimensional, n > 1.
  !! This uses the TRI_AREA_LENGTH function and thus the computed area
  !! is non-negative independent of the orientation of the triangle.

  pure function tri_area_coord (x) result (area)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: area, l(3)
    l(1) = vector_length(x(:,1)-x(:,2))
    l(2) = vector_length(x(:,2)-x(:,3))
    l(3) = vector_length(x(:,3)-x(:,1))
    area = tri_area_length(l)
  end function tri_area_coord

  !! Area of a quadrilateral in the plane.
  pure function quad_area(x) result(area)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: area
    area = 0.5_r8*( (x(1,3)-x(1,1))*(x(2,4)-x(2,2)) - (x(2,3)-x(2,1))*(x(1,4)-x(1,2)) )
  end function quad_area

  !! Computes the area of a triangle given the lengths of its three sides.
  !! Uses Kahan's algorithm from "Miscalculating area and angles of a needle-
  !! like triangle", 1986, http://www.cs.berkeley.edu/~wkahan/Triangle.pdf.
  !! It is essential that the parentheses are respected. The result is accurate
  !! to within a few ulps. The computed area is necessarily non-negative and is
  !! insensitive to the orientation of the triangle.
  !!
  !! NB: The implementation attempts to use the IEEE exception feature of
  !! Fortran to handle the case of invalid inputs.  The ieee_set_flag call
  !! should be enough to trigger an exception, but on some platforms it does
  !! not. Hence the following two lines which ensure an exception is raised.

  pure function tri_area_length (length) result (area)

    use,intrinsic :: ieee_arithmetic, only: ieee_set_flag, ieee_invalid, &
                                            ieee_value, ieee_signaling_nan

    real(r8), intent(in) :: length(:)
    real(r8) :: area, a, b, c, t

    !! Sort lengths so that A >= B >= C
    a = length(1); b = length(2); c = length(3)
    if (b > a) then
      if (c > a) then
        t = a; a = c; c = t   ! swap A and C
        if (b > a) then
          t = a; a = b; b = t ! swap A and B
        end if
      else
        t = a; a = b; b = t   ! swap A and B
      end if
    else
      if (c > b) then
        t = b; b = c; c = t   ! swap B and C
        if (b > a) then
          t = a; a = b; b = t ! swap A and B
        end if
      end if
    end if

    if (c-(a-b) >= 0.0_r8) then
      area = 0.25_r8 * sqrt((a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c)))
    else  ! not the lengths of a real triangle
      call ieee_set_flag (ieee_invalid, .true.)
      area = ieee_value(area, ieee_signaling_nan)
      area = area + area  ! see the NB above
    end if

  end function tri_area_length

  pure function normalized (v)

    real(r8), intent(in) :: v(:)
    real(r8) :: normalized(size(v))
    real(r8) :: mag

    mag = norm2(v)
    if (mag > 0) then
      normalized = v / mag
    else
      normalized = 0
    end if

  end function normalized


end module cell_geometry
