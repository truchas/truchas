!!
!! Cell geometry procedures from the legacy mesh
!!

#include "f90_assert.fpp"

module legacy_geometry

  use kinds, only: r8
  implicit none
  private

  public :: cell_centroid, face_area, face_centroid_logical

contains

  !! This is derived from the original cell_geometry_module::cell_centroid
  !! procedure.  It has been stripped down to operate on a single cell, whose
  !! vertex coordinates are passed, instead of the entire mesh, and code for
  !! handling 2D cells removed.  The algorithm, whose source is unknown, is
  !! itself unaltered.  It was applied to both degenerate and non-degenerate
  !! hexes.

  subroutine cell_centroid (xv, vol, xc)

    real(r8), intent(in)  :: xv(:,:)  ! cell vertex coordinates
    real(r8), intent(in)  :: vol      ! cell volume
    real(r8), intent(out) :: xc(:)    ! cell centroid coordinates

    integer :: i, ip1, ip2
    real(r8) :: L(3), M(3), N(3), LxD3(3), MxD2(3), NxD1(3)
    real(r8) :: tmp(3), D1(3), D2(3), D3(3), Dv(3), D1xDv(3), D2xDv(3), D3xDv(3)

    ASSERT(size(xv,1) == 3)
    ASSERT(size(xv,2) == 8)
    ASSERT(size(xc) == 3)

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

  end subroutine cell_centroid

  !! This is derived from the original cell_geometry_module::face_area
  !! procedure. It has been stripped down to operate on a single cell, whose
  !! vertex coordinates are passed, instead of the entire mesh, and code for
  !! handling 2D cells removed.  The algorithm is itself unaltered.  It was
  !! applied to both degenerate and non-degenerate hexes; area and normal for
  !! degenerate faces is set to 0.

  subroutine face_area (xv, area, normal)

    use cutoffs_module, only: alittle

    real(r8), intent(in)  :: xv(:,:)
    real(r8), intent(out) :: area(:), normal(:,:)

    integer :: f, i, v1, v2, v3, v4
    real(r8) :: x1(3), x2(3)

    ASSERT(size(xv,1) == 3)
    ASSERT(size(xv,2) == 8)
    ASSERT(size(area) == 6)
    ASSERT(size(normal,1) == 3)
    ASSERT(size(normal,2) == 6)

    do f = 1, 6

      select case (f)
      case (1)
        v1 = 8; v2 = 3; v3 = 7; v4 = 4
      case (2)
        v1 = 6; v2 = 1; v3 = 5; v4 = 2
      case (3)
        v1 = 5; v2 = 4; v3 = 8; v4 = 1
      case (4)
        v1 = 7; v2 = 2; v3 = 6; v4 = 3
      case (5)
        v1 = 3; v2 = 1; v3 = 2; v4 = 4
      case (6)
        v1 = 6; v2 = 8; v3 = 7; v4 = 5
      end select

      x1 = xv(:,v1) - xv(:,v2)
      x2 = xv(:,v3) - xv(:,v4)

      do i = 1, 3
        v1 = modulo(i,3) + 1
        v2 = modulo(v1,3) + 1
        normal(i,f) = 0.5_r8*(x1(v1)*x2(v2) - x2(v1)*x1(v2))
      end do

      normal(:,f) = merge(0.0_r8, normal(:,f), mask=(abs(normal(:,f)) < alittle))
      area(f) = norm2(normal(:,f))

      if (area(f) >= alittle) then
        normal(:,f) = normal(:,f) / area(f)
      else ! presumably a degenerate face
        normal(:,f) = 0.0_r8
        area(f) = 0.0_r8
      end if

    end do

  end subroutine face_area

  !! This is derived from the original cell_geometry_module::face_centroid_logical
  !! procedure. It has been stripped down to operate on a single cell, whose
  !! vertex coordinates are passed, instead of the entire mesh, and code for
  !! handling 2D cells removed.  The algorithm is itself unaltered.  It was
  !! applied to both degenerate and non-degenerate hexes.  The source of this
  !! algorithm and precisely what it is supposed to do are unknown.

  subroutine face_centroid_logical (xv, area, normal, centroid_l)

    use cutoffs_module, only: alittle

    real(r8), intent(in)  :: xv(:,:), area(:), normal(:,:)
    real(r8), intent(out) :: centroid_l(:,:)

    integer :: f, i, ip1, ip2
    integer :: v11, v12, v13, v14, v21, v22, v23, v24, v31, v32, v33, v34
    real(r8) :: coef(3), face_area

    ASSERT(size(xv,1) == 3)
    ASSERT(size(xv,2) == 8)
    ASSERT(size(area) == 6)
    ASSERT(size(normal,1) == 3)
    ASSERT(size(normal,2) == 6)
    ASSERT(size(centroid_l,1) == 3)
    ASSERT(size(centroid_l,2) == 6)

    do f = 1, 6
      centroid_l(:,f) = 0.5_r8
      coef = 1.0_r8
      face_area = 1.0_r8 / (12*(area(f) + alittle))

      select case (f)
      case (1) ! side 1 (vertices 4-8-7-3)
        v11 = 1; v12 = 1; v13 = 1; v14 = 1
        v21 = 7; v22 = 8; v23 = 3; v24 = 4
        v31 = 8; v32 = 4; v33 = 7; v34 = 3
        centroid_l(1,f) = 0.0_r8; coef(1) = 0.0_r8
      case (2) ! side 2 (vertices 1-2-6-5)
        v11 = 1; v12 = 1; v13 = 1; v14 = 1
        v21 = 2; v22 = 1; v23 = 6; v24 = 5
        v31 = 6; v32 = 2; v33 = 5; v34 = 1
        centroid_l(1,f) = 1.0_r8; coef(1) = 0.0_r8
      case (3) ! side 3 (vertices 4-1-5-8)
        v11 = 1; v12 = 4; v13 = 5; v14 = 8
        v21 = 1; v22 = 1; v23 = 1; v24 = 1
        v31 = 5; v32 = 1; v33 = 8; v34 = 4
        centroid_l(2,f) = 0.0_r8; coef(2) = 0.0_r8
      case (4) ! side 4 (vertices 3-7-6-2)
        v11 = 6; v12 = 7; v13 = 2; v14 = 3
        v21 = 1; v22 = 1; v23 = 1; v24 = 1
        v31 = 7; v32 = 3; v33 = 6; v34 = 2
        centroid_l(2,f) = 1.0_r8; coef(2) = 0.0_r8
      case (5) ! side 5 (vertices 4-3-2-1)
        v11 = 2; v12 = 3; v13 = 1; v14 = 4
        v21 = 3; v22 = 4; v23 = 2; v24 = 1
        v31 = 1; v32 = 1; v33 = 1; v34 = 1
        centroid_l(3,f) = 0.0_r8; coef(3) = 0.0_r8
      case (6) ! side 6 (vertices 8-5-6-7)
        v11 = 5; v12 = 8; v13 = 6; v14 = 7
        v21 = 6; v22 = 5; v23 = 7; v24 = 8
        v31 = 1; v32 = 1; v33 = 1; v34 = 1
        centroid_l(3,f) = 1.0_r8; coef(3) = 0.0_r8
      end select

      do i = 1, 3
        ip1 = modulo(i,3) + 1
        ip2 = modulo(ip1,3) + 1
        centroid_l(1,f) = centroid_l(1,f) + coef(1)*face_area*normal(i,f)* &
              ((xv(ip1,v11) - xv(ip1,v12))*(xv(ip2,v13) - xv(ip2,v14)) - &
               (xv(ip2,v11) - xv(ip2,v12))*(xv(ip1,v13) - xv(ip1,v14)))
        centroid_l(2,f) = centroid_l(2,f) + coef(2)*face_area*normal(i,f)* &
              ((xv(ip1,v21) - xv(ip1,v22))*(xv(ip2,v23) - xv(ip2,v24)) - &
               (xv(ip2,v21) - xv(ip2,v22))*(xv(ip1,v23) - xv(ip1,v24)))
        centroid_l(3,f) = centroid_l(3,f) + coef(3)*face_area*normal(i,f)* &
              ((xv(ip1,v31) - xv(ip1,v32))*(xv(ip2,v33) - xv(ip2,v34)) - &
               (xv(ip2,v31) - xv(ip2,v32))*(xv(ip1,v33) - xv(ip1,v34)))
      end do
    end do

  end subroutine face_centroid_logical

end module legacy_geometry
