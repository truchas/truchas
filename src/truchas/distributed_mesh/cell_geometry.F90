!!
!! CELL_GEOMETRY
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 9 Feb 2006
!!

#include "f90_assert.fpp"

module cell_geometry

  use kinds
  implicit none
  private

  public :: cell_volume, tet_volume, hex_volume
  public :: face_normal, tri_face_normal, quad_face_normal
  public :: tet_face_normals, hex_face_normals, eval_hex_volumes, wedge_face_normals
  public :: edge_length
  public :: cross_product, triple_product, vector_length
  public :: cell_face_normals, cell_face_centers
  public :: cell_center

contains

  pure function cross_product (a, b) result (axb)
    real(kind=r8), intent(in) :: a(:), b(:)
    real(kind=r8) :: axb(3)
    axb(1) = a(2)*b(3) - a(3)*b(2)
    axb(2) = a(3)*b(1) - a(1)*b(3)
    axb(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_product

  pure function triple_product (a, b, c) result (abc)
    real(kind=r8), intent(in) :: a(:), b(:), c(:)
    real(kind=r8) :: abc
    abc = a(1)*(b(2)*c(3) - b(3)*c(2)) + &
          a(2)*(b(3)*c(1) - b(1)*c(3)) + &
          a(3)*(b(1)*c(2) - b(2)*c(1))
  end function triple_product
  
  pure function cell_volume (x) result (vol)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: vol
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
  end function cell_volume

  pure function tet_volume (x) result (vol)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: vol
    vol = triple_product(x(:,2)-x(:,1), x(:,3)-x(:,1), x(:,4)-x(:,1)) / 6.0_r8
  end function tet_volume
  
  pure function pyramid_volume (x) result (vol)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: vol, cvol1, cvol2, cvol3, cvol4
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

    cvol(1) = tet_volume(x(:,[1,2,4,5]))
    cvol(2) = tet_volume(x(:,[2,3,1,6]))
    cvol(3) = tet_volume(x(:,[3,4,2,7]))
    cvol(4) = tet_volume(x(:,[4,1,3,8]))
    cvol(5) = tet_volume(x(:,[5,8,6,1]))
    cvol(6) = tet_volume(x(:,[6,5,7,2]))
    cvol(7) = tet_volume(x(:,[7,6,8,3]))
    cvol(8) = tet_volume(x(:,[8,7,5,4]))

    hvol = 0.5_r8 * (sum(cvol) + tet_volume(x(:,[1,3,8,6])) + tet_volume(x(:,[2,4,5,7])))

  end function hex_volume

  pure subroutine eval_hex_volumes (x, hvol, cvol)

    real(r8), intent(in)  :: x(:,:)
    real(r8), intent(out) :: hvol, cvol(:)

    !ASSERT(size(x,dim=1) == 3)
    !ASSERT(size(x,dim=2) == 8)
    !ASSERT(size(cvol) == 8)

    !! Corner tet volumes.
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

    ! incompatible with PURE
    !ASSERT(size(x,dim=1) == 3)
    !ASSERT(size(x,dim=2) == 4)

    !! NB: These must be consistent with the TETRA4 vertex and face labelings
    !! defined in the CELL_TOPOLOGY module.  To avoid a layer of indirection,
    !! its TETRA4_FACE_VERT array was not used here.
    a(:,1) = 0.5_r8 * cross_product(x(:,1)-x(:,4), x(:,2)-x(:,4))
    a(:,2) = 0.5_r8 * cross_product(x(:,2)-x(:,4), x(:,3)-x(:,4))
    a(:,3) = 0.5_r8 * cross_product(x(:,3)-x(:,4), x(:,1)-x(:,4))
    a(:,4) = 0.5_r8 * cross_product(x(:,3)-x(:,1), x(:,2)-x(:,1))

  end function tet_face_normals

  pure function pyramid_face_normals (x) result (a)

    real(r8), intent(in) :: x(:,:)
    real(r8) :: a(3,5)

    ! incompatible with PURE
    !ASSERT(size(x,dim=1) == 3)
    !ASSERT(size(x,dim=2) == 5)

    !! NB: These must be consistent with the PYR5 vertex and face labelings
    !! defined in the CELL_TOPOLOGY module.  To avoid a layer of indirection,
    !! its PYR5_FACE_VERT array was not used here.
    a(:,1) = 0.5_r8 * cross_product(x(:,2)-x(:,1), x(:,5)-x(:,1))
    a(:,2) = 0.5_r8 * cross_product(x(:,3)-x(:,2), x(:,5)-x(:,2))
    a(:,3) = 0.5_r8 * cross_product(x(:,4)-x(:,3), x(:,5)-x(:,3))
    a(:,4) = 0.5_r8 * cross_product(x(:,1)-x(:,4), x(:,5)-x(:,4))
    a(:,5) = 0.5_r8 * cross_product(x(:,3)-x(:,1), x(:,2)-x(:,4))

  end function pyramid_face_normals

  pure function wedge_face_normals (x) result (a)

    real(r8), intent(in) :: x(:,:)
    real(r8) :: a(3,5)

    ! incompatible with PURE
    !ASSERT(size(x,dim=1) == 3)
    !ASSERT(size(x,dim=2) == 6)

    !! NB: These must be consistent with the WED6 vertex and face labelings
    !! defined in the CELL_TOPOLOGY module.  To avoid a layer of indirection,
    !! its WED6_FACE_VERT array was not used here.
    a(:,1) = 0.5_r8 * cross_product(x(:,5)-x(:,1), x(:,4)-x(:,2))
    a(:,2) = 0.5_r8 * cross_product(x(:,6)-x(:,2), x(:,5)-x(:,3))
    a(:,3) = 0.5_r8 * cross_product(x(:,4)-x(:,3), x(:,6)-x(:,1))
    a(:,4) = 0.5_r8 * cross_product(x(:,3)-x(:,1), x(:,2)-x(:,1))
    a(:,5) = 0.5_r8 * cross_product(x(:,5)-x(:,4), x(:,6)-x(:,4))

  end function wedge_face_normals
  
  pure function hex_face_normals (x) result (a)

    real(r8), intent(in) :: x(:,:)
    real(r8) :: a(3,6)

    ! incompatible with PURE
    !ASSERT(size(x,dim=1) == 3)
    !ASSERT(size(x,dim=2) == 8)

    !! NB: These must be consistent with the HEX8 vertex and face labelings
    !! defined in the CELL_TOPOLOGY module.  To avoid a layer of indirection,
    !! its HEX8_FACE_VERT array was not used here.
    a(:,1) = 0.5_r8 * cross_product(x(:,6)-x(:,1), x(:,5)-x(:,2))
    a(:,2) = 0.5_r8 * cross_product(x(:,7)-x(:,2), x(:,6)-x(:,3))
    a(:,3) = 0.5_r8 * cross_product(x(:,8)-x(:,3), x(:,7)-x(:,4))
    a(:,4) = 0.5_r8 * cross_product(x(:,5)-x(:,4), x(:,8)-x(:,1))
    a(:,5) = 0.5_r8 * cross_product(x(:,3)-x(:,1), x(:,2)-x(:,4))
    a(:,6) = 0.5_r8 * cross_product(x(:,7)-x(:,5), x(:,8)-x(:,6))

  end function hex_face_normals

  pure function tri_face_normal (x) result (a)
    real(kind=r8), intent(in) :: x(:,:)
    real(kind=r8) :: a(3)
    a = 0.5_r8 * cross_product(x(:,2)-x(:,1),x(:,3)-x(:,2))
  end function tri_face_normal

  pure function quad_face_normal (x) result (a)
    real(kind=r8), intent(in) :: x(:,:)
    real(kind=r8) :: a(3)
    a = 0.5_r8 * cross_product(x(:,3)-x(:,1), x(:,4)-x(:,2))
  end function quad_face_normal

  pure function face_normal (x) result (a)
    real(kind=r8), intent(in) :: x(:,:)
    real(kind=r8) :: a(size(x,dim=1))
    select case (size(x,dim=1))
    case (2)  ! 2-D coordinate space
      select case (size(x,dim=2))
      case (2)  ! just an edge
        a(1) = x(2,2) - x(2,1)
        a(2) = x(1,1) - x(1,2)
      case default
      end select
    case (3)  ! 3-D coordinate space
      select case (size(x,dim=2))
      case (3)  ! triangular face
        a = 0.5_r8 * cross_product(x(:,2)-x(:,1),x(:,3)-x(:,2))
      case (4)  ! quadrilateral face
        a = 0.5_r8 * cross_product(x(:,3)-x(:,1), x(:,4)-x(:,2))
      case default
      end select
    case default
    end select
  end function face_normal
  
  pure function cell_face_centers (x) result (xc)
    real(r8), intent(in) :: x(:,:)
    real(r8), allocatable :: xc(:,:)
    select case (size(x,dim=2))
    case (4)
      allocate(xc(3,4))
      xc(:,1) = (x(:,1) + x(:,2) + x(:,4)) / 3.0_r8
      xc(:,2) = (x(:,2) + x(:,3) + x(:,4)) / 3.0_r8
      xc(:,3) = (x(:,1) + x(:,3) + x(:,4)) / 3.0_r8
      xc(:,4) = (x(:,1) + x(:,2) + x(:,3)) / 3.0_r8
    case (5)
      allocate(xc(3,5))
      xc(:,1) = (x(:,1) + x(:,2) + x(:,5)) / 3.0_r8
      xc(:,2) = (x(:,2) + x(:,3) + x(:,5)) / 3.0_r8
      xc(:,3) = (x(:,3) + x(:,4) + x(:,5)) / 3.0_r8
      xc(:,4) = (x(:,1) + x(:,4) + x(:,5)) / 3.0_r8
      xc(:,5) = polygon_center(x(:,[1,4,3,2]))
    case (6)
      allocate(xc(3,5))
      xc(:,1) = polygon_center(x(:,[1,2,5,4]))
      xc(:,2) = polygon_center(x(:,[2,3,6,5]))
      xc(:,3) = polygon_center(x(:,[1,4,6,3]))
      xc(:,4) = (x(:,1) + x(:,2) + x(:,3)) / 3.0_r8
      xc(:,5) = (x(:,4) + x(:,5) + x(:,6)) / 3.0_r8
    case (8)
      allocate(xc(3,6))
      xc(:,1) = polygon_center(x(:,[1,2,6,5]))
      xc(:,2) = polygon_center(x(:,[2,3,7,6]))
      xc(:,3) = polygon_center(x(:,[3,4,8,7]))
      xc(:,4) = polygon_center(x(:,[1,5,8,4]))
      xc(:,5) = polygon_center(x(:,[1,4,3,2]))
      xc(:,6) = polygon_center(x(:,[5,6,7,8]))
    end select
  end function cell_face_centers
  
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


  pure function vector_length (x) result (l)
  
    real(kind=r8), intent(in) :: x(:)
    real(kind=r8) :: l
    
    real(kind=r8) :: a, b, c, t
    
    select case (size(x))
    case (2)  ! two-vector
      a = abs(x(1))
      b = abs(x(2))
      !! Swap largest value to A.
      if (b > a) then
        t = a
        a = b
        b = t
      end if
      l = a * sqrt(1.0_r8 + (b/a)**2)
    case (3)  ! three-vector
      a = abs(x(1))
      b = abs(x(2))
      c = abs(x(3))
      !! Swap largest value to A.
      if (b > a) then
        if (c > b) then
          t = a
          a = c
          c = t
        else
          t = a
          a = b
          b = t
        end if
      else if (c > a) then
        t = a
        a = c
        c = t
      end if
      if (a == 0.0_r8) then
        l = 0.0_r8
      else
        l = a * sqrt(1.0_r8 + ((b/a)**2 + (c/a)**2))
      end if
    case default  ! FOR ANYTHING ELSE WE RETURN A BOGUS VALUE
      l = -huge(1.0_r8)
    end select
    
  end function vector_length
  
  
  pure function edge_length (x) result (l)
    real(kind=r8), intent(in) :: x(:,:)
    real(kind=r8) :: l
    l = vector_length (x(:,1)-x(:,2))
  end function edge_length


  pure function tri_area (x) result (area)
  
    real(kind=r8), intent(in) :: x(:,:)
    real(kind=r8) :: area, l(3)
    
    l(1) = vector_length(x(:,1)-x(:,2))
    l(2) = vector_length(x(:,2)-x(:,3))
    l(3) = vector_length(x(:,3)-x(:,1))
    area = tri_area_l (l)
    
  end function tri_area
  
  pure function tri_area_l (l) result (area)
  
    real(r8), intent(in) :: l(3)
    real(r8) :: area, a, b, c, t
    
    a = l(1)
    b = l(2)
    c = l(3)
    
    ! Sort so that a >= b >= c
    if (b > a) then
      if (c > a) then
        t = a
        a = c
        c = t
        if (b > a) then
          t = a
          a = b
          b = t
        end if
      else
        t = a
        a = b
        b = t
      end if
    else
      if (c > b) then
        t = b
        b = c
        c = t
        if (b > a) then
          t = a
          a = b
          b = t
        end if
      end if
    end if
    
    if (c-(a-b) < 0.0_r8) then ! not the lengths of a real triangle
      area = 2.0_r8               ! trigger a floating point exception;
      area = sqrt(1.0_r8 - area)  ! must be a bit obtuse about it.
    else
      area = 0.25_r8 * sqrt((a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c)))
    end if

  end function tri_area_l

  pure function cell_center (x) result (c)

    real(r8), intent(in) :: x(:,:)
    real(r8) :: c(3), smom(3), svol

    select case (size(x,dim=2))
    case (4)  ! tet
      c = sum(x,dim=2)/4

    case (5)  ! pyramid
      smom = 0.0_r8; svol = 0.0_r8
      call aux (x(:,[1,2,4,5]), smom, svol)
      call aux (x(:,[2,3,1,5]), smom, svol)
      call aux (x(:,[3,4,2,5]), smom, svol)
      call aux (x(:,[4,1,3,5]), smom, svol)
      c = smom/svol

    case (6)  ! wedge
      smom = 0.0_r8; svol = 0.0_r8
      call aux (x(:,[1,2,3,4]), smom, svol)
      call aux (x(:,[5,4,6,2]), smom, svol)
      call aux (x(:,[2,3,4,6]), smom, svol)
      call aux (x(:,[2,3,1,5]), smom, svol)
      call aux (x(:,[4,6,5,1]), smom, svol)
      call aux (x(:,[1,3,6,5]), smom, svol)
      c = smom/svol

   case (8)  ! hex
      smom = 0.0_r8; svol = 0.0_r8
      call aux (x(:,[1,2,4,5]), smom, svol)
      call aux (x(:,[2,3,1,6]), smom, svol)
      call aux (x(:,[3,4,2,7]), smom, svol)
      call aux (x(:,[4,1,3,8]), smom, svol)
      call aux (x(:,[5,8,6,1]), smom, svol)
      call aux (x(:,[6,5,7,2]), smom, svol)
      call aux (x(:,[7,6,8,3]), smom, svol)
      call aux (x(:,[8,7,5,4]), smom, svol)
      call aux (x(:,[1,3,8,6]), smom, svol)
      call aux (x(:,[2,4,5,7]), smom, svol)
      c = smom/svol

    case default
      c = 0.0_r8
    end select

  contains

    pure subroutine aux (x, smom, svol)
      real(r8), intent(in) :: x(:,:)
      real(r8), intent(inout) :: smom(:), svol
      svol = svol + tet_volume(x)
      smom = smom + tet_volume(x) * sum(x,dim=2)/4
    end subroutine aux

  end function cell_center

end module cell_geometry
