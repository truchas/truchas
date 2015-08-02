module simplex_geometry

  use kinds, only: r8
  implicit none
  private
  
  !public :: vector_length
  public :: edge_length
  public :: tri_area
  public :: tet_volume
  public :: tet_face_normal
  public :: bary_coord
  
  interface tri_area
    procedure tri_area_length, tri_area_coord
  end interface

contains
  
  !! Computes the area-weighted normal vectors to the faces of a tet given the
  !! coordinates of its four vertices. X(:,j) are the coordinates of the jth
  !! vertex, and P(:,j) is the area-weighted normal to the jth face. If the tet
  !! is positively oriented (positive volume) the normals are outward oriented;
  !! otherwise they are inward oriented.  The order of the vertices defines the
  !! orientation of the tet.  NB: the implicit areas of the faces (lengths of
  !! the normals) will not be exactly identical with what would be computed by
  !! TRI_AREA, which uses a different method of computing area.
  pure function tet_face_normal (x) result (p)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: p(3,4)
    p(:,1) = 0.5_r8 * cross_product(x(:,2)-x(:,4), x(:,3)-x(:,4))
    p(:,2) = 0.5_r8 * cross_product(x(:,3)-x(:,4), x(:,1)-x(:,4))
    p(:,3) = 0.5_r8 * cross_product(x(:,1)-x(:,4), x(:,2)-x(:,4))
    p(:,4) = 0.5_r8 * cross_product(x(:,3)-x(:,1), x(:,2)-x(:,1))
  end function tet_face_normal
  
  pure function cross_product (a, b) result (axb)
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: axb(3)
    axb(1) = a(2)*b(3) - a(3)*b(2)
    axb(2) = a(3)*b(1) - a(1)*b(3)
    axb(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_product
    
  
  pure function tet_volume (x) result (vol)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: vol
    vol = triple_product(x(:,2)-x(:,1), x(:,3)-x(:,1), x(:,4)-x(:,1)) / 6.0_r8
  end function tet_volume

  pure function triple_product (a, b, c) result (abc)
    real(r8), intent(in) :: a(:), b(:), c(:)
    real(r8) :: abc
    abc = a(1)*(b(2)*c(3) - b(3)*c(2)) + &
          a(2)*(b(3)*c(1) - b(1)*c(3)) + &
          a(3)*(b(1)*c(2) - b(2)*c(1))
  end function triple_product

  !! Computes the area of a triangle given the lengths of its three sides.
  !! Uses Kahan's algorithm: "Miscalculating area and angles of a needle-like
  !! triangle", 1986, http://www.cs.berkeley.edu/~wkahan/Triangle.pdf. It is
  !! essential that the parentheses are respected. The result is accurate to
  !! within a few ulps.

  pure function tri_area_length (l) result (area)
  
    use,intrinsic :: ieee_arithmetic, only: ieee_invalid, ieee_signaling_nan, ieee_value, ieee_set_flag
  
    real(r8), intent(in) :: l(:)
    real(r8) :: area, a, b, c, t
    
    !! Sort lengths so that A >= B >= C
    a = l(1); b = l(2); c = l(3)
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
      area = area + area  ! ensure an exception (NB: should not be necessary)
    end if

  end function tri_area_length
  
  !! Computes the area of a triangle given the coordinates of its three
  !! vertices: X(:,j) are the coordinates of vertex j. The vertices may
  !! be given in any order and the coordinates n-dimensional, n > 1.
  pure function tri_area_coord (x) result (area)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: area, l(3)
    l(1) = vector_length(x(:,1)-x(:,2))
    l(2) = vector_length(x(:,2)-x(:,3))
    l(3) = vector_length(x(:,3)-x(:,1))
    area = tri_area_length(l)
  end function tri_area_coord

  !! Computes the length of an edge given the coordinates of its two end
  !! points: X(:,j) are the n-dimensional coordinates of the jth end point.
  pure function edge_length (x) result (l)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: l
    l = vector_length (x(:,1)-x(:,2))
  end function edge_length
  
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! BARYCENTRIC COORDINATE PROCEDURES
!!
!! Given the vertex coordinates of a 2- or 3-simplex (i.e., triangle or tetrahedron) and a
!! point, this routine computes the barycentric coordinates of the point.  It is assumed
!! that the point lies in the simplex (if not the coordinates will not lie in [0,1]), and
!! that the 2-simplex is in R^2 or 3-simplex is in R^3.  For a simplex embedded in a higher
!! dimension space, (new) routines could be written that instead solve the normal equations
!! for the underdetermined linear system.
!!
  
  pure function bary_coord (x, p) result (lambda)
  
    real(r8), intent(in) :: x(:,:) ! tetrahedron vertex coordinates
    real(r8), intent(in) :: p(:)   ! point in tetrahedron
    real(r8) :: lambda(size(p)+1)  ! barycentric coordinates of P
    
    integer :: n, j
    real(r8) :: det, y(size(p),size(p)), q(size(p))
    
    n = size(p)
    
    forall (j = 1:n) y(:,j) = x(:,j) - x(:,n+1)
    q = p - x(:,n+1)
    
    !! Mindless application of Cramer's rule...
    select case (n)
    case (2)
      det = y(1,1)*y(2,2) - y(2,1)*y(1,2)
      lambda(1) = (q(1) * y(2,2) - q(2) * y(1,2)) / det
      lambda(2) = (q(2) * y(1,1) - q(1) * y(2,1)) / det
      lambda(3) = 1.0_r8 - lambda(1) - lambda(2)
    
    case (3)
      !det = triple_product(y(:,1), y(:,2), y(:,3))
      det = y(1,1)*(y(2,2)*y(3,3) - y(3,2)*y(2,3)) + &
            y(2,1)*(y(3,2)*y(1,3) - y(1,2)*y(3,3)) + &
            y(3,1)*(y(1,2)*y(2,3) - y(2,2)*y(1,3))
      !lambda(1) = triple_product(q, y(:,2), y(:,3))
      lambda(1) = (q(1) * (y(2,2)*y(3,3) - y(3,2)*y(2,3)) + &
                   q(2) * (y(3,2)*y(1,3) - y(1,2)*y(3,3)) + &
                   q(3) * (y(1,2)*y(2,3) - y(2,2)*y(1,3))) / det

      !lambda(2) = triple_product(q, y(:,3), y(:,1))
      lambda(2) = (q(1) * (y(2,3)*y(3,1) - y(3,3)*y(2,1)) + &
                   q(2) * (y(3,3)*y(1,1) - y(1,3)*y(3,1)) + &
                   q(3) * (y(1,3)*y(2,1) - y(2,3)*y(1,1))) / det

      !lambda(3) = triple_product(q, y(:,1), y(:,2))
      lambda(3) = (q(1) * (y(2,1)*y(3,2) - y(3,1)*y(2,2)) + &
                   q(2) * (y(3,1)*y(1,2) - y(1,1)*y(3,2)) + &
                   q(3) * (y(1,1)*y(2,2) - y(2,1)*y(1,2))) / det

      lambda(4) = 1.0_r8 - lambda(1) - lambda(2) - lambda(3)
    end select
    
  end function bary_coord

end module simplex_geometry
