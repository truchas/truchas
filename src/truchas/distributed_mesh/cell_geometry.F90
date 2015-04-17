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
  public :: tet_face_normals, hex_face_normals, eval_hex_volumes
  public :: edge_length
  public :: cross_product, triple_product, vector_length

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
    case (8)  ! hex
      vol = hex_volume(x)
    case default
      vol = 0.0_r8
    end select
  end function cell_volume

  pure function tet_volume (x) result (tvol)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: tvol
    tvol = triple_product(x(:,2)-x(:,1), x(:,3)-x(:,1), x(:,4)-x(:,1)) / 6.0_r8
  end function tet_volume
  
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

  subroutine eval_hex_volumes (x, hvol, cvol)

    real(kind=r8), intent(in)  :: x(:,:)
    real(kind=r8), intent(out) :: hvol
    real(kind=r8), intent(out) :: cvol(:)

    ASSERT( size(x,dim=1) == 3 )
    ASSERT( size(x,dim=2) == 8 )
    ASSERT( size(cvol) == 8 )

    !! Corner tet volumes.
    cvol(1) = tet_volume(x(:,(/1,2,4,5/)))
    cvol(2) = tet_volume(x(:,(/2,3,1,6/)))
    cvol(3) = tet_volume(x(:,(/3,4,2,7/)))
    cvol(4) = tet_volume(x(:,(/4,1,3,8/)))
    cvol(5) = tet_volume(x(:,(/5,8,6,1/)))
    cvol(6) = tet_volume(x(:,(/6,5,7,2/)))
    cvol(7) = tet_volume(x(:,(/7,6,8,3/)))
    cvol(8) = tet_volume(x(:,(/8,7,5,4/)))

    hvol = 0.5_r8 * (sum(cvol) + tet_volume(x(:,(/1,3,8,6/))) + tet_volume(x(:,(/2,4,5,7/))))

  end subroutine eval_hex_volumes

  pure function tet_face_normals (x) result (a)
    real(r8), intent(in) :: x(:,:)
    real(r8) :: a(3,4)
    a(:,1) = 0.5_r8 * cross_product(x(:,2)-x(:,4), x(:,3)-x(:,4))
    a(:,2) = 0.5_r8 * cross_product(x(:,3)-x(:,4), x(:,1)-x(:,4))
    a(:,3) = 0.5_r8 * cross_product(x(:,1)-x(:,4), x(:,2)-x(:,4))
    a(:,4) = 0.5_r8 * cross_product(x(:,3)-x(:,1), x(:,2)-x(:,1))
  end function tet_face_normals
  
  function hex_face_normals (x) result (a)

    real(kind=r8), intent(in) :: x(:,:)
    real(kind=r8) :: a(3,6)

    ASSERT( size(x,dim=1) == 3 )
    ASSERT( size(x,dim=2) == 8 )

    a(:,1) = 0.5_r8 * cross_product(x(:,8)-x(:,3), x(:,7)-x(:,4))
    a(:,2) = 0.5_r8 * cross_product(x(:,6)-x(:,1), x(:,5)-x(:,2))
    a(:,3) = 0.5_r8 * cross_product(x(:,5)-x(:,4), x(:,8)-x(:,1))
    a(:,4) = 0.5_r8 * cross_product(x(:,7)-x(:,2), x(:,6)-x(:,3))
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
  
end module cell_geometry
