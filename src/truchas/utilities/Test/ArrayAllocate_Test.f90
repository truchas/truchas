PROGRAM ArrayAllocate_Test
  !=======================================================================
  ! Purpose(s):
  !   Test for ArrayAllocate routines.
  !
  ! Author(s): Bryan R. Lally, LANL ESA-EPE (lally@lanl.gov)
  !
  !=======================================================================
  use ArrayAllocate_Module
  use Output_Module

  implicit none

  ! Local Variables
  Real, pointer, dimension(:) :: p1
  Real, pointer, dimension(:,:) :: p2
  Real, pointer, dimension(:,:,:) :: p3

  Integer :: n1 = 2
  Integer :: n2 = 3
  Integer :: n3 = 4
  Integer :: i, j, k, cnt

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! initialize real, rank 1
  Call ArrayCreate(p1,1,n1,'P1n1')
  cnt = 0
  Do i = 1, n1
     cnt = cnt + 1
     p1(i) = cnt
  End Do

  ! initialize real, rank 2
  Call ArrayCreate(p2,1,n1,1,n2,'P2n1n2')
  cnt = 0
  Do i = 1, n1
     Do j = 1, n2
        cnt = cnt + 1
        p2(i,j) = cnt
     End Do
  End Do

  ! initialize real, rank 3
  Call ArrayCreate(p3,1,n1,1,n2,1,n3,'P3n1n2n3')
  cnt = 0
  Do i = 1, n1
     Do j = 1, n2
        Do k = 1, n3
           cnt = cnt + 1
           p3(i,j,k) = cnt
        End Do
     End Do
  End Do

  Write (*,*) 'dump of P1'
  Do i = 1, n1
     Write (*,*) 'p1(',i,') = ',p1(i)
  End Do
  Write (*,*)

  Write (*,*) 'dump of P2'
  Do i = 1, n1
     Do j = 1, n2
        Write (*,*) 'p2(',i,',',j,') = ',p2(i,j)
     End Do
  End Do
  Write (*,*)

  Write (*,*) 'dump of P3'
  Do i = 1, n1
     Do j = 1, n2
        Do k = 1, n3
           Write (*,*) 'p3(',i,',',j,',',k,') = ',p3(i,j,k)
        End Do
     End Do
  End Do

! deallocate arrays
  Call ArrayDestroy(p1,'P1')
  Call ArrayDestroy(p2,'P2')
  Call ArrayDestroy(p3,'P3')

END PROGRAM ArrayAllocate_Test
