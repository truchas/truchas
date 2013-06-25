MODULE DO_SOLVE_MODULE
  use constants_module, only: zero, one, two
  use truchas_logging_services, only: TLS_panic
  use kind_module,      only: int_kind, real_kind

 implicit none
 private

 public :: DO_SV_DECOMP,DO_LU_SOLVE,DO_LU_DECOMP

CONTAINS


 SUBROUTINE DO_SV_DECOMP(a,mrowa,mrow,ncol,w,v)
! ******************************************************************************
! *                                                                            *
! *   SUBROUTINE  :  Singular value decomposition solution                     *
! *                  Derived from algorithm in Numerical Recipes               *
! *                                                                            *
! *                  A ..... Input M by N matrix; Output M by N matrix "U"     *
! *                  MROWa.. Number of rows of A (M.ge.N)                      *
! *                  MROW .. Number of active rows of A                        *
! *                  NCOL .. Number of columns of A (M.ge.N)                   *
! *                  W ..... Output N element weight vector                    *
! *                  V ..... Output N by N matrix (NOT the transpose of "V")   *
! *                  RV1 ... Workspace N element vector                        *
! *                                                                            *
! * Given Ax=B, form U, W and V (not transpose) such that A=(U)(W)(Vtranspose) *
! * U, W and V are MxN matrix, weight vector(N) and NxN matrix, respectively.  *
! * Note that (U)(Ut) = (V)(Vt) = I.                                           *
! *                                                                            *
! * Based on an implementation from Numerical Recipes, 2nd Edition,            *
! * Press, et al, Cambridge University Press                                   *
! *                                                                            *
! ******************************************************************************
    implicit none

      integer(KIND=int_kind),intent(IN) :: mrowa,mrow,ncol
      real(KIND=real_kind),intent(INOUT) :: a(mrowa,ncol)
      real(KIND=real_kind),intent(OUT):: w(ncol),v(ncol,ncol)

      real(KIND=real_kind) f,g,h,scale,anorm,s,c,x,y,z
      integer(KIND=int_kind),parameter :: NMAX=4  ! Maximum number of columns (ndim+1)
      integer(KIND=int_kind),parameter :: MAX_ITERATIONS=30  ! Maximum number of iterations
      real(KIND=real_kind),dimension(NMAX),save :: rv1
      integer(KIND=int_kind) :: i,j,k,l,its,nm

      if(mrow < ncol)then
!       call PUNT((/'Under determined system.'/),'DO_SVDCMP')
      endif

      g = zero
      scale = zero
      anorm = zero

!  Reduce input A to bidiagonal form via Househoulder
      BIDIAG_FORM: do i=1,ncol
        l = i + 1
        rv1(i) = scale*g
        g = zero
        s = zero
        scale = zero
        if(i <= mrow)then
          do k=i,mrow
            scale = scale + abs(a(k,i))
          end do
          if(scale /= zero)then
            do k=i,mrow
              a(k,i) = a(k,i)/scale
              s = s + a(k,i)*a(k,i)
            end do
            f = a(i,i)
            g = -sign(sqrt(s),f)
            h = f*g - s
            a(i,i) = f - g
            if(i /= ncol)then
              do j=l,ncol
                s = zero
                do k=i,mrow
                  s = s + a(k,i)*a(k,j)
                end do
                f = s/h
                do k=i,mrow
                  a(k,j) = a(k,j)+f*a(k,i)
                end do
              end do  ! j=1,ncol
            endif  ! if(i /=
            do k=i,mrow
              a(k,i) = scale*a(k,i)
            end do
          endif  ! if(scale /=
        endif  ! if(i <=
        w(i) = scale*g
        g = zero
        s = zero
        scale = zero
        if (i <= mrow .and. i /= ncol) then
          do k=l,ncol
            scale = scale + abs(a(i,k))
          end do
          if(scale /= zero)then
            do k=l,ncol
              a(i,k) = a(i,k)/scale
              s = s + a(i,k)*a(i,k)
            end do
            f = a(i,l)
            g = -sign(sqrt(s),f)
            h = f*g - s
            a(i,l) = f - g
            do k=l,ncol
              rv1(k) = a(i,k)/h
            end do
            if(i /= mrow)then
              do j=l,mrow
                s = zero
                do k=l,ncol
                  s = s + a(j,k)*a(i,k)
                end do
                do k=l,ncol
                  a(j,k) = a(j,k) + s*rv1(k)
                end do
              end do  ! j=1
            endif  ! if(i /=
            do k=l,ncol
              a(i,k) = scale*a(i,k)
            end do
          endif  ! if(scale /=
        endif  ! if(i <=
        anorm = max( anorm,( abs(w(i)) + abs(rv1(i)) ) )
      end do BIDIAG_FORM

! Form right transformation
      CALC_RT: do i=ncol,1,-1
        if(i < ncol)then
          if(g /= zero)then
            do j=l,ncol
              v(j,i) = (a(i,j)/a(i,l))/g
            end do
            do j=l,ncol
              s= zero
              do k=l,ncol
                s = s + a(i,k)*v(k,j)
              end do
              do k=l,ncol
                v(k,j) = v(k,j) + s*v(k,i)
              end do
            end do  ! j=
          endif  ! if(g /=
          do j=l,ncol
            v(i,j) = zero
            v(j,i) = zero
          end do
        endif  ! if(i <
        v(i,i) = one
        g = rv1(i)
        l = i
      end do CALC_RT

! Form left transformation
      CALC_LT: do i=min(mrow,ncol),1,-1
        l = i + 1
        g = w(i)
        if(i < ncol)then
          do j=l,ncol
            a(i,j) = zero
          end do
        endif
        if(g /= zero)then
          g = one/g
          if(i /= ncol)then
            do j=l,ncol
              s = zero
              do k=l,mrow
                s = s + a(k,i)*a(k,j)
              end do
              f = (s/a(i,i))*g
              do k=i,mrow
                a(k,j) = a(k,j) + f*a(k,i)
              end do
            end do  ! j=1
          endif  ! if(i /=
          do j=i,mrow
            a(j,i) = a(j,i)*g
          end do
        else
          do j= i,mrow
            a(j,i) = zero
          end do
        endif  ! if(g /=
        a(i,i) = a(i,i) + one
      end do CALC_LT

! Diagonlize bidiagonal form
      FINAL_CALC: do k=ncol,1,-1
        IT_LOOP: do its=1,MAX_ITERATIONS

          L_LOOP: do l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm) == anorm)exit L_LOOP  ! Note rv1(1) is always zero
            if((abs(w(nm))+anorm) == anorm)then  ! Cancel out rv1(l) if l>1
              c = zero  
              s = one
              do i=l,k
                f = s*rv1(i)
                if((abs(f)+anorm) == anorm)exit L_LOOP
                g = w(i)
!               h = sqrt(f*f+g*g)
                h = pythag(f,g)
                w(i)  = h
                h = one/h
                c = (g*h)
                s = -(f*h)
                do j=1,mrow
                  y = a(j,nm)
                  z = a(j,i)
                  a(j,nm) = y*c + z*s
                  a(j, i) = z*c - y*s
                end do
              end do  ! i=l
              exit L_LOOP
            end if  ! if(w(nm
          end do L_LOOP
          z = w(k)

          if(l == k)then  ! Solution converged
            if(z < zero)then
              w(k) = -z
              do j=1,ncol
                v(j,k) = -v(j,k)
              end do
            endif
            exit IT_LOOP
          endif

          if (its == MAX_ITERATIONS) &
            call TLS_panic ('DO_SVDCMP: maximum number of iterations exceeded')

!  Shift bottom 2x2 minor
          x = w(l)
          nm = k - 1
          y = w(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ((y-z)*(y+z)+(g-h)*(g+h))/(two*h*y)
!         g = sqrt(f*f+one)
          g = pythag(f,one)
          f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
!  Next QR transformation
          c = one
          s = one
          NM_LOOP: do j=l,nm
            i = j + 1
            g = rv1(i)
            y = w(i)
            h = s*g
            g = c*g
!           z = sqrt(f*f+h*h)
            z = pythag(f,h)
            rv1(j) = z
            c = f/z
            s = h/z
            f = x*c + g*s
            g = g*c - x*s
            h = y*s
            y = y*c
            do nm=1,ncol
              x = v(nm,j)
              z = v(nm,i)
              v(nm,j) = x*c + z*s
              v(nm,i) = z*c - x*s
            end do
!           z = sqrt(f*f+h*h)
            z = pythag(f,h)
            w(j) = z
            if (z /= zero)then  ! Rotation fixed if z /= zero
              z = one/z
              c = f*z
              s = h*z
            endif
            f = c*g + s*y 
            x = c*y - s*g
            do nm=1,mrow
              y = a(nm,j)
              z = a(nm,i)
              a(nm,j) = y*c + z*s 
              a(nm,i) = z*c - y*s
            end do
          end do NM_LOOP
          rv1(l) = zero
          rv1(k) = f
          w(k) = x
        end do IT_LOOP
      end do FINAL_CALC
      END SUBROUTINE DO_SV_DECOMP

      FUNCTION PYTHAG(a,b)
      real(KIND=real_kind), intent(IN) :: a,b
      real(KIND=real_kind)             :: pythag
      real(KIND=real_kind) :: absa, absb
      absa = abs(a)
      absb = abs(b)
      if(absa > absb)then
        pythag=absa*sqrt(one + (absb/absa)**2)
      else
        if(absb == zero)then
          pythag = zero
        else
          pythag=absb*sqrt(one + (absa/absb)**2)
        endif
      endif
      END FUNCTION PYTHAG


  SUBROUTINE DO_LU_DECOMP(A,n,row,col)
    !=======================================================================
    ! Purpose:
    !
    !   Compute LU decomposition of A, with row and column pivoting.
    !   Row and column permutation matricies are stored in arrays row and
    !   row2.  For example, the solution to :
    !
    !      A x = b
    !
    !   is also given by the solution to:
    !
    !      R A C C^-1 x = R b
    !
    !   Where R is a matrix which permutes two rows of A, and C is
    !   a matrix which permutes columns of A.
    !
    !
    !   The action of R on a vector,       (Rx)(i) = x(row(i))
    !   The action of C^-1 on a vector, (C^-1x)(i) = x(col(i))
    !
    !   Using this routine in conjuction with LSLR_SOLVE_pivot(), we
    !   first solve:
    !
    !   R A C y = R b
    !
    !   and then solve  y = C^-1 x, so that x(col(i)) = y(i)
    !
    !   Notes:
    !   1. A must be rank 3 or 4.  Singular rows are solved by
    !   setting te solution on that row to 0
    !   2. A will be overwritten with its LU decomposition.
    !   3. Singular is defined by the normalized relation:
    !   abs(pivot_value(row)) <= alittle
    !
    !
    !=======================================================================
    implicit none
    Integer(KIND=int_kind), intent(IN) :: n

    ! arguments
    integer (int_kind),               intent(OUT)  :: row(n)
    integer (int_kind),               intent(OUT)  :: col(n)
    real (real_kind), dimension(n,n), intent(INOUT)  :: A

    ! local variables
    Integer(KIND=int_kind) :: piv,i,j,rowp,colp
    real (real_kind) :: tmp
    real (real_kind) :: dont_care
  ! dont_care_magnitude sets the relative magnitude below which values
  ! are "zero" wrt a cell face's LU solve component. It use implies that
  ! components with magnitudes less than 'dont_care' are "unphysical" wrt
  ! the solution being sought and can be summarily eliminated. 
  ! The absolute value
  ! for a cell face is the product of the largest magnitude component of
  ! the (ndim+1)x(ndim+1) matrix with the "dont_care_magnitude"
    real (real_kind),parameter :: dont_care_magnitude=1.d-15
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do i=1,n
       row(i)=i
       col(i)=i
    enddo
    dont_care = maxval(abs(A))*dont_care_magnitude
    do piv=1,n
       ! find largest pivot, at location before piviting (rowp,colp)
       tmp=-1
       do i=piv,n
          do j=piv,n
             if ( tmp < abs(A(i,j)) ) then
                tmp = abs(A(i,j))
                rowp=i           
                colp=j
             endif    
          enddo   
       enddo   
       ! check for zero pivots.  
       if (tmp <= dont_care)A(rowp,colp)=2d200  ! component will be 0'd out in solve.
       if (piv .ne. rowp) then
          ! swap rows:  piv & row
          do j=1,n
             tmp=A(piv,j)
             A(piv,j)=A(rowp,j)
             A(rowp,j)=tmp
          enddo
          j = row(piv)
          row(piv)=row(rowp)
          row(rowp)=j
       endif
       if (piv .ne. colp) then
          ! swap columns piv and colp
          do i=1,n
             tmp=A(i,piv)
             A(i,piv)=A(i,colp)
             A(i,colp)=tmp
          enddo
          j = col(piv)
          col(piv)=col(colp)
          col(colp)=j
       endif
       ! row reduce the rest of the matrix:
       do i=piv+1,n
          tmp = A(i,piv)/A(piv,piv)
          do j=piv+1,n
             A(i,j) = A(i,j) - tmp*A(piv,j)
          enddo
          ! store information to act on B:
          A(i,piv)=tmp
       enddo
    enddo
    return
  END SUBROUTINE DO_LU_DECOMP


  SUBROUTINE DO_LU_SOLVE(A,B,n,nm,row,col)
    !=======================================================================
    ! Purpose:
    !
    !   Solve a n x n system of equations, where A is the output
    !   from a previous call to LSLR_SOLVE_LU and B is the RHS vector,
    !   Return the result in ans.  (data in B is trashed)
    !
    !=======================================================================
    implicit none
    Integer(KIND=int_kind), intent(IN)  :: n,nm

    ! arguments
    integer (int_kind),               intent(IN)  :: row(n)
    integer (int_kind),               intent(IN)  :: col(n)
    real (real_kind), dimension(n,n), intent(IN)  :: A
    real (real_kind), dimension(n),   intent(INOUT)  :: B

    ! function return
    real (real_kind)                              :: tmp
    real (real_kind), dimension(n)                :: B2
    ! local variables
    Integer(KIND=int_kind) :: piv,i,j
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! solve the system Ax=B
    ! requres LU decompostion in A, as well as row() and col()
    ! (row and column permutation arrays)
    ! apply row permuations to B
    do i=1,n 
       B2(i)=B(row(i))
    enddo    
    ! apply L to B
    do piv=1,n-1
       do i=piv+1,n
          B2(i) = B2(i) - A(i,piv)*B2(piv)
       enddo 
    enddo    
    ! apply U inverse to B.
    do i=n,1,-1
       if (A(i,i)>1d200) then  ! flag indicating 0 pivot
          B2(i)=0
       else  
          tmp = B2(i)
          do j=i+1,n 
             tmp = tmp  - A(i,j)*B2(j)
          enddo
          B2(i)=tmp/A(i,i)
       endif
    enddo
    ! apply column permutation
    do i=1,n
       B(col(i))=B2(i)
    enddo
    return
  END SUBROUTINE DO_LU_SOLVE


END MODULE DO_SOLVE_MODULE
