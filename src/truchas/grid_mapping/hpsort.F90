!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module hpsort
  !=======================================================================
  ! Purpose(s):
  !
  !    Provide various heapsort sorting routines 
  !
  !    Public Interface(s):
  !
  !      hpsortip
  !      hpsortim
  !      hpsortimp
  !      hpsort1
  ! Author(s): Andrew Kuprat (kuprat@lanl.gov) 
  !      (Note: these are all inspired by the routine 'hpsort' in
  !       Numerical Recipes by Press et al.)
  !=======================================================================

  implicit none 

  PRIVATE

  PUBLIC :: hpsortip
  PUBLIC :: hpsortim
  PUBLIC :: hpsortimp
  PUBLIC :: hpsort1

  integer, parameter :: dp = KIND(1.0d0)

contains

  SUBROUTINE hpsortip(n,ia,ascend,iprm)
    !
    !
    !#######################################################################
    !
    !     PURPOSE -
    !
    !     HPSORTIP ("HeaP SORT an Integer array, producing a
    !     Permutaion") takes an integer N, an integer array IA,
    !     an integer array IPRM of length N, and a real
    !     number ASCEND, and 
    !     reorders IPRM so that IA(IPRM(1)),...,IA(IPRM(N))
    !     is in ascending order if ASCEND is positive and
    !     is in decreasing order if ASCEND is negative.
    !
    !     INPUT ARGUMENTS -
    !
    !        N - number of elements to be sorted.
    !        IA - integer array of values which determine
    !             how IPRM will be reordered.
    !        IPRM - integer array to be reordered.  Strictly
    !                        speaking, IPRM need not be a permutation of
    !                        the integers {1,..,N}, but may be simply
    !                        a mapping from {1,..,N} onto a set of N
    !                        distinct positive integers.
    !        ASCEND - real which controls whether we sort in ascending
    !                 or descending order.
    !
    !     OUTPUT ARGUMENTS -
    !
    !        IPRM - reordered integer array.
    !
    !
    !#######################################################################
    !

    implicit none

    INTEGER n,ia(*)
    REAL(dp) ascend
    INTEGER i,ir,j,l,iprm(n),irra
    if (n.lt.2) return
    l=n/2+1
    ir=n
10  continue
    if(l.gt.1)then
       l=l-1
       irra=iprm(l)
    else
       irra=iprm(ir)
       iprm(ir)=iprm(1)
       ir=ir-1
       if(ir.eq.1)then
          iprm(1)=irra
          return
       endif
    endif
    i=l
    j=l+l
20  if(j.le.ir)then
       if(j.lt.ir)then
          if(ascend*ia(iprm(j)).lt.ascend*ia(iprm(j+1)))j=j+1
       endif
       if(ascend*ia(irra).lt.ascend*ia(iprm(j)))then
          iprm(i)=iprm(j)
          i=j
          j=j+j
       else
          j=ir+1
       endif
       goto 20
    endif
    iprm(i)=irra
    goto 10
  END subroutine hpsortip


  SUBROUTINE hpsortim(n,m,md,ia)
    !
    !
    !#######################################################################
    !
    !     PURPOSE -
    !
    !     HPSORTIM ("HeaP SORT using an Integer M-fold key")
    !     takes an integer N, and an integer array IA, with
    !     first dimension length MD and second dimension length N and
    !     shuffles the columns IA(1:MD,*) so that they are in 
    !     lexicographic order up to the M'th key.  That is, if I<J
    !     and IA(1,I) > IA(1,J), we interchange the I'th and the J'th
    !     columns.  If IA(1,I) = IA(1,J), we then check if 
    !     IA(2,I) > IA(2,J), in which case we again interchange the
    !     columns.  Continuing on in this fashion, we consult the
    !     elements IA(K,I) and IA(K,J) for K up to M if necessary
    !     to break ties.  (M<=MD.) 
    !
    !     INPUT ARGUMENTS -
    !
    !        N - no. of columns to sort into ascending order.
    !        M - maximum number of keys (rows) to consult for comparisons
    !        MD- column length of array IA (M<=MD)
    !        IA - integer array of MD-tuples to be reordered
    !             into lexicographic ascending order up to the M'th key.
    !
    !     OUTPUT ARGUMENTS -
    !
    !        IA - SORTED REAL*8 ARRAY.
    !
    !
    !
    !#######################################################################
    !

    implicit none

    INTEGER n,m,md
    INTEGER ia(md,*)
    INTEGER i,ir,j,l,k,k1
    INTEGER itemp(md)
    if (n.lt.2) return
    l=n/2+1
    ir=n
10  continue
    if(l.gt.1)then
       l=l-1
       do k=1,md
          itemp(k)=ia(k,l)
       enddo
    else
       do k=1,md
          itemp(k)=ia(k,ir)
          ia(k,ir)=ia(k,1)
       enddo
       ir=ir-1
       if(ir.eq.1)then
          do k=1,md
             ia(k,1)=itemp(k)
          enddo
          return
       endif
    endif
    i=l
    j=l+l
20  if(j.le.ir)then
       if(j.lt.ir)then
          !             if(ra(j).lt.ra(j+1))j=j+1
          do k=1,m-1
             if(ia(k,j).lt.ia(k,j+1)) then
                j=j+1
                goto 30
             elseif(ia(k,j).gt.ia(k,j+1)) then
                goto 30
             endif
          enddo
          if(ia(m,j).lt.ia(m,j+1))j=j+1
       endif
30     continue
       do k=1,m-1
          if(itemp(k).lt.ia(k,j)) then
             do k1=1,md
                ia(k1,i)=ia(k1,j)
             enddo
             i=j
             j=j+j
             goto 20
          elseif(itemp(k).gt.ia(k,j)) then
             j=ir+1
             goto 20
          endif
       enddo
       if(itemp(m).lt.ia(m,j)) then
          do k=1,md
             ia(k,i)=ia(k,j)
          enddo
          i=j
          j=j+j
       else
          j=ir+1
       endif
       goto 20
    endif
    do k=1,md
       ia(k,i)=itemp(k)
    enddo
    goto 10
  END subroutine hpsortim


  SUBROUTINE hpsortimp(n,m,md,ia,ascend,iprm)
    !
    !
    !#######################################################################
    !
    !     PURPOSE -
    !
    !     HPSORTIMP ("HeaP SORT using an Integer M-fold key, generating
    !     a Permution array") takes an integer N, an integer array IA, with
    !     first dimension length MD and second dimension length N, an integer
    !     permutation array IPRM of length N, and a real number ASCEND, and 
    !     reorders IPRM so that IA(1,IPRM(1)),...,IA(1,IPRM(N))
    !     is in ascending order if ASCEND is positive and
    !     is in decreasing order if ASCEND is negative.  To 
    !     break ties, we require also that IA(2,IPRM(*)) is
    !     in ascending (descending) order, and so on, until
    !     possibly the M'th key IA(M,*) is consulted
    !     (i.e. ascending lexicographic order of M-tuples).
    !     Of course this requires that M<=MD.
    !
    !     INPUT ARGUMENTS -
    !
    !        N - number of elements to be sorted.
    !        M - we interpret array IA as M-tuples
    !        MD- actual first dimension of array IA (M<=MD)
    !        IA - integer array of values which determine
    !             how IPRM will be reordered, it has 'depth' M for 
    !             purposes of tie-breaking
    !        IPRM - integer array to be reordered.  Strictly
    !             speaking, IPRM need not be a permutation of
    !             the integers {1,..,N}, but may be simply
    !             a mapping from {1,..,N} onto a set of N
    !             distinct positive integers.
    !        ASCEND - real*8 which controls whether we sort in ascending
    !                 or descending order.
    !
    !     OUTPUT ARGUMENTS -
    !
    !        IPRM - reordered integer array.
    !
    !#######################################################################
    !

    implicit none

    INTEGER n,m,md
    INTEGER ia(md,*)
    REAL(dp) ascend
    INTEGER i,ir,j,l,iprm(n),irra,k
    if (n.lt.2) return
    l=n/2+1
    ir=n
10  continue
    if(l.gt.1)then
       l=l-1
       irra=iprm(l)
    else
       irra=iprm(ir)
       iprm(ir)=iprm(1)
       ir=ir-1
       if(ir.eq.1)then
          iprm(1)=irra
          return
       endif
    endif
    i=l
    j=l+l
20  if(j.le.ir)then
       if(j.lt.ir)then
          do k=1,m-1
             if(ascend*ia(k,iprm(j)).lt.ascend*ia(k,iprm(j+1))) then
                j=j+1
                goto 30
             elseif(ascend*ia(k,iprm(j)).gt.ascend*ia(k,iprm(j+1)))then
                goto 30
             endif
          enddo
          if(ascend*ia(m,iprm(j)).lt.ascend*ia(m,iprm(j+1)))j=j+1
       endif
30     continue
       do k=1,m-1
          if(ascend*ia(k,irra).lt.ascend*ia(k,iprm(j))) then
             iprm(i)=iprm(j)
             i=j
             j=j+j
             goto 20
          elseif(ascend*ia(k,irra).gt.ascend*ia(k,iprm(j))) then
             j=ir+1
             goto 20
          endif
       enddo
       if(ascend*ia(m,irra).lt.ascend*ia(m,iprm(j))) then
          iprm(i)=iprm(j)
          i=j
          j=j+j
       else
          j=ir+1
       endif
       goto 20
    endif
    iprm(i)=irra
    goto 10
  END subroutine hpsortimp

  SUBROUTINE hpsort1(n,ra,ascend,iprm)
!
!
!#######################################################################
!
!     PURPOSE -
!
!     HPSORT1 takes an integer N, a real array RA,
!     an integer array IPRM of length N, and a real
!     number ASCEND, and 
!     reorders IPRM so that RA(IPRM(1)),...,RA(IPRM(N))
!     is in ascending order if ASCEND is positive and
!     is in decreasing order if ASCEND is negative.
!
!     INPUT ARGUMENTS -
!
!        N - number of elements to be sorted.
!        RA - real array of values which determine
!             how IPRM will be reordered.
!        IPRM - integer array to be reordered.  Strictly
!                        speaking, IPRM need not be a permutation of
!                        the integers {1,..,N}, but may be simply
!                        a mapping from {1,..,N} onto a set of N
!                        distinct positive integers.
!        ASCEND - real which controls whether we sort in ascending
!                 or descending order.
!
!     OUTPUT ARGUMENTS -
!
!        IPRM - reordered integer array.
!
!#######################################################################
!
 
      implicit none
 
      INTEGER n
      REAL(dp) ra(*),ascend
      INTEGER i,ir,j,l,iprm(n),irra
      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          irra=iprm(l)
        else
          irra=iprm(ir)
          iprm(ir)=iprm(1)
          ir=ir-1
          if(ir.eq.1)then
            iprm(1)=irra
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ascend*ra(iprm(j)).lt.ascend*ra(iprm(j+1)))j=j+1
          endif
          if(ascend*ra(irra).lt.ascend*ra(iprm(j)))then
            iprm(i)=iprm(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        iprm(i)=irra
      goto 10
    END subroutine hpsort1

end module hpsort
