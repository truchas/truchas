MODULE TABULAR_UTILITIES

  !=======================================================================
  ! Purpose(s):
  !
  !  Utility functions for tabular data
  !
  ! Contains: TABLUAR_LINEAR_INTERP
  !
  ! Author(s):  Jim Sicilian (sicilian@lanl.gov)
  ! May 2003
  !
  !=======================================================================
  use kind_module,      only: real_kind, int_kind

  implicit none

  ! Private Module
  private

  ! Public Subroutines
  public :: TABULAR_LINEAR_INTERP


  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

    FUNCTION TABULAR_LINEAR_INTERP(x,xtable,ytable)
    !=======================================================================
    ! Purpose(s):
    !
    !   Perform a linear interpolation within a table.
    !
    !   Public Interface:
    !     * VALUE = TABULAR_LINEAR_INTERP(x,xtable,ytable)
    !                   x is the value at which the dependent value is desired
    !                   xtable is a list of independent variable values
    !                   ytable is a list of dependent variable values
    !
    !               The resulting VALUE is the first dependent variable value if x
    !               is less than the first independent variable; it is the last
    !               dependent variable entry if x is greater than the last 
    !               independent value; otherwise it is linearly interpolated with
    !               the table.
    !
    !               This routine uses a simple binary search through the table until
    !               it brackets x between two values of the independent variable
    !
    !======================================================================

    use constants_module, only: one
    use truchas_logging_services, only: TLS_panic

    implicit none

    real(kind=real_kind)                            :: TABULAR_LINEAR_INTERP


    ! Local Variables
    integer(kind=int_kind)                          :: i, itest, imin, imax, table_length, count
    real(kind=real_kind)                            :: alpha

    ! Argument List
    real(kind=real_kind), INTENT(IN)                :: x
    real(kind=real_kind), dimension(:),INTENT(IN) :: xtable, ytable

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    table_length = SIZE(xtable)

! The independent variable MUST be in ascending order for this to work, so check it
    alpha = xtable(1)
    do i = 2,table_length
       if(xtable(i) <= alpha) then
        call TLS_panic ('TABULAR_LINEAR_INTERP: non-ascending independent variable')
       endif
    end do

! First, check for trivial cases
    if(table_length == 1 .or. x <= xtable(1)) then
        TABULAR_LINEAR_INTERP = ytable(1)
        return
      elseif(x >= xtable(table_length)) then
        TABULAR_LINEAR_INTERP = ytable(table_length)
        return
    endif

! If we get here we actually have to search the table
    imin = 1
    imax = table_length
    count = 0
    do while ( (imax-imin)>1 .and. count < table_length)
        itest = (imin+imax)/2
        count = count +1

!       Special case of x identical to a table entry
        if(x == xtable(itest)) then
            TABULAR_LINEAR_INTERP = ytable(itest)
            return
        endif

!       on which side of xtable(itest) is x?
        if(x > xtable(itest)) then
            imin = itest
        else
            imax = itest
        endif
    end do

!  This is an unlikely event, but check anyway
    if(imax /= imin+1) then
        call TLS_panic ('TABULAR_LINEAR_INTERP: imin /= imax')
    endif

!  x is in the single interval between imin and imax
!     Calculate the averaging factor
    alpha = xtable(imax)-xtable(imin)
    alpha = (x-xtable(imin))/alpha
    TABULAR_LINEAR_INTERP = alpha*ytable(imax) + (one - alpha)*ytable(imin)
    return

    END FUNCTION TABULAR_LINEAR_INTERP

    END MODULE TABULAR_UTILITIES        
