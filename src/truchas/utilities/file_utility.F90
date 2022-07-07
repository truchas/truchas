!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module file_utility

  implicit none
  private
  
  public :: COUNT_TOKENS
  public :: GET_TOKEN

contains

  FUNCTION COUNT_TOKENS (string, separator)
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !    count tokens in a string
    !    
    !    Fixed bug in token count --Sriram
    !---------------------------------------------------------------------------

    ! function result
    Integer :: COUNT_TOKENS

    ! Arguments
    character(*), Intent(IN) :: string
    character(*), Intent(IN) :: separator

    ! Local Variables
    Integer :: i, j

    !---------------------------------------------------------------------------

    j = LEN_TRIM(separator)-1
    
    if ( LEN_TRIM(string) > 0 ) then
       COUNT_TOKENS = 1
    else 
       COUNT_TOKENS = 0
    end if


    do i = 1, LEN_TRIM(string)
       if (string(i:i+j) == TRIM(separator)) then
          COUNT_TOKENS = COUNT_TOKENS + 1
       end if
    end do

  END FUNCTION COUNT_TOKENS


  !-----------------------------------------------------------------------------

  SUBROUTINE GET_TOKEN (token, n, string, separator)
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !    return the n'th token in a string
    !    
    !    Generalized to tokens of arbitrary length --Sriram
    !---------------------------------------------------------------------------

    ! arguments
    integer, intent(IN)  :: n
    character(*), intent(IN)  :: string
    character(*), intent(IN)  :: separator
    character(*), Intent(OUT) :: token

    ! local variables
    integer :: i
    integer :: s ! start
    integer :: e ! end
    integer :: is
    integer :: ip
    integer :: itoken

    !---------------------------------------------------------------------------
    is = len_trim(string)
    ip = len(separator)

    itoken = 1
    s = 1
    e = 1
    token=''
    SEARCH: do i = 1, is
       if ( string(i:i+ip-1) == separator) then
          if ( itoken == n ) then
             exit SEARCH
          else 
             s = i + ip
             itoken = itoken + 1
             if ( s > is ) exit SEARCH
          end if
       end if
    end do SEARCH
    if ( itoken < n-1 ) then
       token=''
       e = 0
    else 
       if ( i <= is) then
          e=i
          token = string(s:e-1)
       else 
          token = string(s:is)
          e=is+1
       end if
    end if

  END SUBROUTINE GET_TOKEN

end module file_utility
