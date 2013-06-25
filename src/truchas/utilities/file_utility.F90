module file_utility
  implicit none
  private
   public :: MAKE_FILE_NAME
   public :: COUNT_TOKENS
   public :: GET_TOKEN

 contains
  FUNCTION MAKE_FILE_NAME (string, number, path, suffix)
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !    prepare a consistent file name from a character string (string), an
    !    optional integer number, and an optional character string "path"
    !
    !    filename will be in the format prefix.id
    !    if number is present, filename will be in the format prefix.id.nnnnn
    !    where nnnnn is a 5 digit integer derived from number
    !    if path is present, use path instead of prefix for the directory component
    !
    !    Added Suffix so that you can send .gmv or whatever --Sriram
    !---------------------------------------------------------------------------
    use kind_module,   only: int_kind
    use output_data_module, only: prefix
    implicit none

    ! argument list
    character (LEN=*),            intent(IN) :: string
    integer (int_kind), optional, intent(IN) :: number
    character (LEN=*),  optional, intent(IN) :: path
    character (LEN=*),  optional, intent(IN) :: suffix

    ! function return
    character (LEN=1024) :: MAKE_FILE_NAME

    ! local variables
    character (LEN=5)    :: number_string
    character (LEN=1024) :: root
    character (LEN=1024) :: buff
    integer              :: start
    integer              :: end

    !---------------------------------------------------------------------------

    ! if the optional path is present
    if (present(path)) then
       ! path specified, find root filename and prefix with path
       start = SCAN(prefix,'/',.true.) + 1
       end   = LEN_TRIM(prefix)
       root  = prefix(start:end)
       buff  = ADJUSTL(path)
       ! if buff is a non-null string, append root name to it
       end = LEN_TRIM(buff)
       if (end > 0) then
          ! make sure buff is terminated with a /
          if (buff(end:end) /= '/') then
#if ( defined(DARWIN_NAG_COMPILER_WORKAROUND) )
             buff = buff(1:len_trim(buff)) // '/'
#else
             buff = TRIM(buff) // '/'
#endif
          end if
#if ( defined(DARWIN_NAG_COMPILER_WORKAROUND) )
          buff = buff(1:len_trim(buff)) // TRIM(root)
#else
          buff = TRIM(buff) // TRIM(root)
#endif
       else
          ! don't append, just copy
          buff = TRIM(root)
       end if

    else
       ! no path specified, use prefix
       buff = TRIM(ADJUSTL(prefix))

    end if

    ! conditionally append string
    if (LEN_TRIM(string) /= 0) then
#if ( defined(DARWIN_NAG_COMPILER_WORKAROUND) )
       buff = buff(1:len_trim(buff)) //  '.' // TRIM(string)
#else
       buff = TRIM(buff) // '.' // TRIM(string)
#endif
    end if

    ! conditionally append number
    if (present(number)) then
       write (number_string,'(i5.5)') number
#if ( defined(DARWIN_NAG_COMPILER_WORKAROUND) )
       buff = buff(1:len_trim(buff)) // '.' // number_string
#else
       buff = TRIM(buff) // '.' // number_string
#endif
    end if

    ! conditionally append suffix
    if ( present(suffix) ) then
#if ( defined(DARWIN_NAG_COMPILER_WORKAROUND) )
       buff = buff(1:len_trim(buff)) // '.' // TRIM(suffix)
#else
       buff = TRIM(buff) // '.' // TRIM(suffix)
#endif
    end if

    ! return file name
    MAKE_FILE_NAME = TRIM(ADJUSTL(buff))

    return

  END FUNCTION MAKE_FILE_NAME


  FUNCTION COUNT_TOKENS (string, separator)
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !    count tokens in a string
    !    
    !    Fixed bug in token count --Sriram
    !---------------------------------------------------------------------------
    implicit none

    ! function result
    Integer :: COUNT_TOKENS

    ! Arguments
    character (LEN=*), Intent(IN) :: string
    character (LEN=*), Intent(IN) :: separator

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

    return

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
    implicit none

    ! arguments
    integer,           intent(IN)  :: n
    character (LEN=*), intent(IN)  :: string
    character (LEN=*), intent(IN)  :: separator
    character (LEN=*), Intent(OUT) :: token

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
    return

  END SUBROUTINE GET_TOKEN

!!$  FUNCTION COUNT_TOKENS (string, seperator)
!!$    !---------------------------------------------------------------------------
!!$    ! Purpose:
!!$    !
!!$    !    count tokens in a string
!!$    !---------------------------------------------------------------------------
!!$    implicit none
!!$
!!$    ! function result
!!$    Integer :: COUNT_TOKENS
!!$
!!$    ! Arguments
!!$    character (LEN=*), Intent(IN) :: string
!!$    character (LEN=*), Intent(IN) :: seperator
!!$
!!$    ! Local Variables
!!$    logical :: lcwt     ! true if Last character Was Token
!!$    Integer :: i
!!$
!!$    !---------------------------------------------------------------------------
!!$
!!$    COUNT_TOKENS = 0
!!$    lcwt = .true.
!!$
!!$    do i = 1, LEN_TRIM(string)
!!$       if (string(i:i) == seperator(1:1)) then
!!$          lcwt = .true.
!!$       else
!!$          if (lcwt) COUNT_TOKENS = COUNT_TOKENS + 1
!!$          lcwt = .false.
!!$       end if
!!$    end do
!!$
!!$    return
!!$
!!$  END FUNCTION COUNT_TOKENS
!!$
!!$  !-----------------------------------------------------------------------------
!!$
!!$  SUBROUTINE GET_TOKEN (token, n, string, seperator)
!!$    !---------------------------------------------------------------------------
!!$    ! Purpose:
!!$    !
!!$    !    return the n'th token in a string
!!$    !---------------------------------------------------------------------------
!!$    implicit none
!!$
!!$    ! arguments
!!$    integer,           intent(IN)  :: n
!!$    character (LEN=*), intent(IN)  :: string
!!$    character (LEN=*), intent(IN)  :: seperator
!!$    character (LEN=*), Intent(OUT) :: token
!!$
!!$    ! local variables
!!$    character (LEN=1024) :: tmpstring
!!$    integer :: i
!!$    integer :: s ! start
!!$    integer :: e ! end
!!$
!!$    !---------------------------------------------------------------------------
!!$
!!$    tmpstring = string
!!$    s = 1
!!$    e = SCAN(tmpstring(1:), seperator)
!!$    if (e == 0) then
!!$       e = LEN_TRIM(tmpstring(1:))
!!$    else
!!$       e = e - 1
!!$    end if
!!$
!!$    do i = 1, n-1
!!$       s = e + 2
!!$       e = SCAN(tmpstring(s:), seperator)
!!$       if (e == 0) then
!!$          e = LEN_TRIM(tmpstring(1:))
!!$       else
!!$          e = e - 1
!!$       end if
!!$    end do
!!$
!!$    token = tmpstring(s:e)
!!$
!!$    return
!!$
!!$  END SUBROUTINE GET_TOKEN

end module file_utility
