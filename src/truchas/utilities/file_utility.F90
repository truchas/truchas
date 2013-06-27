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
    use output_data_module, only: prefix

    ! argument list
    character(*), intent(IN) :: string
    integer, optional, intent(IN) :: number
    character(*), optional, intent(IN) :: path
    character(*), optional, intent(IN) :: suffix

    ! function return
    character(1024) :: MAKE_FILE_NAME

    ! local variables
    character(5)    :: number_string
    character(1024) :: root, buff
    integer :: start, end

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
             buff = TRIM(buff) // '/'
          end if
          buff = TRIM(buff) // TRIM(root)
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
       buff = TRIM(buff) // '.' // TRIM(string)
    end if

    ! conditionally append number
    if (present(number)) then
       write (number_string,'(i5.5)') number
       buff = TRIM(buff) // '.' // number_string
    end if

    ! conditionally append suffix
    if ( present(suffix) ) then
       buff = TRIM(buff) // '.' // TRIM(suffix)
    end if

    ! return file name
    MAKE_FILE_NAME = TRIM(ADJUSTL(buff))

  END FUNCTION MAKE_FILE_NAME


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
