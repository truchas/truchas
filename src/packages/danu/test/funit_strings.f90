!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Fortran Unit Test for string handling
!
module funit_strings_test


    use iso_c_binding
    use fruit
    use funit_utils

    implicit none

    integer, parameter :: STRING_LEN        = 64
    integer, parameter :: STRING_ARRAY_LEN  = 16 
    integer, parameter :: STRING_ARRAY_SIZE = 10

    character(kind=C_CHAR,len=STRING_LEN), save, private :: string
    character(kind=C_CHAR,len=STRING_ARRAY_LEN), dimension(STRING_ARRAY_SIZE), save, private :: string_array

    character(kind=C_CHAR,len=STRING_LEN), save, private :: c_string

    ! C functions to test
    interface
        integer(C_INT) function compute_fortran_trim_len(fort_string, length) &
                                bind(c) 
            use, intrinsic :: iso_c_binding                   
            character(C_CHAR),     intent(in) :: fort_string
            integer(C_INT),        value      :: length
        end function compute_fortran_trim_len
    end interface    
            
    interface 
        subroutine convert_string_f2c(f_str, f_len, c_str, c_len) bind(c)
        use, intrinsic :: iso_c_binding
        character(C_CHAR), intent(in) :: f_str
        integer(C_INT), value :: f_len
        character(C_CHAR), intent(out) :: c_str
        integer(C_INT),value  :: c_len
        end subroutine convert_string_f2c
    end interface    
            
    interface 
        subroutine convert_string_c2f(c_str, f_str, f_len) bind(c)
        use, intrinsic :: iso_c_binding
        character(C_CHAR), intent(in) :: c_str
        character(C_CHAR), intent(out) :: f_str
        integer(C_INT),value  :: f_len
        end subroutine convert_string_c2f
    end interface
  

    contains
    
    
    ! setup_before_all
    subroutine setup_before_all

        string = 'This is a test string #$%@'

        string_array(1)  = 'Hello'
        string_array(2)  = 'World'
        string_array(3)  = 'testme'
        string_array(4)  = '0123456789'
        string_array(5)  = 'Lori'
        string_array(6)  = '  me '
        string_array(7)  = 'A'
        string_array(8)  = 'Bb'
        string_array(9)  = 'ABBA'
        string_array(10) = 'laLLLLL'

        c_string = 'This is a C STRING' // C_NULL_CHAR
        
    end subroutine setup_before_all

    ! setup = setup_before_each
    subroutine setup
    end subroutine setup

    ! teardown_before_all
    ! teardown = teardown_before_each
    subroutine teardown
    end subroutine teardown

    subroutine trim_len_test

        integer        :: f_trim_len,idx
        integer(C_INT) :: c_len

        ! Create a basic output file
        idx = 1
        do while ( idx .le. STRING_ARRAY_SIZE )
            f_trim_len = len_trim(string_array(idx))
            print *, 'STRING=',string_array(idx), 'SIZE=', STRING_ARRAY_LEN
            c_len = compute_fortran_trim_len(string_array(idx), &
                                             STRING_ARRAY_LEN)
            
            if ( f_trim_len .ne. c_len ) then
               write(*,'(a,i3,a)') 'TRIM_LEN case ', idx, ' failed'
               write(*, '(a,a,a,i5,a,i5)') 'FORTRAN STRING=', string_array(idx), &
                                        'TRIM_LEN=', f_trim_len, '  C_TRIM_LEN=', &
                                         c_len               
            else                             
               write(*,'(a,i5,a,i5)') 'F_TRIM_LEN=',f_trim_len,' C_LEN=', c_len 
            endif
            call assertEquals(f_trim_len,c_len)
            idx = idx + 1
        end do    

        ! Single string
        f_trim_len = len_trim(string)
        c_len = compute_fortran_trim_len(string,STRING_LEN)
        call assertEquals(f_trim_len, c_len)

    end subroutine trim_len_test

    subroutine string_convert_test

        character(len=STRING_LEN) :: f_string
        character(len=STRING_LEN+1) :: c_string
        integer :: l, f_len, c_len

        f_len=STRING_LEN
        c_len=STRING_LEN+1
        
        f_string(1:f_len) = ' '
        c_string(1:c_len) = C_NULL_CHAR

        l = len_trim(string)

        call convert_string_f2c(string,STRING_LEN,c_string,c_len)
        call assertEquals(string(1:l),c_string(1:l), &
                          'f2c conversion failed')

        call convert_string_c2f(c_string,f_string,f_len)
        call assertEquals(f_string(1:l), c_string(1:l), &
                         'c2f conversion failed' )        

    end subroutine string_convert_test


end module funit_strings_test


! --- Main test driver program
program strings_test_driver

      use fruit
      use funit_strings_test

      logical :: all_pass

! --- Code

! --- Initialize FRUIT
      call init_fruit

! --- Setup before any test is called
      call setup_before_all

! --- Run tests
      call trim_len_test
      call string_convert_test

! --- Report results      
      call fruit_summary

! --- Exit if any test failed
      call is_all_successful(all_pass)

      if ( .not. all_pass ) then
          print *,'FAILED tests'
          call fail_exit_now
      end if    
      


end program strings_test_driver    
