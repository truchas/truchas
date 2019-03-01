!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Fortran Unit Test: hid_t
!
module funit_hidt_test

    use, intrinsic :: iso_c_binding
    use fruit

    use funit_utils

    implicit none

    contains
    
    subroutine setup_before_all
        call create_test_h5_file
    end subroutine setup_before_all

    subroutine teardown
        call delete_test_file
    end subroutine teardown

    subroutine hid_test

        integer(C_INT64_T) :: hid
        type(C_PTR)        :: hid_ptr

        hid = open_test_file()
        hid_ptr = create_hid_struct(hid)
        call destroy_hid_struct(hid_ptr)
        call close_test_file(hid)
        
    end subroutine hid_test   


end module funit_hidt_test


! --- Main test driver program
program hidt_test_driver

    use fruit
    use funit_hidt_test

    logical :: all_pass

! --- Code

! --- Initialize FRUIT
      call init_fruit

! --- Set up before all
      call setup_before_all

! --- Run tests
      call hid_test

! --- Report results      
      call fruit_summary

! --- Tear down
      call teardown

! --- Exit if any test failed
      call is_all_successful(all_pass)

      if ( .not. all_pass ) then
          print *,'FAILED tests'
          call fail_exit_now
      end if    
      


end program hidt_test_driver    
