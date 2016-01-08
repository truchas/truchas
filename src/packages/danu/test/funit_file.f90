!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Sample Fortran
!
module funit_file_test

    use iso_c_binding
    use fruit

    use funit_utils
    use danu_module

    implicit none

    character(64) :: test_file = "dummy_file.h5"

    contains
    ! setup_before_all
    ! setup = setup_before_each
    subroutine setup
        call create_test_h5_file
    end subroutine setup

    ! teardown_before_all
    ! teardown = teardown_before_each

    subroutine teardown
        call delete_test_file
    end subroutine teardown

    subroutine file_create_test

        type(C_PTR) :: hid_ptr
        integer     :: err

        call danu_file_create(test_file,hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err)
        call danu_file_close(hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err)

    end subroutine file_create_test

    subroutine file_open_test

        type(C_PTR) :: hid_ptr
        integer     :: err

        ! Read only open
        call danu_file_open_rdonly(test_file,hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err)
        call danu_file_close(hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err)

        ! Read Write open
        call danu_file_open_rdwr(test_file,hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err)
        call danu_file_close(hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err)




    end subroutine file_open_test


end module funit_file_test


! --- Main test driver program
program file_test_driver

    use fruit
    use funit_file_test

    logical :: all_pass

! --- Code

! --- Initialize FRUIT
      call init_fruit

! --- Run tests
      call file_create_test
      call file_open_test

! --- Report results      
      call fruit_summary

! --- Exit if any test failed
      call is_all_successful(all_pass)

      if ( .not. all_pass ) then
          print *,'FAILED tests'
          call fail_exit_now
      end if    
      


end program file_test_driver    
