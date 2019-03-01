!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
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

        call danu_file_create(test_file,hid_ptr)
        call danu_file_close(hid_ptr)

    end subroutine file_create_test

!LORI    subroutine output_create_test
!LORI
!LORI        type(C_PTR) :: hid_ptr
!LORI        integer :: err
!LORI
!LORI        ! Create a basic output file
!LORI        call output_file_create(testfile,hid_ptr,stat=err)
!LORI        call assertEquals(DANU_SUCCESS,err)
!LORI
!LORI        call teardown
!LORI
!LORI   end subroutine output_create_test   
!LORI
!LORI   subroutine output_open_rdonly_test
!LORI       
!LORI        type(C_PTR) :: hid_ptr
!LORI        integer :: err
!LORI
!LORI        ! Open output file read only ... this should fail
!LORI        call output_file_open_rdonly(testfile,hid_ptr,stat=err)
!LORI        call assertEquals(DANU_FAILURE,err)
!LORI
!LORI        ! Now create the file
!LORI        call setup
!LORI
!LORI        ! Open read only
!LORI        call output_file_open_rdonly(testfile,hid_ptr,stat=err)
!LORI        call assertEquals(DANU_SUCCESS,err)
!LORI        call output_file_close(hid_ptr)
!LORI
!LORI        call teardown
!LORI
!LORI        
!LORI   end subroutine output_open_rdonly_test

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
      !call output_create_test
      !call output_open_rdonly_test

! --- Report results      
      !call fruit_summary

! --- Exit if any test failed
      call is_all_successful(all_pass)

      if ( .not. all_pass ) then
          print *,'FAILED tests'
          call exit_now(DANU_FAILURE)
      end if    
      


end program file_test_driver    
