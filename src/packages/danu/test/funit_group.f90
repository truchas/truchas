!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Sample Fortran
!
module funit_group_test

    use iso_c_binding
    use fruit

    use funit_utils
    use danu_module

    implicit none

    character(64) :: test_file = "dummy_group.h5"

    contains
    ! setup_before_all
    ! setup = setup_before_each
    subroutine setup
        call create_test_h5_file
    end subroutine setup

    subroutine setup_before_each(fid)

        type(C_PTR), intent(out) :: fid

        call danu_file_create(test_file,fid)
    
    end subroutine setup_before_each  

    subroutine teardown_after_each(hid)
         type(C_PTR), intent(inout) :: hid
         
         call danu_file_close(hid)
         call delete_test_file

    end subroutine teardown_after_each  

    ! teardown_before_all
    ! teardown = teardown_before_each

    subroutine teardown
        call delete_test_file
    end subroutine teardown

    subroutine group_fail_test

        type(C_PTR) :: hid
        integer     :: err

        type(C_PTR) :: gid
        character(len=32) :: grp_name
        logical :: exists


        ! Create test file
        call setup_before_each(hid)

        ! Fail test
        grp_name = 'Group A'
        exists = group_exists(hid,grp_name)
        call assertEquals(.false.,exists, 'Failed to stat group that DNE')

        call group_open(hid,grp_name,gid,err)
        call assertEquals(DANU_FAILURE,err,'Failed to open group DNE')

        call group_close(gid,err)
        call assertEquals(DANU_FAILURE,err,'Failed to close DNE group')

        call teardown_after_each(hid)

    end subroutine group_fail_test

    subroutine group_pass_test

        type(C_PTR) :: hid
        integer     :: err

        character(len=16) :: grp_name 
        type(C_PTR)       :: gid
        logical           :: flag


        ! Open the test file
        call setup_before_each(hid)

        ! Create the group
        grp_name = 'Test Group'
        call group_create(hid,grp_name,gid,err)
        call assertEquals(DANU_SUCCESS,err, 'Faile to create a group')
        call group_close(gid,err)
        call assertEquals(DANU_SUCCESS,err, 'Faile to close a group')

        ! Open the group
        call group_open(hid,grp_name,gid,err)
        call assertEquals(DANU_SUCCESS,err, 'Faile to open a group')
        call group_close(gid,err)
        call assertEquals(DANU_SUCCESS,err, 'Faile to close a group')

        ! Test existence
        flag = group_exists(hid,grp_name)
        call assertEquals(.true.,flag, 'Failed to flag group that exists')

        call teardown_after_each(hid)


    end subroutine group_pass_test


end module funit_group_test


! --- Main test driver program
program group_test_driver

    use fruit
    use funit_group_test

    logical :: all_pass

! --- Code

! --- Initialize FRUIT
      call init_fruit

! --- Run tests
      call group_fail_test
      !call group_pass_test

! --- Report results      
      call fruit_summary

! --- Exit if any test failed
      call is_all_successful(all_pass)

      if ( .not. all_pass ) then
          print *,'FAILED tests'
          call fail_exit_now
      end if    
      


end program group_test_driver    
