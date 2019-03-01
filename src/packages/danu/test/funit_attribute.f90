!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Sample Fortran
!
module funit_attribute_test

    use iso_c_binding
    use fruit

    use funit_utils

    use danu_module

    implicit none

    character(64) :: test_file = "dummy_file.h5"
    integer       :: iseed = 1020304
    integer       :: test_file_num_attr = 0

    contains
    ! setup_before_all
    ! setup = setup_before_each
    subroutine setup
        type(C_PTR) :: hid
        call danu_file_create(test_file,hid)
        call danu_file_close(hid)
    end subroutine setup

    subroutine update_test_file_count
        test_file_num_attr = test_file_num_attr + 1
    end subroutine update_test_file_count


    ! teardown_before_all
    ! teardown = teardown_before_each

    subroutine teardown
        call delete_test_file
    end subroutine teardown

    subroutine attribute_exists_test
    
        type(C_PTR)   :: hid_ptr
        integer       :: err
        logical       :: exists

        character(32) :: attr_exists     = 'Attribute that exists'
        character(32) :: attr_not_exists = 'Attribute that does not exist'

        integer :: attr_value

        call danu_file_open_rdwr(test_file,hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to open test file')

        if ( err .eq. DANU_SUCCESS ) then

            exists = attribute_exists(hid_ptr,attr_not_exists)
            call assertEquals(exists,.false.,'False exist test failed')

            attr_value = generate_random_int(iseed)
            call attribute_write(hid_ptr,attr_exists,attr_value,err)
            call update_test_file_count

            if ( err .eq. DANU_SUCCESS ) then
                exists = attribute_exists(hid_ptr,attr_exists)
                call assertEquals(exists,.true., 'True exist test failed')
            else
                call assertEquals(DANU_SUCCESS,err, 'Failed to write attribute')    
            endif

            call danu_file_close(hid_ptr)

        endif


    end subroutine attribute_exists_test

    subroutine attribute_count_test

        type(C_PTR)   :: hid_ptr
        integer       :: err

        integer       :: num_found

        call danu_file_open_rdonly(test_file,hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to open test file')

        if ( err .eq. DANU_SUCCESS ) then
            call attribute_count(hid_ptr,num_found)
            call assertEquals(test_file_num_attr,num_found, 'Attribute count failed')

            call danu_file_close(hid_ptr,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to close file')

        endif

    end subroutine attribute_count_test

    subroutine attribute_list_test

        type(C_PTR) :: hid_ptr
        integer     :: err

        integer     :: num_found, i

        character(len=64), dimension(:), allocatable :: names

        call danu_file_open_rdonly(test_file,hid_ptr,err)
        if ( err .eq. DANU_SUCCESS ) then
            call attribute_count(hid_ptr,num_found, err)

            allocate(names(num_found))

            names = ' '

            call attribute_list(hid_ptr,names,err)
            write(*,*) 'Number of attributes found', num_found
            i = 1
            do while( i .le. num_found )
              write(*,*) 'Attribute=', trim(names(i))
              i = i + 1
            end do


        endif

        deallocate(names)

    end subroutine attribute_list_test
    
    subroutine attribute_int_test
       
        type(C_PTR)   :: hid_ptr
        integer       :: err

        character(64) :: attr_name_int1 = 'Dummy Integer Attribute 1'
        integer       :: rand_int

        character(64) :: attr_name_int2 = 'Dummy Integer Attribute 2'
        integer*4     :: int4

        integer :: int_save, int_read



        call danu_file_open_rdwr(test_file,hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to open test file')

        if ( err .eq. DANU_SUCCESS ) then

            rand_int = generate_random_int(iseed)
            call attribute_write(hid_ptr,attr_name_int1,rand_int,err)
            call assertEquals(DANU_SUCCESS,err,'Failed to write int attribute')
            call update_test_file_count
           
            int4 = generate_random_int(iseed) 
            call attribute_write(hid_ptr,attr_name_int2,int4,err)
            call assertEquals(DANU_SUCCESS,err,'Failed to write int attribute')
            call update_test_file_count

            call danu_file_close(hid_ptr,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to close test file')

            call danu_file_open_rdonly(test_file,hid_ptr,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to open (rdonly) test file')

            call attribute_read(hid_ptr,attr_name_int1,int_read,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to read int attribute test file')
            call assertEquals(rand_int, int_read,  &
                              'Failed to read coorect int attribute')

            call attribute_read(hid_ptr,attr_name_int2,int_read,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to read int attribute test file')
            call assertEquals(int4, int_read,  &
                              'Failed to read coorect int attribute')

            call danu_file_close(hid_ptr,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to close test file')

        endif

    end subroutine attribute_int_test

    subroutine attribute_real4_test
      
        type(C_PTR)   :: hid_ptr 
        integer       :: err

        character(64) :: attr_name_real1 = 'Dummy REAL4 Attribute 1'
        real          :: rand_real4

        character(64) :: attr_name_real2 = 'Dummy REAL4 Attribute 2'
        real*4        :: real4 = 12345678

        real :: real4_save, real4_read



        call danu_file_open_rdwr(test_file,hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to open test file')

        if ( err .eq. DANU_SUCCESS ) then

            rand_real4 = generate_random_float(iseed)
            call attribute_write(hid_ptr,attr_name_real1,rand_real4,err)
            call assertEquals(DANU_SUCCESS,err,'Failed to write real4 attribute')
            call update_test_file_count
            
            call attribute_write(hid_ptr,attr_name_real2,real4,err)
            call assertEquals(DANU_SUCCESS,err,'Failed to write real4 attribute')
            call update_test_file_count

            call danu_file_close(hid_ptr,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to close test file')

            call danu_file_open_rdonly(test_file,hid_ptr,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to open (rdonly) test file')

            call attribute_read(hid_ptr,attr_name_real1,real4_read,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to read real4 attribute test file')
            call assertEquals(rand_real4, real4_read,  &
                              'Failed to read coorect real4 attribute')

            call attribute_read(hid_ptr,attr_name_real2,real4_read,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to read real4 attribute test file')
            call assertEquals(real4, real4_read,  &
                              'Failed to read coorect real4 attribute')

            call danu_file_close(hid_ptr,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to close test file')

        endif

    end subroutine attribute_real4_test

    subroutine attribute_real8_test
       
        type(C_PTR)   :: hid_ptr
        integer       :: err

        character(64)  :: attr_name_real1 = 'Dummy REAL8 Attribute 1'
        real(C_DOUBLE) :: rand_real8

        character(64) :: attr_name_real2 = 'Dummy REAL8 Attribute 2'
        real*8        :: real8 = 12345678

        real*8 :: real8_save, real8_read



        call danu_file_open_rdwr(test_file,hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to open test file')

        if ( err .eq. DANU_SUCCESS ) then

            rand_real8 = generate_random_double(iseed)
            call attribute_write(hid_ptr,attr_name_real1,rand_real8,err)
            call assertEquals(DANU_SUCCESS,err,'Failed to write real8 attribute')
            call update_test_file_count
            
            call attribute_write(hid_ptr,attr_name_real2,real8,err)
            call assertEquals(DANU_SUCCESS,err,'Failed to write real8 attribute')
            call update_test_file_count

            call danu_file_close(hid_ptr,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to close test file')

            call danu_file_open_rdonly(test_file,hid_ptr,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to open (rdonly) test file')

            call attribute_read(hid_ptr,attr_name_real1,real8_read,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to read real8 attribute test file')
            call assertEquals(rand_real8, real8_read,  &
                              'Failed to read coorect real8 attribute')

            call attribute_read(hid_ptr,attr_name_real2,real8_read,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to read real8 attribute test file')
            call assertEquals(real8, real8_read,  &
                              'Failed to read coorect real8 attribute')

            call danu_file_close(hid_ptr,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to close test file')

        endif

    end subroutine attribute_real8_test

    subroutine attribute_character_test
       
        type(C_PTR)   :: hid_ptr
        integer       :: err

        character(64)    :: attr_name1 = 'Dummy Character Attribute 1'
        character(len=8) :: rand_char1 = '12345678'

        character(64) :: attr_name2 = 'Dummy Character Attribute 2'
        character(len=64)   :: rand_char2 =  'The lazy dog jumpeed'

        character(len=32) :: char_save, char_read



        call danu_file_open_rdwr(test_file,hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to open test file')

        if ( err .eq. DANU_SUCCESS ) then

            call attribute_write(hid_ptr,attr_name1,rand_char1,err)
            call assertEquals(DANU_SUCCESS,err,'Failed to write character attribute')
            call update_test_file_count
            
            call attribute_write(hid_ptr,attr_name2,rand_char2,err)
            call assertEquals(DANU_SUCCESS,err,'Failed to write character attribute')
            call update_test_file_count

            call danu_file_close(hid_ptr,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to close test file')

            call danu_file_open_rdonly(test_file,hid_ptr,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to open (rdonly) test file')

            call attribute_read(hid_ptr,attr_name1,char_read,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to read character attribute test file')
            call assertEquals(rand_char1, char_read,  &
                              'Failed to read correct character attribute')

            call attribute_read(hid_ptr,attr_name2,char_read,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to read character attribute test file')
            call assertEquals(rand_char2, char_read,  &
                              'Failed to read correct character attribute')

            call danu_file_close(hid_ptr,err)
            call assertEquals(DANU_SUCCESS,err, 'Failed to close test file')

        endif

    end subroutine attribute_character_test

    
end module funit_attribute_test


! --- Main test driver program
program attribute_test_driver

    use fruit
    use funit_attribute_test

    logical :: all_pass

! --- Code

! --- Initialize FRUIT
      call init_fruit

! --- Run setup routines
      call setup

! --- Run tests
      call attribute_int_test
      call attribute_real4_test
      call attribute_real8_test
      call attribute_character_test
      call attribute_exists_test
      call attribute_count_test
      call attribute_list_test

! --- Report results      
      call fruit_summary

! --- Tear down
      !call teardown

! --- Exit if any test failed
      call is_all_successful(all_pass)

      if ( .not. all_pass ) then
          print *,'FAILED tests'
          call fail_exit_now
      end if    
      


end program attribute_test_driver    
