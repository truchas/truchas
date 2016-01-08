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
module funit_nsdata_test

    use iso_c_binding
    use fruit

    use funit_utils
    use danu_module

    implicit none

    character(64) :: test_file = "dummy_file.h5"
    character(64) :: test_sim_group = "Test Simulation"

    contains
    ! setup_before_all
    ! setup = setup_before_each
    subroutine setup_output_file(fid,sid)

        type(C_PTR), intent(out) :: fid, sid

        call output_file_create(test_file,fid)
        call simulation_add(fid,test_sim_group,sid)

    end subroutine setup_output_file

    subroutine add_datasets(sid,num,data_names)

        type(C_PTR), intent(in) :: sid
        integer, intent(in)     :: num
        character(len=64),dimension(:) :: data_names

        integer :: i
        integer(C_INT),dimension(2) :: idata
       
        
        idata(1) = 1
        idata(2) = 2
        i = 1
        do while ( i .le. num )
          write(data_names(i), *) 'Dataset', i
          call data_write(sid,data_names(i),idata)
          i = i + 1
        end do  

    end subroutine add_datasets
    
    
    subroutine nsdata_attribute_test
        
        type(C_PTR)   :: fid, sid, hid
        integer       :: err
        character(len=64),dimension(1:1) :: names
        character(len=16) :: attr_name
        character(len=16) :: attr


        call set_unit_name('Non-series Data Attribute')

        call setup_output_file(fid,sid)
        call add_datasets(sid,1,names)
        call output_file_close(fid)


        call output_file_open_rdwr(test_file,fid)
        call simulation_open(fid,test_sim_group,sid)
        call data_open_dataset(sid,names(1),hid,err)
        call assertEquals(err, DANU_SUCCESS,             &
                          'Failed to open non series dataset')
        attr_name='Dummy Char attribute'                 
        attr='CELL DATA'
        call attribute_write(hid,attr_name,attr)



        call output_file_close(fid)

    
    end subroutine nsdata_attribute_test

    subroutine nsdata_exists_test

        type(C_PTR)   :: fid, sid
        integer       :: err
        character(len=64),dimension(1) :: data_name
        logical       :: exists

        call set_unit_name('Non-series Exists')

        call setup_output_file(fid,sid)

        data_name(1) = 'Data Does Not Exist'
        exists = data_exists(sid,data_name(1))
        call assertEquals(.false., exists,'Failed to return correct exists (F)')

        call add_datasets(sid,1,data_name)

        exists =  data_exists(sid,data_name(1))
        call assertEquals(.true., exists,       &
                         'Failed to return correct status (T)')

        call group_close(sid)
        call output_file_close(fid)

    end subroutine nsdata_exists_test

    subroutine nsdata_count_test

        type(C_PTR) :: fid,sid
        integer     :: err
        integer     :: cnt
        integer     :: num
        integer     :: iseed = 3
        character(len=64) , dimension(:), allocatable :: names


        call set_unit_name('Data Count')

        call setup_output_file(fid,sid)
        num = generate_random_bound_int(10,512,iseed)
        allocate(names(num))

        call data_count(sid,cnt,err)
        call assertEquals(DANU_SUCCESS,err,'Failed to count simulations')
        call assertEquals(0,cnt, 'Failed to return 0 nsdata count')
        
        call add_datasets(sid,num,names)
        call data_count(sid,cnt,err)
        call assertEquals(DANU_SUCCESS,err,'Failed to count simulations')
        call assertEquals(num,cnt, 'Failed to return correct nsdata count')

        deallocate(names)
        call group_close(sid)
        call output_file_close(fid)

    end subroutine nsdata_count_test

    subroutine nsdata_list_test

        type(C_PTR) :: fid, sid
        integer :: cnt, num
        integer :: iseed = 12345
        character(len=64), dimension(:), allocatable :: names, rnames
        integer(C_INT) :: err
        
        ! Set up 
        call set_unit_name('Non-series Data List')
        num = generate_random_bound_int(10,128,iseed)
        allocate(names(num))
        
        call setup_output_file(fid,sid)

        ! Query empty file for a list
        call data_list(sid,names,err)
        call assertEquals(DANU_FAILURE,err,              &
                          'Failed to return fail status')

        call add_datasets(sid,num,names)

        ! Send an small array 
        allocate(rnames(num-2))
        call data_list(sid,rnames,err)
        call assertEquals(DANU_FAILURE,err,              &
                          'Failed to return fail status')
        deallocate(rnames)

        ! Send the right size
        allocate(rnames(num))
        call data_list(sid,rnames,err)
        call assertEquals(DANU_SUCCESS,err,              &
                          'Failed to read dataset names')
        
        cnt = 1
        do while ( cnt .le. num )
          call assertEquals(names(cnt),rnames(cnt),          &
                            'Failed to read names correctly')
          cnt = cnt + 1
        end do  


        deallocate(names)                  
        deallocate(rnames)                  
        call group_close(sid) 
        call output_file_close(fid) 
       
    end subroutine nsdata_list_test

    subroutine nsdata_int

        type(C_PTR) :: fid, sid
        integer(C_INT) :: err
        character(len=64) :: data_name 
               
        integer :: rank, i
        integer :: length
        integer :: idata0,iread0, iseed
        integer(C_INT),dimension(:), allocatable :: idata1
        integer(C_INT),dimension(:), allocatable :: iread1
        integer(C_INT),dimension(:,:), allocatable :: idata2
        integer(C_INT),dimension(:,:), allocatable :: iread2
        integer(C_INT),dimension(:,:,:), allocatable :: idata3
        integer(C_INT),dimension(:,:,:), allocatable :: iread3
        integer(C_INT)                    :: flag, num

        ! Set up
        call set_unit_name('Non-series Write INTEGER')
        call setup_output_file(fid,sid)
        iseed = 13579



        ! Tests that should pass

        ! Write/Read single data
        idata0 = generate_random_int(iseed)
        idata0 = 10
        data_name = 'Dummy Integer Data'
        call data_write(sid,data_name,idata0,err)
        call assertEquals(DANU_SUCCESS, err,                       &
                          'Failed to write single integer data')
        call data_read(sid,data_name,iread0,err)
        call assertEquals(DANU_SUCCESS, err,                        &
                         'Failed to read single integer data')
        call assertEquals(idata0,iread0,                             &
                        'Failed to read correct (0) integer data')


        ! Write/Read Rank 1 Array
        rank = 1
        length = 512 
        allocate(idata1(length))
        allocate(iread1(length))
        i = 1
        num = 0
        do while ( i .le. length ) 
           idata1(i) = i
           i =  i + 1
        end do
        num = length

        data_name = 'Dummy Integer Fortran Data Rank 1'
        call data_write(sid,data_name,idata1,err)
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to write 1D integer data')
        call data_read(sid,data_name,iread1,err)                 
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to read 1D integer data')
        flag = int_array_cmp(idata1, iread1,num)
        call assertEquals(0,flag,                                &
                        'Failed to read the correct rank 1 data')

        ! Write/Read Rank 2 Array
        rank = 2
        length = 10
        allocate(idata2(rank,length))
        allocate(iread2(rank,length))
        i = 1
        num = 0
        do while ( i .le. rank ) 
           idata2(i,:) = i
           i =  i + 1
        end do
        num = rank*length

        data_name = 'Dummy Integer Fortran Data Rank 2'
        call data_write(sid,data_name,idata2,err)
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to write 2D integer data')
        call data_read(sid,data_name,iread2,err)                 
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to read 2D integer data')
        flag = int_array_cmp(idata2, iread2,num)
        call assertEquals(0,flag,                                &
                        'Failed to read the correct rank 2 data')
        
        
        ! Write/Read Rank 3 Array
        rank = 3
        length = 8
        allocate(idata3(length,4,rank))
        allocate(iread3(length,4,rank))
        i = 1
        num = 0
        do while ( i .le. length ) 
           idata3(i,:,:) = i
           i =  i + 1
        end do 
        num = length*4*rank

        data_name = 'Dummy Integer Fortran Data Rank 3'
        call data_write(sid,data_name,idata3,err)
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to write 3D integer data')
        call data_read(sid,data_name,iread3,err)                 
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to read 3D integer data')
        flag = int_array_cmp(idata3,iread3,num)
        call assertEquals(0,flag,                                &
                        'Failed to read the correct rank 3 data')
        
       
        ! Clean up

        deallocate(idata2)
        deallocate(idata3)
        deallocate(iread2)
        deallocate(iread3)

        call group_close(sid)
        call output_file_close(fid)

    end subroutine nsdata_int

    subroutine nsdata_real4

        type(C_PTR) :: fid, sid
        integer(C_INT) :: err
        character(len=64) :: data_name 
               
        integer :: rank, i
        integer :: length
        real(C_FLOAT) :: r4data0, r4read0
        real(C_FLOAT),dimension(:), allocatable :: r4data1
        real(C_FLOAT),dimension(:), allocatable :: r4read1
        real(C_FLOAT),dimension(:,:), allocatable :: r4data2
        real(C_FLOAT),dimension(:,:), allocatable :: r4read2
        real(C_FLOAT),dimension(:,:,:), allocatable :: r4data3
        real(C_FLOAT),dimension(:,:,:), allocatable :: r4read3
        integer(C_INT)                    :: flag, num
        integer(C_INT) :: iseed

        call set_unit_name('Non-series Write REAL4')
        call setup_output_file(fid,sid)
        iseed = 123456

        ! Tests that should pass

        ! Write/Read single data
        r4data0 = generate_random_float(iseed)
        data_name = 'Dummy REAL4 Data'
        call data_write(sid,data_name,r4data0,err)
        call assertEquals(DANU_SUCCESS, err,                       &
                          'Failed to write single REAL4 data')
        call data_read(sid,data_name,r4read0,err)
        call assertEquals(DANU_SUCCESS, err,                        &
                         'Failed to read single REAL4 data')
        call assertEquals(r4data0,r4read0,                             &
                        'Failed to read correct (0) REAL4 data')


        ! Write/Read Rank 1 Array
        rank = 1
        length = 512 
        allocate(r4data1(length))
        allocate(r4read1(length))
        i = 1
        num = 0
        do while ( i .le. length ) 
           r4data1(i) = generate_random_float(iseed) 
           i =  i + 1
        end do
        num = length

        data_name = 'Dummy REAL4 Fortran Data Rank 1'
        call data_write(sid,data_name,r4data1,err)
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to write 1D integer data')
        call data_read(sid,data_name,r4read1,err)                 
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to read 1D integer data')
        flag = float_array_cmp(r4data1, r4read1,num)
        call assertEquals(0,flag,                                &
                        'Failed to read the correct rank 1 data')

        ! Write/Read Rank 2 Array
        rank = 2
        length = 10
        allocate(r4data2(rank,length))
        allocate(r4read2(rank,length))
        i = 1
        num = 0
        do while ( i .le. rank ) 
           r4data2(i,:) = generate_random_float(iseed)
           i =  i + 1
        end do
        num = rank*length

        data_name = 'Dummy REAL4 Fortran Data Rank 2'
        call data_write(sid,data_name,r4data2,err)
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to write 2D integer data')
        call data_read(sid,data_name,r4read2,err)                 
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to read 2D integer data')
        flag = float_array_cmp(r4data2, r4read2,num)
        call assertEquals(0,flag,                                &
                        'Failed to read the correct rank 2 data')
        
        
        ! Write/Read Rank 3 Array
        rank = 3
        length = 8
        allocate(r4data3(length,4,rank))
        allocate(r4read3(length,4,rank))
        i = 1
        num = 0
        do while ( i .le. length ) 
           r4data3(i,:,:) = generate_random_float(iseed)
           i =  i + 1
        end do 
        num = length*4*rank

        data_name = 'Dummy REAL4 Fortran Data Rank 3'
        call data_write(sid,data_name,r4data3,err)
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to write 3D integer data')
        call data_read(sid,data_name,r4read3,err)                 
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to read 3D integer data')
        flag = float_array_cmp(r4data3,r4read3,num)
        call assertEquals(0,flag,                                &
                        'Failed to read the correct rank 3 data')
        
       
        ! Clean up

        deallocate(r4data1)
        deallocate(r4data2)
        deallocate(r4data3)
        deallocate(r4read1)
        deallocate(r4read2)
        deallocate(r4read3)

        call group_close(sid)
        call output_file_close(fid)


        

    end subroutine nsdata_real4
    
    subroutine nsdata_real8

        type(C_PTR) :: fid, sid
        integer(C_INT) :: err
        character(len=64) :: data_name 
               
        integer :: rank, i
        integer :: length
        real(C_DOUBLE) :: r8data0, r8read0
        real(C_DOUBLE),dimension(:), allocatable :: r8data1
        real(C_DOUBLE),dimension(:), allocatable :: r8read1
        real(C_DOUBLE),dimension(:,:), allocatable :: r8data2
        real(C_DOUBLE),dimension(:,:), allocatable :: r8read2
        real(C_DOUBLE),dimension(:,:,:), allocatable :: r8data3
        real(C_DOUBLE),dimension(:,:,:), allocatable :: r8read3
        integer(C_INT)                    :: flag, num
        integer(C_INT) :: iseed

        call set_unit_name('Non-series Write INTEGER')
        call setup_output_file(fid,sid)
        iseed = 123456

        ! Tests that should pass

        ! Write/Read single data
        r8data0 = generate_random_double(iseed)
        data_name = 'Dummy REAL8 Data'
        call data_write(sid,data_name,r8data0,err)
        call assertEquals(DANU_SUCCESS, err,                       &
                          'Failed to write single REAL8 data')
        call data_read(sid,data_name,r8read0,err)
        call assertEquals(DANU_SUCCESS, err,                        &
                         'Failed to read single REAL8 data')
        call assertEquals(r8data0,r8read0,                             &
                        'Failed to read correct (0) REAL4 data')


        ! Write/Read Rank 1 Array
        rank = 1
        length = 512 
        allocate(r8data1(length))
        allocate(r8read1(length))
        i = 1
        num = 0
        do while ( i .le. length ) 
           r8data1(i) = generate_random_double(iseed) 
           i =  i + 1
        end do
        num = length

        data_name = 'Dummy REAL8 Fortran Data Rank 1'
        call data_write(sid,data_name,r8data1,err)
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to write 1D integer data')
        call data_read(sid,data_name,r8read1,err)                 
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to read 1D integer data')
        flag = double_array_cmp(r8data1, r8read1,num)
        call assertEquals(0,flag,                                &
                        'Failed to read the correct rank 1 data')

        ! Write/Read Rank 2 Array
        rank = 2
        length = 10
        allocate(r8data2(rank,length))
        allocate(r8read2(rank,length))
        i = 1
        num = 0
        do while ( i .le. rank ) 
           r8data2(i,:) = generate_random_double(iseed)
           i =  i + 1
        end do
        num = rank*length

        data_name = 'Dummy REAL8 Fortran Data Rank 2'
        call data_write(sid,data_name,r8data2,err)
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to write 2D integer data')
        call data_read(sid,data_name,r8read2,err)                 
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to read 2D integer data')
        flag = double_array_cmp(r8data2, r8read2,num)
        call assertEquals(0,flag,                                &
                        'Failed to read the correct rank 2 data')
        
        
        ! Write/Read Rank 3 Array
        rank = 3
        length = 8
        allocate(r8data3(length,4,rank))
        allocate(r8read3(length,4,rank))
        i = 1
        num = 0
        do while ( i .le. length ) 
           r8data3(i,:,:) = generate_random_double(iseed)
           i =  i + 1
        end do 
        num = length*4*rank

        data_name = 'Dummy REAL8 Fortran Data Rank 3'
        call data_write(sid,data_name,r8data3,err)
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to write 3D integer data')
        call data_read(sid,data_name,r8read3,err)                 
        call assertEquals(DANU_SUCCESS,err,                     &
                          'Failed to read 3D integer data')
        flag = double_array_cmp(r8data3,r8read3,num)
        call assertEquals(0,flag,                                &
                        'Failed to read the correct rank 3 data')
        
       
        ! Clean up

        deallocate(r8data1)
        deallocate(r8data2)
        deallocate(r8data3)
        deallocate(r8read1)
        deallocate(r8read2)
        deallocate(r8read3)

        call group_close(sid)
        call output_file_close(fid)

    end subroutine nsdata_real8
 
end module funit_nsdata_test


! --- Main test driver program
program nsdata_test_driver

    use fruit
    use funit_nsdata_test

! --- Local variables
      
      logical :: all_pass

! --- Code

! --- Initialize FRUIT
      call init_fruit

! --- Run tests
      call nsdata_attribute_test
      call nsdata_exists_test
      call nsdata_count_test
      call nsdata_list_test
      call nsdata_int
      call nsdata_real4
      call nsdata_real8

! --- Report results      
      call fruit_summary

! --- Exit if any test failed
      all_pass = .false.
      call is_all_successful(all_pass)

      if ( .not. all_pass ) then
          print *,'FAILED tests'
          call fail_exit_now
      end if    
      


end program nsdata_test_driver    
