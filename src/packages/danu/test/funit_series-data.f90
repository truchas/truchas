!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module funit_data_test

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

    subroutine create_series_data(nid,num,names)

        type(C_PTR), intent(in) :: nid
        integer, intent(in)     :: num 
        character(len=*),dimension(:),intent(out) :: names

        integer        :: i
        real(C_DOUBLE) :: dummy
        integer        :: iseed = 0

        i = 1
        do while ( i .le. num )
          write(names(i), '(a,i3)') 'Dummy Data', i
          dummy = generate_random_double(iseed)
          call simulation_data_write(nid,names(i),dummy)
          i = i + 1
        end do  

          
    end subroutine create_series_data

    

    subroutine data_type_test
        
        type(C_PTR)   :: fid, sid, nid
        integer       :: err
        integer       :: iseed
        character(len=64) :: dname
        logical       :: exists
        integer        :: cyc, code
        real(C_DOUBLE) :: time
        character(len=64),dimension(1:1) :: names

        call set_unit_name('Series Data Type')

        call setup_output_file(fid,sid)
        iseed = 345987
        time = generate_random_double(iseed)
        cyc = generate_random_bound_int(1,1000,iseed)
        call sequence_next_id(sid,cyc,time,nid)
        call create_series_data(nid,1,names)

        call simulation_data_type(nid,names(1),code,err)
        call assertEquals(err, DANU_SUCCESS,                     &
                          'Failed to read data type')

        call assertEquals(DANU_DATASET_DOUBLE, code,             &
                          'Failed to read correct type')        

        call group_close(nid)
        call group_close(sid)
        call output_file_close(fid)



    end subroutine data_type_test

    subroutine data_rank_test
        
        type(C_PTR)   :: fid, sid, nid
        integer       :: err, i, j
        integer       :: iseed
        character(len=64) :: dname
        logical       :: exists
        integer        :: cyc, rank
        integer, dimension(1:2) :: dims
        integer, dimension(:,:), allocatable :: idata 
        real(C_DOUBLE) :: time

        call set_unit_name('Series Data Rank')

        call setup_output_file(fid,sid)
        iseed = 345987
        time = generate_random_double(iseed)
        cyc = generate_random_bound_int(1,1000,iseed)
        call sequence_next_id(sid,cyc,time,nid)

        dname='Sample Dataset'
        call simulation_data_rank(nid,dname,rank,err)
        call assertEquals(err, DANU_FAILURE,                     &
                          'Failed to fail status for dataset DNE')

                          
        dims(1)=generate_random_bound_int(10,20,iseed)
        dims(2)=dims(1)+generate_random_bound_int(1,10,iseed)

        allocate(idata(dims(1),dims(2)) )
        idata=0
        i=1
        do while ( i .le. dims(1) )
          j=1
          do while ( j .le. dims(2) )
            idata(i,j) = generate_random_int(iseed)
            j=j+1
            end do
          i=i+1
        end do  

        call simulation_data_write(nid,dname,idata)

        call simulation_data_rank(nid,dname,rank,err)
        call assertEquals(err, DANU_SUCCESS,                      &
                         'Failed to find dataset rank')
        call assertEquals(rank,2,                                 &
                         'Failed to read correct rank')        


        deallocate(idata)
        call group_close(nid)
        call group_close(sid)
        call output_file_close(fid)



    end subroutine data_rank_test

    subroutine data_dimensions_test
        
        type(C_PTR)   :: fid, sid, nid
        integer       :: err, i, j
        integer       :: iseed
        character(len=64) :: dname
        logical       :: exists
        integer        :: cyc
        integer, dimension(1:2) :: dims, read_dims
        integer, dimension(:,:), allocatable :: idata 
        real(C_DOUBLE) :: time

        call set_unit_name('Series Data Dimensions')

        call setup_output_file(fid,sid)
        iseed = 345987
        time = generate_random_double(iseed)
        cyc = generate_random_bound_int(1,1000,iseed)
        call sequence_next_id(sid,cyc,time,nid)

        dname='Sample Dataset'
        call simulation_data_dimensions(nid,dname,read_dims,err)
        call assertEquals(err, DANU_FAILURE,                     &
                          'Failed to fail status for dataset DNE')

                          
        dims(1)=generate_random_bound_int(10,20,iseed)
        dims(2)=dims(1)+generate_random_bound_int(1,10,iseed)

        allocate(idata(dims(1),dims(2)) )
        idata=0
        i=1
        do while ( i .le. dims(1) )
          j=1
          do while ( j .le. dims(2) )
            idata(i,j) = generate_random_int(iseed)
            j=j+1
            end do
          i=i+1
        end do  

        call simulation_data_write(nid,dname,idata)

        call simulation_data_dimensions(nid,dname,read_dims,err)
        call assertEquals(err, DANU_SUCCESS,                      &
                         'Failed to find dataset dimensions')
        
        i=1
        do while (i .le. 2 )
          call assertEquals(dims(i),read_dims(i),                  &
                            'Failed to read correct dimensions')
          i=i+1
        end do          
                              

        deallocate(idata)
        call group_close(nid)
        call group_close(sid)
        call output_file_close(fid)



    end subroutine data_dimensions_test

    subroutine data_exists_test

        type(C_PTR)   :: fid, sid, nid
        integer       :: err
        integer       :: iseed
        character(len=64) :: dname
        logical       :: exists
        integer        :: cyc
        real(C_DOUBLE) :: time
        character(len=64),dimension(1:1) :: names

        call set_unit_name('Series Data Exists')

        call setup_output_file(fid,sid)
        iseed = 345987
        time = generate_random_double(iseed)
        cyc = generate_random_bound_int(1,1000,iseed)
        call sequence_next_id(sid,cyc,time,nid)

        dname = 'Data that does not exist'
        exists = simulation_data_exists(nid, dname)

        call assertEquals(.false., exists,'Failed to return correct exists (F)')

        call create_series_data(nid,1,names)
        exists = simulation_data_exists(nid,names(1))

        call assertEquals(.true., exists, 'Failed to return correct exists (T)')


        call group_close(nid)
        call group_close(sid)
        call output_file_close(fid)

    end subroutine data_exists_test

    subroutine data_attribute_test

        type(C_PTR)   :: fid, sid, nid, hid
        integer       :: err
        integer       :: iseed
        character(len=64) :: att_name
        integer           :: att_data
        integer        :: cyc
        real(C_DOUBLE) :: time
        character(len=64),dimension(1:1) :: names

        call set_unit_name('Series Data Attribute')

        call setup_output_file(fid,sid)
        iseed = 345987
        time = generate_random_double(iseed)
        cyc = generate_random_bound_int(1,1000,iseed)
        call sequence_next_id(sid,cyc,time,nid)

        call create_series_data(nid,1,names)

        call simulation_open_data(nid,names(1),hid,err)
        call assertEquals(DANU_SUCCESS,err,'Failed to open series data')

        att_name = 'Integer attribute'
        att_data = generate_random_bound_int(1,1000,iseed)
        call attribute_write(hid,att_name,att_data,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to write attribute')

        
        call group_close(nid)
        call group_close(sid)
        call output_file_close(fid)


    end subroutine data_attribute_test

    subroutine data_count_test

        type(C_PTR) :: fid,sid, nid
        integer     :: err
        integer     :: cnt
        integer     :: num
        integer     :: iseed = 12345 
        integer     :: cyc
        real(C_DOUBLE) :: time
        character(len=64), dimension(:), allocatable :: names


        call set_unit_name('Series Data Count')
        call setup_output_file(fid,sid)
        time = generate_random_double(iseed)
        cyc = generate_random_bound_int(1,1000,iseed)
        call sequence_next_id(sid,cyc,time,nid)

        call simulation_data_count(nid,cnt,err)
        call assertEquals(DANU_SUCCESS, err,                    &
                          'Failed to return correctly')
        call assertEquals(0,cnt,                                &
                          'Failed to return correct number 0')

        num = generate_random_bound_int(32,512,iseed)
        allocate(names(num))
        call create_series_data(nid,num,names)
        deallocate(names)

        call simulation_data_count(nid,cnt,err)
        call assertEquals(DANU_SUCCESS, err,                    &
                          'Failed to return correctly')
        call assertEquals(num,cnt,                                &
                          'Failed to return correct number')

        
        call group_close(nid)
        call group_close(sid)
        call output_file_close(fid)

    end subroutine data_count_test
    
    
    subroutine data_list_test
        
        type(C_PTR) :: fid,sid, nid
        integer          :: err
        integer          :: iseed = 678140
        integer          :: num, cnt, cyc
        real(C_DOUBLE)   :: time
        character(len=64), dimension(:), allocatable :: names
        character(len=64), dimension(:), allocatable :: rnames

        call set_unit_name('Series Data List Test')

        ! Set up the file, simulation and sequence
        call setup_output_file(fid,sid)
        time = generate_random_double(iseed)
        num = generate_random_bound_int(10,128,iseed)
        cyc = generate_random_bound_int(1,10000,iseed)
        call sequence_next_id(sid,cyc,time,nid)

        allocate(names(num))
        allocate(rnames(num))


        ! Call without names
        call simulation_data_list(nid,rnames,err)
        call assertEquals(DANU_FAILURE,err,                         &
                          'Failed to return correct status (F)')

                          
        ! Now create the data sets
        call create_series_data(nid,num,names)

        call simulation_data_list(nid,rnames,err)
        call assertEquals(DANU_SUCCESS,err,                         &
                          'Failed to read dataset names')

        cnt = 1
        do while ( cnt .le. num )
          !write(*,*) 'name,rname=', names(cnt), rnames(cnt)
          call assertEquals(names(cnt),rnames(cnt),                 &
                            'Wrong dataset name')
          cnt = cnt + 1
        end do   


        ! Now send array too small
        deallocate(rnames)
        allocate(rnames(num/2))
        call simulation_data_list(nid,rnames,err)
        call assertEquals(DANU_FAILURE,err,                         &
                          'Failed to return correct status (F)')


        deallocate(names)
        deallocate(rnames)
        call group_close(nid)                  
        call group_close(sid)
        call output_file_close(fid)

    end subroutine data_list_test

    subroutine data_integer_test

        type(C_PTR)      :: fid, sid, nid
        integer          :: iseed = 444
        real(C_DOUBLE)   :: time
        integer          :: num, cyc, ierr, i, j, k, d
        integer          :: idata0,rdata0
        integer(C_INT)   :: flag
        integer,dimension(:),allocatable :: dims
        integer,dimension(:),allocatable :: idata1, rdata1
        integer,dimension(:,:),allocatable :: idata2, rdata2
        integer,dimension(:,:,:),allocatable :: idata3, rdata3
        character(len=128) :: data_name

        call set_unit_name('Series Data Integer Test')

        ! Test Setup 
        call setup_output_file(fid,sid)
         time = generate_random_double(iseed)
         cyc = generate_random_bound_int(1,100000,iseed)
        call sequence_next_id(sid,cyc,time,nid)

        ! Write/Read rank 0 Integer data
        data_name = 'Integer Rank 0'
        idata0 = generate_random_int(iseed)
        call simulation_data_write(nid,data_name,idata0,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to write rank 0 integer')

        call simulation_data_read(nid,data_name,rdata0,ierr)
        call assertEquals(DANU_SUCCESS, ierr,                                  &
                          'Failed to read rank 0 integer')
        call assertEquals(idata0,rdata0,                                       &
                          'Failed to read correct data')

        
        ! Write/Read rank 1 Integer data
        num = generate_random_bound_int(128,512,iseed)
        allocate(idata1(num))
        allocate(rdata1(num))
        call generate_random_array(idata1,iseed)
        data_name = 'Integer Rank 1'
        call simulation_data_write(nid,data_name,idata1,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to write rank 1 integer')

        call simulation_data_read(nid,data_name,rdata1,ierr)

        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to read rank 1 integer')
        flag = int_array_cmp(idata1,rdata1,num)                  
        call assertEquals(0,flag,                                              &
                          'Failed to read correct data')
                          
        
        ! Write/Read rank 2 Integer data
        allocate(dims(2))
        call generate_random_bound_array(1,512,dims,iseed)
        allocate(idata2(dims(1),dims(2)))
        allocate(rdata2(dims(1),dims(2)))
        call generate_random_array(idata2,iseed)
        data_name = 'Integer Rank 2'
        call simulation_data_write(nid,data_name,idata2,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to write rank 2 integer')

        call simulation_data_read(nid,data_name,rdata2,ierr)

        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to read rank 2 integer')
        flag = int_array_cmp(idata2,rdata2,num)                  
        call assertEquals(0,flag,                                              &
                          'Failed to read correct data')


        ! Write/Read rank 3 Integer data
        deallocate(dims)
        allocate(dims(3))
        call generate_random_bound_array(1,100,dims,iseed)
        allocate(idata3(dims(1),dims(2),dims(3)))
        allocate(rdata3(dims(1),dims(2),dims(3)))
        call generate_random_array(idata3,iseed)

        data_name = 'Integer Rank 3'
        call simulation_data_write(nid,data_name,idata3,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to write rank 3 integer')

        call simulation_data_read(nid,data_name,rdata3,ierr)

        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to read rank 2 integer')
        flag = int_array_cmp(idata3,rdata3,num)                  
        call assertEquals(0,flag,                                              &
                          'Failed to read correct data')
        
        ! Test Cleanup
        deallocate(dims)
        deallocate(idata1)
        deallocate(rdata1)
        deallocate(idata2)
        deallocate(rdata2)
        deallocate(idata3)
        deallocate(rdata3)
        call group_close(sid)
        call output_file_close(fid)

    end subroutine data_integer_test
  
    subroutine data_real4_test

        type(C_PTR)      :: fid, sid, nid
        integer          :: iseed = 444
        real(C_DOUBLE)   :: time
        integer          :: cyc, num, ierr, i, j, k, d
        integer          :: data0,rdata0
        integer(C_INT)   :: flag
        integer,dimension(:),allocatable :: dims
        real,dimension(:),allocatable :: data1, rdata1
        real,dimension(:,:),allocatable :: data2, rdata2
        real,dimension(:,:,:),allocatable :: data3, rdata3
        character(len=128) :: data_name

        call set_unit_name('Series Data REAL4 Test')

        ! Test Setup 
        call setup_output_file(fid,sid)
         time = generate_random_double(iseed)
         cyc = generate_random_bound_int(1,1000,iseed)
        call sequence_next_id(sid,cyc,time,nid)

        ! Write/Read rank 0 REAL4 data
        data_name = 'REAL4 Rank 0'
        data0 = generate_random_float(iseed)
        call simulation_data_write(nid,data_name,data0,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to write rank 0 real4')

        call simulation_data_read(nid,data_name,rdata0,ierr)
        call assertEquals(DANU_SUCCESS, ierr,                                  &
                          'Failed to read rank 0 real4')
        call assertEquals(data0,rdata0,                                        &
                          'Failed to read correct data')

        
        ! Write/Read rank 1 REAL4 data
        num = generate_random_bound_int(128,512,iseed)
        allocate(data1(num))
        allocate(rdata1(num))
        call generate_random_array(data1,iseed)
        data_name = 'REAL4 Rank 1'
        call simulation_data_write(nid,data_name,data1,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to write rank 1 real4')

        call simulation_data_read(nid,data_name,rdata1,ierr)

        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to read rank 1 real4')
        flag = float_array_cmp(data1,rdata1,num)                  
        call assertEquals(0,flag,                                              &
                          'Failed to read correct data')
                          
        
        ! Write/Read rank 2 REAL4 data
        allocate(dims(2))
        call generate_random_bound_array(1,512,dims,iseed)
        allocate(data2(dims(1),dims(2)))
        allocate(rdata2(dims(1),dims(2)))
        call generate_random_array(data2,iseed)
        data_name = 'REAL4 Rank 2'
        call simulation_data_write(nid,data_name,data2,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to write rank 2 real4')

        call simulation_data_read(nid,data_name,rdata2,ierr)

        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to read rank 2 real4')
        flag = float_array_cmp(data2,rdata2,num)                  
        call assertEquals(0,flag,                                              &
                          'Failed to read correct data')


        ! Write/Read rank 3 REAL4 data
        deallocate(dims)
        allocate(dims(3))
        call generate_random_bound_array(1,100,dims,iseed)
        allocate(data3(dims(1),dims(2),dims(3)))
        allocate(rdata3(dims(1),dims(2),dims(3)))
        call generate_random_array(data3,iseed)

        data_name = 'REAL4 Rank 3'
        call simulation_data_write(nid,data_name,data3,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to write rank 3 real4')

        call simulation_data_read(nid,data_name,rdata3,ierr)

        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to read rank 3 real4')
        flag = float_array_cmp(data3,rdata3,num)                  
        call assertEquals(0,flag,                                              &
                          'Failed to read correct data')
        
        ! Test Cleanup
        deallocate(dims)
        deallocate(data1)
        deallocate(rdata1)
        deallocate(data2)
        deallocate(rdata2)
        deallocate(data3)
        deallocate(rdata3)
        call group_close(sid)
        call output_file_close(fid)

    end subroutine data_real4_test

    subroutine data_real8_test

        type(C_PTR)      :: fid, sid, nid
        integer          :: iseed = 444
        real(C_DOUBLE)   :: time
        integer          :: cyc, num, ierr, i, j, k, d
        integer          :: data0,rdata0
        integer(C_INT)   :: flag
        integer,dimension(:),allocatable :: dims
        real(C_DOUBLE),dimension(:),allocatable :: data1, rdata1
        real(C_DOUBLE),dimension(:,:),allocatable :: data2, rdata2
        real(C_DOUBLE),dimension(:,:,:),allocatable :: data3, rdata3
        character(len=128) :: data_name

        call set_unit_name('Series Data REAL8 Test')

        ! Test Setup 
        call setup_output_file(fid,sid)
         time = generate_random_double(iseed)
         cyc = generate_random_bound_int(1,10000,iseed)
        call sequence_next_id(sid,cyc,time,nid)

        ! Write/Read rank 0 REAL4 data
        data_name = 'REAL8 Rank 0'
        data0 = generate_random_float(iseed)
        call simulation_data_write(nid,data_name,data0,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to write rank 0 real8')

        call simulation_data_read(nid,data_name,rdata0,ierr)
        call assertEquals(DANU_SUCCESS, ierr,                                  &
                          'Failed to read rank 0 real8')
        call assertEquals(data0,rdata0,                                        &
                          'Failed to read correct dat8')

        
        ! Write/Read rank 1 REAL8 data
        num = generate_random_bound_int(128,512,iseed)
        allocate(data1(num))
        allocate(rdata1(num))
        call generate_random_array(data1,iseed)
        data_name = 'REAL8 Rank 1'
        call simulation_data_write(nid,data_name,data1,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to write rank 1 real8')

        call simulation_data_read(nid,data_name,rdata1,ierr)

        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to read rank 1 real8')
        flag = double_array_cmp(data1,rdata1,num)                  
        call assertEquals(0,flag,                                              &
                          'Failed to read correct data')
                          
        
        ! Write/Read rank 2 REAL8 data
        allocate(dims(2))
        call generate_random_bound_array(1,512,dims,iseed)
        allocate(data2(dims(1),dims(2)))
        allocate(rdata2(dims(1),dims(2)))
        call generate_random_array(data2,iseed)
        data_name = 'REAL8 Rank 2'
        call simulation_data_write(nid,data_name,data2,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to write rank 2 real8')

        call simulation_data_read(nid,data_name,rdata2,ierr)

        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to read rank 2 real8')
        flag = double_array_cmp(data2,rdata2,num)                  
        call assertEquals(0,flag,                                              &
                          'Failed to read correct data')


        ! Write/Read rank 3 REAL8 data
        deallocate(dims)
        allocate(dims(3))
        call generate_random_bound_array(1,100,dims,iseed)
        allocate(data3(dims(1),dims(2),dims(3)))
        allocate(rdata3(dims(1),dims(2),dims(3)))
        call generate_random_array(data3,iseed)

        data_name = 'REAL8 Rank 3'
        call simulation_data_write(nid,data_name,data3,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to write rank 3 real8')

        call simulation_data_read(nid,data_name,rdata3,ierr)

        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to read rank 3 real8')
        flag = double_array_cmp(data3,rdata3,num)                  
        call assertEquals(0,flag,                                              &
                          'Failed to read correct data')
        
        ! Test Cleanup
        deallocate(dims)
        deallocate(data1)
        deallocate(rdata1)
        deallocate(data2)
        deallocate(rdata2)
        deallocate(data3)
        deallocate(rdata3)
        call group_close(sid)
        call output_file_close(fid)

    end subroutine data_real8_test
 

end module funit_data_test

program data_test_driver
        
    use fruit
    use funit_data_test

! --- Local variables
      
      logical :: all_pass

! --- Code

! --- Initialize FRUIT
      call init_fruit

! --- Run tests
      call data_exists_test
      call data_count_test
      call data_type_test
      call data_rank_test
      call data_dimensions_test
      !call data_list_test
      call data_attribute_test
      call data_integer_test
      call data_real4_test
      call data_real8_test

      all_pass = .false.
      call is_all_successful(all_pass)

      if ( all_pass ) then
        write(*,*) 'ALL TESTS PASS'
      else
        write(*,*) 'FAIL Tests detected'
        call fail_exit_now
      endif

end program data_test_driver
