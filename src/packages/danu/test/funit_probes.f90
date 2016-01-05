!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module funit_probe_test

    use iso_c_binding
    use fruit

    use funit_utils
    use danu_module

    implicit none

    character(64) :: test_file = "dummy_file-probes.h5"
    character(64) :: test_sim_group = "Test Simulation"

    contains
    ! setup_before_all
    ! setup = setup_before_each
    subroutine setup_output_file(fid,sid)

        type(C_PTR), intent(out) :: fid, sid

        call output_file_create(test_file,fid)
        call simulation_add(fid,test_sim_group,sid)

    end subroutine setup_output_file

    subroutine create_probe_data(sid,num,names)

        type(C_PTR), intent(in) :: sid
        integer, intent(in)     :: num 
        character(len=*),dimension(:),intent(out) :: names

        integer        :: i, n, dlen
        integer(C_INT),dimension(1) :: dummy
        integer        :: iseed = 0
        type(C_PTR)    :: pid

        n = 1
        dlen=1
        i = 1
        do while ( i .le. num )
          write(names(i), '(a,i3)') 'Dummy Probe Data', i
          dummy(1) = generate_random_int(iseed)
          call probe_create_data(sid,names(i),dummy,pid)
          i = i + 1
        end do  

          
    end subroutine create_probe_data

    subroutine probe_exists_test

        type(C_PTR) :: fid,sid
        logical     :: flag
        character(len=64), dimension(1) :: names


        call set_unit_name('Probe Data Exists')
        call setup_output_file(fid,sid)

        ! Call when probe does not exist
        names(1) = 'Probe that is DNE'
        flag = probe_data_exists(sid,names(1))
        call assertEquals(.false.,flag,                         &
                           'Fail to return correct status')

        call create_probe_data(sid,1,names)
        flag = probe_data_exists(sid,names(1))
        call assertEquals(.true.,flag,                          &
                           'Fail to return correct status')

                           
        ! Clean up files                   
        call group_close(sid)
        call output_file_close(fid)
 
    end subroutine probe_exists_test

    subroutine probe_count_test

        type(C_PTR) :: fid,sid
        integer     :: err
        integer     :: cnt
        integer     :: num
        integer     :: iseed = 12345 
        character(len=64), dimension(:), allocatable :: names


        call set_unit_name('Probe Data Count')
        call setup_output_file(fid,sid)

        call probe_data_count(sid,cnt,err)
        call assertEquals(DANU_SUCCESS, err,                    &
                          'Failed to return correctly')
        call assertEquals(0,cnt,                                &
                          'Failed to return correct number 0')

        num = generate_random_bound_int(32,512,iseed)
        allocate(names(num))
        call create_probe_data(sid,num,names)
        deallocate(names)

        call probe_data_count(sid,cnt,err)
        call assertEquals(DANU_SUCCESS, err,                    &
                          'Failed to return correctly')
        call assertEquals(num,cnt,                                &
                          'Failed to return correct number')

        
        call group_close(sid)
        call output_file_close(fid)

    end subroutine probe_count_test
   
    subroutine probe_dimensions_test

        type(C_PTR) :: fid,sid, pid
        integer     :: err
        integer     :: num, dlen, idx
        integer     ::rnum,rdlen
        integer     :: iseed = 12345 
        character(len=32) :: probe_name
        integer, dimension(:,:),allocatable :: idata


        call set_unit_name('Probe Data Dimensions')
        call setup_output_file(fid,sid)

        probe_name='Probe Test Data'
        num = generate_random_bound_int(10,20,iseed)
        dlen = generate_random_bound_int(1,10,iseed)
        allocate(idata(dlen,num))
        call generate_random_array(idata,iseed)
        call probe_create_data(sid,probe_name,idata,pid)


        call probe_data_dimensions(sid,probe_name,rdlen,rnum,err)
        call assertEquals(DANU_SUCCESS, err,                    &
                          'Failed to return correctly')
        call assertEquals(rnum,num,                             &
                          'Failed to read the correct num')
        call assertEquals(rdlen,dlen,                             &
                          'Failed to read the correct length')

                          
        deallocate(idata)
        call group_close(sid)
        call output_file_close(fid)

    end subroutine probe_dimensions_test
    
  
    subroutine probe_list_test
        
        type(C_PTR) :: fid,sid
        integer     :: err
        integer     :: cnt
        integer     :: num
        logical     :: exists
        integer     :: iseed = 0
        character(len=64), dimension(:), allocatable :: names
        character(len=64), dimension(:), allocatable :: rnames

        call set_unit_name('Probe Data List')
        call setup_output_file(fid,sid)

        num = generate_random_bound_int(5,100,iseed)
        allocate(names(num))
        allocate(rnames(num))

        ! Call without probe data
        call probe_data_list(sid,rnames,err)
        call assertEquals(DANU_FAILURE,err,                        &
                          'Failed to return (F) status')
        
        ! Create probes then call again                  
        call create_probe_data(sid,num,names)
        call probe_data_list(sid,rnames,err)
        call assertEquals(DANU_SUCCESS,err,                        &
                          'Failed to return (P) status')

        cnt = 1
        do while ( cnt .le. num )
          exists = probe_data_exists(sid,rnames(cnt))
          call assertEquals(.true.,exists,                         &
                            'Failed to return probe that exists')
          cnt = cnt + 1
        end do   

        deallocate(names)
        deallocate(rnames)
        call group_close(sid)
        call output_file_close(fid)

    end subroutine probe_list_test

    subroutine probe_open_test

        type(C_PTR)       :: fid, sid,pid
        integer           :: ierr
        character(len=32),dimension(1) :: probe_name
        integer           :: dummy_data


        call set_unit_name('Probe Data Open Test')

        ! Test Setup
        call setup_output_file(fid,sid)

        ! Try to open probe that does not exist
        probe_name(1) = 'Probe DNE'
        call probe_data_open(sid,probe_name(1),pid,ierr)

        call assertEquals(DANU_FAILURE, ierr,                          &
                          'Fail to return correct status (F)')
       
        call create_probe_data(sid,1,probe_name)

        call probe_data_open(sid,probe_name(1),pid,ierr)
        call assertEquals(DANU_SUCCESS, ierr,                          &
                          'Failed to return correct status (P)')

        call group_close(sid)
        call output_file_close(fid)
                           
    end subroutine probe_open_test

    subroutine probe_integer_test

        type(C_PTR)      :: fid, sid, pid
        integer          :: iseed = 444
        integer          :: num, tot, dlen, ierr, i, j, k, d
        integer(C_INT)   :: flag
        integer,dimension(:),allocatable   :: idata0,rdata0
        integer,dimension(:,:),allocatable :: idata1, rdata1
        character(len=128) :: data_name

        call set_unit_name('Probe Data Integer Test')

        ! Test Setup 
        call setup_output_file(fid,sid)

        ! Write/Read rank 0 Integer data
        data_name = 'Integer Rank 0'
        num = generate_random_bound_int(1,1024,iseed)
        allocate(idata0(num))
        allocate(rdata0(num))
        call generate_random_array(idata0,iseed)
        call probe_create_data(sid,data_name,idata0,pid,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to create integer probe')
        call probe_data_read(pid,rdata0,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to read integer probe data')

        flag = int_array_cmp(idata0,rdata0,num)
        call assertEquals(flag,0,                       &
                          'Failed to read correct probe integer data')
                          

        ! Now loop through each case
        d = 2
        do while ( d .le. 7 )
          num = generate_random_bound_int(10,128,iseed)
          allocate(idata1(d,num))
          allocate(rdata1(d,num))
          call generate_random_array(idata1,iseed)
          rdata1=0
          write(data_name,'(a,i3)') 'Dummy Probe Integer Data', d
          call probe_create_data(sid,data_name,idata1,pid,ierr)
          call assertEquals(DANU_SUCCESS, ierr,                               &
                            'Failed to write rank > 1 data')
          call probe_data_read(pid,rdata1,ierr)
          call assertEquals(DANU_SUCCESS, ierr,                               &
                            'Failed to read integer data')
          tot=d*num
          flag = int_array_cmp(idata1,rdata1,tot)
          call assertEquals(flag, 0,                                          &
                            'Failed to read correct integer data')
          deallocate(rdata1)
          deallocate(idata1)
          d =d + 1
        end do

        call group_close(sid)
        call output_file_close(fid)

    end subroutine probe_integer_test
 
 
    subroutine probe_real4_test

        type(C_PTR)      :: fid, sid, pid
        integer          :: iseed = 444
        integer          :: num, tot, dlen, ierr, i, j, k, d
        integer(C_INT)   :: flag
        real,dimension(:),allocatable   :: data0,rdata0
        real,dimension(:,:),allocatable :: data1, rdata1
        character(len=128) :: data_name

        call set_unit_name('Probe Data REAL4 Test')

        ! Test Setup 
        call setup_output_file(fid,sid)

        ! Write/Read rank 0 REAL4 data
        data_name = 'REAL Rank 0'
        num = generate_random_bound_int(1,1024,iseed)
        allocate(data0(num))
        allocate(rdata0(num))
        call generate_random_array(data0,iseed)
        call probe_create_data(sid,data_name,data0,pid,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to create real4 probe')
        call probe_data_read(pid,rdata0,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to read real4 probe data')

        flag = float_array_cmp(data0,rdata0,num)
        call assertEquals(flag,0,                       &
                          'Failed to read correct probe real data')
                          

        ! Now loop through each case
        d = 2
        do while ( d .le. 7 )
          num = generate_random_bound_int(10,128,iseed)
          allocate(data1(d,num))
          allocate(rdata1(d,num))
          call generate_random_array(data1,iseed)
          rdata1=0
          write(data_name,'(a,i3)') 'Dummy Probe REAL4 Data', d
          call probe_create_data(sid,data_name,data1,pid,ierr)
          call assertEquals(DANU_SUCCESS, ierr,                               &
                            'Failed to write rank > 1 data')
          call probe_data_read(pid,rdata1,ierr)
          call assertEquals(DANU_SUCCESS, ierr,                               &
                            'Failed to read real4 data')
          tot=d*num
          flag = float_array_cmp(data1,rdata1,tot)
          call assertEquals(flag, 0,                                          &
                            'Failed to read correct real4 data')
          deallocate(rdata1)
          deallocate(data1)
          d =d + 1
        end do

        call group_close(sid)
        call output_file_close(fid)

    end subroutine probe_real4_test
 
    subroutine probe_real8_test

        type(C_PTR)      :: fid, sid, pid
        integer          :: iseed = 444
        integer          :: num, tot, dlen, ierr, i, j, k, d
        integer(C_INT)   :: flag
        real(C_DOUBLE),dimension(:),allocatable   :: data0,rdata0
        real(C_DOUBLE),dimension(:,:),allocatable :: data1, rdata1
        character(len=128) :: data_name

        call set_unit_name('Probe Data REAL8 Test')

        ! Test Setup 
        call setup_output_file(fid,sid)

        ! Write/Read rank 0 REAL4 data
        data_name = 'REAL Rank 0'
        num = generate_random_bound_int(1,1024,iseed)
        allocate(data0(num))
        allocate(rdata0(num))
        call generate_random_array(data0,iseed)
        call probe_create_data(sid,data_name,data0,pid,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to create real8 probe')
        call probe_data_read(pid,rdata0,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to read real8 probe data')

        flag = double_array_cmp(data0,rdata0,num)
        call assertEquals(flag,0,                       &
                          'Failed to read correct probe real data')
                          

        ! Now loop through each case
        d = 2
        do while ( d .le. 7 )
          num = generate_random_bound_int(10,128,iseed)
          allocate(data1(d,num))
          allocate(rdata1(d,num))
          call generate_random_array(data1,iseed)
          rdata1=0
          write(data_name,'(a,i3)') 'Dummy Probe REAL8 Data', d
          call probe_create_data(sid,data_name,data1,pid,ierr)
          call assertEquals(DANU_SUCCESS, ierr,                               &
                            'Failed to write rank > 1 data')
          call probe_data_read(pid,rdata1,ierr)
          call assertEquals(DANU_SUCCESS, ierr,                               &
                            'Failed to read real8 data')
          tot=d*num
          flag = double_array_cmp(data1,rdata1,tot)
          call assertEquals(flag, 0,                                          &
                            'Failed to read correct real8 data')
          deallocate(rdata1)
          deallocate(data1)
          d =d + 1
        end do

        call group_close(sid)
        call output_file_close(fid)

    end subroutine probe_real8_test

    subroutine probe_append_test

        type(C_PTR)      :: fid, sid, pid
        integer          :: iseed = 444
        integer          :: num, tot, dlen, ierr, i, j, k, d, idx, c_size
        integer(C_INT)   :: flag
        real(C_DOUBLE),dimension(:), allocatable   :: chunk0, data0,rdata0
        real(C_DOUBLE),dimension(:,:),allocatable :: chunk, my_data, rdata1
        character(len=128) :: data_name

        call set_unit_name('Probe Data REAL8 Test')

        ! Test Setup 
        call setup_output_file(fid,sid)

        ! Write/Read rank 0 REAL4 data
        data_name = 'REAL Rank 0'
        num = 50
        c_size = 5
        allocate(chunk0(c_size))
        allocate(data0(num))
        allocate(rdata0(num))
        call generate_random_array(chunk0,iseed)
        data0(1:5) = chunk0
        call probe_create_data(sid,data_name,chunk0,pid,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to create chunk real8 probe')
        i = 2
        j=6
        k=10
        do while (i .le. (num/c_size) )
            call generate_random_array(chunk0,iseed)
            call probe_data_write(pid,chunk0,ierr)
            call assertEquals(DANU_SUCCESS, ierr,                                &
                              'Failed to append probe data')
            data0(j:k)=chunk0
            j=k+1
            k=j+4
            i=i+1
        end do

        call probe_data_read(pid,rdata0,ierr)
        call assertEquals(DANU_SUCCESS,ierr,                                   &
                          'Failed to read real8 probe data')

        flag = double_array_cmp(data0,rdata0,num)
        call assertEquals(flag,0,                       &
                          'Failed to read correct probe real data')

        deallocate(chunk0)
        deallocate(data0)
        deallocate(rdata0)
                          

        ! Now loop through each case
        d = 2
        do while ( d .le. 7 )
          allocate(my_data(d,num))
          allocate(rdata1(d,num))
          allocate(chunk(d,c_size))
          call generate_random_array(chunk,iseed)
          rdata1=0
          write(data_name,'(a,i3)') 'Dummy Probe REAL8 Data', d
          call probe_create_data(sid,data_name,chunk,pid,ierr)
          call assertEquals(DANU_SUCCESS, ierr,                               &
                            'Failed to write rank > 1 data')
          my_data(:,1:c_size) = chunk
          j=6
          k=10
          i=2
          do while ( i .le. (num/c_size) )
            call generate_random_array(chunk,iseed)
            call probe_data_write(pid,chunk,ierr)
            call assertEquals(DANU_SUCCESS, ierr,                                &
                              'Failed to append data')
            my_data(:,j:k)=chunk
            j=k+1
            k=j+4
            i=i+1
          end do

          call probe_data_read(pid,rdata1,ierr)
          call assertEquals(DANU_SUCCESS, ierr,                               &
                            'Failed to read real8 data')
          tot=d*num
          flag = double_array_cmp(my_data,rdata1,tot)
          call assertEquals(flag, 0,                                          &
                            'Failed to read correct real8 data')
          deallocate(chunk)
          deallocate(rdata1)
          deallocate(my_data)
          d =d + 1
        end do

        call group_close(sid)
        call output_file_close(fid)

    end subroutine probe_append_test
 
end module funit_probe_test

program probe_test_driver
        
    use fruit
    use funit_probe_test

! --- Local variables
      
      logical :: all_pass

! --- Code

! --- Initialize FRUIT
      call init_fruit

! --- Run tests
      call probe_exists_test
      call probe_count_test
      call probe_dimensions_test
      call probe_open_test
      call probe_list_test
      call probe_integer_test
      call probe_real4_test
      call probe_real8_test
      call probe_append_test

      call fruit_summary

      all_pass = .false.
      call is_all_successful(all_pass)

      if ( all_pass ) then
        write(*,*) 'ALL TESTS PASS'
      else
        write(*,*) 'FAILED TESTS detected'
        call fail_exit_now
      endif

end program probe_test_driver
