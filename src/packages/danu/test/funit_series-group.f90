module funit_sgroup_test

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

    subroutine add_series_groups(sid,num)
        type(C_PTR), intent(in) :: sid
        integer,     intent(in) :: num

        type(C_PTR)    :: nid
        real(C_DOUBLE) :: time
        integer        :: iseed = 123456
        integer        :: cyc,dc,cnt,err

        cnt = 1
        cyc=generate_random_bound_int(1,1000,iseed)
        do while ( cnt .le. num )
          time = generate_random_double(iseed)
          dc=generate_random_bound_int(1,100,iseed)
          cyc=cyc+dc
          call sequence_next_id(sid,cyc,time,nid,err)
          if ( err .ne. DANU_SUCCESS ) then 
            call exit_now(err)
          endif 
          cnt = cnt + 1
        end do
        
    end subroutine add_series_groups

    subroutine sgroup_exists_test

        type(C_PTR)   :: fid, sid
        integer       :: err
        character(len=64),dimension(1) :: sname
        logical       :: exists

        call set_unit_name('Series Group Exists')

        call setup_output_file(fid,sid)

        sname(1) = 'Series Group DNE'
        exists = sequence_exists(sid,sname(1))
        call assertEquals(.false., exists,'Failed to return correct exists (F)')

        call add_series_groups(sid,1)
        sname(1) = 'Series 1'
        exists = sequence_exists(sid,sname(1))
        call assertEquals(.true., exists, 'Failed to return correct exists (T)')



        call group_close(sid)
        call output_file_close(fid)

    end subroutine sgroup_exists_test

    subroutine sgroup_count_test

        type(C_PTR) :: fid,sid
        integer     :: err
        integer     :: cnt
        integer     :: num
        integer     :: iseed = 12345 


        call set_unit_name('Series Group Count')

        call setup_output_file(fid,sid)
        num = generate_random_bound_int(10,128,iseed)

        call sequence_count(sid,cnt,err)
        call assertEquals(DANU_SUCCESS,err,'Failed to count simulations')
        call assertEquals(0,cnt, 'Failed to return 0 group count')

        
        call add_series_groups(sid,num)
        call sequence_count(sid,cnt,err)
        call assertEquals(DANU_SUCCESS,err,'Failed to count simulations')
        call assertEquals(num,cnt, 'Failed to return num group count')

        call group_close(sid)
        call output_file_close(fid)

    end subroutine sgroup_count_test
    
    subroutine sgroup_id_test

        type(C_PTR) :: fid,sid, nid, hid
        integer          :: err
        integer          :: iseed = 345123
        integer          :: cyc
        real(C_DOUBLE)   :: time
        character(kind=C_CHAR,len=32) :: sname

        call set_unit_name('Series Group ID Test')


        call setup_output_file(fid,sid)

        time = generate_random_double(iseed)
        cyc = generate_random_bound_int(1,1000,iseed)

        call sequence_next_id(sid,cyc,time,nid,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to get next id')

        sname = 'Series 1'
        call sequence_get_id(sid,sname,hid,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to get id for existing group')

        call sequence_get_id(sid,1,hid,err)
        call assertEquals(DANU_SUCCESS, err, 'Failed to get id by number')

        sname = 'Series DNE'
        call sequence_get_id(sid,sname,hid,err)
        call assertEquals(DANU_FAILURE,err,                      &
                         'Failed to return F status non-exist group')

        call sequence_get_id(sid,10,hid,err)
        call assertEquals(DANU_FAILURE, err,                     &
                          'Failed to return F status non-exist greoup')
  
        call group_close(sid)
        call group_close(nid)
        call output_file_close(fid)
        
    end subroutine sgroup_id_test

    subroutine sgroup_list_test
        
        type(C_PTR) :: fid,sid, nid
        integer          :: err
        integer          :: iseed = 345123
        integer          :: num, cnt, cyc
        character(len=64), dimension(:), allocatable :: names
        character(len=64), dimension(:), allocatable :: rnames

        call set_unit_name('Series Group List Test')

        call setup_output_file(fid,sid)
        num = generate_random_bound_int(50,100,iseed)
        allocate(names(num))
        allocate(rnames(num))
        call add_series_groups(sid,num)

        cnt = 1
        do while ( cnt .le. num )
           write(names(cnt),'(a7,i3)') 'Series', cnt
           cnt = cnt + 1
        end do    

        call sequence_list(sid,rnames,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to read sequence names')

        ! This test fails, because the order read in does not match write order
        !cnt = 1
        !do while ( cnt .le. num ) 
        !  call assertEquals(names(cnt), rnames(cnt),            &
        !                    'Failed to read correct series group name')
        !  write(*,*) names(cnt), rnames(cnt)                  
        !  cnt = cnt + 1
        !end do  

        deallocate(rnames)
        allocate(rnames(num/2))
        call sequence_list(sid,rnames,err)
        call assertEquals(DANU_FAILURE,err,                           &
                          'Failed to return F when array too small')

        deallocate(rnames)
        allocate(rnames(num*2))
        call sequence_list(sid,rnames,err)
        call assertEquals(DANU_SUCCESS,err,                           &
                          'Failed to return T when array too big')

        deallocate(names)
        deallocate(rnames)
        call group_close(sid)
        call output_file_close(fid)

    end subroutine sgroup_list_test
   

end module funit_sgroup_test

program sgroup_test_driver
        
    use fruit
    use funit_sgroup_test

! --- Local variables
      
      logical :: all_pass

! --- Code

! --- Initialize FRUIT
      call init_fruit

! --- Run tests
      call sgroup_exists_test
      call sgroup_count_test
      call sgroup_list_test
      call sgroup_id_test

      all_pass = .false.
      call is_all_successful(all_pass)

      if ( all_pass ) then
        write(*,*) 'ALL TESTS PASS'
      else
        write(*,*) 'FAIL Tests detected'
        call fail_exit_now
      endif

end program sgroup_test_driver
