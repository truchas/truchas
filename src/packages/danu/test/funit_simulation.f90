!
! Sample Fortran
!
module funit_sim_test

    use iso_c_binding
    use fruit

    use funit_utils
    use danu_module

    implicit none

    character(64) :: test_file = "dummy_simulation.h5"

    integer, parameter :: max_character_len = 64
    integer, parameter :: num_sims = 10

    contains

    subroutine generate_random_simulations(fid,cnt,names)
        type(C_PTR), intent(in) :: fid
        integer, intent(out)    :: cnt
        character(len=max_character_len), dimension(:), intent(out) :: names

        ! local
        integer :: iseed = 2379672
        integer :: idx
        integer :: names_num
        type(C_PTR) :: sid

        names_num = size(names)
        cnt = generate_random_bound_int(1,names_num,iseed)
        idx = 1
        do while (idx .le. cnt )
            write(names(idx),*) 'Simulation ', idx
            call simulation_add(fid,names(idx),sid)
            call group_close(sid)
            idx = 1 + idx
        end do    
            
    end subroutine generate_random_simulations     

    subroutine generate_simulations(fid,cnt,names)
        type(C_PTR), intent(in) :: fid
        integer, intent(in)    :: cnt
        character(len=max_character_len), dimension(:), intent(out) :: names

        ! local
        integer :: iseed = 2379672
        integer :: idx
        integer :: names_num
        type(C_PTR) :: sid

        names_num = size(names)
        idx = 1
        do while (idx .le. cnt )
            write(names(idx),*) 'Simulation ', idx
            call simulation_add(fid,names(idx),sid)
            call group_close(sid)
            idx = 1 + idx
        end do

    end subroutine generate_simulations   


    subroutine sim_exists_test

        type(C_PTR) :: hid_ptr, sid
        integer     :: err
        character(32) :: sim_name
        logical exists

        call set_unit_name('Simulation Exists Test')

        call output_file_create(test_file,hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to create test file')

        sim_name = 'Simulation That Does Not Exist'
        exists = simulation_exists(hid_ptr,sim_name)
        call assertEquals(.false.,exists,'Failed to return correct flag (F)')

        sim_name = 'Test Simulation'
        call simulation_add(hid_ptr,sim_name,sid)

        exists = simulation_exists(hid_ptr,sim_name)
        call assertEquals(.true.,exists, 'Failed to return correct flag (T)')

        call output_file_close(hid_ptr)

    end subroutine sim_exists_test

    subroutine sim_count_test

        type(C_PTR) :: hid_ptr
        integer     :: err
        integer     :: mcount
        integer     :: rcount
        character(len=max_character_len) :: sim_names(32)


        call set_unit_name('Simulation Count Test')

        call output_file_create(test_file,hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err)

        call simulation_count(hid_ptr,mcount,err)
        call assertEquals(DANU_SUCCESS,err,'Failed to count simulations')
        call assertEquals(0,mcount, 'Failed to return 0 sim count')


        ! Gnerate random number of Simulations
        call generate_random_simulations(hid_ptr,rcount,sim_names)

        call simulation_count(hid_ptr,mcount,err)
        call assertEquals(DANU_SUCCESS,err,'Failed to count simulations')
        call assertEquals(rcount,mcount, 'Failed to return (>0) sim count')

        call output_file_close(hid_ptr)

    end subroutine sim_count_test

    subroutine sim_list_test

        type(C_PTR) :: hid_ptr, mid_ptr
        integer     :: err
        integer     :: iseed = 48756
        integer     :: mcount 
        integer     :: read_count
        integer     :: idx
        integer     :: num_fail
        character(kind=C_CHAR, len=max_character_len) :: sim_names(num_sims) 
        character(kind=C_CHAR, len=max_character_len) :: test_names(num_sims) 
        character(kind=C_CHAR, len=max_character_len), dimension(:), allocatable :: gen_sims

        call set_unit_name('Simulation List Test')

        call output_file_create(test_file,hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err)

        call simulation_list(hid_ptr,test_names,err)
        call assertEquals(DANU_FAILURE, err,                                &
                          'Failed to raise error with no sims present')

        ! Create simes
        call generate_simulations(hid_ptr,num_sims,sim_names)
   
        ! Read the sim names
        call simulation_list(hid_ptr,test_names,err)
        call assertEquals(DANU_SUCCESS, err, &
                          "Failed to return successful status")
       
        idx = 1
        do while (idx .le. num_sims)
            call assertEquals(sim_names(idx),test_names(idx))
            idx = idx + 1
        end do 

        call output_file_close(hid_ptr)

        ! Recreate the file with different number of simulations

        ! Case: number sims < than the size of the array
        call output_file_create(test_file,hid_ptr)

        mcount = num_sims + 1
        do while ( mcount .ge. num_sims )
            call generate_random_simulations(hid_ptr,mcount,sim_names)
        end do    

        call simulation_list(hid_ptr,test_names,err)
        call assertEquals(DANU_SUCCESS, err, &
                          "Failed to return successful status")
       
        idx = 1
        do while (idx .le. mcount)
            call assertEquals(sim_names(idx),test_names(idx))
            idx = idx + 1
        end do 
     
        call output_file_close(hid_ptr)                  

        ! Case: number of sims > than the size of the array
        mcount = generate_random_bound_int(num_sims+1, 10*num_sims, iseed)
        allocate(gen_sims(mcount))
        
        call output_file_create(test_file,hid_ptr)

        call generate_simulations(hid_ptr,mcount,gen_sims)

        call simulation_list(hid_ptr,test_names,err)
        call assertEquals(DANU_FAILURE, err, &
                          "Failed to return fail status when array too small")
       
        idx = 1
        do while (idx .le. num_sims)
            call assertEquals(sim_names(idx),test_names(idx))
            idx = idx + 1
        end do 


        call output_file_close(hid_ptr)

       
    end subroutine sim_list_test
    
    subroutine sim_open_test

        type(C_PTR) :: fid
        integer     :: err
        type(C_PTR) :: sid
        character(len=64)   :: sim_name

        call set_unit_name('Simulation Open')
        call output_file_create(test_file,fid,err)
        call assertEquals(DANU_SUCCESS,err)

        ! Try to open simulation that is DNE
        sim_name = 'Test Simulation'
        call simulation_open(fid,sim_name,sid,err)
        call assertEquals(DANU_FAILURE, err, &
                         'Failed to return fail status when sim is DNE')
        
        call output_file_close(fid)

    end subroutine sim_open_test

    subroutine sim_add_test

        type(C_PTR) :: fid
        integer     :: err
        type(C_PTR) :: sid
        character(len=64)   :: sim_name

        call set_unit_name('Simulation Create')

        call output_file_create(test_file,fid,err)
        call assertEquals(DANU_SUCCESS,err)

        ! These should all pass
        sim_name = 'Test Simulation' 
        call simulation_add(fid,sim_name,sid,err) 
        call assertEquals(DANU_SUCCESS,err, 'Failed to add simulation')

        call output_file_close(fid)

    end subroutine sim_add_test

    subroutine sim_link_mesh

        type(C_PTR) :: fid
        integer     :: err
        type(C_PTR) :: sid
        character(len=64)   :: sim_name
        character(len=64)   :: mesh_name

        call set_unit_name('Simulation Create')

        call output_file_create(test_file,fid,err)
        call assertEquals(DANU_SUCCESS,err)

        ! These should all pass
        sim_name = 'Test Simulation' 
        call simulation_add(fid,sim_name,sid,err) 
        call assertEquals(DANU_SUCCESS,err, 'Failed to add simulation')

        mesh_name = 'Test Simulation Mesh'
        call simulation_link_mesh(fid,sid,mesh_name,err)
        call assertEquals(DANU_SUCCESS,err, &
                          'Failed to create soft link to mesh')

        call output_file_close(fid)


    end subroutine sim_link_mesh

   
end module funit_sim_test


! --- Main test driver program
program sim_test_driver

    use fruit
    use funit_sim_test

    logical :: all_pass

! --- Code

! --- Initialize FRUIT
      call init_fruit

! --- Run tests
      call sim_add_test
      call sim_open_test
      call sim_exists_test
      call sim_count_test
      call sim_list_test
      call sim_link_mesh


! --- Report results      
      call fruit_summary

! --- Exit if any test failed
      call is_all_successful(all_pass)

      if ( .not. all_pass ) then
          print *,'FAILED tests'
          call fail_exit_now
      end if    
      


end program sim_test_driver    
