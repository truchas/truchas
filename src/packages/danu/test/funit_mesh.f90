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
module funit_mesh_test

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

    subroutine generate_random_coordinates(x,y,z)
    
        real(C_DOUBLE), dimension(:), intent(inout) :: x
        real(C_DOUBLE), dimension(:), intent(inout), optional :: y
        real(C_DOUBLE), dimension(:), intent(inout), optional :: z

        integer :: nnodes, i
        integer :: seed

        seed = 1234567
        nnodes = size(x)
        i = 1
        do while (i .le. nnodes)
            x(i) = generate_random_double(seed)
            i = i + 1
        end do

        if (present(y)) then
            i = 1
            nnodes = size(y)
            do while ( i .le. nnodes)
                y(i) = generate_random_double(seed)
                i = i + 1
            end do
        end if
        
        if (present(z)) then
            i = 1
            nnodes = size(z)
            do while (i .le. nnodes)
                z(i) = generate_random_double(seed)
                i = i + 1
            end do
        end if

    end subroutine generate_random_coordinates

    subroutine generate_random_connectivity(idata)
    
        integer(C_INT), dimension(:,:), intent(inout) :: idata

        integer :: nelem, elem_order, i, n
        integer :: seed

        seed = 1234567
        elem_order = size(idata,1)
        nelem = size(idata,2)
        n = 1
        do while (n .le. nelem)
            i=1
            do while ( i .le. elem_order )
              idata(i,n) = generate_random_int(seed)
              i = i + 1
            end do  
            n=n+1  
        end do

    end subroutine generate_random_connectivity

    subroutine mesh_exists_test

        type(C_PTR) :: hid_ptr
        integer     :: err
        character(32) :: mesh_name
        logical exists

        call set_unit_name('Mesh Exists Test')

        call output_file_create(test_file,hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to create test file')

        mesh_name = 'Mesh That Does Not Exist'
        call mesh_exists(hid_ptr,mesh_name, exists,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to return successfully')
        call assertEquals(.false.,exists, 'Failed to return correct flag (F)')

        call output_file_close(hid_ptr)

    end subroutine mesh_exists_test

    subroutine mesh_count_test

        type(C_PTR) :: hid_ptr
        integer     :: err
        integer     :: mcount 


        call set_unit_name('Mesh Count Test')

        call output_file_create(test_file,hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err)

        call mesh_count(hid_ptr, mcount,err)
        call assertEquals(DANU_SUCCESS,err,'Failed to count meshes')
        call assertEquals(0,mcount, 'Failed to return 0 mesh count')


        call output_file_close(hid_ptr)

    end subroutine mesh_count_test

    subroutine mesh_list_test

        type(C_PTR) :: hid_ptr, mid_ptr
        integer     :: err
        integer     :: mcount = 10
        integer     :: read_count
        integer     :: idx
        character(kind=C_CHAR, len=32) :: mesh_names(10) 
        character(kind=C_CHAR, len=32) :: test_names(10) 

        character(kind=C_CHAR, len=32), allocatable :: read_mesh_names(:)

        call set_unit_name('Mesh List Test')

        call output_file_create(test_file,hid_ptr,err)
        call assertEquals(DANU_SUCCESS,err)

        call mesh_list(hid_ptr,mesh_names,err)
        call assertEquals(DANU_FAILURE, err,                                &
                          'Failed to raise error with no meshes present')

        ! Create meshes
        idx = 1
        do while ( idx .le. mcount )
          write(test_names(idx), '(a,i4)') 'Mesh ', idx
          call mesh_add_unstructured(hid_ptr,test_names(idx),3,2,mid_ptr,err)
          idx = idx + 1
        end do

        ! Read meshes and verify names are correct
        call mesh_list(hid_ptr,mesh_names,err)
        call assertEquals(DANU_SUCCESS, err,                                &
                          'Failed to read mesh list')

        idx = 1
        do while ( idx .le. mcount )
          write(*,*) test_names(idx), mesh_names(idx)
          call assertEquals(test_names(idx), mesh_names(idx))
          idx = idx + 1
        end do  
          

        call output_file_close(hid_ptr)                  

        ! Reopen file

        call output_file_open_rdonly(test_file,hid_ptr)
        
        call mesh_count(hid_ptr,mcount,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to count meshes')
        read_count = mcount - 2
        allocate(read_mesh_names(read_count))
        call mesh_list(hid_ptr,read_mesh_names,err)
        call assertEquals(DANU_FAILURE,err,                                   &
                          'Failed to raise error with too small array')

        idx = 1
        do while ( idx .le. read_count )
          call assertEquals(test_names(idx), read_mesh_names(idx),            &
                            'Failed to read the correct mesh name')
          !write(*,*) test_names(idx), read_mesh_names(idx)                  
          idx = idx + 1
        end do  

        deallocate(read_mesh_names)

        ! Now send a large array
        read_count = 2*mcount
        allocate(read_mesh_names(read_count))
        call mesh_list(hid_ptr,read_mesh_names,err)
        call assertEquals(DANU_SUCCESS,err,                                   &
                          'Failed to raise error with too small array')

        idx = 1
        do while ( idx .le. mcount )
          call assertEquals(test_names(idx), read_mesh_names(idx),            &
                            'Failed to read the correct mesh name')
          !write(*,*) test_names(idx), read_mesh_names(idx)                  
          idx = idx + 1
        end do 

        idx = mcount + 1
        do while ( idx .le. read_count )
          call assertEquals(0, len_trim(read_mesh_names(idx)),                 &
                            'Failed to return empty string too large array')
          idx = idx + 1
        end do  

        deallocate(read_mesh_names)


        call output_file_close(hid_ptr)                  

    end subroutine mesh_list_test
    
    subroutine mesh_open_test

        type(C_PTR) :: fid
        integer     :: err
        type(C_PTR) :: mid
        integer     :: d
        integer     :: elemorder
        character(len=64)   :: mesh_name

        call set_unit_name('Mesh Open')

        call output_file_create(test_file,fid,err)
        call assertEquals(DANU_SUCCESS,err)
 
        d = 2
        elemorder = 3
        mesh_name = '2D Unstruct Tri'
        call mesh_add_unstructured(fid,mesh_name,elemorder,d,mid,err) 
        call assertEquals(DANU_SUCCESS,err, 'Failed to add 2D tri mesh')

        call output_file_close(fid)

        ! Reopen file and mesh
        call output_file_open_rdonly(test_file,fid,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to open test file (rdonly)')

        ! Open mesh
        call mesh_open(fid,mesh_name,mid,err)
        call assertEquals(DANU_SUCCESS,err, 'Failed to open mesh')

        call output_file_close(fid)

    end subroutine mesh_open_test

    subroutine mesh_add_test

        type(C_PTR) :: fid
        integer     :: err
        type(C_PTR) :: mid
        integer     :: d
        integer     :: elemorder
        character(len=64)   :: mesh_name

        call set_unit_name('Mesh Add Unstruct')

        call output_file_create(test_file,fid,err)
        call assertEquals(DANU_SUCCESS,err)

        ! These should all pass 
        d = 1
        elemorder = LINE_ELEM_ORDER 
        mesh_name = '1D Unstruct'
        call mesh_add_unstructured(fid,mesh_name,elemorder,d,mid,err) 
        call assertEquals(DANU_SUCCESS,err, 'Failed to add 1D mesh')

        d = 2
        elemorder = TRI_ELEM_ORDER
        mesh_name = '2D Unstruct Tri'
        call mesh_add_unstructured(fid,mesh_name,elemorder,d,mid,err) 
        call assertEquals(DANU_SUCCESS,err, 'Failed to add 2D tri mesh')

        d = 2
        elemorder = QUAD_ELEM_ORDER
        mesh_name = '2D Unstruct Quad'
        call mesh_add_unstructured(fid,mesh_name,elemorder,d,mid,err) 
        call assertEquals(DANU_SUCCESS,err, 'Failed to add 2D quad mesh')

        d = 3
        elemorder = TET_ELEM_ORDER
        mesh_name = '3D Unstruct Tet'
        call mesh_add_unstructured(fid,mesh_name,elemorder,d,mid,err) 
        call assertEquals(DANU_SUCCESS,err, 'Failed to add 3D tet mesh')

        d = 3
        elemorder = QUAD_ELEM_ORDER
        mesh_name = '3D Unstruct Quad'
        call mesh_add_unstructured(fid,mesh_name,elemorder,d,mid,err) 
        call assertEquals(DANU_SUCCESS,err, 'Failed to add 3D hex mesh')

        call output_file_close(fid)

    end subroutine mesh_add_test

    subroutine mesh_coordinates
        
        type(C_PTR) :: fid
        integer     :: err
        type(C_PTR) :: mid1,mid2,mid3
        character(kind=C_CHAR,len=32) :: mesh_name
        integer(C_INT) :: n, nnodes, elem_order, mesh_dim
        real(C_DOUBLE), dimension(:), allocatable :: x,y,z
        real(C_DOUBLE), dimension(:), allocatable :: rx,ry,rz, buf
        logical all_match
        integer(C_INT) :: seed

        
        call set_unit_name('Mesh Write Coordinates')

        ! Create the file and mesh
        call output_file_create(test_file,fid,err)
        call assertEquals(DANU_SUCCESS,err)
       
        mesh_name = 'Test Mesh 3D'
        mesh_dim = 3
        elem_order = HEX_ELEM_ORDER 
        call mesh_add_unstructured(fid,mesh_name,elem_order,mesh_dim,mid3,err)
        call assertEquals(DANU_SUCCESS,err,                                    &
                          'Failed to create 3D mesh for coord test')
        
        mesh_name = 'Test Mesh 2D'
        mesh_dim = 2
        elem_order = QUAD_ELEM_ORDER
        call mesh_add_unstructured(fid,mesh_name,elem_order,mesh_dim,mid2,err)
        call assertEquals(DANU_SUCCESS,err,                                    &
                          'Failed to create 2D mesh for coord test')
        
        mesh_name = 'Test Mesh 1D'
        mesh_dim = 1
        elem_order = LINE_ELEM_ORDER
        call mesh_add_unstructured(fid,mesh_name,elem_order,mesh_dim,mid1,err)
        call assertEquals(DANU_SUCCESS,err,                                    &
                          'Failed to create 1D mesh for coord test')
        
        ! Gernerate random node coordinates
        seed = 1234567
        nnodes = generate_random_bound_int(1024,2048,seed)
        allocate(x(nnodes))
        allocate(rx(nnodes))
        allocate(y(nnodes))
        allocate(ry(nnodes))
        allocate(z(nnodes))
        allocate(rz(nnodes))
        allocate(buf(nnodes))
        call generate_random_coordinates(x,y,z)

        ! Calls that should fail
        call set_unit_name('Mesh Write Coordinates Fail Tests')

        call mesh_write_coordinates(mid1,nnodes,x,y,z,err)
        call assertEquals(DANU_FAILURE,err, 'Failed to flag bad mesh id')

        call mesh_write_coordinates(mid3,nnodes,x,y,err)
        call assertEquals(DANU_FAILURE,err, 'Failed to flag bad mesh id')

        ! Calls that should pass
        call set_unit_name('Mesh Write Coordinates Pass Tests')

        call mesh_write_coordinates(mid3,nnodes,x,y,z,err)
        call assertEquals(DANU_SUCCESS,err)

        call mesh_write_coordinates(mid2,nnodes,x,y,err)
        call assertEquals(DANU_SUCCESS,err)

        call mesh_write_coordinates(mid1,nnodes,x,err)
        call assertEquals(DANU_SUCCESS,err)

        call output_file_close(fid)

        ! Now read the coordinates
        call output_file_open_rdonly(test_file,fid);

        ! Calls designed to fail
        call set_unit_name('Mesh Read Coordinate Fail Tests')

        mesh_name = 'Test Mesh 1D'
        call mesh_open(fid,mesh_name,mid1)
        call mesh_read_coordinates(mid1,rx,ry,err)
        call assertEquals(DANU_FAILURE,err, 'Failed to flag dim mismatch')
        
        call mesh_read_coordinates(mid1,rx,ry,rz,err)
        call assertEquals(DANU_FAILURE,err, 'Failed to flag dim mismatch')

        ! Calls designed to pass
        call set_unit_name('Mesh Read Coordinate Pass Tests')

        ! 1D Mesh
        mesh_name = 'Test Mesh 2D'
        call mesh_read_coordinates(mid1,rx,err)
        call assertEquals(DANU_SUCCESS, err, 'Failed to read 1D mesh')

        n = 1
        all_match = .true.
        do while ( n .le. nnodes )
            if ( x(n) .ne. rx(n) ) then
                all_match = .false.
            endif    
            
            n = n + 1
        enddo
        call assertEquals(.true.,all_match, 'Failed to read 1D mesh correctly')

        ! 2D Mesh
        mesh_name = 'Test Mesh 2D'
        call mesh_open(fid,mesh_name,mid2)
        call mesh_read_coordinates(mid2,rx,ry,err)
        call assertEquals(DANU_SUCCESS, err, 'Failed to read 2D mesh')

        n = 1
        all_match = .true.
        do while ( n .le. nnodes )
            if ( x(n) .ne. rx(n) ) then
                all_match = .false.
            endif    
            
            if ( y(n) .ne. ry(n) ) then
                all_match = .false.
            endif

            n = n + 1
        enddo
        call assertEquals(.true.,all_match, 'Failed to read 2D mesh correctly')

        ! Read 2D by index
        call mesh_read_coordinates(mid2,1,buf)
        n = 1
        all_match=.true.
        do while ( n .le. nnodes )
          if ( x(n) .ne. buf(n) ) then
            all_match = .false.
          end if
          n = n + 1
        end do
        call assertEquals(.true.,all_match, 'Failed to read 2D mesh by index')

        call mesh_read_coordinates(mid2,2,buf)
        n = 1
        all_match=.true.
        do while ( n .le. nnodes )
          if ( y(n) .ne. buf(n) ) then
            all_match = .false.
          end if
          n = n + 1
        end do
        call assertEquals(.true.,all_match, 'Failed to read 2D mesh by index')

        
        ! 3D Mesh        
        mesh_name = 'Test Mesh 3D'
        call mesh_open(fid,mesh_name,mid3)
        call mesh_read_coordinates(mid3,rx,ry,rz,err)
        call assertEquals(DANU_SUCCESS, err, 'Failed to read 3D mesh')

        n = 1
        all_match = .true.
        do while ( n .le. nnodes )
            if ( x(n) .ne. rx(n) ) then
                all_match = .false.
            endif    
            
            if ( y(n) .ne. ry(n) ) then
                all_match = .false.
            endif    
 
            if ( z(n) .ne. rz(n) ) then
                all_match = .false.
            endif    
       
            n = n + 1
        end do 
        call assertEquals(.true.,all_match, 'Failed to read 3D mesh correctly')

        ! Read 3D by index
        call mesh_read_coordinates(mid3,1,buf)
        n = 1
        all_match=.true.
        do while ( n .le. nnodes )
          if ( x(n) .ne. buf(n) ) then
            all_match = .false.
          end if
          n = n + 1
        end do
        call assertEquals(.true.,all_match, 'Failed to read 3D mesh by index')

        call mesh_read_coordinates(mid3,2,buf)
        n = 1
        all_match=.true.
        do while ( n .le. nnodes )
          if ( y(n) .ne. buf(n) ) then
            all_match = .false.
          end if
          n = n + 1
        end do
        call assertEquals(.true.,all_match, 'Failed to read 3D mesh by index')
 
        call mesh_read_coordinates(mid3,3,buf)
        n = 1
        all_match=.true.
        do while ( n .le. nnodes )
          if ( z(n) .ne. buf(n) ) then
            all_match = .false.
          end if
          n = n + 1
        end do
        call assertEquals(.true.,all_match, 'Failed to read 3D mesh by index')

        call output_file_close(fid)

        deallocate(x)
        deallocate(rx)
        deallocate(y)
        deallocate(ry)
        deallocate(z)
        deallocate(rz)
        deallocate(buf)

    end subroutine mesh_coordinates

    subroutine mesh_connectivity
   
        type(C_PTR) :: fid
        integer     :: err
        type(C_PTR) :: mid
        character(kind=C_CHAR,len=32) :: mesh_name
        character(kind=C_CHAR,len=48) :: data_name
        integer(C_INT) :: n, i, num, elem_order, mesh_dim
        integer(C_INT), dimension(:,:), allocatable :: idata, ridata
        logical all_match
        integer(C_INT) :: iseed, r_num

        
        call set_unit_name('Mesh Write Connectivity')

        ! Create the file and mesh
        call output_file_create(test_file,fid)
       
        mesh_name = 'Test Mesh 3D'
        mesh_dim = 3
        elem_order = HEX_ELEM_ORDER
        call mesh_add_unstructured(fid,mesh_name,elem_order,mesh_dim,mid)


        ! Write test (PASS)
        iseed = 342567
        num = generate_random_bound_int(2048,5000,iseed)
        allocate(idata(elem_order,num))
        call generate_random_connectivity(idata)

        call set_unit_name('Mesh Write Connectivity (PASS)')

        call mesh_write_connectivity(mid,num,idata,err)
        call assertEquals(DANU_SUCCESS,err,                                    &
                          'Failed to write connectivity data')

       
        call output_file_close(fid)

        ! Reopen for read test
        call set_unit_name('Mesh Read Connectivity (PASS)')

        call output_file_open_rdonly(test_file,fid)

        call mesh_open(fid,mesh_name,mid)

        allocate(ridata(elem_order,num))

        call mesh_read_connectivity(mid,ridata,err)
        call assertEquals(DANU_SUCCESS,err,                                    &
                          'Failed to read connectivity data')

        all_match = .true.
        n = 1
        do while ( n .le. num )
           i=1
           do while ( i .le. elem_order )
             if ( ridata(i,n) .ne. idata(i,n) ) then
                 all_match = .false.
             endif    
             i=i+1
           end do    
           n = n + 1
        end do   

        call assertEquals(.true.,all_match,                                    &
                          'Failed to read connectivity data correctly')
        call output_file_close(fid)

    end subroutine mesh_connectivity

    subroutine mesh_attributes
       
        type(C_PTR) :: fid
        integer     :: err
        type(C_PTR) :: mid
        character(kind=C_CHAR,len=32) :: mesh_name
        character(kind=C_CHAR,len=48) :: data_name
        integer(C_INT) :: mesh_type,  mesh_dim,  elem_order, nnodes, nelem,    &
                          elem_type 
        integer(C_INT) :: rmesh_type, rmesh_dim, relem_order, rnnodes, rnelem, &
                          relem_type
        integer(C_INT), dimension(:,:), allocatable :: idata
        real(C_DOUBLE), dimension(:), allocatable :: x,y,z
        integer(C_INT) :: iseed

        call set_unit_name('Mesh Attributes')

        ! Create a test file without a mesh
         call output_file_create(test_file,fid)

        ! Create a test mesh
        mesh_name = 'Test Mesh Attributes'
        mesh_type = UNSTRUCTURED_MESH
        mesh_dim = 3
        elem_order = TET_ELEM_ORDER
        elem_type = TET_ELEM
        call mesh_add_unstructured(fid,mesh_name,elem_order,mesh_dim,mid, err)
        call assertEquals(DANU_SUCCESS,err,         &
                          'Failed to create a test mesh')

        ! Create garbage data for the mesh                   

        ! Gernerate random node coordinates
        iseed = 1234567
        nnodes = generate_random_bound_int(1000,5000,iseed)
        allocate(x(nnodes))
        allocate(y(nnodes))
        allocate(z(nnodes))
        call generate_random_coordinates(x,y,z)
        call mesh_write_coordinates(mid,nnodes,x,y,z,err)

        ! Generate random connectivity data
        nelem = nnodes / TET_ELEM_ORDER
        allocate(idata(TET_ELEM_ORDER,nelem))
        call generate_random_connectivity(idata)
        call mesh_write_connectivity(mid,nelem,idata)

        ! Close the file
        call output_file_close(fid)

        ! Calls should fail since the mid is no longer vaild after file close
        call mesh_get_type(mid,rmesh_type,err)
        call assertEquals(DANU_FAILURE,err,         &
                          'Failed to raise error when reading mesh type')
      
        call mesh_get_elementtype(mid,relem_type,err)
        call assertEquals(DANU_FAILURE,err,         &
                          'Failed to raise error when reading element type')
       
        call mesh_get_dimension(mid,rmesh_dim,err)
        call assertEquals(DANU_FAILURE,err,         &
                          'Failed to raise error when reading mesh dimension')
       
        call mesh_get_nnodes(mid,rnnodes,err)
        call assertEquals(DANU_FAILURE,err,         &
                          'Failed to raise error when reading number of nodes')
       
        call mesh_get_nelem(mid,rnelem,err)
        call assertEquals(DANU_FAILURE,err,         &
                          'Failed to raise error when reading number of elements')
        
        call mesh_get_elem_order(mid,relem_order,err)
        call assertEquals(DANU_FAILURE,err,         &
                          'Failed to raise error when reading element order')

        ! Re-open file and read attributes
        call output_file_open_rdonly(test_file,fid)

        call mesh_open(fid,mesh_name,mid)

        ! Passing tests
        call mesh_get_type(mid,rmesh_type,err)
        call assertEquals(DANU_SUCCESS,err,         &
                          'Failed to raise error when reading mesh type')
        call assertEquals(mesh_type,rmesh_type,     &
                          'Failed to read correct mesh_type')
       
        call mesh_get_elementtype(mid,relem_type,err)
        call assertEquals(DANU_SUCCESS,err,         &
                          'Failed to raise error when reading element type')
        call assertEquals(elem_type,relem_type,     &
                          'Failed to read correct elem_type')
       
        call mesh_get_dimension(mid,rmesh_dim,err)
        call assertEquals(DANU_SUCCESS,err,         &
                          'Failed to raise error when reading mesh dimension')
        call assertEquals(mesh_dim,rmesh_dim,     &
                          'Failed to read correct mesh_dim')
       
        call mesh_get_nnodes(mid,rnnodes,err)
        call assertEquals(DANU_SUCCESS,err,         &
                          'Failed to raise error when reading number of nodes')
        call assertEquals(nnodes,rnnodes,     &
                          'Failed to read correct nnodes')
       
       
        call mesh_get_nelem(mid,rnelem,err)
        call assertEquals(DANU_SUCCESS,err,         &
                          'Failed to read number of elements')
        call assertEquals(nelem,rnelem,     &
                          'Failed to read correct number of elements nelem')
       
        
        call mesh_get_elem_order(mid,relem_order,err)
        call assertEquals(DANU_SUCCESS,err,         &
                          'Failed to raise error when reading element order')
        call assertEquals(elem_order,relem_order,     &
                          'Failed to read correct elem_order')
       
 
       
        call output_file_close(fid)

    endsubroutine mesh_attributes

    subroutine mesh_wrappers
   
        type(C_PTR) :: fid
        integer     :: err
        type(C_PTR) :: mid
        character(kind=C_CHAR,len=32) :: mesh_name
        integer :: nelem, nnodes
        real(C_DOUBLE), dimension(:), allocatable :: x,y,z
        integer(C_INT), dimension(:,:), allocatable :: idata
        integer(C_INT) :: iseed

        
        call set_unit_name('Mesh Wrappers')

        ! Create the file and mesh
        call output_file_create(test_file,fid)

        ! Write test (PASS) Dummy Data
        iseed = 342567
        nelem = generate_random_bound_int(2048,5000,iseed)
        nnodes=nelem/HEX_ELEM_ORDER
        allocate(idata(HEX_ELEM_ORDER,nelem))
        allocate(x(nnodes))
        allocate(y(nnodes))
        allocate(z(nnodes))
        call generate_random_connectivity(idata)
        call generate_random_coordinates(x,y,z)


        mesh_name = 'Test Mesh 3D'
        call mesh_create_hex_unstruct(fid,mesh_name,nelem,nnodes,x,y,z,idata,mid,err)
        call assertEquals(DANU_SUCCESS,err,         &
                          'Failed to write HEX mesh with a wrapper')
        call output_file_close(fid)

        end subroutine mesh_wrappers   

end module funit_mesh_test


! --- Main test driver program
program mesh_test_driver

    use fruit
    use funit_mesh_test

    logical :: all_pass

! --- Code

! --- Initialize FRUIT
      call init_fruit

! --- Run tests
      !call mesh_add_test
      !call mesh_open_test
      !call mesh_exists_test
      !call mesh_count_test
      !call mesh_coordinates
      !call mesh_connectivity
      !call mesh_list_test
      !call mesh_attributes
      call mesh_wrappers


! --- Report results      
      call fruit_summary

! --- Exit if any test failed
      call is_all_successful(all_pass)

      if ( .not. all_pass ) then
          print *,'FAILED tests'
          call fail_exit_now
      end if    
      


end program mesh_test_driver    
