!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program grid_mapping_test_driver
  !=======================================================================
  ! Purpose(s):
  !
  !    Driver to test GRID_MAPPING_MODULE.  This driver requires
  !    EXODUS for reading exodus binary meshes.  (For a driver that
  !    does not require EXODUS module, use DRIVERA.F90 which uses
  !    ascii mesh files.)
  !
  ! Author(s): Andrew Kuprat (kuprat@lanl.gov)
  !
  !=======================================================================
  use grid_mapping_module, only : read_int_volumes, compute_int_volumes, &
       grid_vols, grid_vol_fracs, write_int_volumes, map_cell_field, &
       right_int_volumes, grid_int_vols, destroy_grid_int_vols, gm_mesh
  use exodus_mesh_type
  use exodus_mesh_io, only: read_exodus_mesh
  use grid_mapping_exodus, only : copy_exodus_mesh_to_gm_mesh
  use grid_mapping_utils, only : writegmv_gm_mesh

  implicit none

  integer, parameter :: dp=kind(1.0d0)

  integer, parameter :: task_lun=8, grid_lun=9, gmv_lun=10
  character(80) :: taskfile,filename_mesh_a,filename_mesh_b,label
  character(80) :: title1,title2,title3,fname
  logical :: reverse_order,exactly_conservative,preserve_constants,strict
  logical :: map_exists
  real(dp) :: defval
  type(gm_mesh), allocatable :: mesh_a, mesh_b
  type(exodus_mesh), allocatable :: exodus_mesh_a, exodus_mesh_b
  type(grid_int_vols) :: int_vols
  real(dp), dimension(:), pointer :: field_a=>null(), field_b=>null(), vol_a=>null(), vol_b=>null()
  real(dp), dimension(:), pointer :: vf_a=>null(), vf_b=>null()
  real(dp) :: minval_a, maxval_a, minval_b, maxval_b
  integer :: ier

  taskfile='taskfile'
  open (task_lun,file=taskfile)

  do 

     read (task_lun,*,end=20) title1
     read (task_lun,*) title2
     read (task_lun,*) title3
     read (task_lun,*) label, filename_mesh_a
     read (task_lun,*) label, filename_mesh_b
     read (task_lun,*) label, reverse_order
     read (task_lun,*) label, exactly_conservative
     read (task_lun,*) label, preserve_constants
     read (task_lun,*) label, strict
     read (task_lun,*) label, defval

     write(*,fmt='(/,a)') trim(title1)
     write(*,fmt='(a)') trim(title2)
     write(*,fmt='(a)') trim(title3)
     write(*,fmt='(2a)') 'filename_mesh_a= ',trim(filename_mesh_a)
     write(*,fmt='(2a)') 'filename_mesh_b= ',trim(filename_mesh_b)
     write(*,fmt='(a,l1)') 'reverse_order= ',reverse_order
     write(*,fmt='(a,l1)') 'exactly_conservative= ',exactly_conservative
     write(*,fmt='(a,l1)') 'preserve_constants= ',preserve_constants
     write(*,fmt='(a,l1)') 'strict= ',strict
     write(*,fmt='(a,es20.12)') 'defval= ',defval

     ! Read two exodus binary meshes
     allocate(exodus_mesh_a, exodus_mesh_b)
     call read_exodus_mesh(filename_mesh_a,exodus_mesh_a)
     call read_exodus_mesh(filename_mesh_b,exodus_mesh_b)

     ! Set up gm_mesh data object for both meshes
     allocate(mesh_a)
     call copy_exodus_mesh_to_gm_mesh(exodus_mesh_a,mesh_a)
     deallocate(exodus_mesh_a)

     allocate(mesh_b)
     call copy_exodus_mesh_to_gm_mesh(exodus_mesh_b,mesh_b)
     deallocate(exodus_mesh_b)

     ! Does the right mapping file already exist on disk?
     map_exists = .false.

     OPEN (UNIT = grid_lun, FILE='int_vols_file', STATUS='old', IOSTAT=ier, POSITION='rewind', &
          FORM = 'unformatted')
     
     if (ier.eq.0) then
        
        call read_int_volumes(int_vols, grid_lun, ier)
        close(grid_lun)
        if (ier.eq.0) then
           map_exists = right_int_volumes(mesh_a, mesh_b, int_vols)
        endif
     
     endif
     
        
     if (.not.map_exists) then

        write(*,fmt='(/,a,/)') '*** Recomputing grid mapping'

!!$        call compute_int_volumes(mesh_a,mesh_b,int_vols,maxwarn=100)
        call compute_int_volumes(mesh_a,mesh_b,int_vols)

        ! Write grid-grid overlap statistics
        if (associated(vol_a)) deallocate(vol_a)
        allocate(vol_a(mesh_a%nelt))
        if (associated(vol_b)) deallocate(vol_b)
        allocate(vol_b(mesh_b%nelt))

        call grid_vols(int_vols,vol_a,vol_b)
        minval_a=minval(vol_a)
        maxval_a=maxval(vol_a)
        minval_b=minval(vol_b)
        maxval_b=maxval(vol_b)

        write(*,fmt='(a,es20.12)') 'A mesh has MIN VOLUME= ',minval_a
        write(*,fmt='(a,es20.12)') 'A mesh has MAX VOLUME= ',maxval_a
        write(*,fmt='(a,es20.12)') 'B mesh has MIN VOLUME= ',minval_b
        write(*,fmt='(a,es20.12)') 'B mesh has MAX VOLUME= ',maxval_b

        deallocate (vol_a)
        deallocate (vol_b)

        if (associated(vf_a)) deallocate(vf_a)
        allocate(vf_a(mesh_a%nelt))
        if (associated(vf_b)) deallocate(vf_b)
        allocate(vf_b(mesh_b%nelt))

        call grid_vol_fracs(int_vols,vf_a,vf_b,strict=.false.)
        minval_a=minval(vf_a)
        maxval_a=maxval(vf_a)
        minval_b=minval(vf_b)
        maxval_b=maxval(vf_b)

        write(*,fmt='(/,a)')         'Disregarding block IDs, we have that'
        write(*,fmt='(a,es20.12)') 'A mesh elts have MIN VF in B mesh= ',minval_a
        write(*,fmt='(a,es20.12)') 'A mesh elts have MAX VF in B mesh= ',maxval_a
        write(*,fmt='(a,es20.12)') 'B mesh elts have MIN VF in A mesh= ',minval_b
        write(*,fmt='(a,es20.12)') 'B mesh elts have MAX VF in A mesh= ',maxval_b

        call grid_vol_fracs(int_vols,vf_a,vf_b,strict=.true.)
        minval_a=minval(vf_a)
        maxval_a=maxval(vf_a)
        minval_b=minval(vf_b)
        maxval_b=maxval(vf_b)

        write(*,fmt='(/,a)')         'Respecting block IDs, we have that'
        write(*,fmt='(a,es20.12)') 'A mesh elts have MIN VF in B mesh= ',minval_a
        write(*,fmt='(a,es20.12)') 'A mesh elts have MAX VF in B mesh= ',maxval_a
        write(*,fmt='(a,es20.12)') 'B mesh elts have MIN VF in A mesh= ',minval_b
        write(*,fmt='(a,es20.12)') 'B mesh elts have MAX VF in A mesh= ',maxval_b
        deallocate (vf_a)
        deallocate (vf_b)

        ! Write to disk for possible future use
        OPEN (UNIT = grid_lun, FILE='int_vols_file', IOSTAT=ier, POSITION='rewind',FORM='unformatted')
     
        if (ier.eq.0) then
           call write_int_volumes(int_vols,grid_lun)
           close(grid_lun)
        endif

     else

        write(*,fmt='(/,a,/)') '*** Reusing grid mapping'

     endif

     ! Allocate fields
     if (associated(field_a)) deallocate(field_a)
     allocate(field_a(mesh_a%nelt))
     if (associated(field_b)) deallocate(field_b)
     allocate(field_b(mesh_b%nelt))

     if (.not.reverse_order) then
        ! Test out mapping field_a to field_b with field_a=1.
        field_a=1.d0

        call map_cell_field(field_a,field_b,int_vols, &
             reverse_order,exactly_conservative,preserve_constants, &
             strict, defval)

     else
        ! Test out mapping field_b to field_a with field_b=1.
        field_b=1.d0
        call map_cell_field(field_b,field_a,int_vols, &
             reverse_order,exactly_conservative,preserve_constants, &
             strict, defval)
     endif

! Dump out mesh_a, mesh_b, and a superposition of mesh_a and mesh_b
     fname=trim(title2)//'_mesh_a.gmv'
     open (unit=gmv_lun, FILE=fname, IOSTAT=ier, POSITION='rewind')
     call writegmv_gm_mesh(gmv_lun,mesh_a,field=field_a)
     close (unit=gmv_lun)

     fname=trim(title2)//'_mesh_b.gmv'
     open (unit=gmv_lun, FILE=fname, IOSTAT=ier, POSITION='rewind')
     call writegmv_gm_mesh(gmv_lun,mesh_b,field=field_b)
     close (unit=gmv_lun)

     fname=trim(title2)//'_meshes_ab.gmv'
     open (unit=gmv_lun, FILE=fname, IOSTAT=ier, POSITION='rewind')
     call writegmv_gm_mesh(gmv_lun,mesh_a,mesh_b)
     close (unit=gmv_lun)

     minval_a=minval(field_a)
     maxval_a=maxval(field_a)
     minval_b=minval(field_b)
     maxval_b=maxval(field_b)

     write(*,fmt='(a,es20.12)') 'field_a has MINVAL= ',minval_a
     write(*,fmt='(a,es20.12)') 'field_a has MAXVAL= ',maxval_a
     write(*,fmt='(a,es20.12)') 'field_b has MINVAL= ',minval_b
     write(*,fmt='(a,es20.12)') 'field_b has MAXVAL= ',maxval_b
  
     ! Allocate Volume fields
     if (associated(vol_a)) deallocate(vol_a)
     allocate(vol_a(mesh_a%nelt))
     if (associated(vol_b)) deallocate(vol_b)
     allocate(vol_b(mesh_b%nelt))

     call grid_vols(int_vols,vol_a,vol_b)
     write(*,fmt='(/,a,es20.12)') 'Integral on A mesh= ',dot_product(vol_a,field_a)
     write(*,fmt='(a,es20.12)') 'Integral on B mesh= ',dot_product(vol_b,field_b)
  
     call destroy_grid_int_vols(int_vols)
     deallocate(mesh_a, mesh_b)
  enddo

20 continue
  close(task_lun)
  stop

end program grid_mapping_test_driver

