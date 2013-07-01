subroutine MAPFIELD(ncellsa, nvca, connectivitya, &
                    nnodesa, coordsa, cblocksa,   &
                    ncellsb, nvcb, connectivityb, & 
                    nnodesb, coordsb, cblocksb,   &
                    ndim, maxwarn,                &
                    maprule,                      &                 
                    field_a, field_b, intgrl_a, intgrl_b, defval, ierr) 

  !map a cell fields from mesh 'a' to mesh 'b'

  use grid_mapping_module
  implicit none

  integer, parameter                                            :: double_kind = KIND(1.0d0)
  integer,                                        INTENT(IN)    :: ncellsa, ncellsb, nnodesa, nnodesb, ndim, nvca, nvcb 
  integer, dimension(nvca,ncellsa),               INTENT(IN)    :: connectivitya
  integer, dimension(nvcb,ncellsb),               INTENT(IN)    :: connectivityb
  integer, dimension(ncellsa),                    INTENT(IN)    :: cblocksa
  integer, dimension(ncellsb),                    INTENT(IN)    :: cblocksb
  real(kind=double_kind), dimension(ndim,nnodesa), INTENT(IN)   :: coordsa
  real(kind=double_kind), dimension(ndim,nnodesb), INTENT(IN)   :: coordsb
  real(kind=double_kind), dimension(ncellsa),     INTENT(IN)    :: field_a
  real(kind=double_kind), dimension(ncellsb),     INTENT(INOUT) :: field_b
  real(kind=double_kind),                         INTENT(OUT)   :: intgrl_a, intgrl_b
  real(kind=double_kind),                         INTENT(IN)    :: defval
  integer,                                        INTENT(OUT)   :: ierr
  integer,                                        INTENT(IN)    :: maxwarn
  character(len=*),                               INTENT(IN)    :: maprule
  integer                                                       :: status
  type(gm_mesh)                                                 :: mesh_a, mesh_b
  logical                                                       :: map_exists,reverse_order, &
                                                                   exactly_conservative,preserve_constants,strict
  type(grid_int_vols)                                           :: int_vols
  real(kind=double_kind)                                        :: minval_a, maxval_a, minval_b, maxval_b
  real(kind=double_kind), dimension(:), pointer                 :: vol_a=>null(),vf_a=>null(),vol_b=>null(),vf_b=>null()
  integer, parameter                                            :: task_lun=8, grid_lun=9

  exactly_conservative = .false.
  preserve_constants   = .false.

  if (maprule == 'conservative') then
     exactly_conservative = .true.
     preserve_constants   = .false.
  end if

  if (maprule == 'constants_preserving') then
     exactly_conservative = .false.
     preserve_constants   = .true.
  end if

  intgrl_a       = 0.0d0
  intgrl_b       = 0.0d0

  reverse_order  = .false.

  !need to swap the connectivity and coords arrays around to match the mapping interface

  !connectivityta = TRANSPOSE(connectivitya)
  !connectivitytb = TRANSPOSE(connectivityb)
  !coordsta       = TRANSPOSE(coordsa)
  !coordstb       = TRANSPOSE(coordsb)

  !create gm mesh 'a' and gm mesh 'b'

  !first create gm mesh 'a' typically from truchas XML output
  mesh_a%nnod     = nnodesa
  mesh_a%nelt     = ncellsa
  allocate(mesh_a%pos_node(ndim,mesh_a%nnod))
  mesh_a%pos_node = coordsa
  allocate(mesh_a%node_elt(nvca,mesh_a%nelt))
  allocate(mesh_a%block_elt(mesh_a%nelt))
  mesh_a%node_elt  = connectivitya
  mesh_a%block_elt = cblocksa

  !now create gm mesh 'b' typically from an exodus file
  mesh_b%nnod     = nnodesb
  mesh_b%nelt     = ncellsb
  allocate(mesh_b%pos_node(ndim,mesh_b%nnod))
  mesh_b%pos_node = coordsb
  allocate(mesh_b%node_elt(nvcb,mesh_b%nelt))
  allocate(mesh_b%block_elt(mesh_b%nelt))
  mesh_b%node_elt  = connectivityb
  mesh_b%block_elt = cblocksb

  ! Does the right mapping file already exist on disk?
  map_exists = .false.

  OPEN (UNIT = grid_lun, FILE='int_vols_file', STATUS='old', IOSTAT=ierr, POSITION='rewind', &
       FORM = 'unformatted')

  OPEN (UNIT = 20, FILE = 'diagnostics_file', FORM = 'formatted')
     
  if (ierr.eq.0) then
        
     call read_int_volumes(int_vols, grid_lun, ierr)

     close(grid_lun)
     if (ierr.eq.0) then
        map_exists = right_int_volumes(mesh_a, mesh_b, int_vols)
     endif
     
  endif

  if (map_exists) then
     write(20,fmt='(/,a,/)') '*** Recomputing grid mapping ***'
     !print '(a)', '         *** Using previous grid mapping matrix *** '
  end if

  if (.not.map_exists) then

     write(20,fmt='(/,a,/)') '*** Recomputing grid mapping ***'
     !print '(a)', '         *** Recomputing grid mapping ***'
     
     call compute_int_volumes(mesh_a,mesh_b,int_vols, maxwarn, ier=ierr)

     if (ierr.lt.0) then
        return
     end if

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

     write(20,fmt='(a,es20.12)') 'A mesh has MIN VOLUME= ',minval_a
     write(20,fmt='(a,es20.12)') 'A mesh has MAX VOLUME= ',maxval_a
     write(20,fmt='(a,es20.12)') 'B mesh has MIN VOLUME= ',minval_b
     write(20,fmt='(a,es20.12)') 'B mesh has MAX VOLUME= ',maxval_b

     deallocate (vol_a)
     deallocate (vol_b)

     if (associated(vf_a)) deallocate(vf_a)
     allocate(vf_a(mesh_a%nelt))
     if (associated(vf_b)) deallocate(vf_b)
     allocate(vf_b(mesh_b%nelt))

     call grid_vol_fracs(int_vols,vf_a,vf_b,strict=.true.)
     minval_a=minval(vf_a)
     maxval_a=maxval(vf_a)
     minval_b=minval(vf_b)
     maxval_b=maxval(vf_b)

     write(20,fmt='(/,a)')         'Disregarding block IDs, we have that'
     write(20,fmt='(a,es20.12)') 'A mesh elts have MIN VF in B mesh= ',minval_a
     write(20,fmt='(a,es20.12)') 'A mesh elts have MAX VF in B mesh= ',maxval_a
     write(20,fmt='(a,es20.12)') 'B mesh elts have MIN VF in A mesh= ',minval_b
     write(20,fmt='(a,es20.12)') 'B mesh elts have MAX VF in A mesh= ',maxval_b

     call grid_vol_fracs(int_vols,vf_a,vf_b,strict=.true.)
     minval_a=minval(vf_a)
     maxval_a=maxval(vf_a)
     minval_b=minval(vf_b)
     maxval_b=maxval(vf_b)

     write(20,fmt='(/,a)')         'Respecting block IDs, we have that'
     write(20,fmt='(a,es20.12)') 'A mesh elts have MIN VF in B mesh= ',minval_a
     write(20,fmt='(a,es20.12)') 'A mesh elts have MAX VF in B mesh= ',maxval_a
     write(20,fmt='(a,es20.12)') 'B mesh elts have MIN VF in A mesh= ',minval_b
     write(20,fmt='(a,es20.12)') 'B mesh elts have MAX VF in A mesh= ',maxval_b
     deallocate (vf_a)
     deallocate (vf_b)

     ! Write to disk for possible future use
     OPEN (UNIT = grid_lun, FILE='int_vols_file', IOSTAT=ierr, POSITION='rewind',FORM='unformatted')
     
     if (ierr.eq.0) then
        call write_int_volumes(int_vols,grid_lun)
        close(grid_lun)
     endif

  else

     write(20,fmt='(/,a,/)') '*** Reusing grid mapping'

  endif
  
  strict = .true.
  call map_cell_field(field_a,field_b,int_vols, &
       reverse_order,exactly_conservative,preserve_constants, &
       strict, defval, ier=ierr)

  minval_a=minval(field_a)
  maxval_a=maxval(field_a)
  minval_b=minval(field_b)
  maxval_b=maxval(field_b)

  write(20,fmt='(a,es20.12)') 'field_a has MINVAL= ',minval_a
  write(20,fmt='(a,es20.12)') 'field_a has MAXVAL= ',maxval_a
  write(20,fmt='(a,es20.12)') 'field_b has MINVAL= ',minval_b
  write(20,fmt='(a,es20.12)') 'field_b has MAXVAL= ',maxval_b
  
  ! Allocate Volume fields
  if (associated(vol_a)) deallocate(vol_a)
  allocate(vol_a(mesh_a%nelt))
  if (associated(vol_b)) deallocate(vol_b)
  allocate(vol_b(mesh_b%nelt))
  
  call grid_vols(int_vols,vol_a,vol_b)

  intgrl_a = dot_product(vol_a,field_a)
  intgrl_b = dot_product(vol_b,field_b)

  write(20,fmt='(/,a,es20.12)') 'Integral on A mesh= ',dot_product(vol_a,field_a)
  write(20,fmt='(a,es20.12)') 'Integral on B mesh= ',dot_product(vol_b,field_b)

  call destroy_grid_int_vols(int_vols)
  !close diagnostics file
  close(20)

  return
end subroutine MAPFIELD


