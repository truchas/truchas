program test_gs

  !======================================================================
  ! PURPOSE
  !   Test the Gather & Scatter routines provided by PGSLib
  !======================================================================

  ! $Id: test-gs.F,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $


  USE PGSLib_MODULE
  USE Test_GS_Parallel_MODULE
  USE Test_Globals_MODULE

  Implicit none

  integer:: ia, nproc
  integer, pointer, dimension(:,:):: ib
  real, pointer, dimension(:,:):: rb
  integer, pointer, dimension(:) :: iseed
  integer :: seed_size, e, n
  real:: ra
  integer :: n_elems, n_nodes, n_nodesperelem, el, nd
  integer :: n_elem_1, n_elem_2
  integer :: n_elems_tot, n_nodes_tot, node_vec_size
  integer :: k1, k2, item, minl(1)
  integer, pointer, dimension(:):: n_elems_list, n_nodes_list
  integer, pointer, dimension(:,:):: elements_nodes, elements_nodes_orig
  integer, pointer, dimension(:,:):: elements_nodes_item
  integer, pointer, dimension(:):: elements_numbers
  integer, pointer, dimension(:,:):: elements_nodes_tot
  integer, pointer, dimension(:,:):: elem_dest
  real(PGSLib_Double_Type),    pointer, dimension(:,:):: elem_src_db
  real,   pointer, dimension(:,:):: elem_dest_real
  logical, pointer, dimension(:,:):: mask
  integer, pointer, dimension(:)  :: node_numbers, node_dest
  real,   pointer, dimension(:)  :: node_real_src, node_dest_real
  real,   pointer, dimension(:)  :: node_dest_real_tot_given, node_dest_real_tot_expected
  real(PGSLib_DOUBLE_TYPE),   pointer, dimension(:)  :: node_src_db
  integer, pointer, dimension(:,:):: node_vec_src
  logical, pointer, dimension(:)  :: node_mask
  logical :: fatal_error, local_error
  integer :: s_degree, g_degree
  real :: bt, srt, scatter_time, gather_time, setup_time
  type (PGSLib_GS_Trace),    POINTER        :: Trace

  character (LEN=256):: s1, out_string(2)

  call pgslib_initialize(0, FILE_PREFIX='TEST_GS')
  ia = 5
  nproc = pgslib_inquire_npe()

  ! Can we deallocate a trace which has been initialized but not setup?
  WRITE(Out_String, *) ' At start of trace test'
  call pgslib_output(Out_string(1))

  trace => pgslib_gs_init_trace()
  WRITE(Out_String, *) ' Trace initialized'
  call pgslib_output(Out_string(1))
  
  call pgslib_gs_release_trace(trace)
  WRITE(Out_String, *) ' Trace released'
  call pgslib_output(Out_string(1))
    

!!$  if (PGSLib_Inquire_IO_P()) then
!!$     print *, '<Enter> to continue'
!!$     read (*,*)
!!$  end if
  call pgslib_barrier()
  ! This is for psuedo random numbers
  k1 = 1000003
  k2 = 26147

  WRITE(s1,'("NOT!  This is just a test on ",i8," processors")') nproc
  call pgslib_error(s1)
  ! Test the gs setup stuff.

  ! Setup a simple test case
  n_nodesperelem = 6
  node_vec_size = 1* n_nodesperelem ! Larger than needed 

  n_elems = MOD(MOD(k1,500)*pgslib_inquire_thispe_actual() + k2, 500)
!  n_elems = MOD(50*MOD(k1,50)*pgslib_inquire_thispe_actual() + k2, 50)
  n_elems_tot = pgslib_global_sum(n_elems)
  write(out_string,'(" n_elems = ",i8," n_elems_tot = ", i8," nPE = ",i8)') &
       &               n_elems, n_elems_tot, nproc
  call pgslib_output(out_string(1))
  allocate(n_elems_list(pgslib_inquire_npe()))
  call pgslib_collate(n_elems_list, n_elems)
  call pgslib_bcast(n_elems_list)

  if (sum(n_elems_list) /= n_elems_tot) then
     write(out_string,*) 'ERROR: SUM(n_elems_list)/= n_elems_tot = ', SUM(n_elems_list), n_elems_tot
     call pgslib_error(out_string(1))
  end if

  n_elem_2 = n_elems
  n_elem_1 = n_nodesperelem

  n_nodes = n_elem_1*n_elem_2
  n_nodes_tot = pgslib_global_sum(n_nodes)
  allocate(n_nodes_list(pgslib_inquire_npe()))
  call pgslib_collate(n_nodes_list, n_nodes)
  call pgslib_bcast(n_nodes_list)

  if (sum(n_nodes_list) /= n_nodes_tot) then
     write(out_string,*) 'ERROR: SUM(n_nodes_list)/= n_nodes_tot = ', SUM(n_nodes_list), n_nodes_tot
     call pgslib_error(out_string(1))
  end if

  allocate(elements_nodes(n_elem_1, n_elem_2))
  allocate(elements_nodes_orig(n_elem_1, n_elem_2))

  allocate(elements_numbers(n_elem_2))

  allocate(elem_dest(n_elem_1, n_elem_2))
  allocate(elem_dest_real(n_elem_1, n_elem_2))

  allocate(mask(n_elem_1, n_elem_2))

  allocate(node_numbers(n_nodes))
  allocate(node_vec_src(node_vec_size, n_nodes))
  allocate(node_src_db(n_nodes))


  node_numbers = (/ (nd, nd=1,n_nodes) /)
  node_numbers = node_numbers + SUM(n_nodes_list(1: pgslib_inquire_thispe_actual()-1))

  elements_numbers = (/ (el, el=1,n_elem_2) /)
  elements_numbers = elements_numbers + SUM(n_elems_list(1: pgslib_inquire_thispe_actual()-1))


  CALL RANDOM_SEED()
  ! Compute a distinct seed for each processor
  CALL RANDOM_SEED(SIZE=seed_size)
  allocate(rb(seed_size,nproc))
  allocate(ib(seed_size,nproc))
  allocate(iseed(seed_size))
  call RANDOM_NUMBER(rb)
  ib = 2**29 * rb
  do e= 1 , SIZE(rb,1)
    call pgslib_dist(iseed(e), ib(e,:))
  end do
!  call RANDOM_SEED(PUT=iseed)
  call RANDOM_NUMBER(HARVEST = elem_dest_real)
  elements_nodes = elem_dest_real*REAL(n_nodes_tot)
  elements_nodes = MOD(elements_nodes, n_nodes_tot) 
  elements_nodes = MERGE(elements_nodes, elements_nodes + n_nodes_tot, elements_nodes >= 1)

  if (any(elements_nodes< 1) ) then
     call pgslib_error('BIG ERROR, SOME ELEMENTS_NODES < 1')
  end if

  if (any(elements_nodes> n_nodes_tot) ) then
     call pgslib_error('BIG ERROR, SOME ELEMENTS_NODES > n_nodes_tot')
  end if

  elements_nodes_orig = elements_nodes

  ! Okay to initialize again.
  call PGSLib_Initialize(0)

  !======================================================================
  !
  !========== Gather Tests ==============================================

  !======================================================================
  ! Test gather and scatter with no mask
  call Test_EN_GS_Setup(elements_nodes, N_Nodes)
  if (any(elements_nodes <0 )) then
     write(out_string, *) ' Some element"s nodes are off-pe.'
     call pgslib_output(out_string(1))
  ENDIF

  call pgslib_output(' Passed TEST_EN_GS_Setup')

  node_src_db = 0.0
  call RANDOM_NUMBER(HARVEST = node_src_db)
      
  call test_gather(node_src_db, elements_nodes, elements_nodes_orig)

  !!!!!!!!!! Scatter Tests !!!!!!!!!!

  allocate(elem_src_db(n_elem_1, n_elem_2))
  call RANDOM_NUMBER(HARVEST = elem_src_db)

  allocate(mask(n_elem_1, n_elem_2))
  MASK = (elements_nodes /= 1)
  call Test_Scatter_Sum(elem_src_db, n_nodes, elements_nodes, elements_nodes_orig, MASK)
  call Test_Scatter_Max(elem_src_db, n_nodes, elements_nodes, elements_nodes_orig, MASK)

  ! First test w/o temporary mask
  call Test_Scatter_Sum(elem_src_db, n_nodes, elements_nodes, elements_nodes_orig)
  call pgslib_output('back from test_scatter_sum')
  call Test_Scatter_Max(elem_src_db, n_nodes, elements_nodes, elements_nodes_orig)
  call pgslib_output('back from test_scatter_max')


  !======================================================================
  ! Test gather and scatter with a temporary requestor mask only.

  !!      call TEST_EN_GS_Setup(elements_nodes, elements_nodes_pe, elements_nodes_pe_flag, node_numbers)
  if (any(elements_nodes < 0)) then
     write(out_string, *) ' Some element"s nodes are off-pe.'
     call pgslib_output(out_string(1))
  ENDIF

  allocate(mask(n_elem_1, n_elem_2))
  MASK = (elements_nodes /= 1)

!!$  elem_dest = 0
!!$  call TEST_EN_Gather(elem_dest, node_numbers, elements_nodes, MASK=mask)
!!$
!!$  IF (ANY((elem_dest .NE. elements_nodes_orig) .AND. mask)) THEN
!!$     write(out_string, *) ' Some elem_dest /= elements_nodes'
!!$     call pgslib_output(out_string(1))
!!$  ELSE
!!$     write(out_string, *) ' REJOICE: Gather passed a temporary requestor mask test.'   
!!$     call pgslib_output(out_string(1))
!!$  ENDIF

  call Test_Scatter_Sum(elem_src_db, n_nodes, elements_nodes, elements_nodes_orig, MASK)
  call Test_Scatter_Max(elem_src_db, n_nodes, elements_nodes, elements_nodes_orig, MASK)

  call PGSLib_Trace_Degree(s_degree, g_degree, EN_Trace)

  call pgslib_gs_release_trace(EN_Trace)

  !======================================================================
  ! Test gather and scatter with a permenant mask.

  elements_nodes = elements_nodes_orig
  call TEST_EN_GS_Setup(elements_nodes, N_nodes, MASK=MASK)
  if (any(elements_nodes <0 )) then
     write(out_string, *) ' Some element"s nodes are off-pe.'
     call pgslib_output(out_string(1))
  ENDIF

  call pgslib_output(' Passed TEST_EN_GS_Setup with MASK')

  call Test_Scatter_Sum(elem_src_db, n_nodes, elements_nodes, elements_nodes_orig, MASK)
  call Test_Scatter_Max(elem_src_db, n_nodes, elements_nodes, elements_nodes_orig, MASK)




  deallocate(elem_src_db)


  ! Print some timing information
  setup_time = PGSLib_Read_Maximum_Time(Setup_Trace_STATISTICS())
  gather_time = PGSLib_Read_Maximum_Time(GATHER_BUFFER_STATISTICS())
  scatter_time = PGSLib_Read_Maximum_Time(SCATTER_BUFFER_STATISTICS())
  bt  = PGSLib_Barrier_Time()
  bt  = pgslib_global_maxval(bt)
  srt = PGSLib_SR_Time()
  srt = pgslib_global_maxval(srt)
  s_degree = pgslib_global_maxval(s_degree)
  g_degree = pgslib_global_maxval(g_degree)
  if (PGSLib_Inquire_IO_P()) then
     PRINT *, 'Time for Setup', setup_time
     PRINT *, 'Time for Gather Buffers', gather_time 
     PRINT *, 'Time for Scatter Buffers', scatter_time 
     print *, 'Graph degree, scatter, gather:', s_degree, g_degree
     PRINT *, 'Time in MPI_Barrier for gs', bt
     PRINT *, 'Time in send-rcv core for gs', srt
  end if

  call pgslib_gs_release_trace(EN_Trace)
  Call pgslib_finalize()
  stop

!!$  node_dest_real = 0
!!$  elem_src_real = 1.0
!!$  call TEST_EN_Scatter_SUM(node_dest_real, elem_src_real, elements_nodes, MASK=mask)
!!$
!!$  ! Need the node mask to check the result.
!!$  allocate(node_mask(n_nodes))
!!$  node_mask = .FALSE.
!!$  call TEST_EN_Scatter_OR(node_mask, mask, elements_nodes, MASK=mask)
!!$
!!$  IF (ANY( (node_dest .NE. ((n_elem_2* n_elem_1)/ n_nodes)) .AND. Node_Mask) ) THEN
!!$     write(out_string, *) ' Some node_dest /= ', (n_elem_2* n_elem_1)/n_nodes
!!$     call pgslib_output(out_string(1))
!!$  ELSE
!!$     write(out_string, *) ' REJOICE: Scatter passed a temporary requestor mask test.'   
!!$     call pgslib_output(out_string(1))
!!$  ENDIF

!$VECTOR$!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$VECTOR$!  !
!$VECTOR$!  ! Test vector gather with no mask
!$VECTOR$!  !
!$VECTOR$!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$VECTOR$!
!$VECTOR$!  do item = 1, node_vec_size
!$VECTOR$!     node_vec_src(item, :) = item + node_numbers
!$VECTOR$!  end do
!$VECTOR$!
!$VECTOR$!  allocate(elements_nodes_item(n_elem_1, n_elem_2))
!$VECTOR$!  do item = 1, n_elem_1
!$VECTOR$!     elements_nodes_item(item,:) = item
!$VECTOR$!  enddo
!$VECTOR$!
!$VECTOR$!  elem_dest = 0
!$VECTOR$!
!$VECTOR$!  call TEST_EN_Gather(elem_dest, node_vec_src, elements_nodes, elements_nodes_item)
!$VECTOR$!  
!$VECTOR$!  IF (ANY((elem_dest) .NE. (elements_nodes_orig + elements_nodes_item) ) ) THEN
!$VECTOR$!     write(out_string, *) 'After TEST_EN_Gather (vector)'
!$VECTOR$!     call pgslib_error(out_string(1))
!$VECTOR$!     write(out_string, *) 'Some elem_dest /= elements_nodes_item'
!$VECTOR$!     call pgslib_error(out_string(1))
!$VECTOR$!  ELSE
!$VECTOR$!     write(out_string, *) ' REJOICE: Gather vector passed a simple test.'   
!$VECTOR$!     call pgslib_output(out_string(1))
!$VECTOR$!  ENDIF
!$VECTOR$!
!$VECTOR$!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test gather and scatter with a temporary requestor mask and a temporary node mask.

!!$      call TEST_EN_GS_Setup(elements_nodes, elements_nodes_pe, node_numbers)
!!$      if (any(.not. elements_nodes_pe_flag)) then
!!$         write(out_string, *) ' Some elements_nodes_pe_flag == .FALSE.'
!!$         call pgslib_output(out_string(1))
!!$      ENDIF
!!$
!!$      node_mask = .FALSE.
!!$      call TEST_EN_Scatter_OR(node_mask, mask, elements_nodes, elements_nodes_pe_flag, MASK=mask)
!!$
!!$      elem_dest = 0
!!$      call TEST_EN_Gather(elem_dest, node_numbers, elements_nodes, elements_nodes_pe_flag, MASK=mask, SRC_MASK = node_mask)
!!$      IF (ANY(((elem_dest+1) .NE. elements_nodes) .AND. mask)) THEN
!!$         write(out_string, *) ' Some elem_dest /= elements_nodes'
!!$         call pgslib_output(out_string(1))
!!$      ELSE
!!$         write(out_string, *) ' REJOICE: Gather passed a temporary requestor mask and temporary node mask test.'   
!!$         call pgslib_output(out_string(1))
!!$      ENDIF
!!$
!!$      node_dest = 0
!!$      elem_src = 1
!!$      call TEST_EN_Scatter_SUM(node_dest, elem_src, elements_nodes, elements_nodes_pe_flag, MASK=mask, DEST_MASK = node_mask)
!!$      
!!$      IF (ANY( (node_dest .NE. ((n_elem_2* n_elem_1)/ n_nodes)) .AND. Node_Mask) ) THEN
!!$         write(out_string, *) ' Some node_dest /= ', (n_elem_2* n_elem_1)/n_nodes
!!$         call pgslib_output(out_string(1))
!!$      ELSE
!!$         write(out_string, *) ' REJOICE: Scatter passed a temporary requestor mask and temporary node mask test.'   
!!$         call pgslib_output(out_string(1))
!!$      ENDIF
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$! Test gather and scatter with a permenant requestor mask only.
!!$
!!$      elements_nodes_pe = -1
!!$      elements_nodes_pe_flag = .false.
!!$
!!$      MASK = (elements_nodes /= 1)
!!$      call TEST_EN_GS_Setup(elements_nodes, elements_nodes_pe, elements_nodes_pe_flag, node_numbers, MASK=mask)
!!$      if (any(.not. elements_nodes_pe_flag)) then
!!$         write(out_string, *) ' Some elements_nodes_pe_flag == .FALSE.'
!!$         call pgslib_output(out_string(1))
!!$      ENDIF
!!$
!!$      elem_dest = 0
!!$      call TEST_EN_Gather(elem_dest, node_numbers, elements_nodes, elements_nodes_pe_flag, MASK=mask)
!!$      IF (ANY(((elem_dest+1) .NE. elements_nodes) .AND. mask)) THEN
!!$         write(out_string, *) ' Some elem_dest /= elements_nodes'
!!$         call pgslib_output(out_string(1))
!!$      ELSE
!!$         write(out_string, *) ' REJOICE: Gather passed a permenant requestor mask test.'   
!!$         call pgslib_output(out_string(1))
!!$      ENDIF
!!$
!!$      node_dest = 0
!!$      elem_src = 1
!!$      call TEST_EN_Scatter_SUM(node_dest, elem_src, elements_nodes, elements_nodes_pe_flag, MASK=mask)
!!$      
!!$      ! Need the node mask to check the result.
!!$      node_mask = .FALSE.
!!$      call TEST_EN_Scatter_OR(node_mask, mask, elements_nodes, elements_nodes_pe_flag, MASK=mask)
!!$      
!!$      IF (ANY( (node_dest .NE. ((n_elem_2* n_elem_1)/ n_nodes) ) .AND. node_mask) )THEN
!!$         write(out_string, *) ' Some node_dest /= ', (n_elem_2* n_elem_1)/n_nodes
!!$         call pgslib_output(out_string(1))
!!$      ELSE
!!$         write(out_string, *) ' REJOICE: Scatter passed a permenant requestor mask test.'   
!!$         call pgslib_output(out_string(1))
!!$      ENDIF
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$! Test gather and scatter with a permenant requestor and temporary node mask
!!$
!!$      elements_nodes_pe = -1
!!$      elements_nodes_pe_flag = .false.
!!$      MASK = (elements_nodes /= 1)
!!$      call TEST_EN_GS_Setup(elements_nodes, elements_nodes_pe, elements_nodes_pe_flag, node_numbers, MASK=mask)
!!$      if (any(.not. elements_nodes_pe_flag)) then
!!$         write(out_string, *) ' Some elements_nodes_pe_flag == .FALSE.'
!!$         call pgslib_output(out_string(1))
!!$      ENDIF
!!$
!!$      node_mask = .FALSE.
!!$      call TEST_EN_Scatter_OR(node_mask, mask, elements_nodes, elements_nodes_pe_flag, MASK=mask)
!!$      
!!$      elem_dest = 0
!!$      call TEST_EN_Gather(elem_dest, node_numbers, elements_nodes, elements_nodes_pe_flag, MASK=mask, SRC_MASK=node_mask)
!!$      IF (ANY(((elem_dest+1) .NE. elements_nodes) .AND. mask)) THEN
!!$         write(out_string, *) ' Some elem_dest /= elements_nodes'
!!$         call pgslib_output(out_string(1))
!!$      ELSE
!!$         write(out_string, *) ' REJOICE: Gather passed a permenant request mask and temporary node mask test.'  
!!$         call pgslib_output(out_string(1))
!!$      ENDIF
!!$
!!$      node_dest = 0
!!$      elem_src = 1
!!$      call TEST_EN_Scatter_SUM(node_dest, elem_src, elements_nodes, elements_nodes_pe_flag, MASK=mask, DEST_MASK=node_mask)
!!$      
!!$      IF (ANY( (node_dest .NE. ((n_elem_2* n_elem_1)/ n_nodes) ) .AND. node_mask) )THEN
!!$         write(out_string, *) ' Some node_dest /= ', (n_elem_2* n_elem_1)/n_nodes
!!$         call pgslib_output(out_string(1))
!!$      ELSE
!!$         write(out_string, *) ' REJOICE: Scatter passed a  permenant request mask and temporary node mask test.'  
!!$         call pgslib_output(out_string(1))
!!$      ENDIF
!!$

end
