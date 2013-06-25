!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Test_Redux(error)
  ! Complete test of reduction routines in PGSLib

  USE PGSLib_MODULE
  implicit none
  logical(PGSLib_Log_Type) :: error
  
  ! local
  logical(PGSLib_Log_Type) :: local_error

  local_error = .false.
  call test_redux_int(local_error)
  error = error .or. local_error

  call test_redux_real(local_error)
  error = error .or. local_error

  call test_redux_double(local_error)
  error = error .or. local_error

  call test_redux_log(local_error)
  ! test log
  ! test dot_product
  return
end subroutine Test_Redux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test reduction of integers
subroutine test_redux_int(error)
  USE PGSLib_MODULE
  USE error_tests
  implicit none
  logical(PGSLib_Log_Type):: error

  !Local
  logical(PGSLib_Log_Type):: local_error

  integer(pgslib_int_type) :: int_given, int_expected, int_src, src_len
  integer(pgslib_int_type), dimension(1) :: int_v_given, int_v_expected
  integer(pgslib_int_type), pointer, dimension(:) :: list_pes, int_vec1d_src, src_lengths
  integer(pgslib_int_type), pointer, dimension(:) :: int_vec1d_src_b
  integer(pgslib_int_type), pointer, dimension(:) :: int_global_vec1d_src
  integer(pgslib_int_type), pointer, dimension(:) :: int_vec1d_given
  logical(pgslib_log_type), pointer, dimension(:) :: vec1d_mask
  integer, pointer, dimension(:,:) :: int_vec2d_src
  logical, pointer, dimension(:,:) :: vec2d_mask
  integer(pgslib_int_type) :: n, npe, thispe, pe, offset, length, int_temp
  logical(PGSLib_Log_Type) :: src_mask

  integer(pgslib_int_type) :: k1, k2, i, j
  CHARACTER (LEN=2048) :: out_string

  k1 = 1000003
  k2 = 262147

  npe = PGSLib_Inquire_nPE()
  thispe = PGSLib_Inquire_thisPE_Actual()
  local_error = .false.

  ! Setup for permutations
!  k1 = MOD(k1, npe)
!  k2 = MOD(k2, npe)
  
  allocate(list_pes (npe))
  list_pes = (/ (pe, pe=1,npe) /)  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! Scalar integer min/max !!!!!!!!!!
  ! Test min scalar, no mask
  ! Global, implicitly
  int_src = MOD(thispe*k1 + k2, npe)
  int_expected = MINVAL( MOD(list_pes*k1+k2, npe) )
  
  call pgslib_output("Calling int global minval")
  call PGSLib_Flush_Output()
  int_given = pgslib_global_minval(int_src)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .TRUE.
     write(out_string, *) 'FAILED integer scalar no mask minval reduction test:'
     call pgslib_error(out_string)
     call pgslib_output(out_string)
     write(out_string, 10) int_src, int_given, int_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif
10   FORMAT(1x, ' int_src = ', i8, ' int_given = ', i8, ' int_expected = ', i8)

  ! Global, explicitly
  int_src = MOD(thispe*k1 + k2, npe)
  int_expected = MINVAL( MOD(list_pes*k1+k2, npe) )
  
  int_given = pgslib_global_minval(int_src, SCOPE=PGSLib_Global)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .TRUE.
     write(out_string, *) 'FAILED GLOBAL integer scalar no mask minval reduction test:'
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, 10) int_src, int_given, int_expected
     call pgslib_error (out_string)
     call pgslib_output(out_string)
  endif

  ! Local
  int_src = MOD(thispe*k1 + k2, npe)
  int_expected = MINVAL((/int_src/))
  
  int_given = pgslib_global_minval(int_src, SCOPE=PGSLib_Local)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .TRUE.
     write(out_string, *) int_src, int_given, int_expected
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, 10) int_src, int_given, int_expected
     call pgslib_error (out_string)
     call pgslib_output(out_string)
  endif


  ! Test min scalar, with mask
  ! Global, implicitly
  int_src = MOD(thispe*k1 + k2, npe)
  src_mask = MOD(thispe,3) == 0
  int_expected = MINVAL( MOD(list_pes*k1+k2, npe) , MASK=(MOD(list_pes,3) == 0) )

  
  int_given = pgslib_global_minval(int_src, MASK=src_mask)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .TRUE.
     write(out_string, *) 'FAILED integer scalar with MASK MINVAL reduction test:'
     call pgslib_error(out_string)
     call pgslib_output(out_string)
     write(out_string, 20) int_src, int_given, int_expected, src_mask
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif
20   FORMAT(1x,'    int_src = ', i8, ' int_given = ', i8, ' int_expected = ', i8, ' src_mask = ', L4)

  ! Global, explicitly
  int_src = MOD(thispe*k1 + k2, npe)
  src_mask = MOD(thispe,3) == 0
  int_expected = MINVAL( MOD(list_pes*k1+k2, npe) , MASK=(MOD(list_pes,3) == 0) )

  
  int_given = pgslib_global_minval(int_src, MASK=src_mask, SCOPE=PGSLib_GLOBAL)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .TRUE.
     write(out_string, *) 'FAILED GLOBAL integer scalar WITH MASK minval reduction test:'
     call pgslib_error(out_string)
     call pgslib_output(out_string)
     write(out_string, 20) int_src, int_given, int_expected, src_mask
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif

  
  ! Local, explicitly
  int_src = MOD(thispe*k1 + k2, npe)
  src_mask = MOD(thispe,3) == 0
  int_expected = MINVAL( (/ int_src /) , MASK= (/ src_mask /) )

  int_given = pgslib_global_minval(int_src, MASK=src_mask, SCOPE=PGSLib_LOCAL)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .true.
     call pgslib_error("Failed LOCAL integer scalar with mask minval reduction test")
     call pgslib_output("Failed LOCAL integer scalar with mask minval reduction test")
  endif

  
  ! Test max scalar no mask
  ! Global, implicitly
  int_src = MOD(thispe*k1 + k2, npe)
  int_expected = MAXVAL( MOD(list_pes*k1+k2, npe) )
  
  int_given = pgslib_global_maxval(int_src)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .true.
     call pgslib_error("Failed integer scalar no mask maxval reduction test")
     call pgslib_output("Failed integer scalar no mask maxval reduction test")
  endif

  ! Test max scalar with mask
  ! Global, explicitly
  int_src = MOD(thispe*k1 + k2, npe)
  int_expected = MAXVAL( MOD(list_pes*k1+k2, npe) )
  
  int_given = pgslib_global_maxval(int_src, SCOPE=PGSLib_Global)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .true.
     call pgslib_error("Failed GLOBAL integer scalar no mask maxval reduction test")
     call pgslib_output("Failed GLOBAL integer scalar no mask maxval reduction test")
  endif

  ! Local, explicitly
  int_src = MOD(thispe*k1 + k2, npe)
  int_expected = MAXVAL( (/ int_src /) )
  
  int_given = pgslib_global_maxval(int_src, SCOPE=PGSLib_Local)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .true.
     call pgslib_error("Failed LOCAL integer scalar no mask maxval reduction test")
     call pgslib_output("Failed LOCAL integer scalar no mask maxval reduction test")
  endif

  ! Test max scalar with mask
  ! Global, implicitly
  int_src = MOD(thispe*k1 + k2, npe)
  src_mask = MOD(thispe, 3) == 0
  int_expected = MAXVAL( MOD(list_pes*k1+k2, npe), MASK=(MOD(list_pes,3) == 0) )
  
  int_given = pgslib_global_maxval(int_src, MASK=src_mask)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .true.
     call pgslib_error("Failed integer scalar with mask maxval reduction test")
     call pgslib_output("Failed integer scalar with mask maxval reduction test")
  endif

  ! Global, explicitly
  int_src = MOD(thispe*k1 + k2, npe)
  src_mask = MOD(thispe, 3) == 0
  int_expected = MAXVAL( MOD(list_pes*k1+k2, npe), MASK=(MOD(list_pes,3) == 0) )
  
  int_given = pgslib_global_maxval(int_src, MASK=src_mask, SCOPE=PGSLib_GLOBAL)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .true.
     call pgslib_error("Failed GLOBAL integer scalar with mask maxval reduction test")
     call pgslib_output("Failed GLOBAL integer scalar with mask maxval reduction test")
  endif

  ! Local, explicitly
  int_src = MOD(thispe*k1 + k2, npe)
  src_mask = MOD(thispe, 3) == 0
  int_expected = MAXVAL( (/ int_src /), MASK=(/ src_mask /) )
  
  int_given = pgslib_global_maxval(int_src, MASK=src_mask, SCOPE=PGSLib_LOCAL)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .true.
     call pgslib_error("Failed LOCAL integer scalar with mask maxval reduction test")
     call pgslib_output("Failed LOCAL integer scalar with mask maxval reduction test")
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! vector 1d integer min/max
  src_len = MOD(2*thispe, 3) + 1
  allocate(int_vec1d_src(src_len))
  
  int_vec1d_src = (/ (n, n=1,src_len) /)
  int_vec1d_src = MOD(int_vec1d_src*k1 + k2, npe)

  ! Test min 1d vector no mask

  int_expected = PGSLib_Global_Minval( MINVAL(int_vec1d_src) )
  int_given    = PGSLib_Global_Minval( int_vec1d_src )

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .true.
     call pgslib_error("Failed integer 1d vector no mask minval reduction test")
     call pgslib_output("Failed integer 1d vector no mask minval reduction test")
  endif

  ! Test max 1d vector no mask

  int_expected = PGSLib_Global_Maxval( MAXVAL(int_vec1d_src) )
  int_given    = PGSLib_Global_Maxval( int_vec1d_src )

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .true.
     call pgslib_error("Failed integer 1d vector no mask maxval reduction test")
     call pgslib_output("Failed integer 1d vector no mask maxval reduction test")
  endif

  ! Test min 1d vector with mask

  allocate(vec1d_mask(src_len))
  vec1d_mask = MOD(int_vec1d_src,2) == 0

  int_expected = PGSLib_Global_Minval( MINVAL(int_vec1d_src, MASK=vec1d_mask) )
  int_given    = PGSLib_Global_Minval( int_vec1d_src, MASK=vec1d_mask)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .true.
     call pgslib_error("Failed integer 1d vector with mask minval reduction test")
     call pgslib_output("Failed integer 1d vector with mask minval reduction test")
  endif


  ! Test max 1d vector with mask
  vec1d_mask = MOD(int_vec1d_src,2) == 0

  int_expected = PGSLib_Global_Maxval( MAXVAL(int_vec1d_src, MASK=vec1d_mask) )
  int_given    = PGSLib_Global_Maxval( int_vec1d_src, MASK=vec1d_mask)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .true.
     call pgslib_error("Failed integer 1d vector with mask maxval reduction test")
     call pgslib_output("Failed integer 1d vector with mask maxval reduction test")
  endif


  deallocate(int_vec1d_src)
  deallocate(vec1d_mask)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! vector 2d integer min/max
  src_len = MOD(2*thispe, 3) + 1
  allocate(int_vec2d_src(src_len, 2*npe))
  
  do i=1,src_len
     do j = 1,2*npe
        int_vec2d_src(i,j) = i+j
     enddo
  enddo
  
  int_vec2d_src = MOD((int_vec2d_src+thispe)*k1 + k2, npe)

  ! Test min 2d vector no mask

  int_expected = PGSLib_Global_Minval( MINVAL(int_vec2d_src) )
  int_given    = PGSLib_Global_Minval( int_vec2d_src )

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .true.
     call pgslib_error("Failed integer 2d vector no mask minval reduction test")
     call pgslib_output("Failed integer 2d vector no mask minval reduction test")
  endif

  ! Test max 2d vector no mask

  int_expected = PGSLib_Global_Maxval( MAXVAL(int_vec2d_src) )
  int_given    = PGSLib_Global_Maxval( int_vec2d_src )

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .true.
     call pgslib_error("Failed integer 2d vector no mask maxval reduction test")
     call pgslib_output("Failed integer 2d vector no mask maxval reduction test")
  endif

  ! Test min 2d vector with mask

  allocate(vec2d_mask(src_len, 2*npe))
  vec2d_mask = MOD(int_vec2d_src,2) == 0

  int_expected = PGSLib_Global_Minval( MINVAL(int_vec2d_src, MASK=vec2d_mask) )
  int_given    = PGSLib_Global_Minval( int_vec2d_src, MASK=vec2d_mask)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .true.
     call pgslib_error("Failed integer 2d vector with mask minval reduction test")
     call pgslib_output("Failed integer 2d vector with mask minval reduction test")
  endif

  ! Test max 2d vector with mask

  vec2d_mask = MOD(int_vec2d_src,2) == 0

  int_expected = PGSLib_Global_Maxval( MAXVAL(int_vec2d_src, MASK=vec2d_mask) )
  int_given    = PGSLib_Global_Maxval( int_vec2d_src, MASK=vec2d_mask)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .true.
     call pgslib_error("Failed integer 2d vector with mask maxval reduction test")
     call pgslib_output("Failed integer 2d vector with mask maxval reduction test")
  endif

  deallocate(int_vec2d_src)
  deallocate(vec2d_mask)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!! TEST SUM !!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test sum scalar no mask integer
  int_src = thispe
  int_expected = SUM(list_pes)
  int_given = pgslib_global_sum( int_src )
  
  local_error = local_error .OR. (int_given .ne. int_expected)
  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     write(out_string, *) "FAILED integer SCALAR NO MASK sum reduction test"
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, 10) int_src, int_given, int_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif
  
  ! Test sum scalar with mask integer
  int_src = thispe
  int_expected = SUM(list_pes, MASK=(MOD(list_pes,2) == 0))
  
  int_given = pgslib_global_sum( int_src , MASK=(MOD(thispe,2) == 0))
  
  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .TRUE.
     write(out_string, *) 'FAILED integer scalar with MASK SUM reduction test:'
     call pgslib_error(out_string)
     call pgslib_output(out_string)
     write(out_string, 20) int_src, int_given, int_expected, (MOD(thispe,2) == 0)
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test sum 1d no mask integer
  src_len = thispe 
  allocate(int_vec1d_src(src_len))
  int_vec1d_src = (/ (thispe+i, i=0,thispe) /)
  int_vec1d_src = MOD(int_vec1d_src*k1 + k2, npe)

  int_expected = PGSLib_Global_Sum( SUM(int_vec1d_src) )
  int_given    = PGSLib_Global_Sum( int_vec1d_src )

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .TRUE.
     write(out_string, *) "FAILED integer 1D vector NO MASK sum reduction test"
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     call PGSLib_Flush_Output()
     write(out_string, * ) '    vec_src'
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     call PGSLib_Flush_Output()
     write(out_string, * ) int_vec1d_src
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     call PGSLib_Flush_Output()
     write(out_string, 36) int_given, int_expected
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     call PGSLib_Flush_Output()
  endif
36   FORMAT(1x, '    int_given = ', i8,' int_expected = ', i8)



  ! Test sum 1d withmask integer
  allocate(vec1d_mask(src_len))
  vec1d_mask = MOD(int_vec1d_src, 3) == 0

  int_expected = PGSLib_Global_Sum( SUM(int_vec1d_src, MASK=vec1d_mask) )
  int_given    = PGSLib_Global_Sum( int_vec1d_src, MASK=vec1d_mask )

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .TRUE.
     call pgslib_error("Failed integer 1d vector with mask sum reduction test")
     call pgslib_output("Failed integer 1d vector with mask sum reduction test")
  endif

  deallocate(vec1d_mask)
  deallocate(int_vec1d_src)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test sum 2d integer no mask
  src_len = thispe
  allocate(int_vec2d_src(src_len, thispe))
  do i = 1, src_len
     do j = 1, thispe 
        int_vec2d_src(i,j) = i*j
     enddo
  enddo
  

  int_vec2d_src = MOD(int_vec2d_src*k1 + k2, npe)

  int_expected = PGSLib_Global_Sum( SUM(int_vec2d_src) )
  int_given    = PGSLib_Global_Sum( int_vec2d_src )

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .TRUE.
     write(out_string, *) 'FAILED integer 2D vector NO MASK sum reduction test'
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, *)  int_vec2d_src
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, 36)  int_given, int_expected
     call pgslib_error (out_string)
     call pgslib_output(out_string)
  endif

  ! Test sum 2d integer with mask
  allocate(vec2d_mask(src_len, thispe))
  
  vec2d_mask = MOD(int_vec2d_src, 3) == 0

  int_expected = PGSLib_Global_Sum( SUM(int_vec2d_src, MASK= vec2d_mask) )
  int_given    = PGSLib_Global_Sum( int_vec2d_src, MASK= vec2d_mask )

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .TRUE.
     call pgslib_error("Failed integer 2d vector with mask sum reduction test")
     call pgslib_output("Failed integer 2d vector with mask sum reduction test")
  endif

  deallocate(vec2d_mask)

  ! Test sum 2d integer with mask, no true
  allocate(vec2d_mask(src_len, thispe))
  
  vec2d_mask = .FALSE.

  int_expected = 0
  int_given    = PGSLib_Global_Sum( int_vec2d_src, MASK= vec2d_mask )

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .TRUE.
     call pgslib_error("Failed integer 2d vector .FALSE. mask sum reduction test")
     call pgslib_output("Failed integer 2d vector .FALSE. mask sum reduction test")
  endif

  deallocate(vec2d_mask)
  
  ! Test minloc 1d no mask
  ! Setup global array as ICs, then dist to each PE
  allocate(src_lengths(npe))
  src_lengths = list_pes + 1
  src_len = SUM(src_lengths)
  allocate(int_global_vec1d_src(src_len))
  offset = 1
  do pe = 1, npe
     length = pe + 1
     int_global_vec1d_src(offset:offset+length-1) = (/ (i*pe, i=1, length) /)
     offset = offset + length
  enddo
  int_global_vec1d_src = MOD(int_global_vec1d_src*k1 + k2, 11*npe)

  ! Dist array to use as source
  allocate(int_vec1d_src(src_lengths(thispe)))

  call pgslib_dist(int_vec1d_src, int_global_vec1d_src, src_lengths)

  int_v_expected = MINLOC(int_global_vec1d_src)
  int_v_given    = PGSLib_Global_MINLOC(int_vec1d_src)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .TRUE.
     call pgslib_error("Failed integer 1d vector MINLOC reduction test")
     call pgslib_output("Failed integer 1d vector MINLOC reduction test")
  endif

  deallocate(int_global_vec1d_src)
  deallocate(int_vec1d_src)
  deallocate(src_lengths)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !          DOT_PRODUCT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test dot_product 1d no mask
  src_len = thispe 
  allocate(int_vec1d_src(src_len))
  allocate(int_vec1d_src_b(src_len))
  int_vec1d_src = (/ (thispe+i, i=0,thispe) /)
  int_vec1d_src = MOD(int_vec1d_src*k1 + k2, npe)
  int_vec1d_src_b = MOD(int_vec1d_src*k1 + k2, npe)

  int_expected = PGSLib_Global_Sum( DOT_PRODUCT(int_vec1d_src, int_vec1d_src_b) )
  int_given    = PGSLib_Global_DOT_PRODUCT(int_vec1d_src, int_vec1d_src_b)

  if (WITHIN_TOLERANCE(int_expected, int_given)) then
     local_error = .TRUE.
     write(out_string, *) 'FAILED integer 1d VECTOR NO MASK dot_product reduction test'
     call pgslib_error(out_string)
     call pgslib_output(out_string)
     write(out_string, *) '   vec_src_1'
     call pgslib_error(out_string)
     call pgslib_output(out_string)
     write(out_string, *) int_vec1d_src
     call pgslib_error(out_string)
     call pgslib_output(out_string)
     write(out_string, *) '   vec_src_2'
     call pgslib_error(out_string)
     call pgslib_output(out_string)
     write(out_string, *) int_vec1d_src_b
     call pgslib_error(out_string)
     call pgslib_output(out_string)
     write(out_string, 36) int_given, int_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif


  if (.not. local_error) then
     call pgslib_output("Passed integer reduction test")
  endif
  error = error .or. local_error

  return
end subroutine test_redux_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test reduction of reals
subroutine test_redux_real(error)
  USE PGSLib_MODULE
  use error_tests
  implicit none
  logical(PGSLib_Log_Type):: error

  !Local
  logical(PGSLib_Log_Type):: local_error

  real(pgslib_real_type) :: real_given, real_expected, real_src
  real(pgslib_real_type), pointer, dimension(:) :: real_vec1d_src, real_temp
  real(pgslib_real_type), pointer, dimension(:) :: real_vec1d_src_b
  real(pgslib_real_type), pointer, dimension(:) :: real_global_vec1d_src
  logical(pgslib_log_type), pointer, dimension(:) :: vec1d_mask
  real, pointer, dimension(:,:) :: real_vec2d_src
  logical, pointer, dimension(:,:) :: vec2d_mask
  integer (pgslib_Int_type), pointer, dimension(:):: list_pes, int_temp, src_lengths
  integer (pgslib_Int_type) :: n, npe, thispe, pe, src_len
  integer (pgslib_Int_type), dimension(1) :: int_v_given, int_v_expected
  integer (pgslib_Int_type) :: offset, length
  logical(PGSLib_Log_Type) :: src_mask

  integer (pgslib_int_type) :: k1, k2, i, j

  CHARACTER (LEN=2048) :: out_string

  k1 = 1000003
  k2 = 262147


  npe = PGSLib_Inquire_nPE()
  thispe = PGSLib_Inquire_thisPE_Actual()

  local_error = .false.

  ! Setup for permutations
!  k1 = MOD(k1, npe)
!  k2 = MOD(k2, npe)
  
  allocate(list_pes (npe))
  list_pes = (/ (pe, pe=1,npe) /)  
  allocate(real_temp(npe))
  allocate(int_temp(npe))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! Scalar real min/max !!!!!!!!!!
  ! Test min scalar, no mask
  real_src = cos(REAL(MOD(thispe*k1 + k2, npe)) )
  real_expected = MINVAL( cos(REAL(MOD(list_pes*k1+k2, npe))) )
  
  real_given = pgslib_global_minval(real_src)

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     Local_error  = .true.
     call pgslib_error("Failed real scalar no mask minval reduction test")
     call pgslib_output("Failed real scalar no mask minval reduction test")
  endif

  
  ! Test min scalar, with mask
  real_src = cos(REAL(MOD(thispe*k1 + k2, npe)) )
  src_mask = MOD(thispe,3) == 0
  real_expected = MINVAL( cos(real(MOD(list_pes*k1+k2, npe))) , MASK=(MOD(list_pes,3) == 0) )
  
  real_given = pgslib_global_minval(real_src, MASK=src_mask)

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .true.
     call pgslib_error("Failed real scalar with mask minval reduction test")
     call pgslib_output("Failed real scalar with mask minval reduction test")
  endif

  
  ! Test max scalar no mask
  real_src = cos(REAL(MOD(thispe*k1 + k2, npe)) )
  real_expected = MAXVAL( cos(REAL(MOD(list_pes*k1+k2, npe))) )
  
  real_given = pgslib_global_maxval(real_src)

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .true.
     call pgslib_error("Failed real scalar no mask maxval reduction test")
     call pgslib_output("Failed real scalar no mask maxval reduction test")
  endif

  ! Test max scalar with mask
  real_src = cos(REAL(MOD(thispe*k1 + k2, npe)) )
  src_mask = MOD(thispe,3) == 0
  real_expected = MAXVAL( cos(real(MOD(list_pes*k1+k2, npe))) , MASK=(MOD(list_pes,3) == 0) )
  
  real_given = pgslib_global_maxval(real_src, MASK=src_mask)

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .true.
     call pgslib_error("Failed real scalar with mask maxval reduction test")
     call pgslib_output("Failed real scalar with mask maxval reduction test")
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! vector 1d real min/max
  src_len = MOD(2*thispe, 3) + 1
  allocate(real_vec1d_src(src_len))
  
  real_vec1d_src = cos(real((/ (n, n=1,src_len) /)))

  ! Test min 1d vector no mask

  real_expected = PGSLib_Global_Minval( MINVAL(real_vec1d_src) )
  real_given    = PGSLib_Global_Minval( real_vec1d_src )

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .TRUE.
     WRITE(out_string, *) "FAILED real 1d vector NO MASK MINVAL reduction test"
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, *) '    vec_src'
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, *) real_vec1d_src
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, 60) real_given, real_expected, ABS(real_given - real_expected)
     call pgslib_error (out_string)
     call pgslib_output(out_string)
  endif
60   FORMAT(1x, '     real_given = ', g12.6, ' real_expected = ', g12.6, ' diff = ', g12.6)


  ! Test max 1d vector no mask

  real_expected = PGSLib_Global_Maxval( MAXVAL(real_vec1d_src) )
  real_given    = PGSLib_Global_Maxval( real_vec1d_src )

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .TRUE.
     call pgslib_error("Failed real 1d vector no mask maxval reduction test")
     call pgslib_output("Failed real 1d vector no mask maxval reduction test")
  endif

  ! Test min 1d vector with mask

  allocate(vec1d_mask(src_len))
  vec1d_mask = real_vec1d_src .ge. 0.0

  real_expected = PGSLib_Global_Minval( MINVAL(real_vec1d_src, MASK=vec1d_mask) )
  real_given    = PGSLib_Global_Minval( real_vec1d_src, MASK=vec1d_mask)

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .TRUE.
     call pgslib_error("Failed real 1d vector with mask minval reduction test")
     call pgslib_output("Failed real 1d vector with mask minval reduction test")
  endif


  ! Test max 1d vector with mask
  real_expected = PGSLib_Global_Maxval( MAXVAL(real_vec1d_src, MASK=vec1d_mask) )
  real_given    = PGSLib_Global_Maxval( real_vec1d_src, MASK=vec1d_mask)

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .true.
     call pgslib_error("Failed real 1d vector with mask maxval reduction test")
     call pgslib_output("Failed real 1d vector with mask maxval reduction test")
  endif


  deallocate(real_vec1d_src)
  deallocate(vec1d_mask)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! vector 2d real min/max
  src_len = MOD(2*thispe, 3) + 1
  allocate(real_vec2d_src(src_len, 2*npe))
  
  do i=1,src_len
     do j = 1,2*npe
        real_vec2d_src(i,j) = cos(real(i+j))
     enddo
  enddo
  
  ! Test min 2d vector no mask

  real_expected = PGSLib_Global_Minval( MINVAL(real_vec2d_src) )
  real_given    = PGSLib_Global_Minval( real_vec2d_src )

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .true.
     call pgslib_error("Failed real 2d vector no mask minval reduction test")
     call pgslib_output("Failed real 2d vector no mask minval reduction test")
  endif

  ! Test max 2d vector no mask

  real_expected = PGSLib_Global_Maxval( MAXVAL(real_vec2d_src) )
  real_given    = PGSLib_Global_Maxval( real_vec2d_src )

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .true.
     call pgslib_error("Failed real 2d vector no mask maxval reduction test")
     call pgslib_output("Failed real 2d vector no mask maxval reduction test")
  endif

  ! Test min 2d vector with mask

  allocate(vec2d_mask(src_len, 2*npe))
  vec2d_mask = real_vec2d_src .ge. 0.0

  real_expected = PGSLib_Global_Minval( MINVAL(real_vec2d_src, MASK=vec2d_mask) )
  real_given    = PGSLib_Global_Minval( real_vec2d_src, MASK=vec2d_mask)

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .true.
     call pgslib_error("Failed real 2d vector with mask minval reduction test")
     call pgslib_output("Failed real 2d vector with mask minval reduction test")
  endif

  ! Test max 2d vector with mask

  real_expected = PGSLib_Global_Maxval( MAXVAL(real_vec2d_src, MASK=vec2d_mask) )
  real_given    = PGSLib_Global_Maxval( real_vec2d_src, MASK=vec2d_mask)

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .true.
     call pgslib_error("Failed real 2d vector with mask maxval reduction test")
     call pgslib_output("Failed real 2d vector with mask maxval reduction test")
  endif

  deallocate(real_vec2d_src)
  deallocate(vec2d_mask)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!! TEST SUM !!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test sum scalar no mask real
  real_src = thispe
  real_src = cos(real_src)
  real_temp = list_pes
  real_expected = SUM( cos( real_temp) )
  real_given = pgslib_global_sum( real_src )
  
  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .TRUE.
     write(out_string,*)'FAILED real scalar NO MASK sum reduction test'
     call pgslib_error(out_string)
     call pgslib_output(out_string)
     write(out_string,50) real_src, real_given, real_expected, ABS(real_given - real_expected)
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif
50   FORMAT(1x, ' real_src = ', g12.6, ' real_given = ', g12.6, ' real_expected =', g12.6, &
          &        ' diff = ', g12.6)
  ! Test sum scalar with mask real
  real_temp = list_pes
  real_temp = cos(real_temp)

  real_expected = SUM(real_temp, MASK=(MOD(list_pes,2) == 0))
  
  real_given = pgslib_global_sum( real_src , MASK=(MOD(thispe,2) == 0))
  
  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .TRUE.
     call pgslib_error("Failed real scalar with mask sum reduction test")
     call pgslib_output("Failed real scalar with mask sum reduction test")
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test sum 1d no mask real
  src_len = thispe 
  allocate(real_vec1d_src(src_len))
  real_vec1d_src = cos(real( (/ (thispe+i, i=0,thispe) /)) )

  real_expected = PGSLib_Global_Sum( SUM(real_vec1d_src) )
  real_expected = PGSLib_Global_Sum( SUM(real_vec1d_src) )
  real_given    = PGSLib_Global_Sum( real_vec1d_src )
  real_given    = PGSLib_Global_Sum( real_vec1d_src )

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .TRUE.
     WRITE(out_string, *) "FAILED real 1d vector NO MASK SUM reduction test"
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, *) '    vec_src'
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, *) real_vec1d_src
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, 60) real_given, real_expected, ABS(real_given - real_expected)
     call pgslib_error (out_string)
     call pgslib_output(out_string)
  endif



  ! Test sum 1d withmask real
  allocate(vec1d_mask(src_len))
  vec1d_mask = real_vec1d_src .ge. 0.0

  real_expected = PGSLib_Global_Sum( SUM(real_vec1d_src, MASK=vec1d_mask) )
  real_given    = PGSLib_Global_Sum( real_vec1d_src, MASK=vec1d_mask )

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .TRUE.
     write(out_string,*) 'FAILED real 1D WITH MASK sum reduction test'
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, *) '    vec_src'
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, *) real_vec1d_src
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, *) '    vec_mask'
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, *) vec1d_mask
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string,66) real_given, real_expected, ABS(real_given - real_expected)
     call pgslib_error (out_string)
     call pgslib_output(out_string)
  endif
66 FORMAT(1x, '    real_given = ', g12.6, ' real_expected = ', g12.6, &
          &        ' diff = ', g12.6)

  deallocate(vec1d_mask)
  deallocate(real_vec1d_src)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test sum 2d real no mask
  src_len = thispe 
  allocate(real_vec2d_src(src_len, thispe))
  do i = 1, thispe 
     do j = 1, thispe 
        real_vec2d_src(i,j) = cos(real(i*j))
     enddo
  enddo
  
  real_expected = PGSLib_Global_Sum( SUM(real_vec2d_src) )
  real_given    = PGSLib_Global_Sum( real_vec2d_src )

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .TRUE.
     write(out_string, *) "FAILED real 2D vector NO MASK sum reduction test"
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, *) real_vec2d_src
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string, 60) real_given, real_expected, ABS(real_given - real_expected)
     call pgslib_error (out_string)
     call pgslib_output(out_string)
  endif

  ! Test sum 2d real with mask
  allocate(vec2d_mask(src_len, thispe))
  
  vec2d_mask = real_vec2d_src .ge. 0.0

  real_expected = PGSLib_Global_Sum( SUM(real_vec2d_src, MASK= vec2d_mask) )
  real_given    = PGSLib_Global_Sum( real_vec2d_src, MASK= vec2d_mask )

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .TRUE.
     write(out_string,*)'FAILED real 2D WITH MASK sum reduction test'
     call pgslib_error(out_string)
     call pgslib_output(out_string)
     write(out_string, *) '    vec_src'
     call pgslib_error(out_string)
     call pgslib_output(out_string)
     write(out_string, *) real_vec2d_src
     call pgslib_error(out_string)
     call pgslib_output(out_string)
     write(out_string, *) '    vec_mask'
     call pgslib_error(out_string)
     call pgslib_output(out_string)
     write(out_string, *) vec2d_mask
     call pgslib_error(out_string)
     call pgslib_output(out_string)
     write(out_string,60) real_given, real_expected, &
                          ABS(real_given - real_expected)
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif
  
  ! Test sum 2d real with mask, no true
  vec2d_mask = .FALSE.

  real_expected = 0.0
  real_given    = PGSLib_Global_Sum( real_vec2d_src, MASK= vec2d_mask )

  if (WITHIN_TOLERANCE(real_expected, real_given)) then
     local_error = .TRUE.
     call pgslib_error("Failed real 2d vector .FALSE. mask sum reduction test")
     call pgslib_output("Failed real 2d vector .FALSE. mask sum reduction test")
  endif


  deallocate(vec2d_mask)
  deallocate(real_vec2d_src)

  ! Test minloc 1d no mask
  
  ! Setup source by making a test vector on All the PE's, the distributing IO_ROOT's
  allocate(src_lengths(npe))
  src_lengths = list_pes + 1
  src_len = SUM(src_lengths)
  allocate(real_global_vec1d_src(src_len))
  offset = 1
  do pe = 0, npe-1
     length = pe+2
     real_global_vec1d_src(offset:offset+length-1) = -(/ (i*(pe+1), i=1, length) /)
     offset = offset + length
  enddo
  real_global_vec1d_src = npe*cos(real_global_vec1d_src)

  ! Now setup local sources by distributing from IO_ROOT
  allocate(real_vec1d_src(thispe+1))
  call pgslib_dist(real_vec1d_src, real_global_vec1d_src, src_lengths)

  int_v_expected = MINLOC(real_global_vec1d_src)
  int_v_given    = PGSLib_Global_MINLOC(real_vec1d_src)

  if (WITHIN_TOLERANCE(int_v_given(1), int_v_expected(1))) then
     local_error = .true.
     call pgslib_error("Failed real 1d vector MINLOC reduction test")
     call pgslib_output("Failed real 1d vector MINLOC reduction test")
  endif

  deallocate(src_lengths)
  deallocate(real_global_vec1d_src)
  deallocate(real_vec1d_src)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !          DOT_PRODUCT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test dot_product 1d no mask
  src_len = thispe 
  allocate(real_vec1d_src(src_len))
  allocate(real_vec1d_src_b(src_len))
  real_vec1d_src = cos(real( (/ (thispe+i, i=0,thispe) /)) )
  real_vec1d_src_b = sin(real( (/ (thispe+i, i=0,thispe) /)) )

  real_expected = PGSLib_Global_SUM(DOT_PRODUCT(real_vec1d_src, real_vec1d_src_b))
  real_given    = PGSLib_Global_DOT_PRODUCT(real_vec1d_src, real_vec1d_src_b)

  if (WITHIN_TOLERANCE( real_expected, real_given )) then
     local_error = .true.
     call pgslib_error("Failed single 1d vector no mask dot_product reduction test")
     call pgslib_output("Failed single 1d vector no mask dot_product reduction test")
     write(out_string,*) ' real_given = ',real_given, ' real_expected = ', real_expected
     Call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif

  deallocate(real_vec1d_src)
  deallocate(real_vec1d_src_b)

  if (.not. local_error) then
     call pgslib_output("Passed real reduction test")
  endif
  error = error .or. local_error
  
  return
end subroutine test_redux_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test reduction of doubles
subroutine test_redux_double(error)
  USE PGSLib_MODULE
  use error_tests
  implicit none
  logical(PGSLib_Log_Type):: error

  !Local
  logical(PGSLib_Log_Type):: local_error

  real (pgslib_double_type) :: double_given, double_expected, double_src
  real (pgslib_double_type), pointer, dimension(:) :: double_vec1d_src, double_temp
  real (pgslib_double_type), pointer, dimension(:) :: double_vec1d_src_b
  logical(pgslib_log_type), pointer, dimension(:) :: vec1d_mask
  real, pointer, dimension(:,:) :: double_vec2d_src , double_vec2d_src_b
  logical, pointer, dimension(:,:) :: vec2d_mask
  integer (pgslib_Int_type), pointer, dimension(:):: list_pes
  integer (pgslib_Int_type) :: n, npe, thispe, pe, src_len
  logical(PGSLib_Log_Type) :: src_mask

  integer (pgslib_int_type) :: k1, k2, i, j

  CHARACTER (LEN=2048) :: out_string

  k1 = 1000003
  k2 = 262147

  npe = PGSLib_Inquire_nPE()
  thispe = PGSLib_Inquire_thisPE_Actual()
  local_error = .false.

  ! Setup for permutations
  k1 = MOD(k1, npe)
  k2 = MOD(k2, npe)
  
  allocate(list_pes (npe))
  allocate(double_temp(npe))
  list_pes = (/ (pe, pe=1,npe) /)  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! Scalar double min/max !!!!!!!!!!
  ! Test min scalar, no mask
  double_src = MOD(thispe*k1 + k2, npe)
  double_src = cos(double_src)
  double_temp = MOD(list_pes*k1+k2, npe)

  double_expected = MINVAL( cos(double_temp) )
  double_given = pgslib_global_minval(double_src)

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double scalar no mask minval reduction test")
     call pgslib_output("Failed double scalar no mask minval reduction test")
     write(out_string,*) ' double_given = ',double_given, '  double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif

  
  ! Test min scalar, with mask
  double_src = MOD(thispe*k1 + k2, npe)
  double_src = cos(double_src)
  src_mask = MOD(thispe+2,3) == 0
  double_temp = MOD(list_pes*k1+k2, npe)
  double_expected = MINVAL( cos(double_temp) , MASK=(MOD(list_pes+2,3) == 0) )
  
  double_given = pgslib_global_minval(double_src, MASK=src_mask)

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double scalar with mask minval reduction test")
     call pgslib_output("Failed double scalar with mask minval reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif

  
  ! Test max scalar no mask
  double_src = MOD(thispe*k1 + k2, npe)
  double_src = cos(double_src)
  double_temp = MOD(list_pes*k1+k2, npe)

  double_expected = MAXVAL( cos(double_temp) )
  double_given = pgslib_global_maxval(double_src)

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double scalar no mask maxval reduction test")
     call pgslib_output("Failed double scalar no mask maxval reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif

  ! Test max scalar with mask
  double_src = MOD(thispe*k1 + k2, npe)
  double_src = cos(double_src)
  src_mask = MOD(thispe+2,3) == 0
  double_temp = MOD(list_pes*k1+k2, npe)
  double_expected = MAXVAL( cos(double_temp) , MASK=(MOD(list_pes+2,3) == 0) )
  
  double_given = pgslib_global_maxval(double_src, MASK=src_mask)

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double scalar with mask maxval reduction test")
     call pgslib_output("Failed double scalar with mask maxval reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! vector 1d double min/max
  src_len = MOD(2*thispe, 3) + 1
  allocate(double_vec1d_src(src_len))
  
  double_vec1d_src = cos(real((/ (n, n=1,src_len) /)))

  ! Test min 1d vector no mask

  double_expected = PGSLib_Global_Minval( MINVAL(double_vec1d_src) )
  double_given    = PGSLib_Global_Minval( double_vec1d_src )

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double 1d vector no mask minval reduction test")
     call pgslib_output("Failed double 1d vector no mask minval reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif

  ! Test max 1d vector no mask

  double_expected = PGSLib_Global_Maxval( MAXVAL(double_vec1d_src) )
  double_given    = PGSLib_Global_Maxval( double_vec1d_src )

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double 1d vector no mask maxval reduction test")
     call pgslib_output("Failed double 1d vector no mask maxval reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif

  ! Test min 1d vector with mask

  allocate(vec1d_mask(src_len))
  vec1d_mask = double_vec1d_src .ge. 0.0

  double_expected = PGSLib_Global_Minval( MINVAL(double_vec1d_src, MASK=vec1d_mask) )
  double_given    = PGSLib_Global_Minval( double_vec1d_src, MASK=vec1d_mask)

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double 1d vector with mask minval reduction test")
     call pgslib_output("Failed double 1d vector with mask minval reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif


  ! Test max 1d vector with mask
  double_expected = PGSLib_Global_Maxval( MAXVAL(double_vec1d_src, MASK=vec1d_mask) )
  double_given    = PGSLib_Global_Maxval( double_vec1d_src, MASK=vec1d_mask)

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double 1d vector with mask maxval reduction test")
     call pgslib_output("Failed double 1d vector with mask maxval reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif


  deallocate(double_vec1d_src)
  deallocate(vec1d_mask)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! vector 2d double min/max
  src_len = MOD(2*thispe, 3) + 1
  allocate(double_vec2d_src(src_len, 2*npe))
  
  do i=1,src_len
     do j = 1,2*npe
        double_vec2d_src(i,j) = cos(real(i+j))
     enddo
  enddo
  
  ! Test min 2d vector no mask

  double_expected = PGSLib_Global_Minval( MINVAL(double_vec2d_src) )
  double_given    = PGSLib_Global_Minval( double_vec2d_src )

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double 2d vector no mask minval reduction test")
     call pgslib_output("Failed double 2d vector no mask minval reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif

  ! Test max 2d vector no mask

  double_expected = PGSLib_Global_Maxval( MAXVAL(double_vec2d_src) )
  double_given    = PGSLib_Global_Maxval( double_vec2d_src )

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double 2d vector no mask maxval reduction test")
     call pgslib_output("Failed double 2d vector no mask maxval reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif

  ! Test min 2d vector with mask

  allocate(vec2d_mask(src_len, 2*npe))
  vec2d_mask = double_vec2d_src .ge. 0.0

  double_expected = PGSLib_Global_Minval( MINVAL(double_vec2d_src, MASK=vec2d_mask) )
  double_given    = PGSLib_Global_Minval( double_vec2d_src, MASK=vec2d_mask)

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double 2d vector with mask minval reduction test")
     call pgslib_output("Failed double 2d vector with mask minval reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif

  ! Test max 2d vector with mask

  double_expected = PGSLib_Global_Maxval( MAXVAL(double_vec2d_src, MASK=vec2d_mask) )
  double_given    = PGSLib_Global_Maxval( double_vec2d_src, MASK=vec2d_mask)

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double 2d vector with mask maxval reduction test")
     call pgslib_output("Failed double 2d vector with mask maxval reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif

  deallocate(double_vec2d_src)
  deallocate(vec2d_mask)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!! TEST SUM !!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test sum scalar no mask double
  double_src = thispe 
  double_src = cos(double_src)
  double_temp = list_pes 

  double_expected = SUM( cos(double_temp))
  double_given = pgslib_global_sum( double_src )
  
  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double scalar no mask sum reduction test")
     call pgslib_output("Failed double scalar no mask sum reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif
  
  ! Test sum scalar with mask double
  double_temp = list_pes + 1
  double_src = thispe + 1
  double_src = cos(double_src)
  double_expected = SUM(cos(double_temp), MASK=(MOD(list_pes,2) == 0))
  double_given = pgslib_global_sum( double_src , MASK=(MOD(thispe,2) == 0))
  
  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double scalar with mask sum reduction test")
     call pgslib_output("Failed double scalar with mask sum reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test sum 1d no mask double
  src_len = thispe + 1
  allocate(double_vec1d_src(src_len))
  double_vec1d_src = cos(real( (/ (thispe+i, i=0,thispe) /)) )

  double_expected = PGSLib_Global_Sum( SUM(double_vec1d_src) )
  double_given    = PGSLib_Global_Sum( double_vec1d_src )

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double 1d vector no mask sum reduction test")
     call pgslib_output("Failed double 1d vector no mask sum reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif


  ! Test sum 1d withmask double
  allocate(vec1d_mask(src_len))
  vec1d_mask = double_vec1d_src .ge. 0.0

  double_expected = PGSLib_Global_Sum( SUM(double_vec1d_src, MASK=vec1d_mask) )
  double_given    = PGSLib_Global_Sum( double_vec1d_src, MASK=vec1d_mask )

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double 1d vector with mask sum reduction test")
     call pgslib_output("Failed double 1d vector with mask sum reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif

  deallocate(vec1d_mask)
  deallocate(double_vec1d_src)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test sum 2d double no mask
  src_len = thispe 
  allocate(double_vec2d_src(src_len, thispe))
  do i = 1, thispe 
     do j = 1, thispe 
        double_vec2d_src(i,j) = cos(real(i*j))
     enddo
  enddo
  
  double_expected = PGSLib_Global_Sum( SUM(double_vec2d_src) )
  double_given    = PGSLib_Global_Sum( double_vec2d_src )

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .true.
     call pgslib_error("Failed double 2d vector no mask sum reduction test")
     call pgslib_output("Failed double 2d vector no mask sum reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif

  ! Test sum 2d double with mask
  allocate(vec2d_mask(src_len, thispe))
  
  vec2d_mask = double_vec2d_src .ge. 0.0

  double_expected = PGSLib_Global_Sum( SUM(double_vec2d_src, MASK= vec2d_mask) )
  double_given    = PGSLib_Global_Sum( double_vec2d_src, MASK= vec2d_mask )

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .TRUE.
     write(out_string, *) "Failed double 2d vector with mask sum reduction test"
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string,*) '   vec_src'
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string,*) double_vec2d_src
     call pgslib_error (out_string)
     call pgslib_output(out_string)
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif

  ! Test sum 2d double with mask, no true
  vec2d_mask = .FALSE.

  double_expected = 0.0
  double_given    = PGSLib_Global_Sum( double_vec2d_src, MASK= vec2d_mask )

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .TRUE.
     call pgslib_error("Failed double 2d vector .FALSE. mask sum reduction test")
     call pgslib_output("Failed double 2d vector .FALSE. mask sum reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif


  deallocate(vec2d_mask)
  deallocate(double_vec2d_src)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !          DOT_PRODUCT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test dot_product 1d no mask
  src_len = thispe + 1
  allocate(double_vec1d_src(src_len))
  allocate(double_vec1d_src_b(src_len))
  double_vec1d_src = cos(real( (/ (thispe+i, i=0,thispe) /)) )
  double_vec1d_src_b = sin(real( (/ (thispe+i, i=0,thispe) /)) )

  double_expected = PGSLib_Global_SUM(DOT_PRODUCT(double_vec1d_src, double_vec1d_src_b))
  double_given    = PGSLib_Global_DOT_PRODUCT(double_vec1d_src, double_vec1d_src_b)

  if (WITHIN_TOLERANCE(double_expected, double_given)) then
     local_error = .TRUE.
     call pgslib_error("Failed double 1d vector no mask dot_product reduction test")
     call pgslib_output("Failed double 1d vector no mask dot_product reduction test")
     write(out_string,*) ' double_given = ',double_given, ' double_expected = ', double_expected
     call pgslib_error(out_string)
     call pgslib_output(out_string)
  endif


  deallocate(double_vec1d_src)
  deallocate(double_vec1d_src_b)

  if (.not. local_error) then
     call pgslib_output("Passed double reduction test")
  endif
  error = error .or. local_error
  
  return
end subroutine test_redux_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test reduction of logicals
subroutine test_redux_log(error)
  USE PGSLib_MODULE
  implicit none
  logical(PGSLib_Log_Type):: error

  !Local
  logical(PGSLib_Log_Type):: local_error

  logical (pgslib_log_type) :: log_src
  integer (pgslib_int_type) :: src_len, int_given, int_expected
  integer(pgslib_int_type), pointer, dimension(:) :: list_pes
  logical (pgslib_log_type), pointer, dimension(:) :: log_vec1d_src
  logical (pgslib_log_type), pointer, dimension(:,:) :: log_vec2d_src
  integer(pgslib_int_type) :: n, npe, thispe, pe

  integer(pgslib_int_type) :: k1, k2, i, j

  k1 = 1000003
  k2 = 262147

  npe = PGSLib_Inquire_nPE()
  thispe = PGSLib_Inquire_thisPE_Actual()
  local_error = .false.

  ! Setup for permutations
  k1 = MOD(k1, npe)
  k2 = MOD(k2, npe)
  
  allocate(list_pes (npe))
  list_pes = (/ (pe, pe=1,npe) /)  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! logical COUNT !!!!!!!!!!
  ! Test count scalar, no mask
  log_src = MOD(thispe*k1 + k2, npe) == 0
  int_expected = COUNT( MOD(list_pes*k1+k2, npe) == 0)
  
  int_given = pgslib_global_count(log_src)

  if (int_given .ne. int_expected) then
     local_error = .true.
     call pgslib_error("Failed logical scalar no mask count reduction test")
     call pgslib_output("Failed logical scalar no mask count reduction test")
  endif

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! vector 1d logical COUNT
  src_len = MOD(2*thispe, 3) + 1
  allocate(log_vec1d_src(src_len))
  
  log_vec1d_src = MOD((/ (n, n=1,src_len) /), npe) == 0

  ! Test count 1d vector no mask

  int_expected = PGSLib_Global_SUM( COUNT(log_vec1d_src) )
  int_given    = PGSLib_Global_COUNT( log_vec1d_src )

  if (int_given .ne. int_expected) then
     local_error = .true.
     call pgslib_error("Failed logical 1d vector no mask count reduction test")
     call pgslib_output("Failed logical 1d vector no mask count reduction test")
  endif

  deallocate(log_vec1d_src)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! vector 2d logical COUNT
  src_len = MOD(2*thispe, 3) + 1
  allocate(log_vec2d_src(src_len, 2*npe))
  
  do i=1,src_len
     do j = 1,2*npe
        log_vec2d_src(i,j) = MOD(MOD((i+j)*k1 + k2, npe),3) == 0
     enddo
  enddo
  
  ! Test count 2d vector no mask

  int_expected = PGSLib_Global_SUM( COUNT(log_vec2d_src) )
  int_given    = PGSLib_Global_COUNT( log_vec2d_src )

  if (int_given .ne. int_expected) then
     local_error = .true.
     call pgslib_error("Failed logical 2d vector no mask count reduction test")
     call pgslib_output("Failed logical 2d vector no mask count reduction test")
  endif

  deallocate(log_vec2d_src)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!! LOGICAL ALL !!!!!!!!!!!!!!!!!

  ! Test ALL
  ! Test ANY

  if (.not. local_error) then
     call pgslib_output("Passed logical reduction test")
  endif
  error = error .or. local_error

  return
end subroutine test_redux_log

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tests enough of redux to proceed with assurance to the bcast tests.
subroutine test_redux_prelim()
  USE PGSLib_MODULE
  implicit none

  logical :: log_given, log_expected, log_src
  logical :: error, abort

  log_expected = .false.
  log_given    = .false.

  error = .false.
  if (pgslib_global_any(log_given)) then
     error = .true.
     call pgslib_output("Found TRUE, PGSLIb_Global_Any FAILED in test_redux_prelim")
  endif

  log_given  = .true.
  error = .false.
  if (.NOT. pgslib_global_any(log_given)) then
     error = .true.
     call pgslib_output("Found FALSE, PGSLIb_Global_Any FAILED in test_redux_prelim")
  endif

  ! Set logical to .TRUE. on odd numbered PEs.
  log_given = MOD(PGSLib_Inquire_thisPE(),2) == 1
  if (.NOT.pgslib_global_any(log_given)) then
     error = .true.
     call pgslib_error("Found a FALSE, PGSLIb_Global_Any FAILED in test_redux_prelim")
     call pgslib_output("Found a FALSE, PGSLIb_Global_Any FAILED in test_redux_prelim")
  endif

  abort = pgslib_global_any(error)
  if (abort) then
     call pgslib_fatal_error("Failed test_redux_prelim")
  else
     call pgslib_output("Passed test_redux_prelim")
  endif

  return
end subroutine test_redux_prelim



