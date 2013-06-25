! -*- Mode: f90 -*-

  ! Args:
  integer(Ubik_int_type), intent(in) :: scale, nsizeb, nrows, plot
  real(Ubik_real_type), intent(in) :: dvalue
  type(Ubik_control_type), intent(inout) :: control
  real(Ubik_real_type), intent(out) :: proc_time, elapsed_time
  real(Ubik_real_type), intent(out), dimension(:) :: b, bsave, ref, x

  ! Local:
  integer(Ubik_int_type) :: incr, status, &
       count_start, count_end, count_rate, count_max

  real(Ubik_real_type) :: time_start, time_end
!!$  real(Ubik_real_type) :: tmp

#ifdef use_PGSLib
  ! Integers for testing with PGSLib.
  integer(Ubik_int_type) :: npe, nrows_pe, rcount=0

  ! PGSLib derived type.
  type(PGSLib_GS_struct) :: trace

  ! Dynamically-allocated arrays for testing with PGSLib.
  real(Ubik_real_type), dimension(:) :: b_pe, x_pe
  real(Ubik_real_type), dimension(:,:) :: a_trans, a_pe, b_tmp
  integer(Ubik_int_type), dimension(:,:) :: ja_trans, ja_pe
  integer(Ubik_int_type), dimension(:) :: lengths, ja_aux
#endif

  select case (prob)
  case (Read_System)

     ! Read in system.
     call Ubik_Load (1, coeff, x, b)
     close (unit=1)

  case default

     ! Set up Saad's problem.
     call SetUpSaad (nsizeb, dvalue, coeff, b, ref, x)
  end select

  ! set trivial initial guess
  x = 0.0

!!$  ! set random initial guess
!!$  call random_number(x)
!!$  tmp = 0.00000001
!!$  x = 1.0 + tmp*(x - 0.5)
!!$  tmp = SUM(x) / REAL(SIZE(x))
!!$  write(6,*) MinVal(x), tmp, MaxVal(x)

  ! Plot structure of coefficient.
  if (Ubik_on_IO_PE()) then
     
     ! Create data file.
     open (unit=1, file='UbikSolve_mat.dat', status='replace')
     call Ubik_Dump (1, coeff, x, b)
     close (unit=1)

     ! Create gnuplot command file.
     open (unit=1, file='UbikSolve_mat.gpl', status='replace')
     write(1,400) 'set title "Coefficient Sparsity Structure"'
     write(1,400) 'set size square'
     write(1,400) 'set nokey'
     write(1,400) 'set xlabel "j"'
     write(1,400) 'set ylabel "i"'
     write(1,400) 'set xrange [1:', nrows, ']'
     write(1,400) 'set yrange [', nrows, ':1]'
     incr = nrows / 10
     write(1,401) 'set xtics (1,', incr, ',', 2*incr, ',', 3*incr, &
          ',', 4*incr, ',', 5*incr, ',', 6*incr, ',', 7*incr, &
          ',', 8*incr, ',', 9*incr, ',', nrows, ')'
     write(1,401) 'set ytics (1,', incr, ',', 2*incr, ',', 3*incr, &
          ',', 4*incr, ',', 5*incr, ',', 6*incr, ',', 7*incr, &
          ',', 8*incr, ',', 9*incr, ',', nrows, ')'
     write(1,400) 'plot "UbikSolve_mat.dat" index 0 using 2:1 with dots linetype 2'
     close (unit=1)

     ! Create gnuplot driver.
     open (unit=1, file='UbikSolve_mat.gpl_driver', status='replace')
     write(1,400) 'load "UbikSolve_mat.gpl"'
     write(1,400) 'pause -1 "Hit return to continue. "' 
     close (unit=1)

     ! Run gnuplot to display the plot.
     if (plot == 1) then
        call SYSTEM ('gnuplot '//'UbikSolve_mat.gpl_driver')
     endif

  endif

  ! Write out system and initial guess prior to scaling.
  call Ubik_Spew (Ubik_luout(control), coeff, x, b, 'SYSTEM PRIOR TO MANIPULATION (coeff, x, b):')

  ! For ELL storage this ensures that main diagonal elements are in the
  ! first column and that elements in each row are left-justified.
  ! For full storage the routine simply returns.
  call Ubik_EnsureFormat (coeff)

  ! Scale system as desired.
  if (scale >= 1) then
     call Ubik_ScaleSystemByRow (coeff, b)
     if (scale == 2) then
        call Ubik_ScaleSystemByDiagonal (coeff, b)
     end if
  end if

  ! Write out system to be solved, along with initial guess.
  call Ubik_Spew (Ubik_luout(control), coeff, x, b, 'SYSTEM TO BE SOLVED (coeff, x, b):')

  call Ubik_FlushUnit (Ubik_luout(control), status)
  if (status /= 0) then
     write(6,*)
     write(6,*) 'Error in Ubik_FlushUnit.'
     stop
  end if

  ! Copy system.  For iterative methods it's used for comparison after solution
  ! just to make sure it doesn't get munged.  For direct methods it's also used
  ! to calculate final residual since the coeff. gets overwritten.
  coeff_save = coeff
  bsave = b

#ifdef use_PGSLib

  ! This is all bogus and out of date.
  a_trans = TRANSPOSE(a)
  ja_trans = TRANSPOSE(ja)

  do i=1,maxnz
     call PGSLib_dist (a_pe(i,:), a_trans(i,:), lengths)
     call PGSLib_dist (ja_pe(i,:), ja_trans(i,:), lengths)
  end do

!!$  do i=1,maxnz
!!$     call PGSLib_dist (a_pe(:,i), a(:,i), lengths)
!!$     call PGSLib_dist (ja_pe(:,i), ja(:,i), lengths)
!!$  end do

  call PGSLib_dist (x_pe, x, lengths)
  call PGSLib_dist (b_pe, b, lengths)

  call PGSLib_GS_setup (ja_pe, ja_aux, nrows_pe, trace, mask=(ja_pe /= 0))
#endif

  ! Start solution.
  if (Ubik_on_IO_PE()) then
     write(6,*)
     write(6,*) 'Set-up complete.  Starting solve...'
  end if

  ! Initialize and start timers.
  call CPU_Time (time_start)
  call System_Clock (count_start, count_rate, count_max)

  ! Compute norm if necessary for stopping test.
  if (Ubik_stopping_axb(control)) then
     call Ubik_set_anorm (control, Ubik_Norm(Ubik_norm_type(control), coeff))
  end if

  if (Ubik_method(control) == Ubik_method_LU) then

     call Ubik_LUdecomp (coeff)
     x = b
     call Ubik_SolveLTriang (unit_diagonal=.TRUE., use_transpose=.FALSE., a=coeff, x=x)
     call Ubik_SolveUTriang (unit_diagonal=.FALSE., use_transpose=.FALSE., a=coeff, x=x)

  else if (Ubik_method(control) == Ubik_method_Jacobi .or. &
           Ubik_method(control) == Ubik_method_SOR .or. &
           Ubik_method(control) == Ubik_method_SSOR) then

     iterate_control%method = Ubik_method(control)
     call Ubik_Stationary (x, b, control, MatVec, Iterate)
  
  else

     ! Need #ifdef use_PGSLib blocks for each solver call that send b_pe and x_pe
     ! instead of b and x.

     if (precond_control%precond /= Ubik_precond_none) then
        call GeneratePrecond (precond_control, coeff, precond)
     end if

     select case (Ubik_method(control))
     case (Ubik_method_CG)
        call Ubik_CG (x, b, control, MatVec, ApplyPrecond)
     case (Ubik_method_GMRES_Left_saad)
        call Ubik_GMRES_Left_Saad (x, b, control, MatVec, ApplyPrecond)
     case (Ubik_method_FGMRES)
        call Ubik_FGMRES (x, b, control, MatVec, ApplyPrecond)
     case (Ubik_method_BCGS, Ubik_method_TFQMR)
        call Ubik_BCGS (x, b, control, MatVec, ApplyPrecond)
     case (Ubik_method_BCGSTAB)
        call Ubik_BCGSTAB (x, b, control, MatVec, ApplyPrecond)
     case (Ubik_method_GMRES_Left_Radicati)
        call Ubik_GMRES_Left_Radicati (x, b, control, MatVec, ApplyPrecond)
     end select
  end if

  ! Stop timers.
  call CPU_Time (time_end)
  proc_time = time_end - time_start

  call System_Clock (count_end, count_rate, count_max)
  elapsed_time = REAL(count_end - count_start) / REAL(count_rate)

#ifdef use_PGSLib
  call PGSLib_collate (x, x_pe)
  call PGSLib_bcast (x)
#endif

  ! ===> This is broken in parallel !!!!!
  if ((Ubik_residual_default(control) .or. Ubik_residual_update(control)) .and. &
       .not.Ubik_stopping_true(control)) then
#ifdef use_PGSLib
     call PGSLib_output( 'true residual calc.' )

     call PGSLib_collate (x, x_pe)
     if (Ubik_on_IO_PE()) then
        write(44,*)
        write(44,*) 'npe: ', npe
        write(44,*)
        write(44,*) 'x:'
        do i=1,nrows
           write(44,'(i4,3x,e22.15)') i, x(i)
        end do
     end if

     b_tmp = 0.0
     call PGSLib_gather (b_tmp, x_pe, ja_pe, ja_aux, trace, mask=(ja_pe /= 0))
     b_pe = b_pe - SUM(a_pe*b_tmp, dim=1)
!!$     b_pe = b_pe - SUM(a_pe*b_tmp, dim=2)

     call PGSLib_collate (b, b_pe)
     if (Ubik_on_IO_PE()) then
        write(44,*)
        write(44,*) 'b:'
        do i=1,nrows
           write(44,'(i4,3x,e22.15)') i, b(i)
        end do
     end if
#else
     b = b - Ubik_MatVec(coeff_save, x)
#endif
     call Ubik_set_rnormt (control, Ubik_Norm(Ubik_norm_type(control), b))
  end if

  return
400 format (a, i7, a, i7, a, f7.4, a)
401 format (a, 10(i7, a), a)
402 format (a, f5.3)
