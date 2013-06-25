! -*- Mode: f90 -*-

!**********************************************************************
!    $Id: UbikTest.F,v 1.32 2003/10/19 22:29:22 turner Exp $
!
!**********************************************************************
program UbikTest
  use Ubik_module

  use Driver_module
  use ELL_data_module
  use Full_data_module
  use Iterate_module
  use MatVec_module
  use Precond_module
  use SetUpSaad_module
  use Usage_module
  use test_data_module
#ifdef use_PGSLib
  use PGSLib_module
#endif
#ifdef NAG_COMPILER
  use f90_unix
  use f90_unix_proc
#endif

  implicit none

  character(len=1) :: backslash, dummy_char
  character(len=80) :: string, sysfile
  character(len=15) :: title2
  character(len=40) :: title1
  character(len=23) :: string0
  character(len=9) :: string1, string2
  character(len=30) :: string3
  character(len=11) :: string4
  character(len=36), dimension(Ubik_method_min:Ubik_method_max) :: method_options
  character(len=21), dimension(Ubik_method_min:Ubik_method_max) :: precond_options

#ifndef NAG_COMPILER
  integer(Ubik_int_type) :: IARGC
#endif
  integer(Ubik_int_type) :: i, ios, iplot, iscale, iter, status, nrows, maxnz

  ! unused
!  integer(Ubik_int_type) :: blocks, blocksize, nelem, lower_bandwidth, upper_bandwidth

  real(Ubik_real_type) :: proc_time, diffmax, elapsed_time, &
       cpu, rnorm, err, rnormt, errt

  ! Derived types.
  type(Ubik_control_type) :: control

  ! Control vars, with defaults as appropriate
  integer(Ubik_int_type) :: nsizea=80, nsizeb=20
  integer(Ubik_int_type) :: method, subspace=0
  integer(Ubik_int_type) :: luerr=98, luout=99, lutty=6, freq=0
  integer(Ubik_int_type) :: itmax=0, norm_type=2, residual=0
  real(Ubik_real_type) :: eps=1.0d-6, omega=1.0, dvalue=0.2

  ! Dynamically-allocated vectors.
  real(Ubik_real_type), allocatable, dimension(:) :: b, bsave, ref, x

  ! Dynamically-allocated matrices (others are in modules).
  real(Ubik_real_type), allocatable, dimension(:,:) :: coeff_full_save
  type(Ubik_ELL_matrix) :: coeff_ell_save

#ifdef use_PGSLib
  ! Integers for testing with PGSLib.
  integer(Ubik_int_type) :: npe, nrows_pe, rcount=0

  ! PGSLib derived type.
  type(PGSLib_GS_struct) :: trace

  ! Dynamically-allocated arrays for testing with PGSLib.
  real(Ubik_real_type), allocatable, dimension(:) :: b_pe, x_pe
  real(Ubik_real_type), allocatable, dimension(:,:) :: a_trans, a_pe, b_tmp
  integer(Ubik_int_type), allocatable, dimension(:,:) :: ja_trans, ja_pe
  integer(Ubik_int_type), allocatable, dimension(:) :: lengths, ja_aux
#endif

  ! Set backslash character.
  backslash = CHAR(92)

  ! Set method strings.
  method_options(Ubik_method_LU) = 'LU Decomposition'
  method_options(Ubik_method_CG) = 'CG'
  method_options(Ubik_method_GMRES_Left_Radicati) = 'left-preconditioned GMRES (Radicati)'
  method_options(Ubik_method_GMRES_Left_Saad) = 'left-preconditioned GMRES (Saad)'
  method_options(Ubik_method_FGMRES) = 'FGMRES'
  method_options(Ubik_method_BCGS) = 'BCGS'
  method_options(Ubik_method_TFQMR) = 'TFQMR'
  method_options(Ubik_method_BCGSTAB) = 'Bi-CGSTAB'
  method_options(Ubik_method_SOR) = 'SOR'
  method_options(Ubik_method_SSOR) = 'SSOR'
  method_options(Ubik_method_Jacobi) = 'Jacobi'

  ! Set preconditioner strings.
  precond_options(Ubik_precond_none) = 'none'
  precond_options(Ubik_method_Jacobi) = 'm-step Jacobi'
  precond_options(Ubik_method_SSOR) = 'm-step SSOR'
  precond_options(Ubik_method_IC) = 'IC'
  precond_options(Ubik_method_ILU) = 'ILU'

  ! Set default values for control structures.
#ifdef use_PGSLib
    call PGSLib_initialize(0)
#endif
  call Ubik_set_defaults (control)
  call SetPrecondDefaults (precond_control)
  call SetIterateDefaults (iterate_control)
  
  ! process input
  if (Ubik_on_IO_PE()) then

     ! can't run with no args
     if (IARGC() == 0) then
        call BasicOptions ()
        stop
     end if

     ! help options
     string = ' '
     call Ubik_ParseArg ('help', string)
     select case (string)
     case ('basic')
        call BasicOptions ()
        stop
     case ('all')
        call AllOptions ()
        stop
     end select

     ! unit numbers for output
     call Ubik_ParseArg ('luout', luout)
     call Ubik_ParseArg ('luerr', luerr)
     call Ubik_ParseArg ('lutty', lutty)

     ! frequency of intermediate status output
     call Ubik_ParseArg ('freq', freq)

     ! which problem to solve...
     call Ubik_ParseArg ('problem', prob)
     if (prob == 0) then
        ! ... one read from a file...
        sysfile = 'system.in'
        call Ubik_ParseArg ('sys', sysfile)
        open (unit=1, file=TRIM(sysfile), status='old', iostat=ios)
        if (ios /= 0) then
           print *, 'error opening ', sysfile
           stop
        end if
        close (unit=1)
     else
        ! ... or a predefined problem (Saad's)
        call Ubik_ParseArg ('asymm', dvalue)
        call Ubik_ParseArg ('unknowns', nsizea)
        call Ubik_ParseArg ('blocksize', nsizeb)
     end if

     ! scaling
     string = ' '
     call Ubik_ParseArg ('scale', string)
     select case (string)
     case ('row', 'Row', 'ROW')
        iscale = 1
     case ('diag', 'Diag', 'DIAG', 'diagonal', 'Diagonal', 'DIAGONAL')
        iscale = 2
     case ('both', 'Both', 'BOTH')
        iscale = 3
     case default
        iscale = 0
     end select

     ! norm type, max iters, convergence criterion, and verbosity
     call Ubik_ParseArg ('norm', norm_type)
     call Ubik_ParseArg ('itmax', itmax)
     call Ubik_ParseArg ('eps', eps)
     if (eps < 0.0d0) then
        print *, 'invalid value for eps'
        call BasicOptions ()
        stop
     end if
     string = 'warnings'
     call Ubik_ParseArg ('verbosity', string)
     select case (string)
     case ('none')
        call Ubik_set_outmode_none (control)
     case ('errors')
        call Ubik_set_outmode_errors (control)
     case ('warnings')
        call Ubik_set_outmode_warnings (control)
     case ('summary')
        call Ubik_set_outmode_summary (control)
     case ('system')
        call Ubik_set_outmode_system (control)
     case ('iterates')
        call Ubik_set_outmode_iterates (control)
     case ('full')
        call Ubik_set_outmode_full (control)
     end select

     ! stopping test
     string = 'b'
     call Ubik_ParseArg ('stop', string)
     select case (string)
     case ('relchg')
        call Ubik_set_stopping_relchg (control)
     case ('b')
        call Ubik_set_stopping_b (control)
     case ('axb')
        call Ubik_set_stopping_axb (control)
     case ('x')
        call Ubik_set_stopping_x (control)
     case ('r0')
        call Ubik_set_stopping_r0 (control)
     case ('r')
        call Ubik_set_stopping_r (control)
     case default
        print *, 'invalid value for stopping test'
        call BasicOptions ()
        stop
     end select
     string = ' '
     call Ubik_ParseArg ('use_true_residual', string)
     select case (string)
     case ('y', 'Y', 'yes', 'YES', 'Yes', '1')
        call Ubik_set_stopping_true (control)
     end select

     ! residual update control
     call Ubik_ParseArg ('residual', residual)

     ! residual plotting
     string = ' '
     call Ubik_ParseArg ('plot', string)
     select case (string)
     case ('y', 'Y', 'yes', 'YES', 'Yes', '1')
        iplot = 1
     case default
        iplot = 0
     end select

     ! method
     string = 'gmres'
     call Ubik_ParseArg ('method', string)
     select case (string)
     case ('BCGS', 'bcgs', 'CGS', 'cgs')
        method = Ubik_method_BCGS
     case ('CG', 'cg')
        method = Ubik_method_CG
     case ('GMRES_Left_Radicati', 'gmres_left_radicati', 'GMRES', 'gmres')
        method = Ubik_method_GMRES_Left_Radicati
     case ('GMRES_Left_Saad', 'gmres_left_saad')
        method = Ubik_method_GMRES_Left_Saad
     case ('FGMRES', 'fgmres')
        method = Ubik_method_FGMRES
     case ('TFQMR', 'tfqmr')
        method = Ubik_method_TFQMR
     case ('BCGSTAB', 'bcgstab', 'BICGSTAB', 'bicgstab', 'Bi-CGstab', 'Bi-CGSTAB', 'bi-cgstab')
        method = Ubik_method_BCGSTAB
     case ('SOR', 'sor')
        method = Ubik_method_SOR
     case ('SSOR', 'ssor')
        method = Ubik_method_SSOR
     case ('jac', 'jacobi', 'JAC', 'Jacobi')
        method = Ubik_method_Jacobi
     case ('GE', 'LU')
        method = Ubik_method_LU
     case default
        print *, 'invalid value for method'
        call BasicOptions ()
        stop
     end select
     call Ubik_ParseArg ('subspace', subspace)

     ! preconditioning
     string = 'none'
     call Ubik_ParseArg ('pre', string)
     select case (string)
     case ('none', 'NONE', 'None')
        precond_control%precond = Ubik_precond_none
     case ('jac', 'jacobi', 'JAC', 'Jacobi')
        precond_control%precond = Ubik_method_Jacobi
     case ('SSOR', 'ssor')
        precond_control%precond = Ubik_method_SSOR
     case ('IC', 'ic')
        precond_control%precond = Ubik_method_IC
     case ('ILU', 'ilu')
        precond_control%precond = Ubik_method_ILU
     case default
        print *, 'invalid value for pre'
        call BasicOptions ()
        stop
     end select
     precond_control%steps = 1
     call Ubik_ParseArg ('steps', precond_control%steps)
     call Ubik_ParseArg ('relax', omega)
     precond_control%omega = omega
     iterate_control%omega = omega
  end if

  ! At this point everything necessary for calc has been set.

  if (subspace == 0) then
     subspace = MAX(20,INT(SQRT(FLOAT(nsizea))))
  end if
  if (itmax == 0) then
     itmax = MAX(20,3*INT(SQRT(FLOAT(nsizea))))
     if (method == Ubik_method_Jacobi .OR. &
         method == Ubik_method_SOR .OR. &
         method == Ubik_method_SSOR) then
        itmax = 10*itmax
     else if (precond_control%precond == Ubik_precond_none) then
        itmax = 2*itmax
     end if
  end if

#ifdef dec
  call Ubik_set_control (control, &
#else
  call Ubik_set (control, &
#endif
       eps = eps, &
       output_frequency = freq, &
       itmax = itmax, &
       luout = luout, &
       luerr = luerr, &
       lutty = lutty, &
       method = method, &
       norm_type = norm_type, &
       omega = omega, &
       residual = residual, &
       subspace = subspace)
  call Ubik_create (control)

#ifdef use_PGSLib
  call PGSLib_bcast (prob)
  call PGSLib_bcast (dvalue)
  call PGSLib_bcast (nsizea)
  call PGSLib_bcast (nsizeb)
  call PGSLib_bcast (iscale)
  call PGSLib_bcast (iplot)
  call PGSLib_bcast (Ubik_subspace(control))
  call PGSLib_bcast (Ubik_luerr(control))
  call PGSLib_bcast (Ubik_luout(control))
  call PGSLib_bcast (Ubik_lutty(control))
  call PGSLib_bcast (Ubik_output_frequency(control))
  call PGSLib_bcast (Ubik_method(control))
  call PGSLib_bcast (Ubik_omega(control))
  call PGSLib_bcast (Ubik_stopping(control))
  call PGSLib_bcast (Ubik_norm_type(control))
  call PGSLib_bcast (Ubik_eps(control))
  call PGSLib_bcast (Ubik_itmax(control))
  call PGSLib_bcast (Ubik_residual(control))
  call PGSLib_bcast (Ubik_outmode(control))
  call PGSLib_bcast (precond_control%precond)
  call PGSLib_bcast (precond_control%steps)
  call PGSLib_bcast (precond_control%omega)
#endif

  ! Open output and error file.
  if (Ubik_on_IO_PE()) then
     open (unit=Ubik_luout(control), file='UbikSolve.out', status='replace')
     open (unit=Ubik_luerr(control), file='UbikSolve.err', status='replace')
  end if

  ! Set other parameters.
  select case (prob)
  case (Read_System)
     if (Ubik_on_IO_PE()) then
        open (unit=1, file=TRIM(sysfile), status='old')
        read(1,*) dummy_char, storage_format, nrows
        if (storage_format == Ubik_ELL) then
           read(1,*) dummy_char, maxnz
        end if
     end if
  case (Saad_Full, Saad_ELL)
     nrows = nsizea
     if (prob == Saad_Full) then
        storage_format = Ubik_Full
        maxnz = nsizea
     else
        storage_format = Ubik_ELL
        select case (nsizea/nsizeb)
        case (1)
           maxnz = 3
        case (2)
           maxnz = 4
        case default
           maxnz = 5
        end select
     end if
  end select

#ifdef use_PGSLib
  npe = PGSLib_inquire_npe()
  nrows_pe = nrows / npe

  call PGSLib_bcast (nrows_pe)

  ! Allocate arrays for PGSLib.

  ! Arrays for ELL transpose format.
  allocate( a_trans(maxnz,nrows), stat=status )
  allocate( ja_trans(maxnz,nrows), stat=status )

  ! Arrays for each PE.
  allocate( b_pe(nrows_pe), stat=status )
  allocate( x_pe(nrows_pe), stat=status )
!!$  allocate( a_pe(nrows_pe,maxnz), stat=status )
!!$  allocate( ja_pe(nrows_pe,maxnz), stat=status )
!!$  allocate( b_tmp(nrows_pe,maxnz), stat=status )
  allocate( a_pe(maxnz,nrows_pe), stat=status )
  allocate( ja_pe(maxnz,nrows_pe), stat=status )
  allocate( b_tmp(maxnz,nrows_pe), stat=status )

  ! Other.
  allocate( ja_aux(nrows_pe), stat=status )
  allocate( lengths(npe), stat=status )

  ! Initialization.
  lengths(1:npe) = nrows_pe
  a_pe = 0.0
  ja_pe = 0
#endif

  ! Allocate standard arrays and call driver routine.

  ! Vectors.
  allocate( b(nrows), stat=status )
  allocate( bsave(nrows), stat=status )
  allocate( ref(nrows), stat=status )
  allocate( x(nrows), stat=status )

  ! Coeff. arrays.

  select case (storage_format)
  case (Ubik_Full)

     allocate( coeff_full(nrows,maxnz), stat=status )
     allocate( precond_full(nrows,maxnz), stat=status )
     allocate( coeff_full_save(nrows,maxnz), stat=status )

     call Driver (coeff_full_save, &
          nsizeb, nrows, dvalue, iscale, iplot, &
          bsave, ref, x, b, control, proc_time, elapsed_time)

  case (Ubik_ELL)

     call Ubik_Create (coeff_ell, nrows, maxnz)
     call Ubik_Create (precond_ell, nrows, maxnz)
     call Ubik_Create (coeff_ell_save, nrows, maxnz)

     call Driver (coeff_ell_save, &
          nsizeb, nrows, dvalue, iscale, iplot, &
          bsave, ref, x, b, control, proc_time, elapsed_time)

  end select

  ! Aliases for final values of cpu time, etc.
  iter = Ubik_iter(control)
  cpu = Ubik_cpu(control)
  rnorm = Ubik_rnorm(control)
  err = Ubik_err(control)
  rnormt = Ubik_rnormt(control)
  errt = Ubik_errt(control)

  call Ubik_Spew (Ubik_luout(control), x, 'FINAL SOLUTION (x):')

  if (Ubik_on_IO_PE()) then

     ! If used iterative method, check to see if system got munged.
!!$     if (method >= 0) then
!!$        do j=1,maxnz
!!$           do i=1,nrows
!!$              if (asave(i,j) /= a(i,j)) then
!!$                 write(lutty,*) 'Coefficient element (', i, j, ') has been changed'
!!$                 write(lutty,*) '  was: ', asave(i,j)
!!$                 write(lutty,*) '  now: ', a(i,j)
!!$              end if
!!$           end do
!!$        end do
!!$        if (prob == 2) then
!!$           do j=1,maxnz
!!$              do i=1,nrows
!!$                 if (jasave(i,j) /= ja(i,j)) then
!!$                    write(lutty,*) 'Map element (', i, j, ') has been changed'
!!$                    write(lutty,*) '  was: ', jasave(i,j)
!!$                    write(lutty,*) '  now: ', ja(i,j)
!!$                 end if
!!$              end do
!!$           end do
!!$        end if
!!$        do i=1,nrows
!!$           if (bsave(i) /= b(i)) then
!!$              write(lutty,*) 'Source element (', i, ') has been changed'
!!$              write(lutty,*) '  was: ', bsave(i)
!!$              write(lutty,*) '  now: ', b(i)
!!$           end if
!!$        end do
!!$     end if

     write(lutty,*)
     select case (prob)
     case (Read_System)
        write(lutty,501) 'System read from input file: ', TRIM(sysfile)
     case (Saad_Full, Saad_ELL)
        write(lutty,*) 'Saad Test Problem'
        write(lutty,*)   'Block size          : ', nsizeb
        write(lutty,500) 'Degree of asymmetry : ', dvalue
     end select
     write(lutty,*)  'Number of unknowns  : ', nrows

     select case (storage_format)
     case (Ubik_Full)
        write(lutty,*) 'Storage format      : Full conventional format'
     case (Ubik_ELL)
        write(lutty,*) 'Storage format      : ELLPACK-ITPACK format'
        write(lutty,*) 'Max. non-zeros      : ', maxnz
     end select

     write(lutty,*)
     write(lutty,*) 'Method used: ', method_options(Ubik_method(control))
     select case (Ubik_method(control))
     case (Ubik_method_GMRES_Left_Saad, Ubik_method_FGMRES, Ubik_method_GMRES_Left_Radicati)
        write(lutty,*) '  Krylov subspace size: ', Ubik_subspace(control)
     case (Ubik_method_SOR, Ubik_method_SSOR, Ubik_method_Jacobi)
        write(lutty,402) '  Relaxation parameter: ', Ubik_omega(control)
     end select
     if (iscale >= 1) then
        write(lutty,*) 'System was row-scaled'
        if (iscale == 2) then
           write(lutty,*) 'Jacobi (diagonal) scaling was applied'
        end if
     end if

     if (Ubik_method(control) == Ubik_method_CG .or. &
         Ubik_method(control) == Ubik_method_GMRES_Left_Saad .or. &
         Ubik_method(control) == Ubik_method_GMRES_Left_Radicati .or. &
         Ubik_method(control) == Ubik_method_FGMRES .or. &
         Ubik_method(control) == Ubik_method_BCGS .or. &
         Ubik_method(control) == Ubik_method_BCGSTAB .or. &
         Ubik_method(control) == Ubik_method_TFQMR .or. &
         Ubik_method(control) == Ubik_method_Jacobi .or. &
         Ubik_method(control) == Ubik_method_SOR .or. &
         Ubik_method(control) == Ubik_method_SSOR) then

        if (Ubik_method(control) == Ubik_method_CG .or. &
            Ubik_method(control) == Ubik_method_GMRES_Left_Saad .or. &
            Ubik_method(control) == Ubik_method_GMRES_Left_Radicati .or. &
            Ubik_method(control) == Ubik_method_FGMRES .or. &
            Ubik_method(control) == Ubik_method_BCGS .or. &
            Ubik_method(control) == Ubik_method_BCGSTAB .or. &
            Ubik_method(control) == Ubik_method_TFQMR) then

           write(lutty,*) 'Preconditioner: ', precond_options(precond_control%precond)
           if (precond_control%precond == Ubik_method_Jacobi .or. &
               precond_control%precond == Ubik_method_SSOR) then
              write(lutty,*) '  Number of steps: ', precond_control%steps
           end if
           if (precond_control%precond > 0) then
              write(lutty,402) '  Relaxation parameter: ', precond_control%omega
           end if

           if (Ubik_residual_default(control)) then
              write(lutty,*) 'b - Ax computed and used as needed by method'
           else if (Ubik_residual_compute_only(control) .or.  &
                    Ubik_residual_compute_update(control)) then
              write(lutty,*) 'b - Ax computed for each new iterate'
              if (Ubik_residual_compute_update(control)) then
                 write(lutty,*) 'b - Ax used every sqrt(n) iterations'
              end if
           else if (Ubik_residual_update(control)) then
              write(lutty,*) 'b - Ax computed and used every sqrt(n) iterations'
           end if
        end if

        write(lutty,*)
        write(lutty,500) 'Convergence criterion: ', Ubik_eps(control)
        if (Ubik_stopping_relchg(control)) then
           write(lutty,*) 'Convergence based on || x - xold || / || x ||'
        else if (Ubik_stopping_axb(control)) then
           write(lutty,*) 'Convergence based on || r || / (|| A ||*|| x || + || b ||)'
        else if (Ubik_stopping_b(control)) then
           write(lutty,*) 'Convergence based on || r || / || b ||'
        else if (Ubik_stopping_x(control)) then
           write(lutty,*) 'Convergence based on || r || / || x ||'
        else if (Ubik_stopping_r0(control)) then
           write(lutty,*) 'Convergence based on || r || / || r0 ||'
        else if (Ubik_stopping_r(control)) then
           write(lutty,*) 'Convergence based on || r ||'
        end if
        if (Ubik_stopping_true(control)) write(lutty,*) ' (true residual used in convergence test)'

        select case (Ubik_norm_type(control))
        case (0)
           write(lutty,*) 'Norm used: infinity norm'
        case (1)
           write(lutty,*) 'Norm used: 1-norm'
        case (2)
           write(lutty,*) 'Norm used: 2-norm for vectors,'
           write(lutty,*) '           Frobenius norm for matrices'
        end select

        write(lutty,*)
        write(lutty,*) 'Iterations completed: ', iter
        if (.not.Ubik_converged(control)) then
           write(lutty,*)
           if (Ubik_invalid_input(control)) then
              write(lutty,*) 'Convergence was not achieved due to invalid input'
           else if (Ubik_alloc_failure(control)) then
              write(lutty,*) 'Convergence was not achieved due to not enough memory'
           else if (Ubik_internal_error(control)) then
              write(lutty,*) 'Convergence was not achieved due to internal error'
           else if (Ubik_itmax_exceeded(control)) then
              write(lutty,*) 'Convergence was not achieved due to iteration limit exceeded'
           else if (Ubik_breakdown(control)) then
              write(lutty,*) 'Convergence was not achieved due to algorithm breakdown'
           else if (Ubik_matmul_error(control)) then
              write(lutty,*) 'Convergence was not achieved due to matrix-vector multiplication error'
           else if (Ubik_precond_error(control)) then
              write(lutty,*) 'Convergence was not achieved due to preconditioning error'
           end if
        end if

        write(lutty,*)
        if (Ubik_method(control) /= Ubik_method_Jacobi .and. &
            Ubik_method(control) /= Ubik_method_SOR .and. &
            Ubik_method(control) /= Ubik_method_SSOR) then
           write(lutty,*) 'Norm of residual at termination     : ', rnorm
        end if
        write(lutty,*) 'Approx. rel. error at termination   : ', err
     end if
  end if

  ! If this was a test problem (rather than a system read in from a
  ! file), compute the absolute difference between the computed and
  ! reference solutions.  
  if (prob /= Read_System) then
     ref = ref - x
     call Ubik_Spew (Ubik_luout(control), ref, 'ABSOLUTE DIFFERENCE (ref):')
     diffmax = Ubik_Norm(0, ref)
  end if

  if (Ubik_on_IO_PE()) then
     write(lutty,*)
     write(lutty,500) 'Norm of true residual at term.      : ', rnormt
     if (Ubik_method(control) >= 0) then
        if (Ubik_residual_compute_only(control) .or. &
            Ubik_residual_compute_update(control) .or. &
            Ubik_stopping_true(control)) then
           write(lutty,500) 'True relative error at termination  : ', errt
        end if
     end if

     ! If this was a test problem (rather than a system read in from a
     ! file), write out the max. absolute difference between the computed and
     ! reference solutions.  
     if (prob /= Read_System) then
        write(lutty,*)
        write(lutty,500) 'Max. abs. diff. in solution    : ', diffmax
     end if

     write(lutty,*)
     write(lutty,500) 'CPU time spent in solver       : ', proc_time
     write(lutty,500) 'Elapsed time in solver         : ', elapsed_time
     if (Ubik_matvecs(control) /= 0) then
        write(lutty,500) 'Number of matrix-vector products computed: ', Ubik_matvecs(control)
     end if
     if (Ubik_preconds(control) /= 0) then
        write(lutty,500) 'Number of preconditioner applications    : ', Ubik_preconds(control)
     end if

     ! If an iterative method was used, write per iteration info to file for
     ! plotting.
     if (Ubik_method(control) >= 0) then
        close (unit=Ubik_luout(control))
        open (unit=Ubik_luout(control), file='UbikSolve.dat', status='replace')

        select case (Ubik_method(control))
        case (Ubik_method_Jacobi, Ubik_method_SOR, Ubik_method_SSOR)

           ! Methods that don't generate rnorm (stationary methods).  Remember that
           ! for these methods control%stop <= 0 since Ubik_Stationary checks for
           ! control%stop > 0, which doesn't make sense, and if true sets
           ! control%stop = -control%stop.
           if (Ubik_stopping_true(control)) then

              ! True residual calculated and used in stopping test.
              do i=1,iter

                 cpu = Ubik_cpu (control, i)
                 errt = Ubik_errt (control, i)

                 write(Ubik_luout(control),102) REAL(i), cpu, errt
              end do
           else

              ! Stopping test 0 used and true residual calculated.
              if (Ubik_residual_compute_only(control)) then
                 do i=1,iter

                    cpu = Ubik_cpu (control, i)
                    rnorm = Ubik_rnorm (control, i)
                    rnormt = Ubik_rnormt (control, i)

                    write(Ubik_luout(control),102) REAL(i), cpu, err, rnormt
                 end do
              else

                 ! Stopping test 0 used and true residual not calculated.
                 do i=1,iter

                    cpu = Ubik_cpu (control, i)
                    rnorm = Ubik_rnorm (control, i)

                    write(Ubik_luout(control),102) REAL(i), cpu, err
                 end do
              end if
           end if

        case (Ubik_method_TFQMR)

           ! TFQMR, which generates two iterates per std. iteration.
           if (Ubik_residual_compute_only(control) .or. &
               Ubik_residual_compute_update(control) .or. &
               Ubik_stopping_true(control)) then
              do i=1,iter

                 cpu = Ubik_cpu (control, i)
                 rnorm = Ubik_rnorm (control, i)
                 rnormt = Ubik_rnormt (control, i)

                 write(Ubik_luout(control),102) 0.5*REAL(i), cpu, rnorm, rnormt
              end do
           else
              do i=1,iter

                 cpu = Ubik_cpu (control, i)
                 rnorm = Ubik_rnorm (control, i)

                 write(Ubik_luout(control),102) 0.5*REAL(i), cpu, rnorm
              end do
           end if

        case default

           ! Other methods that do generate rnorm (non-stationary methods).
           if (Ubik_residual_compute_only(control) .or. &
               Ubik_residual_compute_update(control) .or. &
               Ubik_stopping_true(control)) then
              do i=1,iter

                 cpu = Ubik_cpu (control, i)
                 rnorm = Ubik_rnorm (control, i)
                 rnormt = Ubik_rnormt (control, i)

                 write(Ubik_luout(control),102) REAL(i), cpu, rnorm, rnormt
              end do
           else
              do i=1,iter

                 cpu = Ubik_cpu (control, i)
                 rnorm = Ubik_rnorm (control, i)

                 write(Ubik_luout(control),102) REAL(i), cpu, rnorm
              end do
           end if

        end select

        ! Write out GNUPLOT command file.
        close (unit=Ubik_luout(control))
        open (unit=Ubik_luout(control), file='UbikSolve.gpl', status='replace')

        if (Ubik_method(control) >= 0) then
           write(Ubik_luout(control),400)  &
                'set title "Method: '//TRIM(method_options(Ubik_method(control)))// &
                backslash//'n'// &
                'Preconditioner: '//TRIM(precond_options(precond_control%precond))//'"'
        else
           write(Ubik_luout(control),400) &
                'set title "Method: '//TRIM(method_options(Ubik_method(control)))//'"'
        end if
        write(Ubik_luout(control),400) 'set key below'
        write(Ubik_luout(control),400) 'set zero 1.0e-30'
        write(Ubik_luout(control),400) 'set xlabel "Iteration"'
        write(Ubik_luout(control),400) 'set ylabel "Norm"'
        write(Ubik_luout(control),400) 'set logscale y'
        write(Ubik_luout(control),400) 'set ytics nomirror'
        write(Ubik_luout(control),400) 'set format y "%.0e"'
        write(Ubik_luout(control),400) 'set y2label "CPU time"'
        write(Ubik_luout(control),400) 'set y2range [0:*]'
        write(Ubik_luout(control),400) 'set y2tics'
        write(Ubik_luout(control),400) 'set format y2 "%.2f"'
        string0 = ' "UbikSolve.dat" using '
        string1 = '1:3 title'
        string2 = '1:4 title'
        string3 = '1:2 axes x2y2 title "CPU time"'
        string4 = ' with lines'
        title2 = ' "|| b - Ax ||"' 
        write(Ubik_luout(control),400) 'plot '//backslash

        select case (Ubik_method(control))
        case (Ubik_method_Jacobi, Ubik_method_SOR, Ubik_method_SSOR)
!           title1 = ' "'//TRIM(stopping_options(ABS(Ubik_stopping(control))))//'"'
           write(Ubik_luout(control),400) &
                string0//string1//title1//string4//', '//backslash
           if (Ubik_residual_compute_only(control) .and. Ubik_stopping_relchg(control)) then
              write(Ubik_luout(control),400) &
                   string0//string2//title2//string4//', '//backslash
           end if
        case default
           title1 = ' "|| r ||"'
           write(Ubik_luout(control),400) &
                string0//string1//title1//string4//', '//backslash
           if (Ubik_residual_compute_only(control) .or. &
               Ubik_residual_compute_update(control) .or. &
               Ubik_stopping_true(control)) then
              write(Ubik_luout(control),400) &
                   string0//string2//title2//string4//', '//backslash
           end if
        end select

        write(Ubik_luout(control),400) string0//string3//string4
        close (unit=Ubik_luout(control))
        open (unit=Ubik_luout(control), file='UbikSolve.gpl_driver', status='replace')
        write(Ubik_luout(control),400) 'load "UbikSolve.gpl"'
        write(Ubik_luout(control),400) 'pause -1 "Hit return to continue. "' 
        close (unit=Ubik_luout(control))

     end if

     if (iplot == 1) then
        call SYSTEM ('gnuplot '//'UbikSolve.gpl_driver')
     end if

  end if

!  call Ubik_spew (lutty, 'deallocating coeff and precond...')

  ! Deallocate in reverse order.

  select case (storage_format)
  case (Ubik_Full)

     if (allocated(coeff_full_save)) then
        deallocate( coeff_full_save )
     end if
     if (allocated(precond_full)) then
        deallocate( precond_full )
     end if
     if (allocated(coeff_full)) then
        deallocate( coeff_full )
     end if

  case (Ubik_ELL)

     call Ubik_Destroy (coeff_ell)
     call Ubik_Destroy (precond_ell)
#if !defined(lahey) && \
    !defined(NAG_COMPILER) && \
    !(defined(linux) && (defined(dec) || defined(intel)))
     call Ubik_Destroy (coeff_ell_save)
#endif

  end select

!  call Ubik_spew (lutty, 'deallocating control...')

  ! components of control_type
  call Ubik_Destroy (control)

!  call Ubik_spew (lutty, 'deallocating vectors...')

  ! vectors
  if (allocated(x)) then
     deallocate( x )
  end if
  if (allocated(ref)) then
     deallocate( ref )
  end if
  if (allocated(bsave)) then
     deallocate( bsave )
  end if
  if (allocated(b)) then
     deallocate( b )
  end if

#ifdef use_PGSLib
  if (allocated(lengths)) then
     deallocate( lengths )
  end if
  if (allocated(ja_aux)) then
     deallocate( ja_aux )
  end if
  if (allocated(b_tmp)) then
     deallocate( b_tmp )
  end if
  if (allocated(ja_pe)) then
     deallocate( ja_pe )
  end if
  if (allocated(a_pe)) then
     deallocate( a_pe )
  end if
  if (allocated(x_pe)) then
     deallocate( x_pe )
  end if
  if (allocated(b_pe)) then
     deallocate( b_pe )
  end if
  if (allocated(ja_trans)) then
     deallocate( ja_trans )
  end if
  if (allocated(a_trans)) then
     deallocate( a_trans )
  end if

  ! Clean up PGSLib.
  call PGSLib_Finalize ()
#endif

  stop
102 format (4(2x, e15.8))
400 format (a, i4, a, i4, a, f7.4, a)
401 format (a, 10(i4, a), a)
402 format (a, f6.3)
500 format (1x, a, 1p, g10.3)
501 format (1x, 2a)
end program UbikTest
