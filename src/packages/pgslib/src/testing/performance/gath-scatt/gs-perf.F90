! Measure performance of gathers and scatters
! $Id: gs-perf.F,v 1.1.1.1 2000/10/11 22:44:25 ferrell Exp $
PROGRAM gs_perf
  use pgslib_module
  use perf_types
  use gather_perf_measure

  implicit none

  ! Set up various indices and gather and scatter times with those indices

  ! Local Variables
  integer(PGSLib_Int_Type) :: NSource, NDest, Ndest_tot, NLoop, NI
  integer(PGSLib_Int_Type) :: k1, k2
  integer(PGSLib_Int_Type), dimension(:), POINTER :: Index
  integer(PGSLib_Int_Type), dimension(:), POINTER :: Temp
  integer, parameter :: Long_Int = SELECTED_INT_KIND(16)
  integer(Long_Int), dimension(:), POINTER :: Temp_Long
  integer :: i
  character (LEN=2048), dimension(12) :: out_string

  type (GS_Perf_Timer_Type) :: Timings


  call pgslib_initialize(IO_PE = 0, FILE_PREFIX='GS-PERF')

  nsource = 20000 
  ndest   = 20000
  ndest_tot = PGSLib_Global_SUM(Ndest)
  NLoop   = 8

  ALLOCATE(Index(NSource))
  ALLOCATE(Temp(NSource))

  out_string=''
  write(out_string,'(a)') '======================================================================'
  call std_out_string(out_string)
  out_string = ''
  write(out_string,*) '   Testing on ', pgslib_inquire_npe(),' processors.'
  call std_out_string(out_string)
  write(out_string,35) nsource, ndest, nloop
35 format(1x,'    Source Size = ', I10,'  Dest Size = ', I10, ' Looping ',i4,' times.')
  call std_out_string(out_string)
  out_string=''
  write(out_string,'(a)') '======================================================================'
  call std_out_string(out_string)

  ! First use identity
  Temp = 1
  Index = PGSLib_SUM_PREFIX(Temp)

  Timings = gather_performance(Index, NDest, NLoop)

  call SPrint_Timings(Out_String, 'Performance for IDENTITY Gather ', Timings)
  call std_out_string(out_string)
  call pgslib_output(out_string)

  call Print_Global_Times()

  ! Now use virtual cube, with CShift direction
  Temp = 1
  Index = PGSLib_SUM_PREFIX(Temp)
  Index = MOD(( (Index-1) + 1),NDest_Tot) + 1


  Timings = gather_performance(Index, NDest, NLoop)

  call SPrint_Timings(out_string, 'Performance for CShift Gather ', Timings)
  call std_out_string(out_string)
  call pgslib_output(out_string)

  call Print_Global_Times()

  ! Shift all data over by one processor.  Should get
  ! graph degre == 1, but all data gets moved
  Temp = 1
  Index = PGSLib_SUM_PREFIX(Temp)
  Index = MOD(( (Index-1) + NDest),NDest_Tot) + 1
  ! Now reverse, to test hash insert algorithm
!!$  Temp = Index
!!$  Index = Temp( (/ (NDest - i + 1, i = 1, NDest) /) )

  Timings = gather_performance(Index, NDest, NLoop)

  call SPrint_Timings(out_string, 'Performance for CShift of all data ', Timings)
  call std_out_string(out_string)
  call pgslib_output(out_string)

  call Print_Global_Times()

  ! Now use a random pattern
  ! This is for psuedo random numbers
  k1 = 1000003
  k2 = 26147
  ! Need to use a long integer to avoid overflow
  ALLOCATE(Temp_Long(SIZE(Temp)))
  Temp = 1
  Temp_Long = PGSLib_SUM_PREFIX(Temp) - 1
  Index = MOD(Temp_Long*MOD(K1, NDest_Tot) + INT(k2, Long_Int), INT(NDest_Tot, Long_Int)) + 1

  Timings = gather_performance(Index, NDest, NLoop)

  call SPrint_Timings(out_string, 'Performance for RANDOM Gather ', Timings)
  call std_out_string(out_string)
  call pgslib_output(out_string)

  call Print_Global_Times()


  call pgslib_finalize()

  stop
end PROGRAM gs_perf

subroutine SPrint_Timings(Out_String, Message, Timings)
  use perf_types
  implicit none
  ! Print timing info from the structure into a string
  character (LEN=2048), intent(OUT), dimension(12):: Out_String
  character (LEN=*), intent(IN)                   :: Message
  type (GS_Perf_Timer_Type), intent(IN) :: Timings

  out_string = ''
  write(out_string, 10) Message,                                 &
                        SetupTime(Timings), GatherDegree(Timings),&
                        TotalTime(Timings), TotalData(Timings),  &
                        NetBW(Timings),     NetBWNorm(Timings),  &
                        OffPETime(Timings), offPEData(Timings),  &
                        offPEBW(Timings),   offPEBWNorm(Timings)
10 FORMAT(1x, '***', a, /, '  Times for Gathers',/,&
                                   '   SetupTime = ', g14.4,' Graph Degree = ',  i14, / ,&
                                   '   TotalTime = ', g14.4,'    TotalData = ', g14.4,   &
                                   '       NetBW = ', g14.4,'    NetBWNorm = ', g14.4,/, &
                                   '   OffPETime = ', g14.4,'    OffPEData = ', g14.4,   &
                                   '     OffPEBW = ', g14.4,'  OffPEBWNorm = ', g14.4)

  return
end subroutine SPrint_Timings

subroutine std_out_string(String)
  use pgslib_module
  implicit none
  character (LEN=2048), intent(IN), dimension(12) :: string

  integer :: n

  do n = 1, SIZE(String)
     if (len_trim(string(n)) <= 0) exit
     if (PGSLib_Inquire_IO_P()) write(*,*) TRIM(String(n))
  end do
  return

end subroutine std_out_string

subroutine Print_Global_Times()
  use pgslib_module
  implicit none
  
  character (LEN=2048), dimension(12) :: out_string

  out_string = ''
  write(out_string, 15) 'Total Trace Setup Time', PGSLib_Read_Maximum_Time(SETUP_TRACE_STATISTICS())
  call std_out_string(out_string)
  out_string = ''
  write(out_string, 15) 'Total Supplement Time', PGSLib_Read_Maximum_Time(SETUP_Supplement_STATISTICS())
  call std_out_string(out_string)
  out_string = ''
  write(out_string, 15) 'Supplement G_TO_L Time', PGSLib_Read_Maximum_Time(Supplement_G_TO_L_Statistics())
  call std_out_string(out_string)
  out_string = ''
  write(out_string, 15) 'Add Item To Table Time ', PGSLib_Read_Maximum_Time(Add_Item_To_T_STATISTICS())
  call std_out_string(out_string)
  out_string = ''
  write(out_string, 15) 'Supplement L_From_G Time', PGSLib_Read_Maximum_Time(Supplement_L_From_G_Statistics())
  call std_out_string(out_string)
  out_string = ''
  write(out_string, 15) 'Item Index From Table Time ', PGSLib_Read_Maximum_Time(Item_Index_From_T_STATISTICS())
  call std_out_string(out_string)
  out_string = ''
  write(out_string, 15) 'Total Buffer Setup Time', PGSLib_Read_Maximum_Time(SETUP_Buffers_STATISTICS())
  call std_out_string(out_string)
  out_string = ''
  write(out_string, 15) 'Total Gather Time', PGSLib_Read_Maximum_Time(GATHER_STATISTICS())
  call std_out_string(out_string)
  out_string = ''
  write(out_string, 15) 'Total offPE Gather Time', PGSLib_Read_Maximum_Time(GATHER_Buffer_STATISTICS())
  call std_out_string(out_string)
  out_string = ''
  write(out_string, 15) 'Total Barrier Time', PGSLib_Barrier_Time()
  call std_out_string(out_string)
  out_string = ''
  write(out_string, 15) 'Total Send-Rcv Time', PGSLib_Sr_Time()
  call std_out_string(out_string)

15 FORMAT(1x, '  ', a30,' = ', g20.6)

  return
end subroutine Print_Global_Times


