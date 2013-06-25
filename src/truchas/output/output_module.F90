MODULE OUTPUT_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define quantities and procedures for formatted output.
  !
  !   Public Interface:
  !
  !     * call INITIALIZE_IO ()
  !
  !         Open the input file and all output files.
  !      
  !     string = MAKE_FILE_NAME  (string [,number])
  !        return a file name string in a uniform format.
  !
  !     lun = getlun()
  !        return a logical unit number
  !
  !     status = freelun(lun)
  !        free a logical unit number for reuse
  !
  ! Contains: INITIALIZE_IO
  !           MAKE_FILE_NAME
  !           getlun
  !           freelun
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !=======================================================================
  use output_data_module, only: msh_lun,             &
                                tab_lun,             &
                                input_file,          &
                                prefix,              &
                                outdir,              &
                                indir ,              &
                                title,               &
                                odm_b_out, &
                                odm_b_err, &
                                odm_b_aux, &
                                odm_b_int, &
                                out_lun,  &
                                err_lun,  &
                                tty_lun

  use file_utility,       only: MAKE_FILE_NAME
  use tbrook_module,      only: OUT_TBROOK,         &
                                TTY_TBROOK,         &
                                ERR_TBROOK,         &
                                AUX_TBROOK,         &
                                INT_TBROOK,         &
                                OUT_TTY_TBROOK,     &
                                OUT_ERR_TBROOK,     &
                                OUT_ERR_TTY_TBROOK, &
                                ERR_TTY_TBROOK
  implicit none

  private

  public :: OUT_TBROOK,         &
            TTY_TBROOK,         &
            ERR_TBROOK,         &
            AUX_TBROOK,         &
            INT_TBROOK,         &
            OUT_TTY_TBROOK,     &
            OUT_ERR_TBROOK,     &
            OUT_ERR_TTY_TBROOK, &
            ERR_TTY_TBROOK


  ! Public Procedures
  public :: INITIALIZE_IO,     &
            TERMINATE_IO!,      &
!            MAKE_FILE_NAME,    &
!            getlun,            &
!            freelun

  ! Stuff from output_data_module
!  public :: msh_lun, &
!            tab_lun, &
!            odm_b_out, &
!            odm_b_err, &
!            odm_b_aux, &
!            odm_b_int, &
!            odm_b_tty!, &
!            out_lun,  &
!            err_lun,  &
!            tty_lun

  ! Public Variables
  public :: input_file, &
            prefix,     &
            outdir,     &
            indir,      &
            title

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE OM_Set_Output_Files (TB_ID,  OLDNAME, OLD_STYLE_BROOK, OLD_STYLE)
    !
    ! Sets the output mode of standard output files
    ! Output is always sent to the base brook.
    ! If Old_STYLE is defined, output is sent to the old files
    ! i.e. the <root>.out, <root>.err files, etc
    !
    
    use parallel_info_module, only: p_info
    use tbrook_module, only: ASSIGNMENT(=),   &
                             Operator(+),     &
                             Brook,           &
                             BaseBrookAscii,  &
                             TBrook_Destroy,  &
                             TBrook_Set,      &
                             TB
    implicit none
    integer,          intent(in) :: TB_ID
    character(LEN=*), intent(in) :: OLDNAME
    logical,          intent(in) :: OLD_STYLE
    type(Brook), intent(in), optional, target :: OLD_STYLE_BROOK
    type(brook) :: B
    
    integer  :: i

    i = 0
    if ( p_info%IOP) then
       TB(TB_ID)%B = BaseBrookAscii
    else
       call TBrook_Destroy(B, i)
       TB(TB_ID) = B ! a blank brook
    end if

    TB(TB_ID)%TAG = TRIM(OLDNAME)
    
    if ( present(old_style_brook) .and. OLD_STYLE) then
       TB(TB_ID)%B = TB(TB_ID)%B + OLD_STYLE_BROOK
    end if

    ! Old style brooks are always 'ascii'
    i = 0
    call TBrook_Set(TB_ID, FORM="Ascii", iStatus=i)

    return

  END SUBROUTINE OM_SET_OUTPUT_FILES

  SUBROUTINE INITIALIZE_IO ()
    !=======================================================================
    ! Purpose:
    !   open I/O files and initialize various I/O-related quantities
    !=======================================================================
    use tbrook_module, only: TBROOK_SET,        &
                             Operator(+),       &
                             Assignment(=),     &
                             TB,                &
                             OUT_TBROOK,        &
                             OUT_TTY_TBROOK,    &
                             OUT_ERR_TBROOK,    &
                             ERR_TBROOK,        &
                             ERR_TTY_TBROOK,    &
                             AUX_TBROOK,        &
                             TBrook_SetBasefile,&
                             Tbrook_remove_duplicates

    implicit none

    ! Local variables
    integer :: iStatus

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize the Truchas Brooks
    iStatus = 0
    call TBrook_SetBaseFile(FILE=TRIM(MAKE_FILE_NAME('TBrook.xml')), iStatus=iStatus)
    if ( iStatus /= 0 ) Write(*,*) 'Failed to initialize Base TBrook.xml output file.'

    ! assign logical unit numbers for I/O
    ! we have to assume that these work, as we can't punt from this module (circular dependency)

    ! Assume stdout is unit 6.
    tty_lun = 6
    msh_lun = getlun()
    tab_lun = getlun()
    out_lun = getlun()
    err_lun = getlun()
    
    ! Open the stdout file
    call TBrook_Set(odm_b_out, &
                    FILE=TRIM(MAKE_FILE_NAME('out')), &
                    UNIT=out_lun, &
                    FORM='ascii', &
                    iStatus=iStatus)
    call OM_Set_Output_Files(OUT_TBROOK, 'OUT', odm_b_out, .true.)


    ! Open the stderr file
    call TBrook_Set(odm_b_err, &
                    FILE=TRIM(MAKE_FILE_NAME('err')), &
                    UNIT=err_lun, &
                    FORM='ascii', &
                    iStatus=iStatus)
    call OM_Set_Output_Files(ERR_TBROOK, 'ERR',  odm_b_err, .true.)

    ! Set up combination brooks
    TB(OUT_TTY_TBROOK)%B       = TB(TTY_TBROOK)%B + TB(OUT_TBROOK)%B
    TB(OUT_TTY_TBROOK)%TAG     = 'OUTTTY'

    TB(ERR_TTY_TBROOK)%B       = TB(TTY_TBROOK)%B + TB(ERR_TBROOK)%B
    TB(ERR_TTY_TBROOK)%TAG     = 'ERRTTY'

    TB(OUT_ERR_TBROOK)%B       = TB(OUT_TBROOK)%B + TB(ERR_TBROOK)%B
    TB(OUT_ERR_TBROOK)%TAG     = 'OUTERR'

    TB(OUT_ERR_TTY_TBROOK)%B   = TB(TTY_TBROOK)%B + TB(OUT_TBROOK)%B + TB(ERR_TBROOK)%B
    TB(OUT_ERR_TTY_TBROOK)%TAG = 'OUTERRTTY'

    call TBrook_Remove_Duplicates(OUT_ERR_TBROOK,     iStatus=iStatus)
    call TBrook_Remove_Duplicates(OUT_TTY_TBROOK,     iStatus=iStatus)
    call TBrook_Remove_Duplicates(ERR_TTY_TBROOK,     iStatus=iStatus)
    call TBrook_Remove_Duplicates(OUT_ERR_TTY_TBROOK, iStatus=iStatus)


    ! Open the output auxiliary file
    call TBrook_Set(odm_b_aux, &
                    FILE=TRIM(MAKE_FILE_NAME('aux')), &
                    UNIT=getlun(), &
                    FORM='ascii', &
                    iStatus=iStatus)
    call OM_Set_Output_Files(AUX_TBROOK, 'AUX',  odm_b_aux, .true.)

    ! Open the output interface file
    call TBrook_Set(odm_b_int, &
                    FILE=TRIM(MAKE_FILE_NAME('int')), &
                    StandardHeaders=.false., &
                    UNIT=getlun(), &
                    FORM='binary', &
                    iStatus=iStatus)
    ! There is no need to do this since this will only duplicate the data
    ! call OM_Set_Output_Files(INT_TBROOK, 'INT',  odm_b_int, .true.)

    return

  END SUBROUTINE INITIALIZE_IO

  SUBROUTINE TERMINATE_IO()
    !======================================================================
    ! Purpose(s):
    !   Close all the output files/streams.  If any finalize data
    !   needs to be written, output it here.
    !======================================================================
    
    use tbrook_module, only: TBrook_CleanUP,      &
                             TBrook_Close,        &
                             TBrook_Destroy,      &
                             BaseBrook,           &
                             TBrook_CloseXMLTag
    use output_data_module, only: cycle_tag_open
    ! Local variables
    integer :: out_code


    out_code = 0

    ! Close CYCLE tag if it is open
    if ( cycle_tag_open ) call TBrook_CloseXMLTag(BaseBrook, XMLTag="CYCLE", iStatus=out_code)
    cycle_tag_open = .false.

    ! Close all files we are responsible for

#ifdef NAG_COMPILER_WORKAROUND
    call TBrook_Write(odm_b_aux, Variable='NAG_COMPILER_WORKAROUND',Advance=.true., iStatus=out_code)
    call TBrook_Write(odm_b_int, Variable='NAG_COMPILER_WORKAROUND',Advance=.true., iStatus=out_code)
#endif

    call TBrook_Close(odm_b_out, iStatus=out_code)
    call TBrook_Destroy(odm_b_out, iStatus=out_code)

    call TBrook_Close(odm_b_err, iStatus=out_code)
    call TBrook_Destroy(odm_b_err, iStatus=out_code)

    call TBrook_Close(odm_b_aux, iStatus=out_code)
    call TBrook_Destroy(odm_b_aux, iStatus=out_code)

!!  call TBrook_Close(odm_b_gra, iStatus=out_code)
!!  call TBrook_Destroy(odm_b_gra, iStatus=out_code)

    call TBrook_Close(odm_b_int, iStatus=out_code)
    call TBrook_Destroy(odm_b_int, iStatus=out_code)

    ! clean up the brooks module
    call TBrook_Cleanup(iStatus=out_code)

    return
  END SUBROUTINE TERMINATE_IO

  !-----------------------------------------------------------------------------

   FUNCTION GETLUN ()
     use brook_module, only: brook_AssignUnit
     implicit none
     integer :: getLun
     getLun = Brook_AssignUnit()
   END FUNCTION GETLUN
   !----------------------------------------------------------------------------

   FUNCTION FREELUN (lun)
     use brook_module, only: Brook_ReleaseUnit
     implicit none
     integer :: freelun
     integer :: lun
     FreeLun = 0
     call Brook_releaseUnit(iUnit=lun, iStatus=FreeLUN)
     return
   END FUNCTION FREELUN

END MODULE OUTPUT_MODULE
