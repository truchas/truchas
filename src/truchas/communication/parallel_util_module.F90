MODULE PARALLEL_UTIL_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !    Define parallel utility procedures.
  !
  ! Public Interface(s):
  !
  !   * call PARALLEL_INIT ()
  !
  !      Initialize the p_info derived type.
  !
  ! Contains: PARALLEL_INIT
  !
  ! Author(s): Robert Ferrell (CPCA, Ltd., ferrell@cpca.com)
  !
  !=======================================================================
  use parallel_info_module
  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: PARALLEL_INIT, Is_IO_PE, p_info

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE PARALLEL_INIT (argv)
    !=======================================================================
    ! Purpose(s):
    !
    !      Initialize parallel parameters, such as the number of PEs,
    !      the PE number of this process, and the I/O PE.
    !
    !=======================================================================
    use parallel_info_module, only: p_info
#ifdef USE_PGSLIB
    use pgslib_module,        only: PGSLib_INQUIRE_IO_ROOT_PE,        &
                                    PGSLib_INQUIRE_THISPE_ACTUAL,     &
                                    PGSLib_INQUIRE_NPE,               &
                                    PGSLib_INQUIRE_IO_P,              &
                                    PGSLib_USEGLOBALSERVICES,         &
                                    PGSLib_INITIALIZE,                &
                                    PGSLib_CL_MAX_TOKEN_LENGTH
#endif

    implicit none

#ifdef USE_PGSLIB
    character(len=PGSLib_CL_MAX_TOKEN_LENGTH), dimension(:), pointer :: argv
#endif

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
    ! Preset the PEInfo structure
#ifndef USE_PGSLIB
    p_info%GlobalServicesFlag = -1
    p_info%nPE                =  1
    p_info%thisPE             = p_info%GlobalServicesFlag
    p_info%IO_ROOT_PE         = p_info%GlobalServicesFlag
    p_info%UseGlobalServices  = .true.
    p_info%IOP                = .true.
    p_info%IsParallel         = .false.
#endif

    ! Initialize PGSLIB
#ifdef USE_PGSLIB
    ! PGSLib_UseParallelServices may be called before pgslib_initialize.
    p_info%GlobalServicesFlag = PGSLib_USEGLOBALSERVICES()

    call PGSLib_INITIALIZE (0, argv=argv)

    ! Set various PE numbers
    p_info%nPE               = PGSLib_INQUIRE_NPE()
    p_info%thisPE            = PGSLib_INQUIRE_THISPE_ACTUAL()
    p_info%IO_ROOT_PE        = PGSLib_INQUIRE_IO_ROOT_PE()
    p_info%UseGlobalServices = .false.
    p_info%IOP               = PGSLib_INQUIRE_IO_P()
    p_info%IsParallel        = .true.

#endif

    return

  END SUBROUTINE PARALLEL_INIT
        
  Function Is_IO_PE()
    !=======================================================================
    ! Purpose(s):
    !
    !      Returns .TRUE. for the IO_PE, .FALSE. otherwise.
    !
    !=======================================================================
    use parallel_info_module, only: p_info
    LOGICAL :: Is_IO_PE
    
    Is_IO_PE = p_info%IOP
    return
  end Function Is_IO_PE
  

END MODULE PARALLEL_UTIL_MODULE
