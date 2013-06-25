! -*- Mode: f90 -*-

!=======================================================================
! $Id: ParallelInfo.F,v 1.5 2003/10/19 22:20:32 turner Exp $
!
! Encapsulate the global variables which contain the parallel
! computing information.  This module must be used on ALL
! systems (parallel or not).  This module contains variables
! which indicate whether the system is parallel or not, etc.
!=======================================================================
module Ubik_ParallelInfo_module
  use ubik_types_module
#ifdef use_PGSLib
    use PGSLib_module
#endif

  implicit none
  save

  ! named constants
  integer(int_type), parameter :: scope_global__ = 0
  integer(int_type), parameter :: scope_local__ = 1
  
  ! type containing parallel info
  type PEinfo_type
     private

     integer(kind=int_type) :: nPE
     integer(kind=int_type) :: thisPE
     integer(kind=int_type) :: IO_PE
     integer(kind=int_type) :: GlobalServicesFlag
  end type PEinfo_type

  ! declare a PEinfo type
  type(PEinfo_type) :: PEinfo

#ifdef use_PGSLib
  ! declare a PGSLib scope flag
  type(PGSLib_SCOPE) :: PGSLib_scope_flag
#endif

  private
  public :: PEinfo_type, &
            PEinfo, &
            on_IO_PE, &
            initialize_parallel, &
            initialize_scope

  ! these shouldn't really be public, but control_type_module needs them
  ! they are not made public outside of Ubik
  ! initially they were internal to control_type_module, made public so
  !   this module could access, but that created circular dependency
  !   since control_type_module also uses WriteString, which uses this
  !   module - so this seemed like the best solution
  public :: &
       scope_local__, &
       scope_global__

#ifdef use_PGSLib
#ifdef aix
  public :: PGSLib_scope_flag, PGSLib_SCOPE
#else
  public :: PGSLib_scope_flag
#endif
#endif

contains

  !=======================================================================
  ! <para>
  ! Initialize parallel parameters, such as the number of PEs, the PE
  ! number of this process, and the I/O PE
  ! </para>
  !=======================================================================
  subroutine initialize_parallel ()

#ifdef use_PGSLib
    ! initialize PGSLib

    ! PGSLib_UseParallelServices may be called before PGSLib_initialize

    PEinfo%GlobalServicesFlag = PGSLib_UseGlobalServices()

    PEinfo%nPE = PGSLib_Inquire_nPE()
    PEinfo%thisPE = PGSLib_Inquire_thisPE()
    PEinfo%IO_PE = PGSLib_Inquire_IO_ROOT_PE()
#else
    ! not using PGSLib

    PEinfo%GlobalServicesFlag = -1

    PEinfo%nPE = PEinfo%GlobalServicesFlag
    PEinfo%thisPE = PEinfo%GlobalServicesFlag
    PEinfo%IO_PE = PEinfo%GlobalServicesFlag
#endif

  end subroutine initialize_parallel

  ! <para>
  ! Indicate whether this process is on the I/O processor.
  ! </para>
  function on_IO_PE ()
    logical(logical_type) :: on_IO_PE

    if (PEinfo%thisPE == PEinfo%IO_PE) then
       on_IO_PE = .TRUE.
    else
       on_IO_PE = .FALSE.
    end if

    return
  end function on_IO_PE

  !=======================================================================
  ! <para>
  ! Initialize PGSLib scope flag.
  ! </para>
  !=======================================================================
  subroutine initialize_scope (scope)
    ! <arguments>
    integer(kind=int_type), intent(IN) :: scope
    ! </arguments>

#ifdef use_PGSLib

    if (scope == scope_local__) then
       PGSLib_scope_flag = PGSLib_Local
    else
       PGSLib_scope_flag = PGSLib_Global
    end if
#endif

  end subroutine initialize_scope

end module Ubik_ParallelInfo_module
