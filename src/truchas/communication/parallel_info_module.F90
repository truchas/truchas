MODULE PARALLEL_INFO_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Encapsulate the global variables which contain the parallel
  !   computing information.  This module must be used on ALL
  !   systems (parallel or not).  This module contains variables
  !   which indicate whether the system is parallel or not.
  !
  !   Also contains variables specifying parallel computation
  !   model, how to lay out a structured mesh, how many domains
  !   to use.
  !
  !   type PEInfo contains the following components:
  !
  !      nPE                  total number of PEs (1 for serial runs)
  !      thisPE               integer in range [1,nPE]
  !      IO_ROOT_PE           PE which is allowed to do IO
  !      GlobalServicesFlag   Integer which indicates whether system
  !                           provides global services or not.  If not
  !                           then all global ops must be done with
  !                           PGSLib.
  !      UseGlobalServices    .true. if system provides global services
  !      IOP                  .true. on any PE which is allowed to 
  !                           perform I/O operations.
  !      IsParallel          .true. if run is parallel, .false. if serial
  !
  ! Author(s): Robert Ferrell (CPCA, Ltd., ferrell@cpca.com)
  !
  !=======================================================================
  implicit none
  private 

  ! Public types
  public :: PEINFO

  ! Public Variables
  public :: p_info
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  
  ! Define PEINFO structure
  type PEINFO
     integer :: nPE
     integer :: thisPE
     integer :: IO_ROOT_PE
     integer :: GlobalServicesFlag
     logical :: UseGlobalServices
     logical :: IOP
     logical :: IsParallel
  end type PEINFO

  ! Declare PEINFO structure
  type(PEINFO), save :: p_info

END MODULE PARALLEL_INFO_MODULE
