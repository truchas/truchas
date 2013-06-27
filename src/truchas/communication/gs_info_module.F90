MODULE GS_INFO_MODULE
  !======================================================================
  ! PURPOSE
  !  This module contains variables, parameters and types needed
  !  internally in the gather/scatter routines.
  !
  !======================================================================
  use pgslib_module, ONLY: PGSLib_GS_Trace

  SAVE
  private
  public :: COMM_Type, EN, EE, EN_Trace, EE_Trace, EE_All_Ngbr_Trace
  public :: NN_All_Ngbr_Trace
  public :: PGSLib_GS_Trace, EE_Mask_Initialized, El_Nbr_Mask

  type COMM_TYPE
     integer :: Type
  end type COMM_TYPE

  type (COMM_TYPE), Parameter :: EN = COMM_TYPE(1), EE = COMM_TYPE(2)

  type (PGSLib_GS_Trace), POINTER :: EN_Trace, EE_Trace, EE_All_Ngbr_Trace
  type (PGSLib_GS_Trace), POINTER :: NN_All_Ngbr_Trace

  ! Flag to inidicate that some initializaiton has been done.
  logical, pointer, dimension(:,:) :: El_Nbr_MASK
  logical :: EE_MASK_Initialized = .false.

END MODULE GS_INFO_MODULE
