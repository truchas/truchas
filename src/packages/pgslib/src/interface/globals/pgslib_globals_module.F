MODULE PGSLib_Globals_MODULE
  USE PGSLib_Type_MODULE
  implicit none
  save

  PUBLIC :: PGSLib_PEInfo

  ! Global variables
  ! Flag to determine whether we can use global services or not
  INTEGER (PGSLib_Int_Type), Parameter:: PGSLib_UseGlobalServicesFlag = -1

  ! Global PE information
  TYPE (PGSLib_PEInfo_Struct) PGSLib_PEInfo

  ! Constant to specify whether to use PGSLib data mapping, or a user provided array
  integer (PGslib_Int_Type), Public, Parameter :: PGSLib_User_Map    = 1
  integer (PGslib_Int_Type), Public, Parameter :: PGSLib_Default_Map = 2

  ! Constants to determine scope of some of the operations
  TYPE (PGSLib_SCOPE), Public, Parameter :: PGSLib_GLOBAL = PGSLib_SCOPE(1)
  TYPE (PGSLib_SCOPE), Public, Parameter :: PGSLib_LOCAL  = PGSLib_SCOPE(2)

END MODULE PGSLib_Globals_MODULE
