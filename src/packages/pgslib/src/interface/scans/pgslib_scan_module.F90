MODULE PGSLib_SCAN_MODULE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Provide PREFIX and SUFFIX operations
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! $Id: !

  ! Dummy version, to stand in until get scans complete.
  ! this provides only no-seg version so far.

  use pgslib_type_module
  use pgslib_globals_module
  use pgslib_utility_module, only : PGSLib_Fatal_ERROR, PGSLib_Scope_Check
  use pgslib_scan_no_seg_module
  use pgslib_scan_seg_bit_module
  use pgslib_scan_seg_module

  implicit none
  PRIVATE
  PUBLIC :: PGSLib_SUM_PREFIX, PGSLib_SUM_SUFFIX
  PUBLIC :: PGSLib_PARITY_PREFIX, PGSLib_PARITY_SUFFIX

end MODULE PGSLib_SCAN_MODULE
