MODULE PGSLib_IO_MODULE

  ! $Id: pgslib_io_module.F,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $
  USE PGSLib_IO_BCast_MODULE,   ONLY : PGSLib_Bcast
  USE PGSLib_IO_Dist_MODULE,    ONLY : PGSLib_Dist
  USE PGSLib_IO_Collate_MODULE, ONLY : PGSLib_Collate
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC :: PGSLib_BCast
  PUBLIC :: PGSLib_Dist
  PUBLIC :: PGSLib_Collate
END MODULE PGSLib_IO_MODULE
