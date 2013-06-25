!!CPP!! This file is to be included by pgslib_io_dist_module.F

!!CPPP!! $Id: dist_scalar.fpp,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file"
#endif

    USE PGSLib_Type_MODULE
    USE PGSLib_Utility_MODULE,     ONLY : pgslib_check_error,      &
         &                                PGSLib_Inquire_IO_P,     &
         &                                PGSLib_Inquire_nPE
    implicit none
    _DATA_TYPE_ , intent(OUT):: scalar_out
    _DATA_TYPE_ , intent(IN), dimension(:):: scalarv_in
    ! Local variables
    integer vec_len
    logical error_flag

    error_flag = .FALSE.
    if (PGSLib_Inquire_IO_P()) then
       vec_len = SIZE(scalarv_in)
       error_flag = (vec_len /= PGSLib_Inquire_nPE())
    endif

    call pgslib_check_error(error_flag,"Wrong sized array in PGSLIB_DIST")

    call pgslib_dist_scalar_c(scalar_out, scalarv_in)
    RETURN

#undef _DATA_TYPE_
