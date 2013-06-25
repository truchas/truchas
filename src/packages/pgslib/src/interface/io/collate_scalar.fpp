!!CPP!! This file is to be included by pgslib_io_collate_module.F

!!CPPP!! $Id: collate_scalar.fpp,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file"
#endif

    USE PGSLib_Type_MODULE
    USE PGSLib_Utility_MODULE,  ONLY : pgslib_check_error,     &
         &                             PGSLib_Inquire_IO_P,    &
         &                             PGSLib_Inquire_nPE,     &
         &                             PGSLib_Output,          &
         &                             PGSLib_Flush_Output

    implicit none
    _DATA_TYPE_, intent(OUT), dimension(:):: scalarv_out
    _DATA_TYPE_, intent(IN):: scalar_in
    ! Local variables
    integer vec_len, I
    logical error_flag
    character (LEN=1024) :: out_string

    ! Check that destination is large enough.
    error_flag = .FALSE.
    if (PGSLib_Inquire_IO_P()) then
       vec_len = SIZE(scalarv_out)
       error_flag = (vec_len .LT. PGSLib_Inquire_nPE())
    endif
    call pgslib_check_error(error_flag, "Destination array not large enough in PGSLib_collate for integer scalar.")

#ifdef DEBUG_IO
    write(out_string,*) 'Scalar_in = ', scalar_in
    call pgslib_output(out_string)
#endif

    call pgslib_Collate_scalar_c(scalarv_out, scalar_in)

#ifdef DEBUG_IO
    if (PGSLib_Inquire_IO_P()) then
       DO I = 1, MIN(200, SIZE(Scalarv_out,1))
          write(out_string,*) 'I, Scalarv_out = ', I, Scalarv_out(I)
          call pgslib_output(out_string)
       END DO
    endif
    call pgslib_flush_output()
#endif

    RETURN

#undef _DATA_TYPE_
