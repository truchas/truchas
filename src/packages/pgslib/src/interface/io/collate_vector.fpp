!!CPP!! This file is to be included by pgslib_io_collate_module.F

!!CPPP!! $Id: collate_vector.fpp,v 1.2 2001/03/22 00:26:13 ferrell Exp $

#ifndef _DATA_TYPE_
#error "_DATA_TYPE_ must be defined before including this file"
#endif

    USE PGSLib_Type_MODULE
    USE pgslib_globals_module,  ONLY : PGSLib_PEInfo
    USE PGSLib_Utility_MODULE,  ONLY : pgslib_check_error,     &
         &                             PGSLib_Inquire_IO_P,    &
         &                             PGSLib_Inquire_nPE,     &
         &                             PGSLib_Output,          &
         &                             PGSLib_Flush_Output

    IMPLICIT NONE
    _DATA_TYPE_, INTENT(  OUT), DIMENSION(:):: Vector_Out
    _DATA_TYPE_, INTENT(IN   ), DIMENSION(:):: Vector_In

    ! Local variables
    LOGICAL error_flag
    INTEGER In_Length, I
    INTEGER, DIMENSION(PGSLib_PEInfo%nPE):: Lengths
    character (LEN=1024) :: out_string

    ! First need to collate the lengths from the PE''s onto the root pe
    In_Length = SIZE(Vector_In)
    Call PGSLib_Collate(Lengths, In_Length)

    ! Now check to make sure that Vector_Out is large enough.

    error_flag = .false.
    IF (PGSlib_Inquire_IO_P()) then
       error_flag = SUM(Lengths) > SIZE(Vector_Out)
    ENDIF
    call pgslib_check_error(error_flag, "Destination array not large enough in PGSLib_Collate")

#ifdef DEBUG_IO
    DO I = 1, MIN(200, SIZE(Vector_in,1))
       write(out_string,*) 'I, Vector_in = ', I, Vector_in(I)
       call pgslib_output(out_string)
    END DO
#endif

    call pgslib_Collate_vector_c(Vector_Out, Lengths, Vector_In, In_Length)

#ifdef DEBUG_IO
    if (PGSLib_Inquire_IO_P()) then
       DO I = 1, MIN(200, SIZE(Vector_out,1))
          write(out_string,*) 'I, Vector_out = ', I, Vector_out(I)
          call pgslib_output(out_string)
       END DO
    endif
    call pgslib_flush_output()
#endif


    RETURN

#undef _DATA_TYPE_
