!!CPP!! This file is to be included in pgslib_io_dist_module.F

!!CPPP!! $Id: dist_vector.fpp,v 1.2 2001/03/22 00:26:13 ferrell Exp $

    IMPLICIT NONE
    _DATA_TYPE_              , INTENT(OUT), DIMENSION(:):: Vector_Out
    _DATA_TYPE_              , INTENT(IN),  DIMENSION(:):: Vector_In
    INTEGER (PGSLib_Int_Type), INTENT(IN),  OPTIONAL, DIMENSION(:):: Lengths

    ! Local variables
    LOGICAL error_flag
    INTEGER OutLen
    integer, dimension(PGSLib_PEInfo%nPE) :: lengths_local

    error_flag = .false.
    ! Check to make sure that Lengths has enough entries, if provided
    if (PRESENT(Lengths)) then
       IF (PGSLib_Inquire_IO_P()) then
          error_flag = SIZE(Lengths) < PGSLib_Inquire_nPE()
       ENDIF
       call pgslib_check_error(error_flag, "Wrong size Length array in PGSLib_Dist_Int_Vector_F")
       lengths_local = Lengths(1:SIZE(Lengths_Local))
    else
       IF (PGSLib_Inquire_IO_P()) then
!          ALLOCATE(Lengths_local(PGSLib_Inquire_nPE()))
       endif
       call pgslib_collate(lengths_local, size(vector_out))
    endif

    ! Now check that we aren''t trying to wedge too much data into the destinations
    call pgslib_dist(OutLen, lengths_local)
    error_flag = OutLen > SIZE(Vector_Out)

    call pgslib_check_error(error_flag, "PGSLib_Dist_Int_Vector: Vector_Out not large enough")

    call pgslib_dist_vector_c(Vector_Out, OutLen, Vector_In, Lengths_local)

!    IF (.NOT. PRESENT(Lengths)) DEALLOCATE(Lengths_local)
    RETURN

#undef _DATA_TYPE_
