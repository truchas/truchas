!!CPP!! This file is included by pgslib_shift_module.F
!!CPP!! This file provides global CSHIFT functionality

! $Id: cshift.fpp,v 1.1.1.1 2000/10/11 22:44:28 ferrell Exp $

  function _ROUTINE_NAME_ (Source, SHIFT)
    implicit none
    _DATA_TYPE_,       &
         &   intent(IN   ),          &
         &   dimension(           :) :: Source

    integer (PGSLib_Int_Type),       &
         &   intent(IN   )           :: SHIFT

    _DATA_TYPE_,       &
         &   dimension(SIZE(Source,1)) :: _ROUTINE_NAME_

    ! Local variables
    integer (PGSLib_Int_Type),       &
         &   dimension(SIZE(Source,1)) :: Global_Index, Permute_Index
    integer (PGSLib_Int_Type) :: N_Tot

    ! Global_Index counts globall from 1 to size of array
    Permute_Index = 1
    Global_Index = PGSLib_SUM_PREFIX(Permute_Index)
    N_TOt = PGSLib_Global_SUM(SIZE(Source,1))

    ! Permute_Index is circular shift of Global_Index
    Permute_Index = MODULO((Global_Index - 1 - SHIFT), N_Tot) + 1

    ! Once we do perumation we are done.
    Call PGSLib_Permute(DEST   = _ROUTINE_NAME_, &
                        SOURCE = Source,         &
                        INDEX  = Permute_Index)

    return
  end function _ROUTINE_NAME_

#undef _DATA_TYPE_
#undef _ROUTINE_NAME_
