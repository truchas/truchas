!!CPP!! This file is included by pgslib_shift_module.F
!!CPP!! This file provides global EOSHIFT functionality

! $Id: eoshift.fpp,v 1.1.1.1 2000/10/11 22:44:28 ferrell Exp $

  function _ROUTINE_NAME_(Source, SHIFT, BOUNDARY)
    implicit none
    _DATA_TYPE_,       &
         &   intent(IN   ),          &
         &   dimension(           :) :: Source

    integer (PGSLib_Int_Type),       &
         &   intent(IN   )           :: SHIFT

    _DATA_TYPE_,       &
         &   intent(IN   )           :: BOUNDARY

    _DATA_TYPE_,       &
         &   dimension(SIZE(Source,1)) :: _ROUTINE_NAME_

    ! Local variables
    integer (PGSLib_Int_Type),       &
         &   dimension(SIZE(Source,1)) :: Global_Index, Permute_Index
    logical (PGSLib_Log_Type),       &
         &   dimension(SIZE(Source,1)) :: Permute_Mask, Boundary_Mask
    integer (PGSLib_Int_Type) :: N_Tot, n

    ! The shift causes item at index I to move to I - SHIFT.
    ! This can be phrased as:
    !      Dest(I-SHIFT) = Src(I)         or
    !      Dest(I)       = Src(I+SHIFT)
    ! But, that does not address the boundary issue.
    ! All items which are not in the range 1 <= I + SHIFT <= N_Tot
    ! need to be overwritten by the BOUNDARY value.

    ! We need to know the largest index.
    N_Tot = PGSLib_Global_SUM(SIZE(Source,1))

    ! Global_Index counts globall from 1 to size of array
    Global_Index = PGSLib_SUM_PREFIX( (/ (1, n=1, SIZE(Global_Index)) /) )

    ! Permute_Index is the index to permute (scatter) into
    Permute_Index = Global_Index - SHIFT
    ! But, we do not want to scatter anything for values of the index
    ! which are out of range.
    Permute_Mask = (1 <= Permute_Index) .AND. (Permute_Index <= N_Tot)

    ! Now we can do the permutation, which fills up all of dest except
    ! the boundary part
    Call PGSLib_Permute(DEST   = _ROUTINE_NAME_, &
                        SOURCE = Source,         &
                        INDEX  = Permute_Index,  &
                        MASK   = Permute_Mask)

    ! Now we need to fill up the BOUNDARY values for the Dest.
    ! The mask for that is the "complement" of the Permute_Mask, in the
    ! sense that it indicates which elements did not get shifted into

    Boundary_Mask = ((Global_Index + SHIFT) < 1) .OR. (N_Tot < (Global_Index + SHIFT))
    WHERE(Boundary_Mask)
       _ROUTINE_NAME_ = BOUNDARY
    END WHERE

    return
  end function _ROUTINE_NAME_

#undef _DATA_TYPE_
#undef _ROUTINE_NAME_
