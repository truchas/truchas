!!CPP!! This file is to be included by pgslib_instrument.F
!!CPP!! _FUNC_ must be defined before including this file

!!CPP!! $Id: stats_name.fpp,v 1.1.1.1 2000/10/11 22:44:31 ferrell Exp $

!!CPP!! define the variable names based on the function name
!!CPP!! The contortions are to get _F_NAME_, _SLOT_ & _LOCAL to expand
!!CPP!! _FUNC_ before catentation.
#define _ID_(a) a
#define _CAT_(a,b) _ID_(a)_ID_(b)
#define _F_NAME_NAME_(f) _CAT_(f,_STATISTICS)
#define _F_NAME_ _F_NAME_NAME_(_FUNC_)
#define _SLOT_NAME_(f) _CAT_(f,_Slot)
#define _SLOT_ _SLOT_NAME_(_FUNC_)
#define _LOCAL_NAME_(f) _CAT_(f,_L)
#define _LOCAL_ _LOCAL_NAME_(_FUNC_)


function _F_NAME_()
  implicit none
  type (PGSLib_Instrument_T), POINTER :: _F_NAME_
  ! Local variables
  type (PGSLib_Instrument_Slot), SAVE :: _SLOT_ = PGSLib_Instrument_Slot(-1)
  type (PGSLib_Instrument_T), POINTER,  SAVE :: _LOCAL_

  if (_SLOT_%Slot < 1) then
     _LOCAL_ => Next_Instrument_Item()
     _SLOT_%Slot = _LOCAL_%Slot%SLot
  end if
  _F_NAME_ => _LOCAL_
  return
end function _F_NAME_

#undef _CAT_
#undef _F_NAME_NAME_
#undef _F_NAME_
#undef _SLOT_NAME_
#undef _SLOT_
#undef _LOCAL_
#undef _LOCAL_NAME_
#undef _FUNC_
