! -*- Mode: f90 -*-

module ubik_types_module
  implicit none
  save

  integer, parameter :: logical_type = KIND(.true.)
  integer, parameter :: int_type = KIND(1)
#ifdef unicos
  integer, parameter :: real_type = KIND(1.0)
#else
  integer, parameter :: real_type = KIND(1.0d0)
#endif

end module ubik_types_module
