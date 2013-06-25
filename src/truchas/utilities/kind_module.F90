MODULE KIND_MODULE
  !=============================================================================
  ! Purpose(s):
  !   Define logical, integer, and real kind numbers.
  !
  ! Contains: None
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !
  !=============================================================================
  implicit none
  save

  ! Public Module
  public

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  integer, parameter :: log_kind    = KIND(.true.)
  integer, parameter :: int_kind    = KIND(1)

#if FOUR_BYTES_PER_REAL
  integer, parameter :: real_kind   = KIND(1.00)
#else
  integer, parameter :: real_kind   = KIND(1.0d0)
#endif

  integer, parameter :: single_kind = KIND(1.00)
  integer, parameter :: double_kind = KIND(1.0d0)

END MODULE KIND_MODULE
