MODULE CONSTANTS_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define commonly-used constants.
  !
  ! Contains: None
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use kind_module, only: int_kind, real_kind

  implicit none
  save

  ! Public Module
  public

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Integer preset value
  integer(int_kind), parameter :: ipreset = -1836547290

  ! Real preset value
  real(real_kind), parameter :: preset = -1836547290.e+17

  ! integers
  integer(int_kind), parameter :: KILO = 1024
  integer(int_kind), parameter :: FLOATBYTES = 8

  ! Real whole numbers
  real(real_kind), parameter :: zero          =   0.0
  real(real_kind), parameter :: one           =   1.0
  real(real_kind), parameter :: two           =   2.0
  real(real_kind), parameter :: three         =   3.0
  real(real_kind), parameter :: four          =   4.0
  real(real_kind), parameter :: five          =   5.0
  real(real_kind), parameter :: six           =   6.0
  real(real_kind), parameter :: seven         =   7.0
  real(real_kind), parameter :: eight         =   8.0
  real(real_kind), parameter :: nine          =   9.0
  real(real_kind), parameter :: ten           =  10.0
  real(real_kind), parameter :: twelve        =  12.0
  real(real_kind), parameter :: twenty_four   =  24.0
  real(real_kind), parameter :: one_hundred   = 100.0

  ! Real exponential numbers
  real(real_kind), parameter :: ten_tominus14 = 1.0e-14
  real(real_kind), parameter :: ten_tominus12 = 1.0e-12
  real(real_kind), parameter :: ten_tominus10 = 1.0e-10
  real(real_kind), parameter :: ten_tominus9  = 1.0e-09
  real(real_kind), parameter :: ten_tominus8  = 1.0e-08
  real(real_kind), parameter :: ten_tominus7  = 1.0e-07
  real(real_kind), parameter :: ten_tominus6  = 1.0e-06
  real(real_kind), parameter :: ten_tominus5  = 1.0e-05
  real(real_kind), parameter :: ten_tominus4  = 1.0e-04
  real(real_kind), parameter :: ten_tominus3  = 1.0e-03
  real(real_kind), parameter :: ten_tominus2  = 1.0e-02
  real(real_kind), parameter :: ten_toplus6   = 1.0e+06
  real(real_kind), parameter :: ten_toplus10  = 1.0e+10

  ! Real fractions
  real(real_kind), parameter :: one_168th     = one/168.0
  real(real_kind), parameter :: one_105th     = one/105.0
  real(real_kind), parameter :: one_sixtieth  = one/60.0
  real(real_kind), parameter :: one_42nd      = one/42.0
  real(real_kind), parameter :: one_thirtieth = one/30.0
  real(real_kind), parameter :: one_twentieth = one/20.0
  real(real_kind), parameter :: one_twelfth   = one/12.0
  real(real_kind), parameter :: one_tenth     = one/ten
  real(real_kind), parameter :: one_eighth    = one/eight
  real(real_kind), parameter :: one_seventh   = one/seven
  real(real_kind), parameter :: one_sixth     = one/six
  real(real_kind), parameter :: one_fifth     = one/five
  real(real_kind), parameter :: one_fourth    = one/four
  real(real_kind), parameter :: one_third     = one/three
  real(real_kind), parameter :: one_half      = one/two
  real(real_kind), parameter :: two_thirds    = two/three
  real(real_kind), parameter :: three_fourths = three/four

  ! Real constants
  real(real_kind), parameter :: pi  = 3.1415926535
  real(real_kind), parameter :: big = ten_toplus10

END MODULE CONSTANTS_MODULE
