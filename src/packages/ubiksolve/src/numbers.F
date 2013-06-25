! -*- Mode: f90 -*-

module ubik_numbers_module

  ! Parameters used in the place of numbers.

  use ubik_types_module, only: real_type

  implicit none
  save

  private :: real_type
  
  ! Numbers 0-9.
  real(real_type), parameter ::  &
     zero=0.0, one=1.0, two=2.0,   three=3.0, four=4.0,  &
     five=5.0, six=6.0, seven=7.0, eight=8.0, nine=9.0

  ! Numbers 10-90 by ten.
  real(real_type), parameter ::  &
     ten=10.0, twenty=20.0, thirty=30.0, forty=40.0, fifty=50.0,  &
     sixty=60.0, seventy=70.0, eighty=80.0, ninety=90.0

  ! Fractions.
  real(real_type), parameter ::  &
     half=one/two, third=one/three, fourth=one/four, fifth=one/five,  &
     sixth=one/six, seventh=one/seven, eighth=one/eight, ninth=one/nine

  ! Forms of pi.
  real(real_type), parameter ::  &
     pi=3.141592653589793238462643383279,  &
     sqrtpi=1.77245385091, pisqr=pi*pi, invpi=one/pi,  &
     twopi=two*pi, threepi=three*pi, fourpi=four*pi

  ! Decimal multipliers.
  real(real_type), parameter ::  &
     deca=1.0e1,  hecto=1.0e2,  kilo=1.0e3, mega=1.0e6,   giga=1.0e9,  &
     tera=1.0e12, pecta=1.0e15, exa=1.0e18, zetta=1.0e21, yotta=1.0e24

  real(real_type), parameter ::  &
     deci=1.0e-1,  centi=1.0e-2,  milli=1.0e-3, micro=1.0e-6,  nano=1.0e-9,  &
     pico=1.0e-12, femto=1.0e-15, atto=1.0e-18, zepto=1.0e-21, yocto=1.0e-24

  ! Synonyms for decimal multipliers.
  real(real_type), parameter ::  &
     hundred=hecto, oneh=hecto, thou=kilo, thousand=kilo, onek=kilo, tenth=deci, mev=mega

  ! Common convergence criteria.
  real(real_type), parameter ::  &
     eps1=1.0e-1, eps2=1.0e-2, eps3=1.0e-3,  &
     eps4=1.0e-4, eps5=1.0e-5, eps6=1.0e-6,  &
     eps7=1.0e-7, eps8=1.0e-8, eps9=1.0e-9,  &
     eps10=1.0e-10, eps11=1.0e-11, eps12=1.0e-12,  &
     eps13=1.0e-13, eps14=1.0e-14, eps15=1.0e-15

end module ubik_numbers_module
