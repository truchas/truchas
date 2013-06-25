! -*- Mode: f90 -*-

module ubik_numbers_complete_module

  ! Parameters used in the place of numbers.

  use ubik_types_module, only: real_type
  use ubik_numbers_module

  implicit none
  save

  private :: real_type

  ! Numbers 11-19.
  real(real_type), parameter ::  &
     eleven=11.0,  twelve=12.0,    thirteen=13.0, fourteen=14.0,  &
     fifteen=15.0, sixteen=16.0, seventeen=17.0, eighteen=18.0, nineteen=19.0

  ! Numbers 21-29.
  real(real_type), parameter ::  &
     twenty1=21.0, twenty2=22.0, twenty3=23.0, twenty4=24.0,  &
     twenty5=25.0, twenty6=26.0, twenty7=27.0, twenty8=28.0, twenty9=29.0

  ! Numbers 31-39.
  real(real_type), parameter ::  &
     thirty1=31.0, thirty2=32.0, thirty3=33.0, thirty4=34.0,  &
     thirty5=35.0, thirty6=36.0, thirty7=37.0, thirty8=38.0, thirty9=39.0

  ! Numbers 41-49.
  real(real_type), parameter ::  &
     forty1=41.0, forty2=42.0, forty3=43.0, forty4=44.0,  &
     forty5=45.0, forty6=46.0, forty7=47.0, forty8=48.0, forty9=49.0

  ! Numbers 51-59.
  real(real_type), parameter ::  &
     fifty1=51.0, fifty2=52.0, fifty3=53.0, fifty4=54.0,  &
     fifty5=55.0, fifty6=56.0, fifty7=57.0, fifty8=58.0, fifty9=59.0

  ! Numbers 61-69.
  real(real_type), parameter ::  &
     sixty1=61.0, sixty2=62.0, sixty3=63.0, sixty4=64.0,  &
     sixty5=65.0, sixty6=66.0, sixty7=67.0, sixty8=68.0, sixty9=69.0

  ! Numbers 71-79.
  real(real_type), parameter ::  &
     seventy1=71.0, seventy2=72.0, seventy3=73.0, seventy4=74.0,  &
     seventy5=75.0, seventy6=76.0, seventy7=77.0, seventy8=78.0, seventy9=79.0

  ! Numbers 81-89.
  real(real_type), parameter ::  &
     eighty1=81.0, eighty2=82.0, eighty3=83.0, eighty4=84.0,  &
     eighty5=85.0, eighty6=86.0, eighty7=87.0, eighty8=88.0, eighty9=89.0

  ! Numbers 91-99.
  real(real_type), parameter ::  &
     ninety1=91.0, ninety2=92.0, ninety3=93.0, ninety4=94.0,  &
     ninety5=95.0, ninety6=96.0, ninety7=97.0, ninety8=98.0, ninety9=99.0

  ! Fractions with 2 in denominator.
  real(real_type), parameter ::  &
     three2=three/two, five2=five/two, seven2=seven/two, nine2=nine/two

  ! Fractions with 3 in denominator.
  real(real_type), parameter ::  &
     two3=two/three,   four3=four/three,  &
     five3=five/three, seven3=seven/three, eight3=eight/three

  ! Fractions with 4 in denominator.
  real(real_type), parameter ::  &
     three4=three/four, five4=five/four, seven4=seven/four

  ! Fractions with 5 in denominator.
  real(real_type), parameter ::  &
     two5=two/five, three5=three/five, four5=four/five,  &
     six5=six/five, seven5=seven/five, eight5=eight/five, nine5=nine/five

  ! Fractions with 6 in denominator.
  real(real_type), parameter ::  &
     five6=five/six, seven6=seven/six, nine6=nine/six

  ! Fractions with 7 in denominator.
  real(real_type), parameter ::  &
     two7=two/seven,   three7=three/seven,  &
     four7=four/seven, five7=five/seven, six7=six/seven

  ! Fractions with 8 in denominator.
  real(real_type), parameter ::  &
     three8=three/eight, five8=five/eight,  &
     seven8=seven/eight, nine8=nine/eight

  ! Fractions with 9 in denominator.
  real(real_type), parameter ::  &
     two9=two/nine,   four9=four/nine,  &
     five9=five/nine, seven9=seven/nine, eight9=eight/nine

  ! Forms of pi.
  real(real_type), parameter ::  &
     pi2=pi/two, pi3=pi/three, pi4=pi/four,  pi5=pi/five,  &
     pi6=pi/six, pi7=pi/seven, pi8=pi/eight, pi9=pi/nine

end module ubik_numbers_complete_module
