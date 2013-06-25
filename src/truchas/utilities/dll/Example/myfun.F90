real function myfun1 (x)
  real, intent(in) :: x
  myfun1 = x**2
end function myfun1

real function myfun2 (x)
  real, intent(in) :: x
  myfun2 = cos(x)**2
end function myfun2
