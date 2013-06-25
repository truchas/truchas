real function call_f (f, x)
  real, intent(in) :: x
  interface
    real function f(x)
      real, intent(in) :: x
    end function
  end interface
  call_f = f(x)
end function call_f

