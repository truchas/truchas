function call_dll_scafun (dll_scafun, x, p) result(f)

  use kinds, only: r8

  real(r8), intent(in) :: x(*), p(*)
  real(r8) :: f

  interface
    function dll_scafun(x, p) result(f)
      use kinds, only: r8
      real(r8), intent(in) :: x(*), p(*)
      real(r8) :: f
    end function
  end interface

  f = dll_scafun(x, p)

end function call_dll_scafun
