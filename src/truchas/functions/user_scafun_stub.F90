function user_scafun (index, x, p) result(f)

  use kinds, only: r8

  integer,  intent(in) :: index
  real(r8), intent(in) :: x(:), p(:)
  real(r8) :: f

  write(0,*) 'ERROR: USER_SCAFUN stub procedure was called; you must supply your own.'
  stop

end function user_scafun
