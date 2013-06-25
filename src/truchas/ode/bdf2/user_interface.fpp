  interface
    subroutine user (t, u)
      use bdf2_kinds
      real(kind=rk), intent(in) :: t, u(:)
    end subroutine user
  end interface
