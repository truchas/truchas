    interface
      function enorm (u, du)
        use,intrinsic :: iso_fortran_env, only: r8 => real64
        real(r8), intent(in)  :: u(:), du(:)
        real(r8) :: enorm
      end function enorm
    end interface
