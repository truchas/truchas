    interface
      function enorm (u, du)
        use kinds
        real(r8), intent(in)  :: u(:), du(:)
        real(r8) :: enorm
      end function enorm
    end interface
