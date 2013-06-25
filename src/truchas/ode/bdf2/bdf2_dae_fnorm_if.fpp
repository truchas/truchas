    interface
      subroutine fnorm (t, u, udot, f, error)
        use kinds, only: r8
        real(r8), intent(in) :: t, u(:), udot(:), f(:)
        real(r8), intent(out), optional :: error
      end subroutine fnorm
    end interface
