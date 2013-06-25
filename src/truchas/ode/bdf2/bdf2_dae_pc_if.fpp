    interface
      subroutine pc (t, u, udot, f)
        use kinds, only: r8
        real(r8), intent(in) :: t, u(:), udot(:)
        real(r8), intent(inout) :: f(:)
      end subroutine pc
    end interface
