    interface
      subroutine pcfun (t, u, udot, f)
        use kinds
        real(r8), intent(in)  :: t, u(:), udot(:)
        real(r8), intent(out) :: f(:)
      end subroutine pcfun
    end interface
