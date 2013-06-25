    interface
      subroutine fun (t, u, udot, f)
        use kinds, only: r8
        real(r8), intent(in) :: t, u(:), udot(:)
        real(r8), intent(out) :: f(:)
      end subroutine fun
    end interface
