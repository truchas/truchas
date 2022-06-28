    interface
      subroutine pcfun (t, u, udot, f)
        use,intrinsic :: iso_fortran_env, only: r8 => real64
        real(r8), intent(in)  :: t, u(:), udot(:)
        real(r8), intent(out) :: f(:)
      end subroutine pcfun
    end interface
