    interface
      subroutine pc (t, u, udot, f)
        use,intrinsic :: iso_fortran_env, only: r8 => real64
        real(r8), intent(in) :: t, u(:), udot(:)
        real(r8), intent(inout) :: f(:)
      end subroutine pc
    end interface
