    interface
      subroutine updpc (t, u, h, errc)
        use,intrinsic :: iso_fortran_env, only: r8 => real64
        real(r8), intent(in)  :: t, u(:), h
        integer, intent(out) :: errc
      end subroutine updpc
    end interface
