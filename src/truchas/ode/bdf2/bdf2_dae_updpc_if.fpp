    interface
      subroutine updpc (t, u, h, errc)
        use kinds
        real(r8), intent(in)  :: t, u(:), h
        integer, intent(out) :: errc
      end subroutine updpc
    end interface
