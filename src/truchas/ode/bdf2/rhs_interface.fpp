    interface
      subroutine rhs (t, u, f)
        use bdf2_kinds
        real(kind=rk), intent(in)  :: t, u(:)
        real(kind=rk), intent(out) :: f(:)
      end subroutine rhs
    end interface
