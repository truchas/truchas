    interface
      subroutine jac (t, u, dfdu, errc)
        use bdf2_kinds
        real(kind=rk), intent(in)  :: t, u(:)
        real(kind=rk), intent(out) :: dfdu(:,:)
        integer, intent(out) :: errc
      end subroutine jac
    end interface
