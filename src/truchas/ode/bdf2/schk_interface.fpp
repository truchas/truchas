    interface
      subroutine schk (u, stage, errc)
        use bdf2_kinds
        real(kind=rk), intent(in) :: u(:)
        integer, intent(in)  :: stage
        integer, intent(out) :: errc
      end subroutine schk
    end interface
