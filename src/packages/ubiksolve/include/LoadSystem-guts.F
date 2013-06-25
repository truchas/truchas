! -*- Mode: f90 -*-

    integer(int_type), intent(in) :: lu
    real(real_type), intent(out), dimension(:) :: x, b
    ! </arguments>

    ! Local:
    character(len=1) :: discard

    ! Coefficient.
    call Load (lu, a)

    ! Divider.
    read(lu,'(a)') discard

    ! Unknown and source vectors.
    call Load (lu, x, b)

    return
