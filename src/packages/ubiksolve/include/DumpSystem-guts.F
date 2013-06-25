! -*- Mode: f90 -*-

    integer(int_type), intent(in) :: lu
    real(real_type), intent(in), dimension(:) :: x, b
    ! </arguments>

    ! Local:
    character(len=200) :: out_string

    ! Coefficient.
    call Dump (lu, a)

    ! Divider.
    write(out_string,101) 0, 0, 0.0, 0, 0
    call Spew (lu, out_string)

    ! Unknown and source vectors.
    call Dump (lu, x, b)

    return
101 format (2(i8, 1x), 1pg15.8, 2(1x, i8))
