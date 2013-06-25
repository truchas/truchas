    ! -*- Mode: f90 -*-

    real(Ubik_real_type), intent(in), dimension(:) :: b
    real(Ubik_real_type), intent(inout), dimension(:) :: x
    integer(Ubik_int_type), intent(out) :: status

    ! Initialize return status.
    status = 0

    select case (iterate_control%method)
    case (Ubik_method_Jacobi)

       ! Jacobi with relaxation
       call Ubik_Jacobi (iterate_control%omega, b, coeff, x)
    case (Ubik_method_SOR)

       ! GS with relaxation
       call Ubik_SOR (Ubik_method_SOR, iterate_control%omega, b, coeff, x)
    case (Ubik_method_SSOR)

       ! symmetric GS with relaxation
       call Ubik_SOR (Ubik_method_SSOR, iterate_control%omega, b, coeff, x)
    end select

    return
