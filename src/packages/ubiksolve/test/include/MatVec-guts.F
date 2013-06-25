! -*- Mode: f90 -*-
    
    ! Args:
    type(Ubik_vector_type), intent(in) :: x_vec
    integer(Ubik_int_type), intent(out) :: status
    real(Ubik_real_type), intent(out), dimension(:) :: y

    ! Local:
    real(Ubik_real_type), dimension(:), pointer :: x

    x => Ubik_values_ptr(x_vec)
    y = Ubik_MatVec(coeff, x)

    status = 0
