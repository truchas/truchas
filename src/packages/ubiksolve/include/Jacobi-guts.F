! -*- Mode: f90 -*-

    real(real_type), intent(in) :: omega
    real(real_type), intent(in), dimension(:) :: b
    real(real_type), intent(inout), dimension(:) :: x
    ! </arguments>

    ! Local:
    real(real_type), dimension(SIZE(x)) :: diagonal
#ifdef UBIK_DBG
    real(real_type), dimension(SIZE(x)) :: dot
#endif

    diagonal = Diag(a)

#ifdef UBIK_DBG

    ! Form matrix-vector product.
    dot = MatVec(a, x)

    ! Subtract out diagonal contribution.
    dot = dot - diagonal*x

    ! Construct new iterate.
    dot = (b - dot) / diagonal
    x = x + omega*( dot - x )

#else

    x = x + omega*( ((b - (MatVec(a,x) - diagonal*x)) / diagonal) - x )

#endif

    return
