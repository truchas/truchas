module gaussian_quadrature_vofinit

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  implicit none

contains

  subroutine quadrature_qua4(np, gp_coord, gp_weight)

    integer, intent(in) :: np
    real(r8), intent(out) :: gp_coord(:,:), gp_weight(:)

    real(r8) :: gp, gp0, gp1, gp2, gp3, w0, w1, w2, w3

    gp_coord = 0.0_r8
    gp_weight = 0.0_r8

    select case (np)
    case (1)
      gp_coord(1,1) = 0.0_r8
      gp_coord(2,1) = 0.0_r8

      gp_weight(1) = 4.0_r8

    case (4)
      gp = 1.0_r8/dsqrt(3.0_r8)

      gp_coord(1,1) = gp
      gp_coord(2,1) = gp
      gp_coord(1,2) = gp
      gp_coord(2,2) = -gp
      gp_coord(1,3) = -gp
      gp_coord(2,3) = gp
      gp_coord(1,4) = -gp
      gp_coord(2,4) = -gp

      gp_weight(1) = 1.0_r8
      gp_weight(2) = 1.0_r8
      gp_weight(3) = 1.0_r8
      gp_weight(4) = 1.0_r8

    case (9)
      gp0 = dsqrt(3.0_r8/5.0_r8)
      gp1 = 0.0_r8
      gp2 = -dsqrt(3.0_r8/5.0_r8)

      w0 = 5.0_r8/9.0_r8
      w1 = 8.0_r8/9.0_r8
      w2 = 5.0_r8/9.0_r8

      gp_coord(1,1) = gp0
      gp_coord(2,1) = gp0
      gp_coord(1,2) = gp0
      gp_coord(2,2) = gp1
      gp_coord(1,3) = gp0
      gp_coord(2,3) = gp2
      gp_coord(1,4) = gp1
      gp_coord(2,4) = gp0
      gp_coord(1,5) = gp1
      gp_coord(2,5) = gp1
      gp_coord(1,6) = gp1
      gp_coord(2,6) = gp2
      gp_coord(1,7) = gp2
      gp_coord(2,7) = gp0
      gp_coord(1,8) = gp2
      gp_coord(2,8) = gp1
      gp_coord(1,9) = gp2
      gp_coord(2,9) = gp2

      gp_weight(1) = w0*w0
      gp_weight(2) = w0*w1
      gp_weight(3) = w0*w2
      gp_weight(4) = w1*w0
      gp_weight(5) = w1*w1
      gp_weight(6) = w1*w2
      gp_weight(7) = w2*w0
      gp_weight(8) = w2*w1
      gp_weight(9) = w2*w2

    case (16)
      gp0 =  0.8611363115940526_r8
      gp1 =  0.3399810435848563_r8
      gp2 = -0.3399810435848563_r8
      gp3 = -0.8611363115940526_r8

      w0 = 0.3478548451374538_r8
      w1 = 0.6521451548625461_r8
      w2 = 0.6521451548625461_r8
      w3 = 0.3478548451374538_r8

      gp_coord(1,1) = gp0
      gp_coord(2,1) = gp0
      gp_coord(1,2) = gp0
      gp_coord(2,2) = gp1
      gp_coord(1,3) = gp0
      gp_coord(2,3) = gp2
      gp_coord(1,4) = gp0
      gp_coord(2,4) = gp3
      gp_coord(1,5) = gp1
      gp_coord(2,5) = gp0
      gp_coord(1,6) = gp1
      gp_coord(2,6) = gp1
      gp_coord(1,7) = gp1
      gp_coord(2,7) = gp2
      gp_coord(1,8) = gp1
      gp_coord(2,8) = gp3
      gp_coord(1,9) = gp2
      gp_coord(2,9) = gp0
      gp_coord(1,10) = gp2
      gp_coord(2,10) = gp1
      gp_coord(1,11) = gp2
      gp_coord(2,11) = gp2
      gp_coord(1,12) = gp2
      gp_coord(2,12) = gp3
      gp_coord(1,13) = gp3
      gp_coord(2,13) = gp0
      gp_coord(1,14) = gp3
      gp_coord(2,14) = gp1
      gp_coord(1,15) = gp3
      gp_coord(2,15) = gp2
      gp_coord(1,16) = gp3
      gp_coord(2,16) = gp3

      gp_weight(1) = w0*w0
      gp_weight(2) = w0*w1
      gp_weight(3) = w0*w2
      gp_weight(4) = w0*w3
      gp_weight(5) = w1*w0
      gp_weight(6) = w1*w1
      gp_weight(7) = w1*w2
      gp_weight(8) = w1*w3
      gp_weight(9) = w2*w0
      gp_weight(10) = w2*w1
      gp_weight(11) = w2*w2
      gp_weight(12) = w2*w3
      gp_weight(13) = w3*w0
      gp_weight(14) = w3*w1
      gp_weight(15) = w3*w2
      gp_weight(16) = w3*w3

    case default
      call TLS_panic('incorrect number of quadrature points for quad4')
      stop

    end select

    gp_weight(:) = gp_weight(:) / 4.0_r8

  end subroutine quadrature_qua4

  subroutine transform_qua4(nodes, gp_coord, coord)

    real(r8), intent(in) :: nodes(2,4), gp_coord(2)
    real(r8), intent(out) :: coord(2)

    coord = 0.0_r8

    coord(1) = nodes(1,1)*(1.0_r8-gp_coord(1))*(1.0_r8-gp_coord(2)) &
      + nodes(1,2)*(1.0_r8+gp_coord(1))*(1.0_r8-gp_coord(2)) &
      + nodes(1,3)*(1.0_r8+gp_coord(1))*(1.0_r8+gp_coord(2)) &
      + nodes(1,4)*(1.0_r8-gp_coord(1))*(1.0_r8+gp_coord(2))

    coord(2) = nodes(2,1)*(1.0_r8-gp_coord(1))*(1.0_r8-gp_coord(2)) &
      + nodes(2,2)*(1.0_r8+gp_coord(1))*(1.0_r8-gp_coord(2)) &
      + nodes(2,3)*(1.0_r8+gp_coord(1))*(1.0_r8+gp_coord(2)) &
      + nodes(2,4)*(1.0_r8-gp_coord(1))*(1.0_r8+gp_coord(2))

    coord = coord * 0.25_r8

  end subroutine transform_qua4

end module gaussian_quadrature_vofinit
