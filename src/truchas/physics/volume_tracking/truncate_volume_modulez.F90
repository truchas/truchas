!=======================================================================
! Purpose:
!
!   Define procedures necessary to compute hexahedral volumes
!   truncated by a planar interface
!
!   Public Interface:
!
!     * call FACE_PARAM (option, face, K, Lambda, MUa, MUi, MUp, Nu, V1234, X)
!         Compute and store various face parameters needed for a volume
!         truncation calculation.
!
!     * call TRUNCATE_VOLUME (K, Lambda, MUa, MUi, MUp, Nu, Vf, Vol, V1234, X)
!         Compute the volume truncated by a plane.
!
!     * call TRUNCATE_FACE (f, K, Lambda, MUa, MUi, MUp, Nu, Vf, V1234, X)
!         Compute the volume truncated at the current hex face by a plane.
!
! Author(s): Douglas B. Kothe, LANL (dbk@lanl.gov)
!            S. Jay Mosso, LANL (sjm@lanl.gov)
!
!=======================================================================

module truncate_volume_modulez
  use kinds,  only: r8
  use truchas_logging_services
  implicit none
  private

  public :: truncate_volume, truncate_face, face_param, truncvol_data

  integer, parameter :: nvf=4, nfc=6
  real(r8), parameter :: alittle=epsilon(1.0_r8)
  real(r8), parameter :: eps(4) = (/ 1.0_r8, -1.0_r8, 1.0_r8, -1.0_r8 /)

  type truncvol_data
    real(r8) :: K(3)
    real(r8) :: Lambda, Nu, V1234
    real(r8) :: MUa(nvf), MUi(nvf)
    real(r8) :: X(3,nvf)
    integer  :: MUp(nvf)
  end type truncvol_data

contains

  ! Compute and store various face parameters needed for the volume truncation calculation
  type(truncvol_data) function face_param (cell, option, face, int_flux_vol_coord)
    use hex_types,     only: reconstruction_hex
    use cell_topology
    use cell_geometry, only: cross_product

    class(reconstruction_hex), intent(in) :: cell
    character(9),       intent(in) :: option
    integer,            intent(in) :: face
    real(r8), optional, intent(in) :: int_flux_vol_coord(:,:)

    integer  :: i, j, n
    real(r8) :: Tmp1(3)

    ! Store face vertices in an order that is
    !    counterclockwise relative to a direction
    !    that looks from outside the face into
    !    the cell.
    ! For option = 'full_cell', use the donor cell vertices (Cell_Coord)
    ! For option = 'flux_cell', use the vertices of the flux volume (Xv_flux)
    if (option == 'full_cell') then
      face_param%X = cell%node(:,HEX8_FACES(HEX8_XFACE(face):HEX8_XFACE(face+1)-1))
    else if (option == 'flux_cell') then
      if (.not.present(int_flux_vol_coord)) call TLS_fatal('when calculating flux_cell face_param, int_flux_vol_coord is required')
      face_param%X = int_flux_vol_coord(:,HEX8_FACES(HEX8_XFACE(face):HEX8_XFACE(face+1)-1))
    end if

    ! Compute K   (the Area Vector of the ruled surface)
    face_param%K = cross_product (face_param%X(:,3)-face_param%X(:,1), face_param%X(:,4)-face_param%X(:,2))

    ! Compute V1234
    Tmp1(:) = face_param%X(:,1) - face_param%X(:,2) &
         +    face_param%X(:,3) - face_param%X(:,4)
    face_param%V1234 = 0.5_r8 * sum(Tmp1*face_param%K)

    ! Compute the Mu-i-s.  This is the normal of the interface dotted
    ! with the coordinates of each faces vertex.  Mu-p is the vertex number
    ! of the vertex for each face (1 <= Mu-p <= nvf).  When the Mu-i-s are
    ! in ascending order, it denotes which vertex the interface will
    ! pass through first, second, etc.  The variable Mu-p is the vertex
    ! number of the reordered distances.

    ! do j = 1,nvf
    !   face_param%MUi(j) = sum(cell%P%normal(:)*face_param%X(:,j))
    ! end do

    ! !$omp simd
    ! do j = 1,nvf
    !   face_param%MUp(j) = j
    ! end do
    ! !$omp end simd

    do j = 1,nvf
      face_param%MUi(j) = sum(cell%P%normal(:)*face_param%X(:,j))
      face_param%MUp(j) = j
    end do

    face_param%MUa(:) = face_param%MUi(:)

    ! Here Nu and Lambda are temporaries used to facilitate ordering the
    ! Mu-i into the Mu-a.  Put the minimum distance (the first vertex
    ! that the interface will pass through) into Mu-p(1) and put its
    ! facial vertex number into Mu-p(1).
    face_param%Nu = 0.0_r8
    face_param%Lambda = minval(face_param%MUi(1:4))
    do j = 1, nvf
      if (face_param%MUi(j) == face_param%Lambda .and. face_param%Nu == 0.0_r8) then
        face_param%MUp(j) = 1
        face_param%MUp(1) = j
        face_param%MUa(j) = face_param%MUa(1)
        face_param%MUa(1) = face_param%MUi(j)
        face_param%Nu = 1.0_r8
      end if
    end do

    ! Now that the minimum distance is in element 1, order the
    ! other vertices in ascending order using a bubble sort.
    if (face_param%MUa(3) > face_param%MUa(4)) then
      face_param%Lambda    = face_param%MUp(4)
      face_param%MUp   (4) = face_param%MUp(3)
      face_param%MUp   (3) = face_param%Lambda

      face_param%Lambda    = face_param%MUa(4)
      face_param%MUa   (4) = face_param%MUa(3)
      face_param%MUa   (3) = face_param%Lambda
    end if

    if (face_param%MUa(2) > face_param%MUa(3)) then
      face_param%Lambda    = face_param%MUp(3)
      face_param%MUp   (3) = face_param%MUp(2)
      face_param%MUp   (2) = face_param%Lambda

      face_param%Lambda    = face_param%MUa(3)
      face_param%MUa   (3) = face_param%MUa(2)
      face_param%MUa   (2) = face_param%Lambda
    end if

    if (face_param%MUa(3) > face_param%MUa(4)) then
      face_param%Lambda    = face_param%MUp   (4)
      face_param%MUp   (4) = face_param%MUp   (3)
      face_param%MUp   (3) = face_param%Lambda

      face_param%Lambda    = face_param%MUa   (4)
      face_param%MUa   (4) = face_param%MUa   (3)
      face_param%MUa   (3) = face_param%Lambda
    end if

    ! The definition of Lambda is given in Eqn. 11.5 of Zemach-s notes.
    face_param%Lambda = face_param%MUi(2)*face_param%MUi(4) &
         -              face_param%MUi(1)*face_param%MUi(3)

    ! Nu is the face deviation vector dotted with the interface normal.
    ! If a face is a parallelogram, Nu (and B) will be zero.
    ! If a face is not a parallelogram, the magnitude of B measures
    ! the deviation from the parallelogram.
    face_param%Nu = face_param%MUi(1) - face_param%MUi(2) &
         +          face_param%MUi(3) - face_param%MUi(4)
  end function face_param

  ! Compute the truncation volume.
  real(r8) function truncate_volume (cell, trunc_vol)
    use hex_types, only: reconstruction_hex

    class(reconstruction_hex), intent(in) :: cell
    type(truncvol_data),       intent(in) :: trunc_vol(:)

    integer :: f

    ! Initialize relevant quantities
    !if (count_cases) Cases = 0

    ! Loop over faces, accumulating the truncated volume
    Truncate_Volume = 0.0_r8
    do f = 1,nfc
      Truncate_Volume = Truncate_Volume + truncate_face (cell, trunc_vol(f))
    end do

  end function truncate_volume

  ! Compute the volume truncated at the current
  ! hex face by the plane given by X*Normal - Ro = 0
  function truncate_face (cell, trunc_vol_face) result(Vf)
    use hex_types, only: reconstruction_hex

    class(reconstruction_hex), intent(in) :: cell
    type(truncvol_data),       intent(in) :: trunc_vol_face
    real(r8)                              :: Vf

    real(r8) :: Y(nvf)

    Vf = 0.0_r8
    Y = Y_function (cell, trunc_vol_face%X)

    ! get intersection case
    if ( trunc_vol_face%MUa(1) < cell%P%rho .and. cell%P%rho <= trunc_vol_face%MUa(2) .and. &
         trunc_vol_face%MUa(1) /= trunc_vol_face%MUa(2)) then                                ! case 1
      Vf = truncate_face_N (1, cell%P%rho, Y, trunc_vol_face)
      !write(*,*) 'case 1'
    else if ( trunc_vol_face%MUa(2) < cell%P%rho .and. cell%P%rho <= trunc_vol_face%MUa(3)) then ! case 2 & 5
      if (trunc_vol_face%MUp(2) == (1 + mod(trunc_vol_face%MUp(1)+1, nvf))) then ! case 5
        Vf = truncate_face_N (1, cell%P%rho, Y, trunc_vol_face) + &
             truncate_face_N (2, cell%P%rho, Y, trunc_vol_face)
        !write(*,*) 'case 5'
      else                                                                       ! case 2
        Vf = truncate_face_2 (cell%P%rho, Y, trunc_vol_face)
        !write(*,*) 'case 2'
      end if
    else if ( cell%P%rho > trunc_vol_face%MUa(3)) then                                         ! case 3 & 4
      if (cell%P%rho < trunc_vol_face%MUa(4)) then          ! case 3
        Vf = - truncate_face_N (4, cell%P%rho, Y, trunc_vol_face)
        !write(*,*) 'case 3'
      end if

      Vf = Vf + truncate_face_4 (cell, trunc_vol_face) ! case 4
      !write(*,*) 'case 4'
    end if

  end function truncate_face

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><>

  !=======================================================================
  ! PURPOSE -
  !   Compute an expression (for case "2") that is part of
  !   the volume truncated along a hex face by the plane
  !   X*Normal - Ro = 0
  !=======================================================================
  function truncate_face_2 (cell_rho, Y, trunc_vol_face) result(Q)
    real(r8),            intent(in) :: cell_rho, Y(:)
    type(truncvol_data), intent(in) :: trunc_vol_face
    real(r8)                        :: Q

    integer :: i, k
    real(r8) :: J1, J2, J3, W1, W2, W3, W4, Z1, Z2, Z3

    ! 2nd predecessor of MU(b).
    k = trunc_vol_face%MUp(2)
    Z3 = Y(k) ! Z3 is the Y function of the successor of the A vertex
    ! W1 = the difference between the Mu of the 2nd successor of vertex B and the Mu of vertex A
    i = 1 + mod(k + 1,nvf)
    W1 = trunc_vol_face%MUi(i) - trunc_vol_face%MUa(1)

    if (ABS(W1) > alittle) then
      W1 = 1.0_r8 / W1
    else
      W1 = 0.0_r8
    end if

    ! 2nd predecessor of MU(a).
    k = trunc_vol_face%MUp(1)
    i = 1 + mod(k + 1,nvf)
    Z1 = Y(k) ! Z1 is the Y function of the A vertex
    Z2 = Y(i) ! Z2 is the Y function of the 2nd successor of the A vertex
    W3 = eps(k)*trunc_vol_face%Nu*W1
    Q  = eps(k)*trunc_vol_face%V1234*0.5_r8*W1*W1
    ! W2 is the difference between:
    !   the Mu of the 2nd successor of vertex A and
    !   the Mu of vertex B.
    W2 = trunc_vol_face%MUi(i) - trunc_vol_face%MUa(2)

    if (abs(W2) > alittle) then
      W2 = 1.0_r8 / W2
    else
      W2 = 0.0_r8
    end if

    if ((abs(W3) > 1.0e-2) .and. (W3 /= -1.0_r8)) then
      W4 = 1.0_r8/W3
      J1 = (1.0_r8 - log(abs(1.0_r8 + W3))*W4)*W4
      J2 = (0.5_r8 - J1)*W4
      J3 = (1.0_r8/3.0_r8 - J2)*W4
    else
      J3 = 0.25_r8 - W3/5.0_r8 + W3*W3/6.0_r8 - W3*W3*W3/7.0_r8 + 0.125_r8*W3**4
      J2 = 1.0_r8/3.0_r8 - W3*J3
      J1 = 0.5_r8 - W3*J2
    end if

    Q = Q*(J1*(Cell_Rho - trunc_vol_face%MUa(1))**2 -  &
         2.0_r8*(trunc_vol_face%MUa(2) - trunc_vol_face%MUa(1))*&
         (Cell_Rho - trunc_vol_face%MUa(1))*J2 &
         + J3*(trunc_vol_face%MUa(2) - trunc_vol_face%MUa(1))**2) &

         + W1*Z1*(2.0_r8*Cell_Rho - trunc_vol_face%MUa(1) - trunc_vol_face%MUa(2))/6.0_r8 &
         + (W1*W2*(Z2 - Z3)*(Cell_Rho - trunc_vol_face%MUa(2))**2)/6.0_r8
  end function truncate_face_2

  !=======================================================================
  ! PURPOSE -
  !   Compute an expression (for case "4") that is part of
  !   the volume truncated along a hex face by the plane
  !   X*Normal - Ro = 0
  !=======================================================================
  function truncate_face_4 (cell, trunc_vol_face) result(Q)
    use hex_types, only: reconstruction_hex

    class(reconstruction_hex), intent(in) :: cell
    type(truncvol_data),       intent(in) :: trunc_vol_face
    real(r8)                              :: Q

    real(r8) :: Tmp(3)

    tmp(:) = 0.25_r8 * sum(trunc_vol_face%X(:,:), dim=2) - cell%P%rho*cell%P%normal(:)
    Q = sum( tmp * trunc_vol_face%K ) / 6.0_r8

  end function truncate_face_4

  !=======================================================================
  ! PURPOSE -
  !   Compute an expression (for case "n") that is part of
  !   the volume truncated along a hex face by the plane
  !   X*Normal - Ro = 0
  !=======================================================================
  function truncate_face_n (n, plane_const, Y, trunc_vol_face) result(Q)
    integer,  intent(in) :: n
    real(r8), intent(in) :: plane_const
    real(r8), intent(in) :: Y(:)
    type(truncvol_data), intent(in) :: trunc_vol_face
    real(r8)                        :: Q

    integer  :: k
    real(r8) :: J1, S, T, W1, W3

    T = plane_const - trunc_vol_face%MUa(n)
    W3 = trunc_vol_face%Lambda + trunc_vol_face%Nu*trunc_vol_face%MUa(n)

    if (abs(W3) > alittle) then
      W3 = 1.0_r8/W3

      W1 = T*T*W3
      T = T * trunc_vol_face%Nu * W3

      if (abs(T) > alittle) then
        S = 1.0_r8/T
      else
        S = 0.0_r8
      end if

      if (T /= -1.0_r8) then
        J1 = (1.0_r8-log(abs(1.0_r8+T))*S)*S
      else
        J1 = 0.0_r8
      end if

      if (abs(T) > 1.0e-2) then
        W3 = J1 + (-2.0_r8/3.0_r8 + 2.0_r8*J1 + (J1 - 0.5_r8) * S) * S
      else
        W3 = 1.0_r8/12.0_r8 + T*(-1.0_r8/30.0_r8 + T*(1.0_r8/60.0_r8 + T*(-1.0_r8/105.0_r8 + T*1.0_r8/168.0_r8)))
      end if

      k = trunc_vol_face%MUp(n)
      Q = eps(k)*(Y(k)*W1/6.0_r8 + 0.5_r8*trunc_vol_face%V1234*W3*W1*W1)
    else
      Q = 0.0_r8
    end if

  end function truncate_face_n

  !=======================================================================
  ! PURPOSE -
  !   THIS ROUTINE COMPUTES THE VARIABLE :
  !   Yi = (Xi'''-Normal.Ro) . (Xi-Normal.Ro) X (Xi'-Normal.Ro)
  !
  ! Note:
  !   Xi     - is the i-th vertex of face "FACE".
  !   Xi'''  - is the predecessor vertex (third successor vertex) to
  !            vertex Xi
  !   Xi'    - is the successor vertex to vertex Xi
  !   Normal - is the normal vector to the material interface
  !   Ro     - is the distance of vertex Xi from the origin along
  !            the interface normal "NORMAL"
  !   X      - the coordinates of the four vertices that define
  !            call face "FACE".
  !   Y      - the value of the Y-FUNCTION for this value of RO.
  !
  ! This routine is called from truncate_face.
  !=======================================================================
  function Y_function (cell, x)
    use hex_types, only: reconstruction_hex
    use cell_geometry, only: cross_product

    class(reconstruction_hex), intent(in) :: cell
    real(r8),                  intent(in) :: x(:,:)
    real(r8)                              :: Y_FUNCTION(nvf)

    integer  :: i
    real(r8) :: R(3,4), S1(3), S3(3), tmp(3)

    tmp = cell%P%rho*cell%P%normal
    ! !$omp simd
    do i = 1,4
      R(:,i) = x(:,i) - tmp
    end do
    ! !$omp end simd

    S1 = cross_product (R(:,1),R(:,2)) ! S1 = (X1 - Normal.Ro) X (X2 - Normal.Ro)
    S3 = cross_product (R(:,3),R(:,4)) ! S3 = (X3 - Normal.Ro) X (X4 - Normal.Ro)

    !  Y1 = (X4 - Normal.Ro) . (X1 - Normal.Ro) X (X2 - Normal.Ro)
    Y_function(1) = sum( R(:,4)*S1 )

    !  Y2 = (X1 - Normal.Ro) . (X2 - Normal.Ro) X (X3 - Normal.Ro)
    !     However, we take advantage of the fact the vector identity:
    !     A . B X C = C . A X B  which means that we don-t need to
    !     compute the S2 vector as above.  We can use the S1 vector:
    !     Y2 = (X3 - Normal.Ro) . (X1 - Normal.Ro) X (X2 - Normal.Ro)
    Y_function(2) = sum( R(:,3)*S1 )

    !  Y3 = (X2 - Normal.Ro) . (X3 - Normal.Ro) X (X4 - Normal.Ro)
    Y_function(3) = sum( R(:,2)*S3 )

    !  Y4 = (X3 - Normal.Ro) . (X4 - Normal.Ro) X (X1 - Normal.Ro)
    !     However, we take advantage of the fact the vector identity:
    !     A . B X C = C . A X B  which means that we don-t need to
    !     compute the S4 vector as above.  We can use the S3 vector:
    !     Y4 = (X1 - Normal.Ro) . (X3 - Normal.Ro) X (X4 - Normal.Ro)
    Y_function(4) = sum( R(:,1)*S3 )

  end function Y_function

end module truncate_volume_modulez
