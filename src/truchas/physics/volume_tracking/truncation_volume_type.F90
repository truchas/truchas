!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module truncation_volume_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  implicit none
  private

  integer, parameter :: nvf = 4
  real(r8), parameter :: alittle = epsilon(1.0_r8)
  real(r8), parameter :: eps(4) = [1.0_r8, -1.0_r8, 1.0_r8, -1.0_r8]

  type :: truncation_volume_face
    real(r8) :: k(3), lambda, nu, V1234, MUa(nvf), MUi(nvf), x(3,nvf), plane_normal(3)
    integer  :: MUp(nvf)
  contains
    procedure :: init => init_truncation_volume_face
    procedure :: truncate_face
    procedure, private :: truncate_face_2
    procedure, private :: truncate_face_4
    procedure, private :: truncate_face_n
    procedure, private :: Y_function
  end type truncation_volume_face

  type, public :: truncation_volume
    private
    type(truncation_volume_face), allocatable :: face_params(:)
  contains
    procedure :: init => init_truncation_volume
    procedure :: volume
  end type truncation_volume

contains

  subroutine init_truncation_volume(this, nodex, plane_normal)

    use cell_topology

    class(truncation_volume), intent(out) :: this
    real(r8), intent(in) :: nodex(:,:), plane_normal(:)

    integer :: f, nfc
    real(r8) :: fnodex(3,4)

    ! get the number of faces for this cell type
    select case (size(nodex, dim=2))
    case (4) ! tet
      nfc = 4
    case (5,6) ! pyramid, wedge
      nfc = 5
    case (8) ! hex
      nfc = 6
    case default
      call TLS_fatal('unaccounted topology in truncation_volume_type')
    end select

    allocate(this%face_params(nfc))
    do f = 1, size(this%face_params)
      ! Store face vertices in an order that is counterclockwise
      ! relative to a direction that looks from outside the face
      ! into the cell. Treat triangle faces as quadrilaterals
      ! with two nodes in the same location.
      select case (size(nodex, dim=2))
      case (4) ! tet
        fnodex(:,:3) = nodex(:,TET4_FACES(TET4_XFACE(f):TET4_XFACE(f+1)-1))
        fnodex(:,4) = fnodex(:,3)
      case (5) ! pyramid
        fnodex(:,:PYR5_FSIZE(f)) = nodex(:,PYR5_FACES(PYR5_XFACE(f):PYR5_XFACE(f+1)-1))
        if (PYR5_FSIZE(f) < 4) fnodex(:,4) = fnodex(:,3)
      case (6) ! wedge
        fnodex(:,:WED6_FSIZE(f)) = nodex(:,WED6_FACES(WED6_XFACE(f):WED6_XFACE(f+1)-1))
        if (WED6_FSIZE(f) < 4) fnodex(:,4) = fnodex(:,3)
      case (8) ! hex
        fnodex = nodex(:,HEX8_FACES(HEX8_XFACE(f):HEX8_XFACE(f+1)-1))
      end select

      call this%face_params(f)%init(fnodex, plane_normal)
    end do

  end subroutine init_truncation_volume

  ! calculates the truncation volume
  real(r8) function volume(this, plane_rho)

    class(truncation_volume), intent(in) :: this
    real(r8), intent(in) :: plane_rho

    integer :: f

    volume = 0
    do f = 1,size(this%face_params)
      volume = volume + this%face_params(f)%truncate_face(plane_rho)
    end do

  end function volume

  subroutine init_truncation_volume_face(this, nodex, plane_normal)

    use cell_geometry, only: cross_product

    class(truncation_volume_face), intent(out) :: this
    real(r8), intent(in) :: nodex(:,:), plane_normal(:)

    integer  :: j
    real(r8) :: tmp1(3)

    this%plane_normal = plane_normal
    this%X = nodex

    ! Compute K (the Area Vector of the ruled surface)
    this%k = cross_product(this%X(:,3)-this%X(:,1), this%X(:,4)-this%X(:,2))

    ! Compute V1234
    tmp1 = this%X(:,1) - this%X(:,2) + this%X(:,3) - this%X(:,4)
    this%V1234 = sum(tmp1*this%k) / 2

    ! Compute the Mu-i-s.  This is the normal of the interface dotted
    ! with the coordinates of each faces vertex.  Mu-p is the vertex number
    ! of the vertex for each face (1 <= Mu-p <= nvf).  When the Mu-i-s are
    ! in ascending order, it denotes which vertex the interface will
    ! pass through first, second, etc.  The variable Mu-p is the vertex
    ! number of the reordered distances.
    do j = 1,nvf
      this%MUi(j) = dot_product(this%plane_normal, this%X(:,j))
      this%MUp(j) = j
    end do

    this%MUa = this%MUi

    ! Here Nu and Lambda are temporaries used to facilitate ordering the
    ! Mu-i into the Mu-a.  Put the minimum distance (the first vertex
    ! that the interface will pass through) into Mu-p(1) and put its
    ! facial vertex number into Mu-p(1).
    this%nu = 0
    this%lambda = minval(this%MUi)
    do j = 1, nvf
      if (this%MUi(j) == this%lambda .and. this%nu == 0) then
        this%MUp(j) = 1
        this%MUp(1) = j
        this%MUa(j) = this%MUa(1)
        this%MUa(1) = this%MUi(j)
        this%nu = 1
      end if
    end do

    ! Now that the minimum distance is in element 1, order the
    ! other vertices in ascending order using a bubble sort.
    if (this%MUa(3) > this%MUa(4)) then
      this%Lambda    = this%MUp(4)
      this%MUp   (4) = this%MUp(3)
      this%MUp   (3) = this%Lambda

      this%Lambda    = this%MUa(4)
      this%MUa   (4) = this%MUa(3)
      this%MUa   (3) = this%Lambda
    end if

    if (this%MUa(2) > this%MUa(3)) then
      this%Lambda    = this%MUp(3)
      this%MUp   (3) = this%MUp(2)
      this%MUp   (2) = this%Lambda

      this%Lambda    = this%MUa(3)
      this%MUa   (3) = this%MUa(2)
      this%MUa   (2) = this%Lambda
    end if

    if (this%MUa(3) > this%MUa(4)) then
      this%Lambda    = this%MUp   (4)
      this%MUp   (4) = this%MUp   (3)
      this%MUp   (3) = this%Lambda

      this%Lambda    = this%MUa   (4)
      this%MUa   (4) = this%MUa   (3)
      this%MUa   (3) = this%Lambda
    end if

    ! The definition of Lambda is given in Eqn. 11.5 of Zemach-s notes.
    this%lambda = this%MUi(2)*this%MUi(4) - this%MUi(1)*this%MUi(3)

    ! Nu is the face deviation vector dotted with the interface normal.
    ! If a face is a parallelogram, Nu (and B) will be zero.
    ! If a face is not a parallelogram, the magnitude of B measures
    ! the deviation from the parallelogram.
    this%Nu = this%MUi(1) - this%MUi(2) + this%MUi(3) - this%MUi(4)

  end subroutine init_truncation_volume_face

  ! Compute the volume truncated at the current
  ! hex face by the plane given by X*Normal - Ro = 0
  function truncate_face(this, plane_rho) result(Vf)

    class(truncation_volume_face), intent(in) :: this
    real(r8), intent(in) :: plane_rho
    real(r8) :: Vf

    real(r8) :: Y(nvf)

    Vf = 0
    Y = this%Y_function(plane_rho)

    ! get intersection case
    if (this%MUa(1) < plane_rho .and. plane_rho <= this%MUa(2) .and. &
        this%MUa(1) /= this%MUa(2)) then ! case 1
      Vf = this%truncate_face_N(1, plane_rho, Y)
    else if (this%MUa(2) < plane_rho .and. plane_rho <= this%MUa(3)) then ! case 2 & 5
      if (this%MUp(2) == (1 + mod(this%MUp(1)+1, nvf))) then ! case 5
        Vf = this%truncate_face_N(1, plane_rho, Y) + this%truncate_face_N(2, plane_rho, Y)
      else ! case 2
        Vf = this%truncate_face_2(plane_rho, Y)
      end if
    else if (plane_rho > this%MUa(3)) then ! case 3 & 4
      if (plane_rho < this%MUa(4)) then ! case 3
        Vf = - this%truncate_face_N(4, plane_rho, Y)
      end if

      Vf = Vf + this%truncate_face_4(plane_rho) ! case 4
    end if

  end function truncate_face

  ! Compute an expression (for case "2") that is part of
  ! the volume truncated along a hex face by the plane
  ! X*Normal - Ro = 0
  function truncate_face_2(this, plane_rho, Y) result(Q)

    class(truncation_volume_face), intent(in) :: this
    real(r8), intent(in) :: plane_rho, Y(:)
    real(r8) :: Q

    integer :: i, k
    real(r8) :: J1, J2, J3, W1, W2, W3, W4, Z1, Z2, Z3

    ! 2nd predecessor of MU(b).
    k = this%MUp(2)
    Z3 = Y(k) ! Z3 is the Y function of the successor of the A vertex
    ! W1 = the difference between the Mu of the 2nd successor of vertex B and the Mu of vertex A
    i = 1 + mod(k + 1,nvf)
    W1 = this%MUi(i) - this%MUa(1)

    if (abs(W1) > alittle) then
      W1 = 1 / W1
    else
      W1 = 0
    end if

    ! 2nd predecessor of MU(a).
    k = this%MUp(1)
    i = 1 + mod(k + 1,nvf)
    Z1 = Y(k) ! Z1 is the Y function of the A vertex
    Z2 = Y(i) ! Z2 is the Y function of the 2nd successor of the A vertex
    W3 = eps(k)*this%Nu*W1
    Q  = eps(k)*this%V1234*W1*W1 / 2
    ! W2 is the difference between:
    !   the Mu of the 2nd successor of vertex A and
    !   the Mu of vertex B.
    W2 = this%MUi(i) - this%MUa(2)

    if (abs(W2) > alittle) then
      W2 = 1 / W2
    else
      W2 = 0
    end if

    if ((abs(W3) > 1.0e-2) .and. (W3 /= -1.0_r8)) then
      W4 = 1 / W3
      J1 = (1 - log(abs(1 + W3))*W4)*W4
      J2 = (0.5_r8 - J1)*W4
      J3 = (1.0_r8/3.0_r8 - J2)*W4
    else
      J3 = 0.25_r8 - W3/5 + W3*W3/6 - W3*W3*W3/7 + 0.125_r8*W3**4
      J2 = 1.0_r8/3.0_r8 - W3*J3
      J1 = 0.5_r8 - W3*J2
    end if

    Q = Q*(J1*(plane_rho - this%MUa(1))**2 -  &
         2*(this%MUa(2) - this%MUa(1))*&
         (Plane_Rho - this%MUa(1))*J2 &
         + J3*(this%MUa(2) - this%MUa(1))**2) &

         + W1*Z1*(2*plane_rho - this%MUa(1) - this%MUa(2))/6 &
         + (W1*W2*(Z2 - Z3)*(plane_rho - this%MUa(2))**2)/6

  end function truncate_face_2

  ! Compute an expression (for case "4") that is part of
  ! the volume truncated along a hex face by the plane
  ! X*Normal - Ro = 0
  real(r8) function truncate_face_4(this, plane_rho)

    class(truncation_volume_face), intent(in) :: this
    real(r8), intent(in) :: plane_rho

    real(r8) :: tmp(3)

    tmp = 0.25_r8 * sum(this%X, dim=2) - plane_rho*this%plane_normal
    truncate_face_4 = sum(tmp * this%K) / 6

  end function truncate_face_4

  ! Compute an expression (for case "n") that is part of
  ! the volume truncated along a hex face by the plane
  ! X*Normal - Ro = 0
  function truncate_face_n(this, n, plane_rho, Y) result(Q)

    class(truncation_volume_face), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(in) :: plane_rho, Y(:)
    real(r8) :: Q

    integer :: k
    real(r8) :: J1, S, T, W1, W3

    T = plane_rho - this%MUa(n)
    W3 = this%lambda + this%Nu*this%MUa(n)

    if (abs(W3) > alittle) then
      W3 = 1 / W3

      W1 = T*T*W3
      T = T * this%Nu * W3

      if (abs(T) > alittle) then
        S = 1 / T
      else
        S = 0
      end if

      if (T /= -1) then
        J1 = (1 - log(abs(1+T))*S)*S
      else
        J1 = 0
      end if

      if (abs(T) > 1.0e-2) then
        W3 = J1 + (-2.0_r8/3.0_r8 + 2*J1 + (J1 - 0.5_r8) * S) * S
      else
        W3 = 1.0_r8/12.0_r8 + T*(-1.0_r8/30.0_r8 + &
            T*(1.0_r8/60.0_r8 + T*(-1.0_r8/105.0_r8 + T / 168)))
      end if

      k = this%MUp(n)
      Q = eps(k)*(Y(k)*W1/6 + this%V1234*W3*W1*W1 / 2)
    else
      Q = 0
    end if

  end function truncate_face_n

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
  function Y_function(this, plane_rho)

    use cell_geometry, only: cross_product

    class(truncation_volume_face), intent(in) :: this
    real(r8), intent(in) :: plane_rho
    real(r8) :: Y_function(nvf)

    integer  :: i
    real(r8) :: R(3,4), S1(3), S3(3), tmp(3)

    tmp = plane_rho*this%plane_normal
    do i = 1,4
      R(:,i) = this%x(:,i) - tmp
    end do

    S1 = cross_product(R(:,1),R(:,2)) ! S1 = (X1 - Normal.Ro) X (X2 - Normal.Ro)
    S3 = cross_product(R(:,3),R(:,4)) ! S3 = (X3 - Normal.Ro) X (X4 - Normal.Ro)

    !  Y1 = (X4 - Normal.Ro) . (X1 - Normal.Ro) X (X2 - Normal.Ro)
    Y_function(1) = sum(R(:,4)*S1)

    !  Y2 = (X1 - Normal.Ro) . (X2 - Normal.Ro) X (X3 - Normal.Ro)
    !     However, we take advantage of the fact the vector identity:
    !     A . B X C = C . A X B  which means that we don-t need to
    !     compute the S2 vector as above.  We can use the S1 vector:
    !     Y2 = (X3 - Normal.Ro) . (X1 - Normal.Ro) X (X2 - Normal.Ro)
    Y_function(2) = sum(R(:,3)*S1)

    !  Y3 = (X2 - Normal.Ro) . (X3 - Normal.Ro) X (X4 - Normal.Ro)
    Y_function(3) = sum(R(:,2)*S3)

    !  Y4 = (X3 - Normal.Ro) . (X4 - Normal.Ro) X (X1 - Normal.Ro)
    !     However, we take advantage of the fact the vector identity:
    !     A . B X C = C . A X B  which means that we don-t need to
    !     compute the S4 vector as above.  We can use the S3 vector:
    !     Y4 = (X1 - Normal.Ro) . (X3 - Normal.Ro) X (X4 - Normal.Ro)
    Y_function(4) = sum(R(:,1)*S3)

  end function Y_function

end module truncation_volume_type
