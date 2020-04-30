!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE TRUNCATE_VOLUME_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures necessary to compute hexahedral volumes
  !   truncated by a planar interface
  !
  !   Public Interface:
  !
  !     * call FACE_PARAM (option, face, K, Lambda, MUa, MUi, MUp, Nu, V1234, X)
  !
  !         Compute and store various face parameters needed for a volume
  !         truncation calculation.
  !
  !     * call TRUNCATE_VOLUME (K, Lambda, MUa, MUi, MUp, Nu, Vf, Vol, V1234, X)
  !
  !         Compute the volume truncated by a plane.
  !
  !     * call TRUNCATE_FACE (f, K, Lambda, MUa, MUi, MUp, Nu, Vf, V1234, X)
  !
  !         Compute the volume truncated at the current hex face by a plane.
  !
  ! Contains: FACE_PARAM
  !           TRUNCATE_VOLUME
  !           TRUNCATE_FACE
  !           TRUNCATE_FACE_2
  !           TRUNCATE_FACE_2
  !           TRUNCATE_FACE_4
  !           TRUNCATE_FACE_N
  !           Y_FUNCTION
  !
  ! Author(s): Douglas B. Kothe, LANL (dbk@lanl.gov)
  !            S. Jay Mosso, LANL (sjm@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use legacy_mesh_api, only: ndim, nvf
  implicit none
  private

  public :: TRUNCATE_VOLUME, TRUNCATE_FACE, FACE_PARAM, Trunc_Vol, TRUNCVOL_DATA

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  type TRUNCVOL_DATA

    real(r8) :: K(ndim)
    real(r8) :: Lambda
    real(r8) :: MUa(nvf)
    real(r8) :: MUi(nvf)
    real(r8) :: Nu
    real(r8) :: V1234
    real(r8) :: X(nvf,ndim)
    integer  :: MUp(nvf)

  end type TRUNCVOL_DATA

  ! Declare an array of TRUNCVOL_DATA types
  type(TRUNCVOL_DATA), dimension(:,:), pointer :: Trunc_Vol

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE FACE_PARAM (option, face)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute and store various face parameters needed
    !   for the volume truncation calculation
    !
    !=======================================================================
    use interface_module, only: Int_Geom, Int_Flux
    use parameter_module, only: nicells

    ! Arguments
    character(9), intent(IN) :: option
    integer, intent(IN) :: face

    ! Local Variables
    integer :: i, j, n, v1, v2, v3, v4
    real(r8), dimension(nicells,ndim) :: Tmp1, Tmp2

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Store face vertices in an order that is
    !    counterclockwise relative to a direction
    !    that looks from outside the face into
    !    the cell.
    ! For option = 'full_cell', use the donor cell vertices (Cell_Coord)
    ! For option = 'flux_cell', use the vertices of the flux volume (Xv_flux)
    select case (face)
    case(1)
       v1 = 4; v2 = 8; v3 = 7; v4 = 3      ! Left face
    case(2)
       v1 = 5; v2 = 1; v3 = 2; v4 = 6      ! Right face
    case(3)
       v1 = 5; v2 = 8; v3 = 4; v4 = 1      ! Front face
    case(4)
       v1 = 6; v2 = 2; v3 = 3; v4 = 7      ! Back face
    case(5)
       v1 = 3; v2 = 2; v3 = 1; v4 = 4      ! Bottom face
    case(6)
       v1 = 7; v2 = 8; v3 = 5; v4 = 6      ! Top face
    end select

    if (option == 'full_cell') then
       do n=1,ndim
          Trunc_Vol(:,face)%X(1,n) = Int_Geom%Cell_Coord(n,v1)
          Trunc_Vol(:,face)%X(2,n) = Int_Geom%Cell_Coord(n,v2)
          Trunc_Vol(:,face)%X(3,n) = Int_Geom%Cell_Coord(n,v3)
          Trunc_Vol(:,face)%X(4,n) = Int_Geom%Cell_Coord(n,v4)
       end do
    else if (option == 'flux_cell') then
       do n=1,ndim
          Trunc_Vol(:,face)%X(1,n) = Int_Flux%Flux_Vol_Coord(n,v1)
          Trunc_Vol(:,face)%X(2,n) = Int_Flux%Flux_Vol_Coord(n,v2)
          Trunc_Vol(:,face)%X(3,n) = Int_Flux%Flux_Vol_Coord(n,v3)
          Trunc_Vol(:,face)%X(4,n) = Int_Flux%Flux_Vol_Coord(n,v4)
       end do
    endif

    ! Compute K   (the Area Vector of the ruled surface)
    do i = 1,ndim
       Tmp1(:,i) = Trunc_Vol(:,face)%X(3,i) - Trunc_Vol(:,face)%X(1,i)
       Tmp2(:,i) = Trunc_Vol(:,face)%X(4,i) - Trunc_Vol(:,face)%X(2,i)
    end do

    Trunc_Vol(:,face)%K(1) = Tmp1(:,2)*Tmp2(:,3) - Tmp1(:,3)*Tmp2(:,2)
    Trunc_Vol(:,face)%K(2) = Tmp1(:,3)*Tmp2(:,1) - Tmp1(:,1)*Tmp2(:,3)
    Trunc_Vol(:,face)%K(3) = Tmp1(:,1)*Tmp2(:,2) - Tmp1(:,2)*Tmp2(:,1)

    ! Compute V1234
    Trunc_Vol(:,face)%V1234 = 0.0_r8
    do i = 1,ndim
       Tmp1(:,i) = Trunc_Vol(:,face)%X(1,i) - Trunc_Vol(:,face)%X(2,i) + &
            Trunc_Vol(:,face)%X(3,i) - Trunc_Vol(:,face)%X(4,i)
       Trunc_Vol(:,face)%V1234 = Trunc_Vol(:,face)%V1234 +       &
                                       Tmp1(:,i)*Trunc_Vol(:,face)%K(i)
    end do
    Trunc_Vol(:,face)%V1234 = 0.5_r8 * Trunc_Vol(:,face)%V1234

    ! Compute the Mu-i-s.  This is the normal of the interface dotted
    ! with the coordinates of each faces vertex.  Mu-p is the vertex number
    ! of the vertex for each face (1 <= Mu-p <= nvf).  When the Mu-i-s are
    ! in ascending order, it denotes which vertex the interface will
    ! pass through first, second, etc.  The variable Mu-p is the vertex
    ! number of the reordered distances.
    do j = 1, nvf
       Trunc_Vol(:,face)%MUi(j) = 0.0_r8
       do i = 1,ndim
          Trunc_Vol(:,face)%MUi(j) = Trunc_Vol(:,face)%MUi(j) +     &
                 Int_Geom%Normal(i)*Trunc_Vol(:,face)%X(j,i)
       end do
       Trunc_Vol(:,face)%MUp(j) = j
       Trunc_Vol(:,face)%MUa(j) = Trunc_Vol(:,face)%MUi(j)
    end do

    ! Here Nu and Lambda are temporaries used to facilitate ordering the
    ! Mu-i into the Mu-a.  Put the minimum distance (the first vertex
    ! that the interface will pass through) into Mu-p(1) and put its
    ! facial vertex number into Mu-p(1).
    Trunc_Vol(:,face)%Nu = 0.0_r8
    Trunc_Vol(:,face)%Lambda = MIN(Trunc_Vol(:,face)%MUi(1),Trunc_Vol(:,face)%MUi(2), &
                                   Trunc_Vol(:,face)%MUi(3),Trunc_Vol(:,face)%MUi(4))
    do j = 1, nvf
       where (Trunc_Vol(:,face)%MUi(j) == Trunc_Vol(:,face)%Lambda .and.   &
               Trunc_Vol(:,face)%Nu == 0.0_r8)
          Trunc_Vol(:,face)%MUp(j) = 1
          Trunc_Vol(:,face)%MUp(1) = j
          Trunc_Vol(:,face)%MUa(j) = Trunc_Vol(:,face)%MUa(1)
          Trunc_Vol(:,face)%MUa(1) = Trunc_Vol(:,face)%MUi(j)
          Trunc_Vol(:,face)%Nu = 1.0_r8
       end where
    end do

    ! Now that the minimum distance is in element 1, order the
    ! other vertices in ascending order using a bubble sort.
    where (Trunc_Vol(:,face)%MUa(3) > Trunc_Vol(:,face)%MUa(4))

       Trunc_Vol(:,face)%Lambda   = Trunc_Vol(:,face)%MUp(4)
       Trunc_Vol(:,face)%MUp   (4) = Trunc_Vol(:,face)%MUp(3)
       Trunc_Vol(:,face)%MUp   (3) = Trunc_Vol(:,face)%Lambda

       Trunc_Vol(:,face)%Lambda   = Trunc_Vol(:,face)%MUa(4)
       Trunc_Vol(:,face)%MUa   (4) = Trunc_Vol(:,face)%MUa(3)
       Trunc_Vol(:,face)%MUa   (3) = Trunc_Vol(:,face)%Lambda

    end where

    where (Trunc_Vol(:,face)%MUa(2) > Trunc_Vol(:,face)%MUa(3))

       Trunc_Vol(:,face)%Lambda   = Trunc_Vol(:,face)%MUp(3)
       Trunc_Vol(:,face)%MUp   (3) = Trunc_Vol(:,face)%MUp(2)
       Trunc_Vol(:,face)%MUp   (2) = Trunc_Vol(:,face)%Lambda

       Trunc_Vol(:,face)%Lambda   = Trunc_Vol(:,face)%MUa(3)
       Trunc_Vol(:,face)%MUa   (3) = Trunc_Vol(:,face)%MUa(2)
       Trunc_Vol(:,face)%MUa   (2) = Trunc_Vol(:,face)%Lambda

    end where

    where (Trunc_Vol(:,face)%MUa(3) > Trunc_Vol(:,face)%MUa(4))

       Trunc_Vol(:,face)%Lambda   = Trunc_Vol(:,face)%MUp   (4)
       Trunc_Vol(:,face)%MUp   (4) = Trunc_Vol(:,face)%MUp   (3)
       Trunc_Vol(:,face)%MUp   (3) = Trunc_Vol(:,face)%Lambda

       Trunc_Vol(:,face)%Lambda   = Trunc_Vol(:,face)%MUa   (4)
       Trunc_Vol(:,face)%MUa   (4) = Trunc_Vol(:,face)%MUa   (3)
       Trunc_Vol(:,face)%MUa   (3) = Trunc_Vol(:,face)%Lambda

    end where

    ! The definition of Lambda is given in Eqn. 11.5 of Zemach-s notes.
    Trunc_Vol(:,face)%Lambda = Trunc_Vol(:,face)%MUi(2)*Trunc_Vol(:,face)%MUi(4) &
                             - Trunc_Vol(:,face)%MUi(1)*Trunc_Vol(:,face)%MUi(3)

    ! Nu is the face deviation vector dotted with the interface normal.
    ! If a face is a parallelogram, Nu (and B) will be zero.
    ! If a face is not a parallelogram, the magnitude of B measures
    ! the deviation from the parallelogram.
    Trunc_Vol(:,face)%Nu = Trunc_Vol(:,face)%MUi(1) - Trunc_Vol(:,face)%MUi(2) + &
                             Trunc_Vol(:,face)%MUi(3) - Trunc_Vol(:,face)%MUi(4)

  END SUBROUTINE FACE_PARAM

  SUBROUTINE TRUNCATE_VOLUME (Volume_Trunc_Total)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the truncation volume.
    !
    !=======================================================================
    use parameter_module, only: nicells
    use legacy_mesh_api, only: nfc
    use vof_data_module,  only: Cases, count_cases

    ! Arguments
    real(r8), dimension(nicells), intent(OUT) :: Volume_Trunc_Total

    ! Local Variables
    integer :: f
    real(r8), dimension(nicells) :: Volume_Trunc_Face

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize relevant quantities
    Volume_Trunc_Total     = 0.0_r8
    if (count_cases) Cases = 0

    ! Loop over faces, accumulating the truncated volume
    do f = 1, nfc

       Volume_Trunc_Face = 0.0_r8
       call TRUNCATE_FACE (f, Volume_Trunc_Face)

       Volume_Trunc_Total = Volume_Trunc_Total + Volume_Trunc_Face

    end do

  END SUBROUTINE TRUNCATE_VOLUME

  SUBROUTINE TRUNCATE_FACE (f, Vf)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the volume truncated at the current
    !   hex face by the plane given by X*Normal - Ro = 0
    !
    !=======================================================================
    use interface_module, only: Int_Geom
    use parameter_module, only: nicells
    use legacy_mesh_api,  only: nfc
    use vof_data_module,  only: Cases, count_cases
#ifdef TRUCHAS_ALLOW_UNSAFE_VECTORIZATION
    use ieee_exceptions
#endif

    ! Arguments
    integer, intent(IN)  :: f
    real(r8), dimension(nicells), intent(OUT) :: Vf

    ! Local Variables
    integer :: i
    integer, dimension(nicells,nfc) :: Truncation_Case
    logical, dimension(nicells,5) :: Face_Trun
    real(r8), dimension(nicells) :: Q
    real(r8), dimension(nicells,nvf) :: Y

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

#ifdef TRUCHAS_ALLOW_UNSAFE_VECTORIZATION
    logical :: old_halting_mode
    call ieee_get_halting_mode(ieee_divide_by_zero, old_halting_mode)
    call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
#endif


    Y = 0.0_r8
    call Y_FUNCTION (f, Y)

    Vf = 0.0_r8

    ! Case #1 intersection
    Q = 0.0_r8
    call TRUNCATE_FACE_N (1, f, Q, Y)
    where (Int_Geom%Rho >  Trunc_Vol(:,f)%MUa(1) .and. &
           Int_Geom%Rho <= Trunc_Vol(:,f)%MUa(2) .and. &
         Trunc_Vol(:,f)%MUa(1) /= Trunc_Vol(:,f)%MUa(2))
       Vf = Vf + Q
    end where

    ! Count the number of case #1 intersections
    if (count_cases) then
       Face_Trun(:,1) = Int_Geom%Rho <= Trunc_Vol(:,f)%MUa(2) .and. &
                        Int_Geom%Rho >  Trunc_Vol(:,f)%MUa(1) &
            .and. Trunc_Vol(:,f)%MUa(2) /= Trunc_Vol(:,f)%MUa(1)
       cases(1) = cases(1) + COUNT(Face_Trun(:,1))
    end if

    ! Case #5 intersection
    where (Int_Geom%Rho >  Trunc_Vol(:,f)%MUa(2) .and. &
           Int_Geom%Rho <= Trunc_Vol(:,f)%MUa(3) .and. &
         Trunc_Vol(:,f)%MUp(2) == (1 + MOD(Trunc_Vol(:,f)%MUp(1)+1, nvf)))
       Vf = Vf + Q
    end where

    Q = 0.0_r8
    call TRUNCATE_FACE_N (2, f, Q, Y)
    where (Int_Geom%Rho >  Trunc_Vol(:,f)%MUa(2) .and. &
           Int_Geom%Rho <= Trunc_Vol(:,f)%MUa(3) .and. &
         Trunc_Vol(:,f)%MUp(2) == (1 + MOD(Trunc_Vol(:,f)%MUp(1)+1, nvf)))
       Vf = Vf + Q
    end where

    ! Count the number of case #5 intersections
    if (count_cases) then
       Face_Trun(:,5) = Int_Geom%Rho >  Trunc_Vol(:,f)%MUa(2) .and. &
                        Int_Geom%Rho <= Trunc_Vol(:,f)%MUa(3) &
            .and. Trunc_Vol(:,f)%MUp(2) == (1+MOD(Trunc_Vol(:,f)%MUp(1)+1,nvf))
       cases(5) = cases(5) + COUNT(Face_Trun(:,5))
    end if

    ! Case #3 intersection
    Q = 0.0_r8
    call TRUNCATE_FACE_N (4, f, Q, Y)
    where (Int_Geom%Rho > Trunc_Vol(:,f)%MUa(3) .and. &
           Int_Geom%Rho < Trunc_Vol(:,f)%MUa(4)) Vf = Vf - Q

    ! Count the number of case #3 intersections
    if (count_cases) then
       Face_Trun(:,3) = Int_Geom%Rho > Trunc_Vol(:,f)%MUa(3) .and. &
                        Int_Geom%Rho < Trunc_Vol(:,f)%MUa(4)
       cases(3) = cases(3) + COUNT(Face_Trun(:,3))
    end if

    ! Case #4 intersection.  Note, this term is in common with
    !    a term in the calculation of case #3.  That is why the
    !    conditional doesn-t exclude case #3.
    Q = 0.0_r8
    call TRUNCATE_FACE_4 (f, Q)
    where (Int_Geom%Rho > Trunc_Vol(:,f)%MUa(3)) Vf = Vf + Q

    ! Count the number of case #4 intersections
    if (count_cases) then
       Face_Trun(:,4) = Int_Geom%Rho > Trunc_Vol(:,f)%MUa(3)
       cases(4) = cases(4) + COUNT(Face_Trun(:,4))
    end if

    ! Case #2 intersection
    Q = 0.0_r8
    call TRUNCATE_FACE_2 (f, Q, Y)
    where (Int_Geom%Rho >  Trunc_Vol(:,f)%MUa(2) .and. &
           Int_Geom%Rho <= Trunc_Vol(:,f)%MUa(3) .and. &
         Trunc_Vol(:,f)%MUp(2) /= (1+MOD(Trunc_Vol(:,f)%MUp(1)+1,nvf)))
       Vf = Vf + Q
    end where

    ! Count the number of case #2 intersections
    if (count_cases) then
       Face_Trun(:,2) = Int_Geom%Rho >  Trunc_Vol(:,f)%MUa(2) .and. &
                        Int_Geom%Rho <= Trunc_Vol(:,f)%MUa(3) &
            .and. Trunc_Vol(:,f)%MUp(2) /= (1+MOD(Trunc_Vol(:,f)%MUp(1)+1,nvf))
       cases(2) = cases(2) + COUNT(Face_Trun(:,2))
       ! Store the case types in each face
       Truncation_Case(:,f) = 0
       do i = 1,5
          where (Face_Trun(:,i)) Truncation_Case(:,f) = i
       end do

    end if
#ifdef TRUCHAS_ALLOW_UNSAFE_VECTORIZATION
    call ieee_set_halting_mode(ieee_divide_by_zero, old_halting_mode)
#endif

  END SUBROUTINE TRUNCATE_FACE

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><>

  SUBROUTINE TRUNCATE_FACE_2 (face, Q, Y)
    !=======================================================================
    ! PURPOSE -
    !   Compute an expression (for case "2") that is part of
    !   the volume truncated along a hex face by the plane
    !   X*Normal - Ro = 0
    !=======================================================================
    use cutoffs_module,   only: alittle
    use interface_module, only: Int_Geom
    use parameter_module, only: nicells
    use vof_data_module,  only: Eps

    ! Arguments
    integer, intent(IN)  :: face
    real(r8), dimension(nicells,nvf), intent(IN)  :: Y
    real(r8), dimension(nicells), intent(OUT) :: Q

    ! Local Variables
    integer :: i, k
    real(r8), dimension(nicells) :: J1, J2, J3, W1, W2, W3, W4, Z1, Z2, Z3

    real(r8), parameter :: one_third     = 1.0_r8 / 3.0_r8
    real(r8), parameter :: one_fourth    = 0.25_r8
    real(r8), parameter :: one_fifth     = 0.2_r8
    real(r8), parameter :: one_sixth     = 1.0_r8 / 6.0_r8
    real(r8), parameter :: one_seventh   = 1.0_r8 / 7.0_r8
    real(r8), parameter :: one_eighth    = 0.125_r8

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Q  = 0.0_r8
    J1 = 0.0_r8
    J2 = 0.0_r8
    J3 = 0.0_r8
    W1 = 0.0_r8
    W2 = 0.0_r8
    W3 = 0.0_r8
    W4 = 0.0_r8
    Z1 = 0.0_r8
    Z2 = 0.0_r8
    Z3 = 0.0_r8

    do k = 1, nvf
       i = 1 + mod(k + 1,nvf)

       ! Scan the four face indicies until we get to the
       ! 2nd predecessor of MU(b).
       where (Trunc_Vol(:,face)%MUp(2) == k)

          ! Z3 is the Y function of the successor of the A vertex
          Z3 = Y(:,k)

          ! W1 is the difference between:
          !   the Mu of the 2nd successor of vertex B and the Mu of vertex A.
          W1 = (Trunc_Vol(:,face)%MUi(i) - Trunc_Vol(:,face)%MUa(1))
       endwhere
    end do

    where (ABS(W1) > alittle)
       W1 = 1.0_r8 / W1
    elsewhere
       W1 = 0.0_r8
    endwhere

    do k = 1, nvf
       i = 1 + mod(k + 1,nvf)

       ! Scan the four face indicies until we get to the
       ! 2nd predecessor of MU(a).
       where (Trunc_Vol(:,face)%MUp(1) == k)

          ! Z1 is the Y function of the A vertex
          Z1 = Y(:,k)

          ! Z2 is the Y function of the 2nd successor of the A vertex
          Z2 = Y(:,i)
          W3 = eps(k)*Trunc_Vol(:,face)%Nu*W1
          Q  = eps(k)*Trunc_Vol(:,face)%V1234*0.5_r8*W1*W1

          ! W2 is the difference between:
          !   the Mu of the 2nd successor of vertex A and
          !   the Mu of vertex B.
          W2 = Trunc_Vol(:,face)%MUi(i) - Trunc_Vol(:,face)%MUa(2)
       endwhere
    end do

    where (ABS(W2) > alittle)
       W2 = 1.0_r8 / W2
       elsewhere
       W2 = 0.0_r8
    endwhere

    where ((ABS(W3) > 1.0e-2) .and. (W3 /= -1.0_r8))
       W4 = 1.0_r8/W3
       J1 = (1.0_r8 - LOG(ABS(1.0_r8 + W3))*W4)*W4
       J2 = (0.5_r8 - J1)*W4
       J3 = (one_third - J2)*W4
    elsewhere
       J3 = one_fourth - one_fifth*W3 + one_sixth*W3*W3 - &
            one_seventh*W3*W3*W3 + one_eighth*W3**4
       J2 = one_third - W3*J3
       J1 = 0.5_r8 - W3*J2
    end where

    Q = Q*(J1*(Int_Geom%Rho - Trunc_Vol(:,face)%MUa(1))**2 -  &
         2.0_r8*(Trunc_Vol(:,face)%MUa(2) - Trunc_Vol(:,face)%MUa(1))*&
         (Int_Geom%Rho - Trunc_Vol(:,face)%MUa(1))*J2 &
         + J3*(Trunc_Vol(:,face)%MUa(2) - Trunc_Vol(:,face)%MUa(1))**2)

    Q = Q + W1*Z1*(2.0_r8*Int_Geom%Rho - Trunc_Vol(:,face)%MUa(1) - &
            Trunc_Vol(:,face)%MUa(2))/6.0_r8

    Q = Q + (W1*W2*(Z2 - Z3)*(Int_Geom%Rho - Trunc_Vol(:,face)%MUa(2))**2)/6.0_r8

  END SUBROUTINE TRUNCATE_FACE_2

  SUBROUTINE TRUNCATE_FACE_4 (face, Q)
    !=======================================================================
    ! PURPOSE -
    !   Compute an expression (for case "4") that is part of
    !   the volume truncated along a hex face by the plane
    !   X*Normal - Ro = 0
    !=======================================================================
    use interface_module, only: Int_Geom
    use parameter_module, only: nicells

    ! Arguments
    integer, intent(IN) :: face
    real(r8), dimension(nicells), intent(INOUT) :: Q

    ! Local Variables
    integer :: i
    real(r8), dimension(nicells,ndim) :: Tmp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do i = 1,ndim
       Tmp(:,i) = 0.25_r8*(Trunc_Vol(:,face)%X(1,i) + Trunc_Vol(:,face)%X(2,i) +  &
                              Trunc_Vol(:,face)%X(3,i) + Trunc_Vol(:,face)%X(4,i)) - &
                              Int_Geom%Rho*Int_Geom%Normal(i)
    end do

    Q = 0.0_r8
    do i = 1,ndim
       Q = Q + Tmp(:,i)*Trunc_Vol(:,face)%K(i)
    end do
    Q = Q/6.0_r8

  END SUBROUTINE TRUNCATE_FACE_4

  SUBROUTINE TRUNCATE_FACE_N (n, face, Q, Y)
    !=======================================================================
    ! PURPOSE -
    !   Compute an expression (for case "n") that is part of
    !   the volume truncated along a hex face by the plane
    !   X*Normal - Ro = 0
    !=======================================================================
    use cutoffs_module,   only: alittle
    use interface_module, only: Int_Geom
    use parameter_module, only: nicells
    use vof_data_module,  only: Eps
    ! Arguments
    integer, intent(IN) :: n, face
    real(r8), dimension(nicells,nvf), intent(IN) :: Y
    real(r8), dimension(nicells), intent(INOUT) :: Q

    ! Local Variables
    integer :: k
    real(r8), dimension(nicells) :: J1, S, T, W1, W3

    real(r8), parameter :: two_thirds    = 2.0_r8 / 3.0_r8
    real(r8), parameter :: one_twelfth   = 1.0_r8 / 12.0_r8
    real(r8), parameter :: one_thirtieth = 1.0_r8 / 30.0_r8
    real(r8), parameter :: one_sixtieth  = 1.0_r8 / 60.0_r8
    real(r8), parameter :: one_168th     = 1.0_r8 / 168.0_r8
    real(r8), parameter :: one_105th     = 1.0_r8 / 105.0_r8

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Q = 0.0_r8

    T = Int_Geom%Rho - Trunc_Vol(:,face)%MUa(n)
    W3 = Trunc_Vol(:,face)%Lambda + Trunc_Vol(:,face)%Nu*Trunc_Vol(:,face)%MUa(n)

    where (ABS(W3) > alittle)
       W3 = 1.0_r8/W3
    elsewhere
       W3 = 0.0_r8
    endwhere

    W1 = T*T*W3
    T = Trunc_Vol(:,face)%Nu*T*W3

    S = 0.0_r8
    where (ABS(T) > alittle) S = 1.0_r8/T

    J1 = 0.0_r8
    where (T /= -1.0_r8) J1 = (1.0_r8-LOG(ABS(1.0_r8+T))*S)*S

    where (ABS(T) > 1.0e-2)
       W3 = J1 + (-two_thirds + 2.0_r8*J1 + (J1 - 0.5_r8) * S) * S
    elsewhere
       W3 = one_twelfth + T*(-one_thirtieth + T*(one_sixtieth + T*(-one_105th + T*one_168th)))
    endwhere

    do k = 1, nvf
       where (Trunc_Vol(:,face)%MUp(n) == k)
          Q = eps(k)*(Y(:,k)*W1/6.0_r8 + 0.5_r8*Trunc_Vol(:,face)%V1234*W3*W1*W1)
       end where
    end do

  END SUBROUTINE TRUNCATE_FACE_N

  SUBROUTINE Y_FUNCTION (face, Y)
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
    use interface_module, only: Int_Geom
    use parameter_module, only: nicells

    ! Arguments
    integer, intent(IN) :: face
    real(r8), dimension(nicells,nvf), intent(INOUT) :: Y

    ! Local Variables
    integer :: i
    integer, parameter :: ICMP = 1, JCMP = 2, KCMP = 3
    real(r8), dimension(nicells,ndim) :: R1, R2, R3, R4, S1, S3

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do i = 1,ndim
       R4(:,i) = Int_Geom%Rho * Int_Geom%Normal(i)
    end do

     do i = 1,nicells
    R1(i,:) = Trunc_Vol(i,face)%X(1,:) - R4(i,:)
    R2(i,:) = Trunc_Vol(i,face)%X(2,:) - R4(i,:)
    R3(i,:) = Trunc_Vol(i,face)%X(3,:) - R4(i,:)
    R4(i,:) = Trunc_Vol(i,face)%X(4,:) - R4(i,:)
     end do

    !  S1 = (X1 - Normal.Ro) X (X2 - Normal.Ro)
    S1(:,ICMP) = R1(:,JCMP)*R2(:,KCMP) - R1(:,KCMP)*R2(:,JCMP)
    S1(:,JCMP) = R1(:,KCMP)*R2(:,ICMP) - R1(:,ICMP)*R2(:,KCMP)
    S1(:,KCMP) = R1(:,ICMP)*R2(:,JCMP) - R1(:,JCMP)*R2(:,ICMP)

    !  S3 = (X3 - Normal.Ro) X (X4 - Normal.Ro)
    S3(:,ICMP) = R3(:,JCMP)*R4(:,KCMP) - R3(:,KCMP)*R4(:,JCMP)
    S3(:,JCMP) = R3(:,KCMP)*R4(:,ICMP) - R3(:,ICMP)*R4(:,KCMP)
    S3(:,KCMP) = R3(:,ICMP)*R4(:,JCMP) - R3(:,JCMP)*R4(:,ICMP)

    !  Y1 = (X4 - Normal.Ro) . (X1 - Normal.Ro) X (X2 - Normal.Ro)
    Y(:,1) = R4(:,ICMP)*S1(:,ICMP) + R4(:,JCMP)*S1(:,JCMP) + R4(:,KCMP)*S1(:,KCMP)

    !  Y2 = (X1 - Normal.Ro) . (X2 - Normal.Ro) X (X3 - Normal.Ro)
    !     However, we take advantage of the fact the vector identity:
    !     A . B X C = C . A X B  which means that we don-t need to
    !     compute the S2 vector as above.  We can use the S1 vector:
    !     Y2 = (X3 - Normal.Ro) . (X1 - Normal.Ro) X (X2 - Normal.Ro)
    Y(:,2) = R3(:,ICMP)*S1(:,ICMP) + R3(:,JCMP)*S1(:,JCMP) + R3(:,KCMP)*S1(:,KCMP)

    !  Y3 = (X2 - Normal.Ro) . (X3 - Normal.Ro) X (X4 - Normal.Ro)
    Y(:,3) = R2(:,ICMP)*S3(:,ICMP) + R2(:,JCMP)*S3(:,JCMP) + R2(:,KCMP)*S3(:,KCMP)

    !  Y4 = (X3 - Normal.Ro) . (X4 - Normal.Ro) X (X1 - Normal.Ro)
    !     However, we take advantage of the fact the vector identity:
    !     A . B X C = C . A X B  which means that we don-t need to
    !     compute the S4 vector as above.  We can use the S3 vector:
    !     Y4 = (X1 - Normal.Ro) . (X3 - Normal.Ro) X (X4 - Normal.Ro)
    Y(:,4) = R1(:,ICMP)*S3(:,ICMP) + R1(:,JCMP)*S3(:,JCMP) + R1(:,KCMP)*S3(:,KCMP)

  END SUBROUTINE Y_FUNCTION

END MODULE TRUNCATE_VOLUME_MODULE
