MODULE LOCATE_PLANE_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures necessary to locate the planar interfaces
  !   in all interface cells.
  !
  !   Public Interface:
  !
  !     * call LOCATE_PLANE ()
  !
  !       Compute the constant Rho in the interface plane equation
  !       X*Normal - Rho = 0 (given the plane Normal) subject to a
  !       volume conservation constraint. The interface constant
  !       Rho effectively locates the plane in the cell.
  !
  ! Contains: LOCATE_PLANE
  !           RHO_BRACKET
  !           RHO_BRENT
  !           RHO_NEWTON
  !           DVOL_DRO
  !           DVOL_DRO_FACE
  !           DVOL_DRO_FACE
  !           DVOL_DRO_FACE_2
  !           DVOL_DRO_FACE_N
  !
  ! Author(s): Stewart J. Mosso, LANL X-HM (sjm@lanl.gov)
  !            Douglas B. Kothe, LANL (dbk@lanl.gov)
  !
  !=======================================================================
  implicit none
 
  ! Private Module
  private
 
  ! Public Subroutines
  public :: LOCATE_PLANE
 
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
CONTAINS
 
  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE LOCATE_PLANE ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Given the equation for interface planes: X*Normal - Rho = 0,
    !   compute the parameter Rho (given the plane Normal) subject
    !   to a volume conservation constraint.  This is done in
    !   an iterative fashion, using either Newton's or Brent's
    !   method.  The iteration is first bracketed by finding the
    !   interval [Rho_Min,Rho_Max] bracketing Rho.  Given Rho and the
    !   volume flux geometry, the plane equation parameters
    !   (Rho, Normal) are then used to compute the truncated
    !   flux volume (Vp).
    !
    !=======================================================================
    use kind_module,      only: real_kind
    use parameter_module, only: nicells, nvc
    use timing_tree
    use vof_data_module,  only: volume_track_brents_method
 
    implicit none
 
    ! Local Variables
    real(real_kind), dimension(nicells,nvc) :: V_v
    real(real_kind), dimension(nicells)     :: Rho_Min, Rho_Max, V_Min, V_Max
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the locate plane timer.
    call start_timer ("Locate Plane")

    ! Bracket the correct value of Rho [Rho_Min,Rho_Max] to insure
    ! the subsequent iteration will converge efficiently
    call RHO_BRACKET (Rho_Min, Rho_Max, V_Min, V_Max, V_v)

    ! Compute Rho, which parameterizes the plane
    ! characterized by the equation Normal*X - Rho = 0

    if (volume_track_brents_method) then 
       ! Brents method iteration for Rho
       call RHO_BRENT (Rho_Min, Rho_Max, V_Min, V_Max)
    else
       ! Newton-Raphson iteration for Rho
       call RHO_NEWTON (Rho_Min, Rho_Max, V_Min, V_Max, V_v)
    end if
 
    ! Stop the locate plane timer.
    call stop_timer("Locate Plane")

    return
 
  END SUBROUTINE LOCATE_PLANE
 
  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE RHO_BRACKET (Rho_Min, Rho_Max, V_Min, V_Max, V_v)
    !=======================================================================
    ! PURPOSE -
    !   Compute the plane parameter Rho in the plane equation:
    !                    X*Normal - Rho = 0
    !   iteratively using Brents method as documented in Numerical
    !   Recipes (Press, Flannery, Teukolsky and Vetterling; Cambridge
    !   University Press, 1986), p. 253
    !=======================================================================
    use constants_module,       only: zero
    use interface_module,       only: Int_Geom
    use kind_module,            only: int_kind, real_kind
    use parameter_module,       only: ndim, nfc, nicells, nvc
    use truncate_volume_module, only: TRUNCATE_FACE, FACE_PARAM
    use truchas_logging_services

    implicit none
 
    ! Arguments
    real(real_kind), dimension(nicells),     intent(OUT) :: Rho_Min, Rho_Max, V_min, V_max
    real(real_kind), dimension(nicells,nvc), intent(OUT) :: V_v
 
    ! Local Variables
    integer(int_kind)                     :: f, n, v
    real(real_kind),   dimension(nicells) :: Volume, Rho
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize the bracketed values of truncation volume [V_Min,V_Max]
    ! and the associated plane constant [Rho_Min,Rho_Max]. 
    V_v = zero
 
    ! Compute volumes truncated by the interface passing through each of the vertices.
    FACE_LOOP: do f = 1,nfc
 
       ! Compute and store face parameters.
       call FACE_PARAM ('full_cell', f)
 
       VERTEX_LOOP: do v = 1,nvc

          ! Get the value of Rho for this vertex.
          Int_Geom%Rho = zero
          do n = 1,ndim
             Int_Geom%Rho = Int_Geom%Rho + Int_Geom%Cell_Coord(n,v) * Int_Geom%Normal(n)
          end do

          ! Get this face's contribution to the truncation volume.
          Volume = zero
          call TRUNCATE_FACE (f, Volume)
 
          V_v(:,v) = V_v(:,v) + Volume
 
       end do VERTEX_LOOP
 
    end do FACE_LOOP

    ! Set up the limits on the volume bracketing. Make sure that
    ! V_v is bounded zero = V_Min <= V_v <= V_Max = Cell Volume.
    ! If V_v > Cell Volume (because of roundoff), set V_v = Cell Volume

    V_Min = zero
    V_Max = Int_Geom%Cell_Volume
    Rho_Max = -HUGE(zero)
    Rho_Min =  HUGE(zero)

    ! If any V_v is outside of allowed bounds, force V_min <= V_v <= V_max
    do v = 1,nvc
       where (V_v(:,v) > V_Max) V_v(:,v) = V_Max
       where (V_v(:,v) < zero)  V_v(:,v) = V_Min
       Rho_Max = MAX(Rho_Max,V_v(:,v))
       Rho_Min = MIN(Rho_Min,V_v(:,v))
    end do

    ! Make sure that V_v = V_max for at least one vertex.
    ! If not, then set the maximum V_v to V_max.
    if (ANY(Rho_Max < V_Max)) then
      do v = 1,nvc
         where (V_v(:,v) == Rho_Max) V_v(:,v) = V_Max
      end do
    end if

    ! Make sure that V_v = V_min for at least one vertex.
    ! If not, then set the minimum V_v to V_min.
    if (ANY(Rho_Min > V_Min)) then
      do v = 1,nvc
         where (V_v(:,v) == Rho_Min) V_v(:,v) = V_Min
      end do
    end if

    ! Now bracket Rho.
    Rho_Max  =  HUGE(zero)
    Rho_Min  = -HUGE(zero)

    do v = 1,nvc

       ! Get the value of Rho for this vertex.
       Int_Geom%Rho = zero
       do n = 1,ndim
          Int_Geom%Rho = Int_Geom%Rho + Int_Geom%Cell_Coord(n,v) * Int_Geom%Normal(n)
       end do

       where ((V_v(:,v) <= Int_Geom%Vof*Int_Geom%Cell_Volume .and. V_v(:,v) > V_Min) &
              .or. V_v(:,v) == V_Min)
          Rho_Min = Int_Geom%Rho
          V_Min = V_v(:,v)
       end where
 
       where ((V_v(:,v) > Int_Geom%Vof*Int_Geom%Cell_Volume .and. V_v(:,v) < V_Max) &
              .or. V_v(:,v) == V_Max)
          Rho_Max = Int_Geom%Rho
          V_Max = V_v(:,v)
       end where

    end do
 
    ! At this point, [V_Min,V_Max] brackets V_v, and [Rho_Min,Rho_Max]
    ! are those plane constants corresponding to [V_Min,V_Max]. The
    ! bracketed volumes are monotonic, but the Rho's aren't necessarily
    ! monotonic, so make sure that Rho_Min < R_Max. Be sure to interchange
    ! V_Min and V_Max if R_Min and R_Max have to be interchanged.

    where (Rho_Max < Rho_Min)
       Rho = Rho_Max
       Volume = V_Max
       Rho_Max = Rho_Min
       V_Max = V_Min
       Rho_Min = Rho
       V_Min = Volume
    end where

    ! If we haven't properly set Rho everywhere, punt.
    if (ANY(Rho_Min == -HUGE(zero)) .or. ANY(Rho_Max == HUGE(zero))) then
       call TLS_panic ('RHO_BRACKET: unable to bracket plane constant Rho')
    end if

    return

  END SUBROUTINE RHO_BRACKET

  SUBROUTINE RHO_BRENT (Rho_Min, Rho_Max, V_Min, V_Max)
    !=======================================================================
    ! PURPOSE -
    !   Compute the plane parameter Rho in the plane equation:
    !                    X*Normal - Rho = 0
    !   iteratively using Brents method as documented in Numerical
    !   Recipes (Press, Flannery, Teukolsky and Vetterling; Cambridge
    !   University Press, 1986), p. 253
    !=======================================================================
    use cutoffs_module,         only: alittle
    use interface_module,       only: Int_Geom, Int_Flux
    use kind_module,            only: int_kind, log_kind, real_kind
    use parameter_module,       only: nicells
    use scalars_module,         only: zero, one_half, one, two, three
    use truncate_volume_module, only: TRUNCATE_VOLUME
    use vof_data_module,        only: volume_track_iter_max, volume_track_iter_tol
 
    implicit none
 
    ! Arguments
    real(real_kind), dimension(nicells), intent(IN) :: Rho_Min, Rho_Max, V_Min, V_Max
 
    ! Local Variables
    integer(int_kind)                     :: i
    logical(log_kind), dimension(nicells) :: Quad_intrp
    real(real_kind),   dimension(nicells) :: Rho_a, Rho_b, Rho_c, Rho_d,   &
                                             Rho_e, Rho_mid, Rho_tol, V_a, &
                                             V_b, V_c, P, Q, R, S
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
    ! Preset values before the iteration loop
    i = 0
    Int_Flux%Iter = 0
    Rho_a = Rho_Min
    Rho_b = Rho_Max
    Rho_c = zero
    Rho_d = zero
    Rho_e = zero
    P = zero
    Q = zero
    R = zero
    S = zero
    V_a = V_Min - Int_Geom%Vof*Int_Geom%Cell_Volume
    V_b = V_Max - Int_Geom%Vof*Int_Geom%Cell_Volume
    V_c = V_b
    Quad_intrp = .false.
 
    ! Main iteration loop
    do while (MAXVAL(ABS(V_b/MAX(Int_Geom%Cell_Volume,alittle))) > &
       volume_track_iter_tol .and. i < volume_track_iter_Max)
 
       i = i + 1
 
       where (V_b*V_c > zero)
          Rho_c = Rho_a
          V_c = V_a
          Rho_d = Rho_b - Rho_a
          Rho_e = Rho_d
       end where
 
       where (ABS(V_c) < ABS(V_b))
          Rho_a = Rho_b
          Rho_b = Rho_c
          Rho_c = Rho_a
          V_a = V_b
          V_b = V_c
          V_c = V_a
       end where
 
       Rho_tol = two*alittle*ABS(Rho_b) + one_half*volume_track_iter_tol
       Rho_mid = one_half*(Rho_c - Rho_b)
 
       where (ABS(Rho_mid) > Rho_tol .and. V_b /= zero) Int_Flux%Iter = i
 
       Quad_intrp = ABS(Rho_e) >= Rho_tol .and. ABS(V_a) > ABS(V_b)
 
       where (Quad_intrp) S = V_b/V_a

       where (Quad_intrp .and. Rho_a == Rho_c)
          P = two*Rho_mid*S
          Q = one - S
       end where

       where (Quad_intrp .and. Rho_a /= Rho_c)
          Q = V_a/V_c
          R = V_b/V_c
          P = S*(two*Rho_mid*Q*(Q - R) - (Rho_b - Rho_a)*(R - one))
          Q = (Q - one)*(R - one)*(S - one)
       end where

       where (Quad_intrp .and. P > zero) Q = -Q
       where (Quad_intrp) P = ABS(P)
 
       where (Quad_intrp)
          R = MIN(three*Rho_mid*Q - ABS(Rho_tol*Q), ABS(Rho_e*Q))
       end where
 
       where (Quad_intrp .and. two*P < R)
          Rho_e = Rho_d
          Rho_d = P/Q
       end where

       where (Quad_intrp .and. two*P >= R)
          Rho_d = Rho_mid
          Rho_e = Rho_d
       end where
 
       where (.not. Quad_intrp)
          Rho_d = Rho_mid
          Rho_e = Rho_d
       end where
 
       Rho_a = Rho_b ! Move last best guess to A
       V_a  = V_b
 
       where (ABS(Rho_d) > Rho_tol) ! Evaluate new trial root
          Rho_b = Rho_b + Rho_d
       elsewhere
          Rho_b = Rho_b + SIGN(Rho_tol,Rho_mid)
       end where
       Int_Geom%Rho = Rho_b
 
       call TRUNCATE_VOLUME (V_b)
 
       V_b = V_b - Int_Geom%Vof*Int_Geom%Cell_Volume
 
    end do
 
    Int_Geom%Rho = Rho_b
 
    return
 
  END SUBROUTINE RHO_BRENT
 
  SUBROUTINE RHO_NEWTON (Rho_Min, Rho_Max, V_Min, V_Max, V_v)
    !=======================================================================
    ! PURPOSE -
    !   Compute the plane parameter Rho in the plane equation:
    !                    X*Normal - Rho = 0
    !   iteratively using Newtons method.  The derivative d(V)/d(Rho)
    !   is computed analytically via a call to another routine.
    !=======================================================================
    use constants_module,       only: one_half, zero
    use interface_module,       only: Int_Geom, Int_Flux
    use kind_module,            only: int_kind, real_kind
    use parameter_module,       only: ndim, nicells, nvc
    use truncate_volume_module, only: TRUNCATE_VOLUME
    use vof_data_module,        only: volume_track_iter_max, volume_track_iter_tol
 
    implicit none
 
    ! Arguments
    real(real_kind), dimension(nicells)     :: Rho_Min, Rho_Max, V_Min, V_Max
    real(real_kind), dimension(nicells,nvc) :: V_v
 
    ! Local Variables
    integer(int_kind)                        :: i
    real(real_kind), dimension(nicells)      :: Rho, V_guess, mK, dV_guess
    real(real_kind), dimension(nicells,ndim) :: Temp
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
    ! Preset values before the iteration loop
 
    i = 0
    Int_Flux%Iter = 0
    Rho = -10000.0
 
    ! Main iteration loop
 
    where (V_Max > V_Min)
       Int_Geom%Rho = Rho_Min+(Int_Geom%Vof*Int_Geom%Cell_Volume - V_Min)* &
                       (Rho_Max-Rho_Min)/(V_Max-V_Min)
    elsewhere
       Int_Geom%Rho = Rho_Min
    endwhere
 
    where (Int_Geom%Vof > 0.99) Int_Geom%Rho = Rho_Min
 
    where (Int_Geom%Vof < 0.01) Int_Geom%Rho = Rho_Max
 
    Int_Geom%Rho = MAX(MIN(Int_Geom%Rho, Rho_Max), Rho_Min)
 
    call TRUNCATE_VOLUME (V_guess)
 
    do while (MAXVAL(ABS(Int_Geom%Vof-V_guess/Int_Geom%Cell_Volume)) >  &
         volume_track_iter_tol .and. i < volume_track_iter_Max)
 
       where (ABS(Int_Geom%Vof-V_guess/Int_Geom%Cell_Volume) > &
            volume_track_iter_tol) Int_Flux%Iter = Int_Flux%Iter + 1
 
       i = i + 1
 
       ! Where Iter < i, we have converged and found a value
       ! of Rho, so set the limits of the search interval to
       ! to the value of Int_Geom%Rho ( = Rho)
 
       where (Int_Flux%Iter < i)
          Rho_Min = Int_Geom%Rho
          Rho_Max = Int_Geom%Rho
          V_Min  = V_guess
          V_Max  = V_guess
       end where
 
       Rho = Int_Geom%Rho
 
       call DVOL_DRO (dV_guess, mK, Temp, V_v)
 
       dV_guess = MAX(dV_guess,zero)
 
       where (dV_guess /= zero)
          Int_Geom%Rho = Rho + (Int_Geom%Vof*Int_Geom%Cell_Volume - V_guess) / dV_guess
       elsewhere
          Int_Geom%Rho = Rho
       endwhere
 
       where ((Int_Geom%Rho<Rho_Min.or.Int_Geom%Rho > Rho_Max).and. &
            V_guess > Int_Geom%Vof*Int_Geom%Cell_Volume)
          Rho_Max = Rho
          V_Max = V_guess
          Int_Geom%Rho = (Rho_Min + Rho_Max)*one_half
       end where
 
       where ((Int_Geom%Rho < Rho_Min.or.Int_Geom%Rho > Rho_Max).and. &
            V_guess < Int_Geom%Vof*Int_Geom%Cell_Volume)
          Rho_Min = Rho
          V_Min = V_guess
          Int_Geom%Rho = (Rho_Min + Rho_Max)*one_half
       end where
 
       where (Int_Geom%Rho == Rho .and. V_guess > Int_Geom%Vof*Int_Geom%Cell_Volume)
          Rho_Max = Int_Geom%Rho
          V_Max = V_guess
          Int_Geom%Rho = (Rho_Min + Rho_Max)*one_half
       end where
 
       where (Int_Geom%Rho == Rho .and. V_guess < Int_Geom%Vof*Int_Geom%Cell_Volume)
          Rho_Min = Int_Geom%Rho
          V_Min = V_guess
          Int_Geom%Rho = (Rho_Min + Rho_Max)*one_half
       end where
 
       call TRUNCATE_VOLUME (V_guess)
 
    end do
 
    return
 
  END SUBROUTINE RHO_NEWTON

  SUBROUTINE DVOL_DRO (dV, mK, Temp, Vi)
    !=======================================================================
    ! PURPOSE - 
    !   Compute the derivative d(Vol)/d(Ro), where Vol is the
    !   volume truncated in a hex by a plane described by the
    !   equation X*Normal - Ro = 0
    !=======================================================================
    use constants_module,       only: zero, one_half
    use interface_module,       only: Int_Geom
    use kind_module,            only: int_kind, real_kind
    use parameter_module,       only: ndim, nfc, nicells, nvc
    use truncate_volume_module, only: Trunc_Vol

    implicit none

    ! Arguments
    real(real_kind), dimension(nicells)      :: dV, mK
    real(real_kind), dimension(nicells,ndim) :: Temp
    real(real_kind), dimension(nicells,nvc)  :: Vi

    ! Local Variables
    integer(int_kind)                        :: face, i
    real(real_kind), dimension(nicells)      :: Vf
    real(real_kind), dimension(nicells,ndim) :: T1, T2

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    dV = zero

    do face = 1, nfc

       do i = 1,ndim
          T1(:,i) = Trunc_Vol(:,face)%X(1,i) - Trunc_Vol(:,face)%X(2,i)
          T2(:,i) = Trunc_Vol(:,face)%X(4,i) - Trunc_Vol(:,face)%X(1,i)
       end do
       Temp(:,1)=T2(:,2)*T1(:,3)-T2(:,3)*T1(:,2)
       Temp(:,2)=T2(:,3)*T1(:,1)-T2(:,1)*T1(:,3)
       Temp(:,3)=T2(:,1)*T1(:,2)-T2(:,2)*T1(:,1)

       Vi(:,1) = zero
       do i = 1,ndim
          Vi(:,1) = Vi(:,1) + Temp(:,i) * Int_Geom%Normal(i)
       end do

       do i = 1,ndim
          T2(:,i) = Trunc_Vol(:,face)%X(2,i) - Trunc_Vol(:,face)%X(3,i)
       end do
       Temp(:,1)=T1(:,2)*T2(:,3)-T1(:,3)*T2(:,2)
       Temp(:,2)=T1(:,3)*T2(:,1)-T1(:,1)*T2(:,3)
       Temp(:,3)=T1(:,1)*T2(:,2)-T1(:,2)*T2(:,1)

       Vi(:,2) = zero
       do i = 1,ndim
          Vi(:,2) = Vi(:,2) + Temp(:,i) * Int_Geom%Normal(i)
       end do
        do i = 1,nicells
       T1(i,:) = Trunc_Vol(i,face)%K(:)
        end do
       mK = zero
       do i = 1,ndim
          mK = mK + T1(:,i) * Int_Geom%Normal(i)
       end do

       Vi(:,3) = mK - Vi(:,1)
       Vi(:,4) = mK - Vi(:,2)
       mK = -one_half*mK

       call dvol_dro_face (face, mK, Vf, Vi)

       dV = dV + Vf

    end do

    return

  END SUBROUTINE DVOL_DRO

  SUBROUTINE DVOL_DRO_FACE (face, mK, Vf, uX)
    !=======================================================================
    ! PURPOSE - 
    !   Compute the derivative d(Vol)/d(Ro) at the current hex
    !   face, where Vol is the volume truncated in the hex by
    !   a plane described by the equation X*Normal - Ro = 0
    !=======================================================================
    use constants_module,       only: zero
    use interface_module,       only: Int_Geom
    use kind_module,            only: int_kind, real_kind
    use parameter_module,       only: nicells, nvc, nvf
    use truncate_volume_module, only: Trunc_Vol

    implicit none

    ! Arguments
    integer(int_kind)                       :: face
    real(real_kind), dimension(nicells)     :: mK, Vf
    real(real_kind), dimension(nicells,nvc) :: uX

    ! Local Variables
    real(real_kind), dimension(nicells) :: Q

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Vf = zero

    Q = zero
    call dvol_dro_face_n (1, face, uX, Q)

    where (Int_Geom%Rho <= Trunc_Vol(:,face)%MUa(2) .and. &
           Int_Geom%Rho >  Trunc_Vol(:,face)%MUa(1)) Vf = Vf + Q

    where (Int_Geom%Rho >  Trunc_Vol(:,face)%MUa(2) .and. &
           Int_Geom%Rho <= Trunc_Vol(:,face)%MUa(3) .and. &
         Trunc_Vol(:,face)%MUp(2) == (1+MOD(Trunc_Vol(:,face)%MUp(1)+1,nvf)) ) &
                       Vf = Vf + Q 

    Q = zero
    call dvol_dro_face_n (2, face, uX, Q)

    where (Int_Geom%Rho >  Trunc_Vol(:,face)%MUa(2) .and.  &
           Int_Geom%Rho <= Trunc_Vol(:,face)%MUa(3) .and.  &
         Trunc_Vol(:,face)%MUp(2) == (1+MOD(Trunc_Vol(:,face)%MUp(1)+1,nvf)) ) &
                   Vf = Vf + Q

    Q = zero
    call dvol_dro_face_n (4, face, uX, Q)

    where (Int_Geom%Rho > Trunc_Vol(:,face)%MUa(3) .and. &
           Int_Geom%Rho < Trunc_Vol(:,face)%MUa(4)) Vf = Vf - Q

    where (Int_Geom%Rho > Trunc_Vol(:,face)%MUa(3)) Vf = Vf + mK

    Q = zero
    call dvol_dro_face_2 (face, uX, Q)

    where (Int_Geom%Rho >  Trunc_Vol(:,face)%MUa(2) .and.  &
           Int_Geom%Rho <= Trunc_Vol(:,face)%MUa(3) .and.  &
         Trunc_Vol(:,face)%MUp(2) /= (1+MOD(Trunc_Vol(:,face)%MUp(1)+1,nvf)) ) &
                  Vf = Vf + Q 

    return

  END SUBROUTINE DVOL_DRO_FACE

  SUBROUTINE DVOL_DRO_FACE_2 (face, UX, Q)
    !=======================================================================
    ! PURPOSE - 
    !   Compute the derivative d(Vol)/d(Ro) at the current hex
    !   face (for case "2"), where Vol is the volume truncated
    !   in the hex by a plane described by the equation 
    !   X*Normal - Ro = 0
    !=======================================================================
    use constants_module
    use cutoffs_module,         only: alittle
    use interface_module,       only: Int_Geom
    use kind_module,            only: int_kind, real_kind
    use parameter_module,       only: nicells, nvc, nvf
    use truncate_volume_module, only: Trunc_Vol
    use vof_data_module,        only: Eps

    implicit none

    ! Arguments
    integer(int_kind)                         :: face
    real(real_kind),   dimension(nicells,nvc) :: UX
    real(real_kind),   dimension(nicells)     :: Q

    ! Local Variables
    integer(int_kind)                   :: k
    real(real_kind), dimension(nicells) :: J1, J2, N, P, R, S, T, U

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    P = zero
    Q = zero
    R = zero
    S = zero
    T = zero
    U = zero

    do k = 1, nvf
       where (Trunc_Vol(:,face)%MUp(2) == k)
          R = Trunc_Vol(:,face)%MUi(1+mod(k+1,nvf)) - Trunc_Vol(:,face)%MUa(1)
          U = -UX(:,k)
       endwhere
    end do

    where (ABS(R) > alittle)
       R = one/R
       elsewhere
       R = zero
    endwhere

    do k = 1, nvf
       where (Trunc_Vol(:,face)%MUp(1) == k)
          S = Trunc_Vol(:,face)%MUi(1+mod(k+1,nvf)) - Trunc_Vol(:,face)%MUa(2)
          U = U+UX(:,1+mod(k+1,nvf))
          T = eps(k)*Trunc_Vol(:,face)%Nu*R
          Q = eps(k)*Trunc_Vol(:,face)%V1234*R*R
          P = UX(:,k)
       endwhere
    end do

    where (ABS(S) > alittle)
       S = one/S
    elsewhere
       S = zero
    endwhere

    where (ABS(T) > alittle)
       N =one/T
       elsewhere
       N = zero
    endwhere

    where ((ABS(T) > 1.0e-2) .and. (T /= -one))
       J1 = (one - LOG(ABS(one+T))*N)*N
       J2 = (one_half - J1)*N
    elsewhere
       J1 = one_fourth + T*(-one_fifth + &
            T*(one_sixth + T*(-one_seventh + T*one_eighth)))
       J2 = one_third - T*J1 
       J1 = one_half - T*J2
    endwhere

    Q = Q*((Int_Geom%Rho-Trunc_Vol(:,face)%MUa(1))*J1 - &
           (Trunc_Vol(:,face)%MUa(2) - Trunc_Vol(:,face)%MUa(1))*J2)

    Q = Q - one_half*R*P*(two*Int_Geom%Rho - Trunc_Vol(:,face)%MUa(1) - &
                  Trunc_Vol(:,face)%MUa(2))

    Q = Q - one_half*R*S*U*(Int_Geom%Rho - Trunc_Vol(:,face)%MUa(2))**2

    return

  END SUBROUTINE DVOL_DRO_FACE_2

  SUBROUTINE DVOL_DRO_FACE_N (n, face, UX, Q)
    !=======================================================================
    ! PURPOSE - 
    !   Compute the derivative d(Vol)/d(Ro) at the current hex
    !   face (for case "n"), where Vol is the volume truncated
    !   in the hex by a plane described by the equation 
    !   X*Normal - Ro = 0
    !=======================================================================
    use constants_module
    use cutoffs_module,         only: alittle
    use interface_module,       only: Int_Geom
    use kind_module,            only: int_kind, real_kind
    use parameter_module,       only: nicells, nvc, nvf
    use truncate_volume_module, only: Trunc_Vol
    use vof_data_module,        only: Eps

    implicit none

    ! Arguments
    integer(int_kind)                       :: n, face
    real(real_kind), dimension(nicells,nvc) :: UX
    real(real_kind), dimension(nicells)     :: Q

    ! Local Variables
    integer(int_kind)                   :: k
    real(real_kind), dimension(nicells) :: P, R, S, T, J1

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Q = zero

    T = Int_Geom%Rho - Trunc_Vol(:,face)%MUa(n)
    R = Trunc_Vol(:,face)%Lambda + Trunc_Vol(:,face)%Nu*Trunc_Vol(:,face)%MUa(n)
    where (ABS(R) > alittle)
       R = one/R
    elsewhere
       R = zero
    endwhere

    S = T*R*Trunc_Vol(:,face)%Nu
    where (ABS(S) > alittle)
       P = one/S
    elsewhere
       P = zero
    endwhere

    where (ABS(S) > 1.0e-2 .and. S /= -one)
       J1 = (one-LOG(ABS(one+S))*P)*P 
       J1 = J1 + (J1 - one_half)*P
    elsewhere
       J1 = one_sixth + S*(-one_twelfth + S*(one_twentieth + S*(-one_thirtieth + S*one_42nd)))
    endwhere

    do k = 1, nvf
       where (Trunc_Vol(:,face)%MUp(n) == k)
          Q = Trunc_Vol(:,face)%V1234*J1*T*R - one_half*UX(:,k)
          Q = Q*eps(k)*T*T*R
       endwhere
    end do

    return

  END SUBROUTINE DVOL_DRO_FACE_N

END MODULE LOCATE_PLANE_MODULE
