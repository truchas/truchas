!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private
 
  public :: LOCATE_PLANE
 
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
    use parameter_module, only: nicells
    use legacy_mesh_api,  only: nvc
    use vof_data_module,  only: volume_track_brents_method
    use truchas_timers
 
    ! Local Variables
    real(r8), dimension(nicells,nvc) :: V_v
    real(r8), dimension(nicells)     :: Rho_Min, Rho_Max, V_Min, V_Max
 
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
    use interface_module,       only: Int_Geom
    use parameter_module,       only: nicells
    use legacy_mesh_api,        only: ndim, nfc, nvc
    use truncate_volume_module, only: TRUNCATE_FACE, FACE_PARAM
 
    ! Arguments
    real(r8), dimension(nicells),     intent(OUT) :: Rho_Min, Rho_Max, V_min, V_max
    real(r8), dimension(nicells,nvc), intent(OUT) :: V_v
 
    ! Local Variables
    integer :: f, n, v
    real(r8), dimension(nicells) :: Volume, Rho
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize the bracketed values of truncation volume [V_Min,V_Max]
    ! and the associated plane constant [Rho_Min,Rho_Max]. 
    V_v = 0.0_r8
 
    ! Compute volumes truncated by the interface passing through each of the vertices.
    FACE_LOOP: do f = 1,nfc
 
       ! Compute and store face parameters.
       call FACE_PARAM ('full_cell', f)
 
       VERTEX_LOOP: do v = 1,nvc

          ! Get the value of Rho for this vertex.
          Int_Geom%Rho = 0.0_r8
          do n = 1,ndim
             Int_Geom%Rho = Int_Geom%Rho + Int_Geom%Cell_Coord(n,v) * Int_Geom%Normal(n)
          end do

          ! Get this face's contribution to the truncation volume.
          Volume = 0.0_r8
          call TRUNCATE_FACE (f, Volume)
 
          V_v(:,v) = V_v(:,v) + Volume
 
       end do VERTEX_LOOP
 
    end do FACE_LOOP

    ! Set up the limits on the volume bracketing. Make sure that
    ! V_v is bounded zero = V_Min <= V_v <= V_Max = Cell Volume.
    ! If V_v > Cell Volume (because of roundoff), set V_v = Cell Volume

    V_Min = 0.0_r8
    V_Max = Int_Geom%Cell_Volume
    Rho_Max = -HUGE(0.0_r8)
    Rho_Min =  HUGE(0.0_r8)

    ! If any V_v is outside of allowed bounds, force V_min <= V_v <= V_max
    do v = 1,nvc
       where (V_v(:,v) > V_Max) V_v(:,v) = V_Max
       where (V_v(:,v) < 0.0_r8)  V_v(:,v) = V_Min
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
    Rho_Max  =  HUGE(0.0_r8)
    Rho_Min  = -HUGE(0.0_r8)

    do v = 1,nvc

       ! Get the value of Rho for this vertex.
       Int_Geom%Rho = 0.0_r8
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
    if (ANY(Rho_Min == -HUGE(0.0_r8)) .or. ANY(Rho_Max == HUGE(0.0_r8))) then
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
    use parameter_module,       only: nicells
    use truncate_volume_module, only: TRUNCATE_VOLUME
    use vof_data_module,        only: volume_track_iter_max, volume_track_iter_tol
 
    ! Arguments
    real(r8), dimension(nicells), intent(IN) :: Rho_Min, Rho_Max, V_Min, V_Max
 
    ! Local Variables
    integer :: i
    logical, dimension(nicells) :: Quad_intrp
    real(r8), dimension(nicells) :: Rho_a, Rho_b, Rho_c, Rho_d,   &
                                    Rho_e, Rho_mid, Rho_tol, V_a, &
                                    V_b, V_c, P, Q, R, S
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
    ! Initialize values before the iteration loop
    i = 0
    Int_Flux%Iter = 0
    Rho_a = Rho_Min
    Rho_b = Rho_Max
    Rho_c = 0.0_r8
    Rho_d = 0.0_r8
    Rho_e = 0.0_r8
    P = 0.0_r8
    Q = 0.0_r8
    R = 0.0_r8
    S = 0.0_r8
    V_a = V_Min - Int_Geom%Vof*Int_Geom%Cell_Volume
    V_b = V_Max - Int_Geom%Vof*Int_Geom%Cell_Volume
    V_c = V_b
    Quad_intrp = .false.
 
    ! Main iteration loop
    do while (MAXVAL(ABS(V_b/MAX(Int_Geom%Cell_Volume,alittle))) > &
       volume_track_iter_tol .and. i < volume_track_iter_Max)
 
       i = i + 1
 
       where (V_b*V_c > 0.0_r8)
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
 
       Rho_tol = 2.0_r8*alittle*ABS(Rho_b) + 0.5_r8*volume_track_iter_tol
       Rho_mid = 0.5_r8*(Rho_c - Rho_b)
 
       where (ABS(Rho_mid) > Rho_tol .and. V_b /= 0.0_r8) Int_Flux%Iter = i
 
       Quad_intrp = ABS(Rho_e) >= Rho_tol .and. ABS(V_a) > ABS(V_b)
 
       where (Quad_intrp) S = V_b/V_a

       where (Quad_intrp .and. Rho_a == Rho_c)
          P = 2.0_r8*Rho_mid*S
          Q = 1.0_r8 - S
       end where

       where (Quad_intrp .and. Rho_a /= Rho_c)
          Q = V_a/V_c
          R = V_b/V_c
          P = S*(2.0_r8*Rho_mid*Q*(Q - R) - (Rho_b - Rho_a)*(R - 1.0_r8))
          Q = (Q - 1.0_r8)*(R - 1.0_r8)*(S - 1.0_r8)
       end where

       where (Quad_intrp .and. P > 0.0_r8) Q = -Q
       where (Quad_intrp) P = ABS(P)
 
       where (Quad_intrp)
          R = MIN(3.0_r8*Rho_mid*Q - ABS(Rho_tol*Q), ABS(Rho_e*Q))
       end where
 
       where (Quad_intrp .and. 2.0_r8*P < R)
          Rho_e = Rho_d
          Rho_d = P/Q
       end where

       where (Quad_intrp .and. 2.0_r8*P >= R)
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
 
  END SUBROUTINE RHO_BRENT
 
  SUBROUTINE RHO_NEWTON (Rho_Min, Rho_Max, V_Min, V_Max, V_v)
    !=======================================================================
    ! PURPOSE -
    !   Compute the plane parameter Rho in the plane equation:
    !                    X*Normal - Rho = 0
    !   iteratively using Newtons method.  The derivative d(V)/d(Rho)
    !   is computed analytically via a call to another routine.
    !=======================================================================
    use interface_module,       only: Int_Geom, Int_Flux
    use parameter_module,       only: nicells
    use legacy_mesh_api,        only: ndim, nvc
    use truncate_volume_module, only: TRUNCATE_VOLUME
    use vof_data_module,        only: volume_track_iter_max, volume_track_iter_tol
 
    ! Arguments
    real(r8), dimension(nicells)     :: Rho_Min, Rho_Max, V_Min, V_Max
    real(r8), dimension(nicells,nvc) :: V_v
 
    ! Local Variables
    integer :: i
    real(r8), dimension(nicells) :: Rho, V_guess, mK, dV_guess
    real(r8), dimension(nicells,ndim) :: Temp
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
    ! Initialize values before the iteration loop
 
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
 
       dV_guess = MAX(dV_guess,0.0_r8)
 
       where (dV_guess /= 0.0_r8)
          Int_Geom%Rho = Rho + (Int_Geom%Vof*Int_Geom%Cell_Volume - V_guess) / dV_guess
       elsewhere
          Int_Geom%Rho = Rho
       endwhere
 
       where ((Int_Geom%Rho<Rho_Min.or.Int_Geom%Rho > Rho_Max).and. &
            V_guess > Int_Geom%Vof*Int_Geom%Cell_Volume)
          Rho_Max = Rho
          V_Max = V_guess
          Int_Geom%Rho = (Rho_Min + Rho_Max)*0.5_r8
       end where
 
       where ((Int_Geom%Rho < Rho_Min.or.Int_Geom%Rho > Rho_Max).and. &
            V_guess < Int_Geom%Vof*Int_Geom%Cell_Volume)
          Rho_Min = Rho
          V_Min = V_guess
          Int_Geom%Rho = (Rho_Min + Rho_Max)*0.5_r8
       end where
 
       where (Int_Geom%Rho == Rho .and. V_guess > Int_Geom%Vof*Int_Geom%Cell_Volume)
          Rho_Max = Int_Geom%Rho
          V_Max = V_guess
          Int_Geom%Rho = (Rho_Min + Rho_Max)*0.5_r8
       end where
 
       where (Int_Geom%Rho == Rho .and. V_guess < Int_Geom%Vof*Int_Geom%Cell_Volume)
          Rho_Min = Int_Geom%Rho
          V_Min = V_guess
          Int_Geom%Rho = (Rho_Min + Rho_Max)*0.5_r8
       end where
 
       call TRUNCATE_VOLUME (V_guess)
 
    end do
 
  END SUBROUTINE RHO_NEWTON

  SUBROUTINE DVOL_DRO (dV, mK, Temp, Vi)
    !=======================================================================
    ! PURPOSE - 
    !   Compute the derivative d(Vol)/d(Ro), where Vol is the
    !   volume truncated in a hex by a plane described by the
    !   equation X*Normal - Ro = 0
    !=======================================================================
    use interface_module,       only: Int_Geom
    use parameter_module,       only: nicells
    use legacy_mesh_api,        only: ndim, nfc, nvc
    use truncate_volume_module, only: Trunc_Vol

    ! Arguments
    real(r8), dimension(nicells)      :: dV, mK
    real(r8), dimension(nicells,ndim) :: Temp
    real(r8), dimension(nicells,nvc)  :: Vi

    ! Local Variables
    integer :: face, i
    real(r8), dimension(nicells) :: Vf
    real(r8), dimension(nicells,ndim) :: T1, T2

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    dV = 0.0_r8

    do face = 1, nfc

       do i = 1,ndim
          T1(:,i) = Trunc_Vol(:,face)%X(1,i) - Trunc_Vol(:,face)%X(2,i)
          T2(:,i) = Trunc_Vol(:,face)%X(4,i) - Trunc_Vol(:,face)%X(1,i)
       end do
       Temp(:,1)=T2(:,2)*T1(:,3)-T2(:,3)*T1(:,2)
       Temp(:,2)=T2(:,3)*T1(:,1)-T2(:,1)*T1(:,3)
       Temp(:,3)=T2(:,1)*T1(:,2)-T2(:,2)*T1(:,1)

       Vi(:,1) = 0.0_r8
       do i = 1,ndim
          Vi(:,1) = Vi(:,1) + Temp(:,i) * Int_Geom%Normal(i)
       end do

       do i = 1,ndim
          T2(:,i) = Trunc_Vol(:,face)%X(2,i) - Trunc_Vol(:,face)%X(3,i)
       end do
       Temp(:,1)=T1(:,2)*T2(:,3)-T1(:,3)*T2(:,2)
       Temp(:,2)=T1(:,3)*T2(:,1)-T1(:,1)*T2(:,3)
       Temp(:,3)=T1(:,1)*T2(:,2)-T1(:,2)*T2(:,1)

       Vi(:,2) = 0.0_r8
       do i = 1,ndim
          Vi(:,2) = Vi(:,2) + Temp(:,i) * Int_Geom%Normal(i)
       end do
        do i = 1,nicells
       T1(i,:) = Trunc_Vol(i,face)%K(:)
        end do
       mK = 0.0_r8
       do i = 1,ndim
          mK = mK + T1(:,i) * Int_Geom%Normal(i)
       end do

       Vi(:,3) = mK - Vi(:,1)
       Vi(:,4) = mK - Vi(:,2)
       mK = -0.5_r8*mK

       call dvol_dro_face (face, mK, Vf, Vi)

       dV = dV + Vf

    end do

  END SUBROUTINE DVOL_DRO

  SUBROUTINE DVOL_DRO_FACE (face, mK, Vf, uX)
    !=======================================================================
    ! PURPOSE - 
    !   Compute the derivative d(Vol)/d(Ro) at the current hex
    !   face, where Vol is the volume truncated in the hex by
    !   a plane described by the equation X*Normal - Ro = 0
    !=======================================================================
    use interface_module,       only: Int_Geom
    use parameter_module,       only: nicells
    use legacy_mesh_api,        only: nvc, nvf
    use truncate_volume_module, only: Trunc_Vol

    ! Arguments
    integer :: face
    real(r8), dimension(nicells) :: mK, Vf
    real(r8), dimension(nicells,nvc) :: uX

    ! Local Variables
    real(r8), dimension(nicells) :: Q

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Vf = 0.0_r8

    Q = 0.0_r8
    call dvol_dro_face_n (1, face, uX, Q)

    where (Int_Geom%Rho <= Trunc_Vol(:,face)%MUa(2) .and. &
           Int_Geom%Rho >  Trunc_Vol(:,face)%MUa(1)) Vf = Vf + Q

    where (Int_Geom%Rho >  Trunc_Vol(:,face)%MUa(2) .and. &
           Int_Geom%Rho <= Trunc_Vol(:,face)%MUa(3) .and. &
         Trunc_Vol(:,face)%MUp(2) == (1+MOD(Trunc_Vol(:,face)%MUp(1)+1,nvf)) ) &
                       Vf = Vf + Q 

    Q = 0.0_r8
    call dvol_dro_face_n (2, face, uX, Q)

    where (Int_Geom%Rho >  Trunc_Vol(:,face)%MUa(2) .and.  &
           Int_Geom%Rho <= Trunc_Vol(:,face)%MUa(3) .and.  &
         Trunc_Vol(:,face)%MUp(2) == (1+MOD(Trunc_Vol(:,face)%MUp(1)+1,nvf)) ) &
                   Vf = Vf + Q

    Q = 0.0_r8
    call dvol_dro_face_n (4, face, uX, Q)

    where (Int_Geom%Rho > Trunc_Vol(:,face)%MUa(3) .and. &
           Int_Geom%Rho < Trunc_Vol(:,face)%MUa(4)) Vf = Vf - Q

    where (Int_Geom%Rho > Trunc_Vol(:,face)%MUa(3)) Vf = Vf + mK

    Q = 0.0_r8
    call dvol_dro_face_2 (face, uX, Q)

    where (Int_Geom%Rho >  Trunc_Vol(:,face)%MUa(2) .and.  &
           Int_Geom%Rho <= Trunc_Vol(:,face)%MUa(3) .and.  &
         Trunc_Vol(:,face)%MUp(2) /= (1+MOD(Trunc_Vol(:,face)%MUp(1)+1,nvf)) ) &
                  Vf = Vf + Q 

  END SUBROUTINE DVOL_DRO_FACE

  SUBROUTINE DVOL_DRO_FACE_2 (face, UX, Q)
    !=======================================================================
    ! PURPOSE - 
    !   Compute the derivative d(Vol)/d(Ro) at the current hex
    !   face (for case "2"), where Vol is the volume truncated
    !   in the hex by a plane described by the equation 
    !   X*Normal - Ro = 0
    !=======================================================================
    use cutoffs_module,         only: alittle
    use interface_module,       only: Int_Geom
    use parameter_module,       only: nicells
    use legacy_mesh_api,        only: nvc, nvf
    use truncate_volume_module, only: Trunc_Vol
    use vof_data_module,        only: Eps

    ! Arguments
    integer :: face
    real(r8), dimension(nicells,nvc) :: UX
    real(r8), dimension(nicells) :: Q

    ! Local Variables
    integer :: k
    real(r8), dimension(nicells) :: J1, J2, N, P, R, S, T, U

    real(r8), parameter :: one_third     = 1.0_r8 / 3.0_r8
    real(r8), parameter :: one_fourth    = 0.25_r8
    real(r8), parameter :: one_fifth     = 0.2_r8
    real(r8), parameter :: one_sixth     = 1.0_r8 / 6.0_r8
    real(r8), parameter :: one_seventh   = 1.0_r8 / 7.0_r8
    real(r8), parameter :: one_eighth    = 0.125_r8

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    P = 0.0_r8
    Q = 0.0_r8
    R = 0.0_r8
    S = 0.0_r8
    T = 0.0_r8
    U = 0.0_r8

    do k = 1, nvf
       where (Trunc_Vol(:,face)%MUp(2) == k)
          R = Trunc_Vol(:,face)%MUi(1+mod(k+1,nvf)) - Trunc_Vol(:,face)%MUa(1)
          U = -UX(:,k)
       endwhere
    end do

    where (ABS(R) > alittle)
       R = 1.0_r8/R
       elsewhere
       R = 0.0_r8
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
       S = 1.0_r8/S
    elsewhere
       S = 0.0_r8
    endwhere

    where (ABS(T) > alittle)
       N =1.0_r8/T
       elsewhere
       N = 0.0_r8
    endwhere

    where ((ABS(T) > 1.0e-2) .and. (T /= -1.0_r8))
       J1 = (1.0_r8 - LOG(ABS(1.0_r8+T))*N)*N
       J2 = (0.5_r8 - J1)*N
    elsewhere
       J1 = one_fourth + T*(-one_fifth + &
            T*(one_sixth + T*(-one_seventh + T*one_eighth)))
       J2 = one_third - T*J1 
       J1 = 0.5_r8 - T*J2
    endwhere

    Q = Q*((Int_Geom%Rho-Trunc_Vol(:,face)%MUa(1))*J1 - &
           (Trunc_Vol(:,face)%MUa(2) - Trunc_Vol(:,face)%MUa(1))*J2)

    Q = Q - 0.5_r8*R*P*(2.0_r8*Int_Geom%Rho - Trunc_Vol(:,face)%MUa(1) - &
                  Trunc_Vol(:,face)%MUa(2))

    Q = Q - 0.5_r8*R*S*U*(Int_Geom%Rho - Trunc_Vol(:,face)%MUa(2))**2

  END SUBROUTINE DVOL_DRO_FACE_2

  SUBROUTINE DVOL_DRO_FACE_N (n, face, UX, Q)
    !=======================================================================
    ! PURPOSE - 
    !   Compute the derivative d(Vol)/d(Ro) at the current hex
    !   face (for case "n"), where Vol is the volume truncated
    !   in the hex by a plane described by the equation 
    !   X*Normal - Ro = 0
    !=======================================================================
    use cutoffs_module,         only: alittle
    use interface_module,       only: Int_Geom
    use parameter_module,       only: nicells
    use legacy_mesh_api,        only: nvc, nvf
    use truncate_volume_module, only: Trunc_Vol
    use vof_data_module,        only: Eps

    ! Arguments
    integer :: n, face
    real(r8), dimension(nicells,nvc) :: UX
    real(r8), dimension(nicells) :: Q

    ! Local Variables
    integer                   :: k
    real(r8), dimension(nicells) :: P, R, S, T, J1

    real(r8), parameter :: one_sixth     = 1.0_r8 / 6.0_r8
    real(r8), parameter :: one_twelfth   = 1.0_r8 / 12.0_r8
    real(r8), parameter :: one_twentieth = 1.0_r8 / 20.0_r8
    real(r8), parameter :: one_thirtieth = 1.0_r8 / 30.0_r8
    real(r8), parameter :: one_42nd      = 1.0_r8 / 42.0_r8

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Q = 0.0_r8

    T = Int_Geom%Rho - Trunc_Vol(:,face)%MUa(n)
    R = Trunc_Vol(:,face)%Lambda + Trunc_Vol(:,face)%Nu*Trunc_Vol(:,face)%MUa(n)
    where (ABS(R) > alittle)
       R = 1.0_r8/R
    elsewhere
       R = 0.0_r8
    endwhere

    S = T*R*Trunc_Vol(:,face)%Nu
    where (ABS(S) > alittle)
       P = 1.0_r8/S
    elsewhere
       P = 0.0_r8
    endwhere

    where (ABS(S) > 1.0e-2 .and. S /= -1.0_r8)
       J1 = (1.0_r8-LOG(ABS(1.0_r8+S))*P)*P 
       J1 = J1 + (J1 - 0.5_r8)*P
    elsewhere
       J1 = one_sixth + S*(-one_twelfth + S*(one_twentieth + S*(-one_thirtieth + S*one_42nd)))
    endwhere

    do k = 1, nvf
       where (Trunc_Vol(:,face)%MUp(n) == k)
          Q = Trunc_Vol(:,face)%V1234*J1*T*R - 0.5_r8*UX(:,k)
          Q = Q*eps(k)*T*T*R
       endwhere
    end do

  END SUBROUTINE DVOL_DRO_FACE_N

END MODULE LOCATE_PLANE_MODULE
