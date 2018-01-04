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
!
! Author(s): Stewart J. Mosso, LANL X-HM (sjm@lanl.gov)
!            Douglas B. Kothe, LANL (dbk@lanl.gov)
!
!=======================================================================

#include "f90_assert.fpp"

module locate_plane_modulez

  use kinds,     only: r8
  use hex_types, only: reconstruction_hex
  use truchas_logging_services
  implicit none
  private

  real(r8), parameter :: volume_track_iter_tol = 1.0d-8

  type, extends(reconstruction_hex), public :: locate_plane_hex
  contains
    procedure :: init => init_locate_plane_hex
    procedure :: locate_plane
    procedure, private :: rho_bracket
    procedure, private :: rho_brent
  end type locate_plane_hex

contains

  subroutine init_locate_plane_hex (this, norm, vof, volume, node)
    class(locate_plane_hex), intent(out) :: this
    real(r8),                intent(in)  :: norm(:), vof, volume, node(:,:)

    this%P%normal = norm
    this%vof      = vof
    this%volume   = volume
    this%node     = node

  end subroutine init_locate_plane_hex

  !=======================================================================
  ! Purpose(s):
  !
  !   Given the equation for interface planes: X*Normal - Rho = 0,
  !   compute the parameter Rho (given the plane Normal) subject
  !   to a volume conservation constraint.  This is done in
  !   an iterative fashion, using either Brent's
  !   method.  The iteration is first bracketed by finding the
  !   interval [Rho_Min,Rho_Max] bracketing Rho.  Given Rho and the
  !   volume flux geometry, the plane equation parameters
  !   (Rho, Normal) are then used to compute the truncated
  !   flux volume (Vp).
  !
  !=======================================================================
  subroutine locate_plane (this, iter)
    use truncate_volume_modulez, only: truncvol_data

    class(locate_plane_hex), intent(inout) :: this
    integer,                 intent(out)   :: iter

    real(r8)            :: Rho_Min, Rho_Max, V_Min, V_Max
    type(truncvol_data) :: trunc_vol(6)

    ! WARNING: Need to figure out how to make the timer run in parallel with OpenMP.
    !          Currently, it is a global variable, and not thread-safe.
    !call start_timer ("Locate Plane")

    ! Bracket the correct value of Rho [Rho_Min,Rho_Max] to insure
    ! the subsequent iteration will converge efficiently
    call this%rho_bracket (Rho_Min, Rho_Max, V_Min, V_Max, trunc_vol)

    ! Compute Rho, which parameterizes the plane
    ! characterized by the equation Normal*X - Rho = 0
    ! using Brents method iteration
    call this%rho_brent (iter, Rho_Min, Rho_Max, V_Min, V_Max, trunc_vol, 100)

    !call stop_timer("Locate Plane")

  end subroutine locate_plane

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  !=======================================================================
  ! PURPOSE -
  !   Compute the plane parameter Rho in the plane equation:
  !                    X*Normal - Rho = 0
  !   iteratively using Brents method as documented in Numerical
  !   Recipes (Press, Flannery, Teukolsky and Vetterling; Cambridge
  !   University Press, 1986), p. 253
  !=======================================================================
  subroutine rho_bracket (this, Rho_Min, Rho_Max, V_Min, V_Max, trunc_vol)
    use truncate_volume_modulez, only: truncate_volume,face_param,truncvol_data

    class(locate_plane_hex), intent(inout) :: this
    real(r8),                intent(out) :: Rho_Min, Rho_Max, V_min, V_max
    type(truncvol_data),     intent(out) :: trunc_vol(:)

    integer  :: f,v,i
    real(r8) :: V_v(8), tmp

    ! Initialize the bracketed values of truncation volume [V_Min,V_Max]
    ! and the associated plane constant [Rho_Min,Rho_Max]

    ! Set up the limits on the volume bracketing. Make sure that
    ! V_v is bounded zero = V_Min <= V_v <= V_Max = Cell Volume.
    ! If V_v > Cell Volume (because of roundoff), set V_v = Cell Volume
    V_Min = 0.0_r8
    V_Max = This%Volume

    ! Compute volumes truncated by the interface passing through each of the vertices.
    do f = 1,6
      ! Compute and store face parameters.
      trunc_vol(f) = face_param (this, 'full_cell', f)
    end do

    do v = 1,8
      ! Get the value of Rho for this vertex.
      this%P%rho = sum(this%node(:,v) * this%P%normal(:))

      ! get this vertex's truncation volume, and force it to lie between V_min and V_max
      V_v(v) = min(max(truncate_volume (this, trunc_vol), V_min), V_max)
    end do

    ! Make sure that V_v = V_max for at least one vertex,
    ! and            V_v = V_min for at least one vertex
    ! If not, then set the maximum V_v to V_max
    ! and              the minimum V_v to V_min
    if (maxval(V_v) < V_max) V_v( maxloc(V_v, dim=1) ) = V_max
    if (minval(V_v) > V_min) V_v( minloc(V_v, dim=1) ) = V_min

    ! Now bracket Rho
    Rho_Max =  huge(0.0_r8)
    Rho_Min = -huge(0.0_r8)

    do v = 1,8
      if ((V_min < V_v(v) .and. V_v(v) <= this%vof*this%volume) .or. V_v(v) == V_min) then
        Rho_Min = sum(this%node(:,v)*this%P%normal) ! the value of rho for this vertex
        V_Min = V_v(v)
      else if ((this%vof*this%volume < V_v(v) .and. V_v(v) < V_max) .or. V_v(v) == V_max) then
        Rho_Max = sum(this%node(:,v)*this%P%normal) ! the value of rho for this vertex
        V_Max = V_v(v)
      end if
    end do

    ! At this point, [V_Min,V_Max] brackets V_v, and [Rho_Min,Rho_Max]
    ! are those plane constants corresponding to [V_Min,V_Max]. The
    ! bracketed volumes are monotonic, but the Rho's aren't necessarily
    ! monotonic, so make sure that Rho_Min < R_Max. Be sure to interchange
    ! V_Min and V_Max if R_Min and R_Max have to be interchanged.
    if (Rho_max < Rho_min) then
      tmp = Rho_max
      Rho_max = Rho_min
      Rho_min = tmp

      tmp = V_max
      V_max = V_min
      V_min = tmp
    end if

    ! if we haven't properly set rho everywhere, punt
    if (Rho_min == -huge(0.0_r8) .or. Rho_max == huge(0.0_r8)) &
         call TLS_fatal ('RHO_BRACKET: unable to bracket plane constant Rho')

  end subroutine rho_bracket

  !=======================================================================
  ! PURPOSE -
  !   Compute the plane parameter Rho in the plane equation:
  !                    X*Normal - Rho = 0
  !   iteratively using Brents method as documented in Numerical
  !   Recipes (Press, Flannery, Teukolsky and Vetterling; Cambridge
  !   University Press, 1986), p. 253
  !=======================================================================
  subroutine rho_brent (this, iter, Rho_Min, Rho_Max, V_Min, V_Max, trunc_vol, iter_max)
    use truncate_volume_modulez, only: truncate_volume, truncvol_data

    class(locate_plane_hex), intent(inout) :: this
    integer,                 intent(out)   :: iter
    real(r8),                intent(in)    :: Rho_Min, Rho_Max, V_Min, V_Max
    type(truncvol_data),     intent(in)    :: trunc_vol(:)
    integer,                 intent(in)    :: iter_max

    integer :: i
    real(r8) :: Rho_a, Rho_b, Rho_c, Rho_d, Rho_e, Rho_mid, Rho_tol, V_a, V_b, V_c, P, Q, R, S

    ! Initialize values before the iteration loop
    i = 0; iter = 0
    Rho_a = Rho_Min; Rho_b = Rho_Max
    Rho_c = 0.0_r8; Rho_d = 0.0_r8; Rho_e = 0.0_r8
    P = 0.0_r8; Q = 0.0_r8; R = 0.0_r8; S = 0.0_r8
    V_a = V_Min - this%Vof*this%Volume
    V_b = V_Max - this%Vof*this%Volume
    V_c = V_b

    ! Main iteration loop
    do while (abs(V_b/max(this%Volume,epsilon(1.0_r8))) > volume_track_iter_tol &
        .and. i < iter_Max)
      if (V_b*V_c > 0.0_r8) then
        Rho_c = Rho_a
        V_c = V_a
        Rho_d = Rho_b - Rho_a
        Rho_e = Rho_d
      end if

      if (abs(V_c) < abs(V_b)) then
        Rho_a = Rho_b
        Rho_b = Rho_c
        Rho_c = Rho_a
        V_a = V_b
        V_b = V_c
        V_c = V_a
      end if

      Rho_tol = 2.0_r8*epsilon(1.0_r8)*abs(Rho_b) + 0.5_r8*volume_track_iter_tol
      Rho_mid = 0.5_r8*(Rho_c - Rho_b)

      if (abs(Rho_mid) > Rho_tol .and. V_b /= 0.0_r8) iter = i

      if (abs(Rho_e) >= Rho_tol .and. abs(V_a) > abs(V_b)) then
        S = V_b/V_a

        if (Rho_a == Rho_c) then
          P = 2.0_r8*Rho_mid*S
          Q = 1.0_r8 - S
        else
          Q = V_a/V_c
          R = V_b/V_c
          P = S*(2.0_r8*Rho_mid*Q*(Q - R) - (Rho_b - Rho_a)*(R - 1.0_r8))
          Q = (Q - 1.0_r8)*(R - 1.0_r8)*(S - 1.0_r8)
        end if

        if (P > 0.0_r8)  Q = -Q

        P = abs(P)
        R = min(3.0_r8*Rho_mid*Q - abs(Rho_tol*Q), abs(Rho_e*Q))

        if (2.0_r8*P < R) then
          Rho_e = Rho_d
          Rho_d = P/Q
        else
          Rho_d = Rho_mid
          Rho_e = Rho_d
        end if
      else
        Rho_d = Rho_mid
        Rho_e = Rho_d
      end if

      Rho_a = Rho_b ! Move last best guess to A
      V_a  = V_b

      if (abs(Rho_d) > Rho_tol) then ! Evaluate new trial root
        Rho_b = Rho_b + Rho_d
      else
        Rho_b = Rho_b + sign(Rho_tol,Rho_mid)
      end if
      this%P%Rho = Rho_b

      V_b = truncate_volume (this, trunc_vol) - this%Vof*this%Volume
      i = i + 1
    end do

    this%P%Rho = Rho_b

    ! if (i==iter_max) then
    !   ! too many iterations if we get here
    !   ! write(*,'(a,2es14.4)') 'error:  ',fx,x
    !   ! write(*,'(a,3es14.4)') 'bounds: ',x_min,x_mid,x_max
    !   call LS_fatal('too many brent iterations!')
    ! end if
  end subroutine rho_brent

end module locate_plane_modulez
