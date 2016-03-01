!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PREDICTOR_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures associated with advancing the velocity field
  !   to a predicted value that is without regard for solenoidality.
  !
  ! Contains: PREDICTOR          (PUBLIC)
  !           PRESSUREEXPLICIT   (PUBLIC)
  !           SOLVE_FOR_VELOCITY (PRIVATE)
  !           PREDICTOR_SETUP    (PRIVATE)
  !           PREDICTOR_CLEANPUP (PRIVATE)
  !
  ! Author(s): Jerry S. Brock (jsbrock@lanl.gov)
  !            Douglas B. Kothe (dbk@lanl.gov)
  !            M. A. Christon (christon@lanl.gov)
  !
  !=======================================================================

  use kinds, only: r8
  implicit none
  private

  public :: PREDICTOR

CONTAINS

  SUBROUTINE PREDICTOR ()
!===============================================================================
! Purpose(s):
!
!   Predictor phase of the incompressible Navier-Stokes solution
!   algorithm. This routine is the main driver for advancing the
!   cell-centered velocity from time "n" to time "*". The "*" time
!   level velocity field is an estimate of the "n+1" velocity field
!   without regard for the solenoidal constraint.
!
!   This routine basically completes the entire RHS for the momentum equations,
!   then calls SOLVE_FOR_VELOCITY to get the predicited velocity, i.e., the *
!   state.
!
!===============================================================================
    use advection_module,        only: ADVECT_MOMENTUM
    use body_data_module,        only: body_force_face_method
    use body_force_module,       only: add_cell_body_force
    use fluid_data_module,       only: fluidRho, fluidRho_n,           &
                                       fluidVof, fluidVof_n,           &
                                       Drag_Coefficient,               &
                                       Mom_Delta,                      &
                                       momentum_solidify_implicitness
    use flow_phase_change,       only: have_solidifying_flow, solidified_rho
    use legacy_mesh_api,         only: ncells, ndim
    use porous_drag_data,        only: porous_flow
    use porous_drag_module,      only: POROUS_DRAG
    use time_step_module,        only: dt
    use timing_tree
    use turbulence_module,       only: turbulence_model, TURBULENCE
    use viscous_data_module,     only: inviscid, stokes, viscous_implicitness
    use viscous_module,          only: viscousExplicit
    use zone_module,             only: Zone
    use surface_tension_module,  only: CSF, surface_tension, csf_tangential, &
                                       csf_boundary

    ! Local Variables
    integer :: i, n
    real(r8) :: tweight
    real(r8), allocatable :: Csftang(:,:)

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start Timer
    call start_timer ("Predictor")

    ! Initialize
    Mom_Delta        = 0
    Drag_Coefficient = 0

    call predictorSetup(dt, Mom_Delta)

    ! Explicit Pressure Gradient and Buoyant Force
    call PressureExplicit(dt, Mom_Delta)

    ! Explicit Viscous Stress.
    if (.not. inviscid) then
       call TURBULENCE(turbulence_model)

       ! Add in the explicit evaluation of the Standard Newtonian stress tensor.
       if(viscous_implicitness < 1) then
          call viscousExplicit(dt, Mom_Delta)
       end if
    end if

    ! Multiply Pressure Gradient and Viscous Stress by FluidVof
    ! (This is a crude way to account for solid material within the cell.
    ! In SOLVE_FOR_VELOCITY we also divide by FluidVof to account for
    ! the mass of fluid in the cell, more or less canceling out this
    ! term.  Momentum advection is specifically excluded because the VOF
    ! is already accounting for the solid material.)
    do i = 1, ncells
       do n = 1, ndim
          Mom_Delta(n,i) = FluidVof(i)*Mom_Delta(n,i)
       end do
    end do

    ! Advect Momentum... No real advection is done here, instead, the increment from
    ! momentum advection is just accumulated into Mom_Delta.
    if (.not. stokes) then
       call ADVECT_MOMENTUM(Mom_Delta)
    end if

    ! Calculate the RHS momentum change due to the solidification of fluid.
    if (HAVE_SOLIDIFYING_FLOW()) then
       ! MAC note: Here, there can be problems with the explicit term of the
       ! momentum  sink due to solidification.  If the fluid velocity at time
       ! level n is zero, there is no contribution to the sink term.  For this
       ! reason, the recommended practice, for now, is to use
       ! momentum_solidify_implicitness = 1.0.  This is the default.
       tweight = 1.0 - momentum_solidify_implicitness
       do i = 1, ncells
          Mom_Delta(:,i) = Mom_Delta(:,i) - tweight*solidified_rho(i)*Zone(i)%Vc
       end do
    end if

    ! Tangential surface tension 
    if (surface_tension .and. (csf_tangential .or. csf_boundary)) then
      ALLOCATE (Csftang(ndim,ncells))

      ! initialize Csftang...
      Csftang = 0.0_r8

      call CSF (dt, Csftang)

      ! Experimental: boundary tangential surface tension (csf_boundary)
      do n=1,ndim
        where (Zone%Rho /= 0.0_r8)
          Csftang(n,:)=Csftang(n,:)*FluidRho*FluidVof/Zone%Rho
        end where
      enddo

      Mom_Delta=Mom_Delta+Csftang
      DEALLOCATE(Csftang)
    end if

    ! Porous drag.
    if (porous_flow) call POROUS_DRAG(dt, Mom_Delta)

    ! Add in the simple cell-centered body force
    if (.not. body_force_face_method) then
       call add_cell_body_force(dt, Mom_Delta)
    endif

    ! Complete the RHS terms in Mom_Delta by including the momentum at time 'n'
    do i = 1, ncells
       do n = 1, ndim
          Mom_Delta(n,i) = Mom_Delta(n,i) + &
               fluidRho_n(i)*fluidVof_n(i)*Zone(i)%Vc_old(n)
       end do
    end do

    ! Solve for Star time velocity 
    call SOLVE_FOR_VELOCITY(Mom_Delta)

    call predictorCleanup()
     
    ! Stop Timer
    call stop_timer("Predictor")

  END SUBROUTINE PREDICTOR

  SUBROUTINE PressureExplicit (dt, Mom_Delta)
!===============================================================================
! Purpose(s):
! 
! Compute the pressure gradient at time-level 'n' for the RHS of the 
! momentum equations.
!
!===============================================================================
    use legacy_mesh_api, only: ncells, ndim
    use fluid_data_module, only: Centered_GradP_Dynamic, fluidRho, fluidRho_n
    
    ! Arguments...
    real(r8), intent(IN) :: dt
    real(r8), dimension(ndim,ncells), intent(INOUT) :: Mom_Delta

    ! Local variables
    integer :: n

    do n = 1,ndim
       where (fluidRho_n == 0) Centered_GradP_Dynamic(n,:) = 0
       where (fluidRho(:) > 0)
          Mom_Delta(n,:) = Mom_Delta(n,:) - &
                           dt*Centered_GradP_Dynamic(n,:)*fluidRho(:)
       endwhere
    end do

  END SUBROUTINE PressureExplicit

  SUBROUTINE SOLVE_FOR_VELOCITY(Mom_Delta)
!===============================================================================
! Purpose(s):
!
!   Solve for the star time velocity either by simply dividing the star
!   time momentum by the new fluid density, or by an implicit solve
!   for the viscous stress
!
!===============================================================================
    use cutoffs_module,          only: cutvof
    use fluid_data_module,       only: fluidRho, fluidVof, &
                                       realfluidVof, cutRho, &
                                       Drag_Coefficient, isPureImmobile,  &
                                       momentum_solidify_implicitness,    &
                                       mass_limiter, mass_limiter_cutoff
    use flow_phase_change,       only: have_solidifying_flow, solidified_rho
    use legacy_mesh_api,         only: ndim, ncells
    use pgslib_module,           only: PGSLIB_GLOBAL_MAXVAL, &
                                       PGSLIB_GLOBAL_MINVAL
    use porous_drag_data,        only: porous_implicitness
    use time_step_module,        only: dt
    use viscous_data_module,     only: inviscid, viscous_implicitness
    use zone_module,             only: zone

    use linear_solution,         only: LINEAR_SOLVER, Ubik_user, &
                                       PRECOND_DIAGONAL
    use preconditioners,         only: PRECONDITION, DIAG_P
    use viscous_data_module,     only: ubik_viscous, viscous_iterations
    use y_eq_Ax_vel,             only: Y_EQ_AX_VELOCITY

    use UbikSolve_module

    ! Argument List
    real(r8), dimension(:,:) :: Mom_Delta

    ! Local Variables
    integer :: n, i, j
    real(r8) :: rhsmin, rhsmax
    real(r8), dimension(ndim*ncells)  :: Solution
    real(r8), dimension(ndim*ncells)  :: RHS
    real(r8), dimension(ndim,ncells)  :: mass
    real(r8), dimension(ncells)       :: rho
    real(r8), dimension(ncells)       :: mscale
    real(r8), dimension(:), pointer   :: diag

    ! Parameters used for the mass limiter. At some point, it may be 
    ! necessary to provide a user interface to adjust the exponential width.  
    ! For now, an on-off toggle has been provided to switch off the 
    ! mass limiter when desired.
    real(r8) :: cmass, width, cutoff, offset, xi

    ! The cutoff should really be cast in terms of a mass rather than volume 
    ! fractions ... but make do with what is here.  The width of the exponential
    ! function is based on 6-decades from the cutvof level of mass by default.
    ! For cutvof values larger than the default, then the width shrinks while
    ! holing the limiter cutoff at 1.0e-2 by default.
    width = mass_limiter_cutoff/cutvof

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Setup the preconditioner as required
    if (.not. inviscid .and. &
        viscous_implicitness/=0 .and. Ubik_User(UBIK_VISCOUS)%precond /= 0) then

       select case(Ubik_User(UBIK_VISCOUS)%precond)

       ! Other cases left for future preconditioners
       case(PRECOND_DIAGONAL)
          allocate(diag(ndim*ncells))
          DIAG_P => diag
          call setupDiagonal(diag)

       end select

    endif

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Setup the volume-fraction averaged fluid density
    rho(:) = fluidRho(:)*fluidVof(:)

    ! Setup the mass limiter if it is enabled (the default).
    if (have_solidifying_flow() .and. mass_limiter) then
       do i = 1, ncells
          ! Setup the mass limiter 
          mscale(i) = 1
          cutoff = width*cutRho(i)
          if(solidified_rho(i) > 0 .and. rho(i) > cutRho(i) .and. rho(i) < cutoff) then
             ! Here, define xi = amp*cutRho(i)/(cutoff - cutRho(i)), but
             ! scale by something roughly proportional to 1/cutRho(i) 
             ! for robust behavior.
             xi = 1/(cutoff - cutRho(i))
             offset = xi + xi*xi/2
             xi = 1/(rho(i) - cutRho(i))
             mscale(i) = 1 + xi + xi*xi/2 - offset
          endif
       end do
    endif

    ! Treat the explicit cases first
    if(viscous_implicitness == 0 .or. inviscid) then

       ! Setup the base 'mass matrix'.  Here, the mass matrix is 
       ! just a generic array used to hold the left-hand-side terms
       ! for the momentum.  These terms can include the average density,
       ! the porous drag, and momentum sink terms due to solidification.
       do i = 1, ncells
          do n = 1, ndim
             mass(n,i) = rho(i)
          end do
       end do

       if (porous_implicitness > 0.0) then
          do i = 1, ncells
             cmass = dt*porous_implicitness*fluidVof(i)
             do n = 1, ndim
                mass(n,i) = mass(n,i) + cmass*Drag_Coefficient(n,i)
             end do
          end do
       endif

       ! Handle the phase change terms -- need time-weighting done correctly
       if (have_solidifying_flow()) then
          do i = 1, ncells
             cmass = momentum_solidify_implicitness*solidified_rho(i)
             do n = 1, ndim
                mass(n,i) = mass(n,i) + cmass
             end do
          end do
       end if

       ! Now calculate the time '*' velocity by ! dividing by the 
       ! mass matrix with all the appropriate terms in it.
       ! By default, the mass limiter version of the predictor is used.
       if (have_solidifying_flow() .and. mass_limiter) then
          do i = 1, ncells
             ! Update the velocities
             do n = 1, Ndim
                ! Overly restrictive conditions, but use for now (MAC)
                if(mass(n,i) > 0 .and. &
                   rho(i) > cutRho(i) .and. realfluidVof(i) > cutvof) then
                   Zone(i)%Vc(n) = Mom_Delta(n,i)/(mass(n,i)*mscale(i))
                else
                   Zone(i)%Vc(n) = 0
                endif
             end do
          end do ! End of Ncells loop

       else
          ! Original update -- can lead to spurious velocities with void present
          do n = 1, ndim
          ! Old treatment of momentum
          ! where (rho > 0) OR generally where (mass(n,:) > 0)  
             where (mass(n,:) > 0)
                Zone%Vc(n) = Mom_Delta(n,:)/mass(n,:)
             elsewhere 
                Zone%Vc(n) = 0
             endwhere
          end do ! End of Ndim loop
       endif

    else

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

       ! Do the linear solve with viscous/drag/solidification effects treated 
       ! implicitly. 
       ! Note: Things are rescaled by 1/FluidVof(:) here.  This affects the
       ! calculation of the terms in setupDiagonal below and in y_eq_A_x.F90.
       !
       ! Why divide Mom_Delta by fluidVof? Does this skew the metrics? (MAC)
       do i = 1, ncells
          do n = 1, ndim
              if(.not. isPureImmobile(i) .and. realFluidVof(i) > cutvof) then
                 Mom_Delta(n,i) = Mom_Delta(n,i)/FluidVof(i)
              else
                 Mom_Delta(n,i) = 0
              endif
          end do
       end do

       ! copy Mom_Delta into the 1D vector RHS for use in linear solve...
       do i = 1, ncells
          do n = 1, ndim
              ! where in X is vel component n for cell i
              j = (i-1)*(ndim)+n
              rhs(j) = Mom_Delta(n,i)
          end do
       end do

       ! Finalize the right-hand-side using either the mass limiter, or
       ! the original enforcement of zero velocity in immobile/void cells
       if (have_solidifying_flow() .and. mass_limiter) then
          ! Update the right-hand-side to acount for the mass scaling
          do i = 1, ncells
             do n = 1, Ndim
                j = (i-1)*ndim+n
                rhs(j) = rhs(j)/mscale(i)
             end do
          end do ! End of Ncells loop
       endif ! End of mass_limiter test

       ! Setup the initial solution for the iterative solver
       do i = 1, ncells
          do n = 1, ndim
             j = (i-1)*ndim+n
             solution(j) = Zone(i)%Vc(n)
          end do
       end do

       ! Skip solution if the RHS is already essentially zero.
       ! BAD Construction here -- partial fix is to the test both the 
       ! min/max values accounting for the possibility that negative 
       ! values are present. 
       !
       ! The previous hard-coded value of 1.0e-10 was reduced to 1.0e-12
       ! to force the viscous solution to iterate.  In general, this
       ! construction is poor, but will have to fix later. (MAC)
       ! Converted hard-coded value to epsilon(0.D0) as the correct way to 
       ! test on a zero rhs.
       rhsmax = pgslib_global_maxval(rhs)
       rhsmin = pgslib_global_minval(rhs)
       rhsmax = max(rhsmax, abs(rhsmin))
       if(rhsmax > epsilon(0.D0)) then
          ! call the linear solver with the y_eq_ax_velocity matvec...
          call LINEAR_SOLVER(Solution, RHS, Ubik_user(UBIK_VISCOUS), &
                             Y_EQ_AX_VELOCITY, PRECONDITION)    
       endif

       ! Store the number of iterations.
       viscous_iterations = Ubik_iter(Ubik_user(ubik_viscous)%control)

       ! Now lets take the 1D solution vector and copy into Zone%Vc...
       do i = 1, ncells
          do n = 1, ndim
              ! where in X is vel component n for cell i
              j = (i-1)*(ndim)+n
              Zone(i)%Vc(n) = Solution(j)
          end do
       end do

       ! Blowup the storage for the preconditioner -- a horrible practice here
       if(.not. inviscid .and. &
          viscous_implicitness/=0 .and. Ubik_User(UBIK_VISCOUS)%precond/=0) then
          select case(Ubik_User(UBIK_VISCOUS)%precond)

          ! Other cases left for future preconditioners
          case(PRECOND_DIAGONAL)
             NULLIFY(DIAG_P)
             deallocate(diag)

          end select
       endif

    endif

  END SUBROUTINE SOLVE_FOR_VELOCITY

  SUBROUTINE predictorSetup(dt, Mom_Delta)
!===============================================================================
!
! Purpose: setup the for the predictor phase -- mostly call viscousSetup
!
!===============================================================================

    use viscous_module, only: viscousSetup
    use legacy_mesh_api, only: ncells, ndim

    ! Arguments...
    real(r8), intent(IN) :: dt
    real(r8), dimension(ndim,ncells), intent(INOUT) :: Mom_Delta

     ! call the setup for viscous physics...
     call viscousSetup(dt, Mom_Delta)

  END SUBROUTINE predictorSetup

  SUBROUTINE predictorCleanup()
!===============================================================================
!
! Purpose: cleanup after the predictor phase
!
!===============================================================================

     use viscous_module, only: viscousCleanup

     ! call the cleanup for viscous physics...
     call viscousCleanup()

  END SUBROUTINE predictorCleanup

  subroutine setupDiagonal(diag)
!===============================================================================
!
! Purpose: setup the diagonal of the Helmholtz operator including porous drag
! terms for the preconditioner
!
!===============================================================================
    use bc_module,            only: BC, FREE_SLIP, Vel
    use bc_operations
    use fluid_data_module,    only: fluidVof, fluidRho, &
                                    Drag_Coefficient, isPureImmobile, &
                                    momentum_solidify_implicitness
    use flow_phase_change,    only: have_solidifying_flow, solidified_rho
    use legacy_mesh_api,      only: ncells, ndim, nfc, Cell
    use time_step_module,     only: dt
    use viscous_data_module,  only: viscous_implicitness, Mu_Face, Mask
    use porous_drag_data,     only: porous_flow, porous_implicitness
 
    ! Arguments
    real(r8), intent(inout) :: diag(:)

    ! Local Variables

    integer :: i, j, f, n
    real(r8), dimension(nfc,ncells) :: beta

    ! Initialize the diagonal
    diag(:) = 0.0

    ! Calculate the geometrical part of the viscous diagonal
    call calcBeta(beta)

    ! Follow the existing implementation for the diagonal terms, even though
    ! this is sort of bogus...  No need for the fluidVof terms here at all!
    ! This is due to the bug fix for porous drag that include the fluidVof term.
    do n = 1, ndim
       ! Handle the viscous terms
       do f = 1, nfc
          mask(:) = .false.
          mask    = (FREE_SLIP(BC%Flag, Vel%Face_bit(f)))
          do i = 1, ncells
             j = (i-1)*ndim + n
             if (mask(i)) then 
               diag(j) = diag(j) - &
               Mu_Face(f,i)*beta(f,i)*cell(i)%Face_Normal(n,f)*cell(i)%Face_Normal(n,f)
             else
               diag(j) = diag(j) - Mu_Face(f,i)*beta(f,i)
             endif
          end do
       end do
    end do

    do i = 1, ncells
       do n = 1, ndim
          j = (i-1)*ndim + n
          diag(j) = diag(j)*dt*viscous_implicitness/Cell(i)%Volume
       end do
    end do

    ! Inertial, porous drag, and solidification terms
    ! Handle the inertial terms first
    do i = 1, ncells
       do n = 1, ndim
          j = (i-1)*ndim+n
          diag(j) = diag(j) + fluidRho(i)
       end do
    end do
    ! Porous drag terms
    if (porous_flow) then
       do i = 1, ncells
          do n = 1, ndim
             j = (i-1)*ndim+n
             diag(j) = diag(j) + dt*porous_implicitness*Drag_Coefficient(n,i)
          end do
       end do
    endif
    ! Momentum solidification terms
    if (have_solidifying_flow()) then
       do i = 1, ncells
          do n = 1, ndim
             j = (i-1)*ndim+n
             if (fluidVof(i) > 0) then
                diag(j) = diag(j) + &
                     momentum_solidify_implicitness*solidified_rho(i)/fluidVof(i)
             endif
          end do
       end do
    end if

    ! For an immobile material where the viscosity/drag is zero,
    ! this provides forces the diagonal and assumes that the RHS entry 
    ! has been set to zero.   In Solve_For_Velocity, we enforce this
    ! constraint on the RHS for consistency.  Note an alternative would
    ! be to use a penalty enforcement of the condition on velocity.
    !
    ! Note: Also need to trap on the case of all-void cells.  This is a 
    ! little bogus, but essentially forces a sane diagonal (unity) for
    ! void cells.  There is a patch-up that occurs the in mat-vec call back
    ! y_eq_Ax_vel.F90
    do i = 1, ncells
       do n = 1, ndim
          j = (i-1)*ndim+n
          if (isPureImmobile(i) .or. diag(j) == 0) diag(j) = 1 
       end do
    end do

  end subroutine setupDiagonal

  subroutine calcBeta(beta)
!===============================================================================
!
! Purpose: precompute the constant terms for the diagonal preconditioner
!
!===============================================================================

  use legacy_mesh_api,      only: ncells, nfc, Cell, EE_GATHER
  use bc_module,            only: BC, FREE_SLIP, DIRICHLET_VEL, DIRICHLET, &
                                  Vel, Prs
  use bc_operations
  use viscous_data_module,  only: mask

  ! Arguments
  real(r8), intent(inout) :: beta(:,:)

  ! Local variables
  logical, dimension(ncells) :: mask1
  integer :: i, f
  real(r8) :: dx, dy, dz, delta_mag

  ! This could be pre-allocated for the diagonal preconditioner a-priori, 
  ! but for now, leave it local to the routine as an automatic (MAC)
  real(r8), dimension(nfc,ncells) :: xnbr, ynbr, znbr

  beta(:,:) = 0.0

  ! Do the silly gather for face-neigbor data for the coordinates
  call EE_GATHER(xnbr, Cell%Centroid(1))
  call EE_GATHER(ynbr, Cell%Centroid(2))
  call EE_GATHER(znbr, Cell%Centroid(3))

  do f = 1, nfc

     ! For free-slip, Face_Velocity is approximated by the component of
     ! the neighbouring Zone%Vc_Old tangential to the face;
     ! Calculate the normal component first
     mask1 = FREE_SLIP(BC%Flag, Vel%Face_bit(f)) .or. &
             DIRICHLET(BC%Flag, Prs%Face_bit(f))

     ! At Dirichlet pressure boundary faces use the tangential component of the
     ! cell centered velocity and the Fluxing_Velocity as the normal component
     mask = DIRICHLET(BC%Flag, Prs%Face_bit(f))

     ! At Dirichlet velocity boundary faces use specified boundary velocity
     mask = DIRICHLET_VEL(BC%Flag, Vel%Face_bit(f))

     mask = mask1 .or. mask

     ! Use a do-loop to avoid promoting scalars to vectors for now
     do i = 1, ncells
        if (mask(i)) then
           ! Compute beta at a boundary -- solid or otherwise
           dx = Cell(i)%Centroid(1) - Cell(i)%Face_Centroid(1,f)
           dy = Cell(i)%Centroid(2) - Cell(i)%Face_Centroid(2,f)
           dz = Cell(i)%Centroid(3) - Cell(i)%Face_Centroid(3,f)
           delta_mag = dx*dx + dy*dy + dz*dz
        else
           ! Compute beta at internal faces
           dx = Cell(i)%Centroid(1) - xnbr(f,i)
           dy = Cell(i)%Centroid(2) - ynbr(f,i)
           dz = Cell(i)%Centroid(3) - znbr(f,i)
           delta_mag = dx*dx + dy*dy + dz*dz
        endif
        beta(f,i) = (dx*Cell(i)%Face_Normal(1,f) +  &
                     dy*Cell(i)%Face_Normal(2,f) +  &
                     dz*Cell(i)%Face_Normal(3,f))*Cell(i)%Face_Area(f)/delta_mag
     end do
     
  end do ! End face loop

  end subroutine calcBeta

END MODULE PREDICTOR_MODULE
