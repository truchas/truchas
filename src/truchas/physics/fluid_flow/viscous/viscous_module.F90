!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE VISCOUS_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define stress related variables and routines.
  !
  ! Contains: VISCOUS
  !           STRESS_GRADIENT
  !
  ! Author(s): Jerry S. Brock (jsbrock@lanl.gov)
  !            Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  ! Public Subroutines
  public :: viscousExplicit
  public :: viscousSetup
  public :: viscousCleanup
  public :: STRESS_GRADIENT

CONTAINS

  SUBROUTINE viscousSetup (dt, Mom_Delta)
    !=======================================================================
    ! Purpose(s):
    !   Evaluate a cell-center momentum delta due to viscous stress.
    !
    !   Mom_Delta = dt * Stress Gradient
    !======================================================================= 
    use legacy_mesh_api,      only: ncells, ndim, nfc
    use viscous_data_module,  only: Stress_Grad_BC, viscous_implicitness

    use viscous_data_module,  only: Mask, Grad, Mu_Face, &
                                    Normal, Face_Velocity, inviscid

    ! Argument List
    real(r8), intent(IN) :: dt
    real(r8), dimension(ndim,ncells), intent(INOUT) :: Mom_Delta

    ! Local Variables
    integer :: status
    integer :: n
    real(r8) :: tweight

    real(r8), dimension(ndim,ncells) :: Vc

    ! Don't do anything if this isn't a viscous simulation
    !=====================================================
        IF(INVISCID) RETURN
    !=====================================================

    ALLOCATE (Mask(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('VISCOUS: Mask(ncells) allocation failed')
    ALLOCATE (Grad(ndim,nfc,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('VISCOUS: Grad(ndim,nfc,ncells) allocation failed')
    ALLOCATE (Mu_Face(nfc,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('VISCOUS: Mu_Face(nfc,ncells) allocation failed')
    ALLOCATE (Normal(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('VISCOUS: Normal(ncells) allocation failed')
    ALLOCATE (Face_Velocity(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('VISCOUS: Face_Velocity(ncells) allocation failed')
    ALLOCATE (Stress_Grad_BC(ndim,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('VISCOUS: Stress_Grad_BC(ndim,ncells) allocation failed')

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Calculate cell face viscosity using harmonic mean...
    call calculateCellFaceViscosity ()

    ! Don't need the prepass if the viscous terms are explicit
    if(viscous_implicitness > 0) then

       Vc = 0
       call STRESS_GRADIENT(Stress_Grad_BC, Vc)

       ! Increment Momentum Delta
       tweight = (1 - viscous_implicitness)
       do n = 1, ndim
          Mom_Delta(n,:) = Mom_Delta(n,:) + dt*tweight*Stress_Grad_BC(n,:)
       end do

    endif

  END SUBROUTINE viscousSetup

  SUBROUTINE viscousExplicit (dt, Mom_Delta)
    !=======================================================================
    ! Purpose(s):
    !   Evaluate a cell-center momentum delta due to viscous stress.
    !
    !   Mom_Delta = dt * Stress Gradient
    !======================================================================= 
    use legacy_mesh_api,      only: ncells, ndim
    use zone_module,          only: Zone
    use viscous_data_module,  only: viscous_implicitness, Stress_Grad_BC

    ! Argument List
    real(r8), intent(IN) :: dt
    real(r8), dimension(ndim,ncells), intent(INOUT) :: Mom_Delta

    ! Local Variables
    integer :: status
    integer :: n
    real(r8) :: tweight
    real(r8), dimension(:,:), allocatable :: Stress_Grad

    real(r8), dimension(ndim,ncells) :: Vc

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start Timer
    call start_timer ("Viscous")

    ALLOCATE (Stress_Grad(ndim,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('VISCOUS: Stress_Grad(ndim,ncells) allocation failed')

    ! Copy a velocity field from Zone into a ncells by ncells array...
    do n = 1, ndim
       Vc(n,:) = Zone(:)%Vc_Old(n)
    end do

    ! Stress_Grad = Sum(Tau*Normal*Area)/Volume
    call STRESS_GRADIENT(Stress_Grad, Vc)

    ! Increment Momentum Delta
    tweight = (1 - viscous_implicitness)
    do n = 1, ndim
       Mom_Delta(n,:) = Mom_Delta(n,:) + dt*tweight*Stress_Grad(n,:)
       ! See Viscous StressDescription in the flow modeling docs on source forge
       if(viscous_implicitness > 0) then
          Mom_Delta(n,:) = Mom_Delta(n,:) - dt*tweight*Stress_Grad_BC(n,:)
       endif
    end do

    DEALLOCATE (Stress_Grad)

    ! Stop Timer
    call stop_timer ("Viscous")

  END SUBROUTINE viscousExplicit

  SUBROUTINE viscousCleanup ()

    !=======================================================================
    ! Purpose(s):
    !   Cleanup viscous physics, deallocat memory etc.
    !
    !   
    !======================================================================= 

    use viscous_data_module,    only: Mask, Grad, Mu_Face, inviscid,        &
                                      Normal, Face_Velocity, Stress_Grad_BC

    ! Don't do anything if this isn't a viscous simulation
    !=====================================================
        IF(INVISCID) RETURN
    !=====================================================

    DEALLOCATE (Mask)
    DEALLOCATE (Grad)
    DEALLOCATE (Mu_Face)
    DEALLOCATE (Normal)
    DEALLOCATE (Face_Velocity)
    DEALLOCATE (Stress_Grad_BC)
     
  END SUBROUTINE viscousCleanup

  SUBROUTINE calculateCellFaceViscosity ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Evaluate the cell-center stress gradient.
    !
    !                    ---
    !                    \
    !                     |  tau_f * normal_f * area_f
    !                    /
    !                    ---
    !   Stress_Grad  =  -------------------------------
    !
    !                               Volume_c
    !
    !                    ---
    !                    \
    !                     |  mu_f * grad_f * normal_f * area_f
    !                    /
    !                    ---
    !                =  ---------------------------------------
    !
    !                               Volume_c
    !
    !=======================================================================
    use fluid_data_module,    only: Solid_Face
    use legacy_mesh_api,      only: ncells, nfc, Mesh, DEGENERATE_FACE, EE_GATHER
    use flow_property_module,      only: get_viscosity
    use truchas_logging_services
    use viscous_data_module,  only: Mu_Face

    ! Local Variables
    integer :: status
    integer :: i, f

    real(r8), dimension(:,:), allocatable :: Mu_Cell_Ngbr
    real(r8), dimension(:),   allocatable :: Mu_Cell

    ! Allocate the memory...
    ALLOCATE (Mu_Cell_Ngbr(nfc,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('VISCOUS: Mu_Cell_Ngbr(nfc,ncells) allocation failed')
    ALLOCATE (Mu_Cell(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('VISCOUS: Mu_Cell(ncells) allocation failed')

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !
    ! Cell-Face Viscosity
    ! Markus - Turns out that FACE_CENTER is too large a stencil for viscosities near
    !          an interface; in particular, a cell with very little real fluid (unlike
    !          void) will have a small volume-averaged cell viscosity, but will have
    !          huge face viscosities based on neighboring cell-centered values.  For 
    !          now, I'll comment out FACE_CENTER and do a simple ortho-like harmonic 
    !          mean calculation.
    !
    call get_viscosity (Mu_Cell)
    ! call FACE_CENTER (Mu_Cell, Mu_Face)
    call EE_GATHER (Mu_Cell_Ngbr, Mu_Cell)
    do i = 1,ncells
       if (Mu_Cell(i) == 0) then
          Mu_Face(:,i) = 0
       else 
          do f = 1,nfc
             if (Mesh(i)%Ngbr_Cell(f) == DEGENERATE_FACE) then 
                Mu_Face(f,i) = Mu_Cell(i)
             else if (Solid_Face(f,i)) then 
                Mu_Face(f,i) = Mu_Cell(i)
             else if (Mesh(i)%Ngbr_Cell(f) == 0) then
                Mu_Face(f,i) = Mu_Cell(i)
             else if (Mu_Cell(i) == 0 .OR. Mu_Cell_Ngbr(f,i) == 0) then
                Mu_Face(f,i) = 0
             else
                Mu_Face(f,i) = 2 / (1/Mu_Cell(i) + 1/Mu_Cell_Ngbr(f,i))
             endif
          end do
       endif
    end do

    ! deallocate memory...
    DEALLOCATE (Mu_Cell_Ngbr)
    DEALLOCATE (Mu_Cell)

  END SUBROUTINE calculateCellFaceViscosity

  SUBROUTINE STRESS_GRADIENT (Stress_Grad, Vc)
    !=======================================================================
    ! Purpose(s):
    !
    !   Evaluate the cell-center stress gradient.
    !
    !                    ---
    !                    \
    !                     |  tau_f * normal_f * area_f
    !                    /
    !                    ---
    !   Stress_Grad  =  -------------------------------
    !
    !                               Volume_c
    !
    !                    ---
    !                    \
    !                     |  mu_f * grad_f * normal_f * area_f
    !                    /
    !                    ---
    !                =  ---------------------------------------
    !
    !                               Volume_c
    !
    !=======================================================================
    use bc_module,            only: BC, bndry_vel, FREE_SLIP, DIRICHLET_VEL, & ! BC_Vel
                                    DIRICHLET, Vel, Prs
    use bc_operations
    use do_interface,         only: DO_Specifier, do_init_ss, do_update_weights, &
                                    do_gradient_face, DO_SOLVE_ORTHO, &
                                    DO_SOLVE_LU_LSLR
    use discrete_ops_data,    only: use_ortho_face_gradient
    use fluid_data_module,    only: fluidVof
    use legacy_mesh_api,      only: ncells, ndim, nfc, Cell
    use tensor_module,        only: Tensor
    use viscous_data_module,  only: Mask, Grad, Mu_Face, &
                                    Normal, Face_Velocity
    use time_step_module,     only: t

    ! Argument List
    real(r8), dimension(ndim,ncells), intent(OUT) :: Stress_Grad
    real(r8), dimension(:,:),         intent(IN)  :: Vc

    ! Local Variables
    type(DO_Specifier),pointer,save :: SG_SS =>NULL()
    integer :: j, f, m, n, SGSolveTech
    logical, dimension(ncells) :: Mask1

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize Variables
    if(.not. ASSOCIATED(SG_SS))then
      SGSolveTech=DO_SOLVE_LU_LSLR; if(use_ortho_face_gradient)SGSolveTech=DO_SOLVE_ORTHO
      call do_init_ss(SG_SS,SOLVETECH=SGSolveTech,BC_SPEC=Pressure_BC,WEIGHTS=fluidVof)
    else
      call do_update_weights(SOLVESPEC=SG_SS,WEIGHTS=fluidVof)
    endif
    Stress_Grad = 0

    ! Stress Gradient
    NDIM_LOOP : do n = 1, ndim

       call DO_GRADIENT_FACE(PHI=Vc(n,:), SOLVESPEC=SG_SS, GRAD=Grad)

       ! Cell-Face Stress Gradients
       FACE_LOOP : do f = 1, nfc

          ! Apply BCs

          ! Define a Face_Velocity at BC faces
          Face_Velocity = 0

          ! For free-slip, Face_Velocity is approximated by the component of
          ! the neighbouring Zone%Vc_Old tangential to the face;
          ! Calculate the normal component first
          Mask1 = FREE_SLIP (BC%Flag, Vel%Face_bit(f))
          Normal = 0
          do m = 1,ndim
             where (Mask1) Normal = Normal + Vc(m,:)*Cell%Face_Normal(m,f)
          end do
          ! Calculate the tangential component as (Zone%Vc_Old - Normal)
          where (Mask1) Face_Velocity = Vc(n,:) - Normal*Cell%Face_Normal(n,f)

          ! At Dirichlet pressure boundary faces set the velocity gradient to zero
          Mask = DIRICHLET (BC%Flag, Prs%Face_bit(f))
          where (Mask) Face_Velocity = Vc(n,:)
          Mask1 = Mask1 .or. Mask

          ! At Dirichlet velocity boundary faces simply use the specified boundary velocity
          Mask = DIRICHLET_VEL (BC%Flag, Vel%Face_bit(f))
          !! NNC, Jan 2014.  Time-dependent dirichlet velocity.
          !ORIG: where (Mask) Face_Velocity = BC_Vel(n,f,:)
          do j = 1, ncells
             if (Mask(j)) Face_Velocity(j) = bndry_vel%get(n,f,j,t)
          end do

          ! Calculate a face gradient based on Face_Velocity and the neighbouring Zone%Vc_Old.
          Mask = Mask1 .or. Mask
          do m = 1,ndim
             where (Mask) Grad(m,f,:) = (Face_Velocity-Vc(n,:))/Cell%Halfwidth(f)*Cell%Face_Normal(m,f)
          end do

          ! Stress Gradient Tensor
          if (use_ortho_face_gradient) then
             do m = 1, ndim
                ! Stress Gradient Components
                Stress_Grad(n,:) = Stress_Grad(n,:) + Mu_Face(f,:) * Grad(m,f,:) &
                                   * Cell%Face_Normal(m,f) * Cell%Face_Area(f)
             end do
          else
             TENSOR_LOOP : do m = 1, ndim
                ! Tensor Matrix Element
                j = Tensor(n,m)

                ! Stress Gradient Components
                Stress_Grad(n,:) = Stress_Grad(n,:) + Mu_Face(f,:) * Grad(m,f,:) &
                                      * Cell%Face_Normal(m,f) * Cell%Face_Area(f)
                Stress_Grad(j,:) = Stress_Grad(j,:) + Mu_Face(f,:) * Grad(j,f,:) &
                                      * Cell%Face_Normal(n,f) * Cell%Face_Area(f)
             end do TENSOR_LOOP
          endif

       end do FACE_LOOP
    end do NDIM_LOOP

    ! Volume Normalization
    VOLUME_LOOP : do n = 1, ndim
       Stress_Grad(n,:) = Stress_Grad(n,:) / Cell%Volume
    end do VOLUME_LOOP

  END SUBROUTINE STRESS_GRADIENT

END MODULE VISCOUS_MODULE
