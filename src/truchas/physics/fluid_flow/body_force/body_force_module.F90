!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE BODY_FORCE_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Body force modeling procedures.
  !
  !   Public Interface:
  !     * call ADD_BODY_FORCE_FACE (dt, Mom_Delta)
  !       Increment momentum delta array by the body force.
  !
  ! Contains: add_cell_body_force
  !           ADD_FACE_BODY_FORCE
  !           COMPUTE_GRAVITYHEAD
  !           gravity_limiter
  !
  ! Author(s): Jerry S. Brock (jsbrock@lanl.gov)
  !            Douglas B. Kothe (dbk@lanl.gov)
  !            Jim Sicilian (sicilian@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use legacy_mesh_api, only: ndim
  use truchas_logging_services
  implicit none
  private

  public :: add_cell_body_force, ADD_FACE_BODY_FORCE, COMPUTE_GRAVITYHEAD

CONTAINS

  SUBROUTINE add_cell_body_force(dt, Mom_Delta)
    !=======================================================================
    ! Purpose(s):
    !
    !   Increment momentum delta array by the body force term rho*g.  This 
    !   is the one routine where it doesn't matter whether the Boussinesq
    !   approximation has been invoked.  If (boussinesq), then we use a 
    !   reference density everywhere else in the calculation.
    !=======================================================================
    use body_data_module,  only: Body_Force, body_force_implicitness
    use fluid_data_module, only: fluidRho, fluidDeltaRho, &
                                 fluidvof, avgRho, avgRho_n
    use legacy_mesh_api,   only: ncells

    ! Argument List
    real(r8), intent(IN) :: dt
    real(r8), dimension(ndim,ncells), intent(INOUT) :: Mom_Delta

    ! Local Variables
    integer :: n
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    avgRho = fluidVof*(fluidRho+fluidDeltaRho)

    ! Center the body-force terms based on body_force_implicitness.
    ! Default is 1/2.
    do n = 1,ndim
       Mom_Delta(n,:) = Mom_Delta(n,:) +  &
            ((1.0_r8 - body_force_implicitness)*avgRho_n(:) + &
              body_force_implicitness*avgRho(:))*dt*Body_Force(n)
    end do

    ! Store off the current average density for the next time-step
    avgRho_n = avgRho

    ! Gravity limiter to keep bound on the mechanical energy goes here if it's
    ! necessary.  Right now, it seems not to be necessary.

  end subroutine add_cell_body_force


  SUBROUTINE ADD_FACE_BODY_FORCE(dt, Mom_Delta_Face)
    !=======================================================================
    ! Purpose(s):
    !
    !   Increment momentum delta array by the body force term rho*g.  This 
    !   is the one routine where it doesn't matter whether the Boussinesq
    !   approximation has been invoked.  If (boussinesq), then we use a 
    !   reference density everywhere else in the calculation.
    !=======================================================================
    use body_data_module,  only: Body_Force
    use fluid_data_module, only: fluidRho, fluidDeltaRho, fluidvof, Rho_Face
    use legacy_mesh_api,   only: nfc, ncells, Cell, Mesh, DEGENERATE_FACE, EE_GATHER

    ! Argument List
    real(r8), intent(IN) :: dt
    real(r8), dimension(ndim, nfc, ncells), intent(INOUT) :: Mom_Delta_Face

    ! Local Variables
    integer :: status
    integer :: f, n, nc
    real(r8) :: weight, weight_ngbr, wt
    real(r8), dimension(:,:), allocatable :: DensityRatio, fluidRho_Ngbr,      &
                                              fluidDeltaRho_Ngbr, Volume_Ngbr,  &
                                              fluidvof_Ngbr,Face_Rho_New
    real(r8), dimension(:), allocatable :: Tmp
    logical, dimension(nfc,ncells) :: gravity_off
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE (DensityRatio(nfc,ncells),       &
              fluidRho_Ngbr(nfc,ncells),      &
              fluidDeltaRho_Ngbr(nfc,ncells), &
              Volume_Ngbr(nfc,ncells),        &
              fluidvof_Ngbr(nfc,ncells),      &
              Tmp(ncells),                    &
              Face_Rho_New(nfc,ncells),           STAT = status)
         if (status /= 0) call TLS_panic ('ADD_FACE_BODY_FORCE: allocation failed')

    Tmp = Cell%Volume

    call EE_GATHER (fluidRho_Ngbr, fluidRho)
    call EE_GATHER (fluidDeltaRho_Ngbr, fluidDeltaRho)
    call EE_GATHER (Volume_Ngbr, Tmp)
    call EE_GATHER (fluidvof_Ngbr, fluidvof)

    ! Calculate a face value of delta_rho/rho via a cell volume-weighted average.
    DensityRatio = 1.0_r8
    Face_Rho_New = 0.0_r8
    do n = 1,ncells
       do f = 1,nfc
          if (Mesh(n)%Ngbr_Face(f) == DEGENERATE_FACE) cycle
          if (Mesh(n)%Ngbr_Face(f) == 0) then
             Face_Rho_New(f,n) = fluidRho(n) + fluidDeltaRho(n)
             if (fluidRho(n) > 0.0_r8) then
                 DensityRatio(f,n) = Face_Rho_New(f,n) / fluidRho(n)
             endif
          else
             weight_ngbr = Volume_Ngbr(f,n)*fluidvof_Ngbr(f,n)
             weight      = Tmp(n)*fluidvof(n)
             wt = weight+weight_Ngbr
             if(wt > 0.0_r8) then
                Face_Rho_New(f,n) = ( (fluidRho(n) + fluidDeltaRho(n))*weight  +      &
                        (fluidRho_Ngbr(f,n) + fluidDeltaRho_Ngbr(f,n))*weight_Ngbr ) / wt
                if (Rho_Face(f,n) > 0.0_r8) then
                   DensityRatio(f,n) = Face_Rho_New(f,n) / Rho_Face(f,n)
                end if
             endif
          end if
       end do
    end do


    ! zero out DensityRatio at faces between two void cells, or at void mesh boundary faces
    do f = 1,nfc
       where (fluidRho + fluidRho_Ngbr(f,:) == 0.0_r8) DensityRatio(f,:) = 0.0_r8
    end do

    do f = 1,nfc
       do n = 1,ndim
          Mom_Delta_Face(n,f,:) = dt * DensityRatio(f,:) * Body_Force(n) 
       end do
    end do

    ! where the energy bounds is violated set gravity_off = .true.
    call gravity_limiter(gravity_off)
    ! where gravity_off = .true. kill gravity term
    do nc=1,ncells
       do f=1,nfc  
          if (gravity_off(f,nc)) then
             Mom_Delta_Face(:,f,nc) = 0.0_r8
          end if
       end do
    end do

    DEALLOCATE (DensityRatio)
    DEALLOCATE (fluidRho_Ngbr)
    DEALLOCATE (fluidDeltaRho_Ngbr)
    DEALLOCATE (Volume_Ngbr)
    DEALLOCATE (fluidvof_Ngbr)
    DEALLOCATE (Tmp)
    DEALLOCATE (Face_Rho_New)

  END SUBROUTINE ADD_FACE_BODY_FORCE


  SUBROUTINE COMPUTE_GRAVITYHEAD()

    use fluid_data_module,      only: fluiddeltarho, fluidRho, fluidVof
    use legacy_mesh_api,        only: ncells, ndim, nfc, Cell, Mesh, EE_GATHER
    use projection_data_module, only: Boundary_Flag, dirichlet_pressure, ghc, ghn
    use body_data_module,       only: body_force

    ! Local Variables
    integer :: i
    integer :: f, n
    real(r8), dimension(ncells)          :: rhoc
    real(r8), dimension(nfc,ncells)      :: rhon
    real(r8), dimension(ndim,nfc,ncells) :: Cell_Ngbr_Coord
    logical,  dimension(nfc,ncells)      :: gravity_off
    
    ! Gather cell centroid coordinates for all neighbors
    do i=1,ndim
      call EE_GATHER(Cell_Ngbr_Coord(i,:,:), Cell%Centroid(i))
    end do

    ! Compute the volume-fraction averaged fluid density for the body force
    ! including the change in density due to temperature
    !rhoc = fluidVof*(fluidRho + fluiddeltarho)
    rhoc = (fluidRho + fluiddeltarho)

    ! Gather the density for all neighbors
    call EE_GATHER(rhon, rhoc)

    ! Compute the cell-buoyancy forces to be includes with the face pressure-gradient
    do n=1,ncells
      do f=1,nfc
        ghc(f,n) = 0.0_r8
        ghn(f,n) = 0.0_r8
        do i=1,ndim
           ghc(f,n) = ghc(f,n) - &
              rhoc(n)*body_force(i)*(Cell(n)%Centroid(i)   -Cell(n)%Face_Centroid(i,f))
        end do

        if (Mesh(n)%Ngbr_cell(f) .ne. 0.0_r8) then
           do i=1,ndim
            ghn(f,n) = ghn(f,n) - &
            rhon(f,n)*body_force(i)*(Cell_Ngbr_Coord(i,f,n)-Cell(n)%Face_Centroid(i,f))
           end do
        end if

        ! At pressure boundary conditions, 
        if (dirichlet_pressure) then
           if (Boundary_Flag(f,n) == 1) then
              ghn(f,n) = 0.0_r8
           end if
        end if

      end do
    end do

    ! Now to turn off gravity where necessary.
    call gravity_limiter(gravity_off)

    do n=1,ncells
       do f=1,nfc  
          if (gravity_off(f,n)) then
             ghc(f,n) = 0.0_r8
             ghn(f,n) = 0.0_r8
          end if
       end do
    end do
  
  END SUBROUTINE COMPUTE_GRAVITYHEAD

  SUBROUTINE gravity_limiter(gravity_off)

    use input_utilities,        only: NULL_R
    use legacy_mesh_api,        only: ncells, ndim, nfc, Cell, EE_GATHER
    use zone_module,            only: Zone
    use body_data_module,       only: body_force, mechanical_energy_bound

    ! Argument List
    logical, dimension(:,:), intent(OUT) :: gravity_off

    ! Local Variables
    integer :: i
    integer :: f, n
    real(r8) :: me_ngbr
    real(r8) :: dot, factor
    real(r8), dimension(ncells)     :: ke, pe, me
    real(r8), dimension(nfc,ncells) :: ke_Ngbr, pe_ngbr

    gravity_off = .false.

    if (mechanical_energy_bound /= NULL_R) then

       factor = 0.0_r8
       do n=1,ndim
          factor = factor + Body_Force(n)**2
       enddo
       factor = 0.5_r8 / factor

       ke(:) = 0.0_r8
       do i = 1, ncells
          do n=1,ndim 
              ke(i) = ke(i) + max(0.0_r8,Zone(i)%Vc(n)*Body_Force(n))
          enddo
          ke(i) = factor * ke(i)**2 
       end do

       call EE_GATHER (ke_Ngbr, ke)

       do i=1,ncells
          dot = 0.0_r8
          do n=1,ndim
             dot = dot + body_force(n)*Cell(i)%Centroid(n)
          end do
          ! The potential energy increases in the direction opposite to the body_force
          pe(i) = -dot
       end do

       call EE_GATHER (pe_ngbr, pe)

       me(:) = ke(:) + pe(:)

       do n=1,ncells
          do f=1,nfc
             me_ngbr = ke_ngbr(f,n) + pe_ngbr(f,n)
             if (me(n) > mechanical_energy_bound .or. me_ngbr > mechanical_energy_bound) then
                gravity_off(f,n) = .true.
             end if
          end do
       end do

    end if
  
  END SUBROUTINE gravity_limiter

END MODULE BODY_FORCE_MODULE
