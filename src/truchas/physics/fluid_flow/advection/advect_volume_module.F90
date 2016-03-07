!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ADVECT_VOLUME_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures that advance material volumes and volume
  !   fractions according to a standard advection equation.
  !
  ! Public Interface(s):
  !
  !   * call ADVECT_VOLUME (Vof, Vof_n, Fluxing_Velocity, Volume_Flux_Tot)
  !
  !     Advect material volumes with a standard advection equation.
  !     Total volume fluxes at faces are computed with the fluxing
  !     velocity provided, and material volume fluxes at faces are
  !     computed with a piecewise-planar volume tracking algorithm.
  !
  !   * call VOF_BOUNDS (Vof, Volume_Flux_Tot)
  !
  ! Contains: ADVECT_VOLUME
  !           ADVECT_CONTINUUM
  !           CONTINUUM_ADVANCE
  !           FLUX_ACCEPTOR
  !           FLUX_BC
  !           FLUX_RENORM
  !           VOLUME_ADVANCE
  !           VOF_BOUNDS
  !           ADJUST_VOFS
  !           ADJUST_FLUX_MATL
  !           ADJUST_FLUX_TOTAL
  !
  ! Author(s): Douglas B. Kothe (LANL, dbk@lanl.gov)
  !            S. Jay Mosso (LANL, sjm@lanl.gov)
  !            Markus Bussmann (bussmann@lanl.gov)
  !            Jim Sicilian (CCS-2, sicilian@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services  ! entities prefixed with TLS_
  implicit none
  private

  public :: ADVECT_VOLUME

  ! Private data
  integer,  save :: WLimit     = 10
  integer,  save :: WCountTot  = 0
  real(r8), save :: WMaxTot    = 0.0_r8
  integer,  save :: WCountMat  = 0
  real(r8), save :: WMaxMat    = 0.0_r8
  integer,  save :: WCountTotU = 0
  real(r8), save :: WMaxTotU   = 0.0_r8
  integer,  save :: WCountMatU = 0
  real(r8), save :: WMaxMatU   = 0.0_r8

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><>

  SUBROUTINE ADVECT_VOLUME (Vof, Vof_n, Fluxing_Velocity, Volume_Flux_Tot)
    !=======================================================================
    ! Purpose(s):
    !
    !   Integrate the volume advection equation in time by accumulating
    !   the volume fluxes.
    !
    !=======================================================================
    use parameter_module,    only: nmat
    use legacy_mesh_api,     only: nfc, ncells
    use time_step_module,    only: dt
    use timing_tree
    use vof_data_module,     only: adv_dt, volume_track_subcycles, volume_track_interfaces
    use volume_track_module, only: VOLUME_TRACK

    ! Arguments
    real(r8), dimension(nmat,ncells),     intent(INOUT) :: Vof
    real(r8), dimension(nmat,ncells),     intent(IN)    :: Vof_n
    real(r8), dimension(nfc,ncells),      intent(IN)    :: Fluxing_Velocity
    real(r8), dimension(nmat,nfc,ncells), intent(OUT)   :: Volume_Flux_Tot

    ! Local Variables
    integer :: status
    integer :: p, vps
    real(r8), dimension(:,:,:), allocatable :: Volume_Flux_Sub

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the volume advection timer.
    call start_timer("Volume Tracking")

    ! Zero the total flux array.
    Volume_Flux_Tot = 0.0_r8

    ! No subcycling if we're not volume tracking.
    if (volume_track_interfaces) then
       vps = volume_track_subcycles
    else
       vps = 1
    end if

    ! Set the advection timestep.
    adv_dt = dt/vps

    ! If we're volume tracking, we'll need an array to keep track of volume
    ! changes in every subcycle.
    if (volume_track_interfaces) then
       ALLOCATE (Volume_Flux_Sub(nmat,nfc,ncells), STAT = status)
       if (status /= 0) call TLS_panic ('ADVECT_VOLUME: Volume_Flux_Sub(nmat,nfc,ncells) allocation failed')
    end if

    FLUXING_PASSES: do p = 1,vps

       ! Evaluate volume changes using the interface tracker.
       if (volume_track_interfaces) then

          ! Initialize the array that'll keep track of volume fluxes in this subcycle.
          Volume_Flux_Sub = 0.0_r8

          ! Get the donor fluxes.
          call VOLUME_TRACK (Vof, Fluxing_Velocity, Volume_Flux_Sub)

          ! Normalize the donor fluxes.
          call FLUX_RENORM (Fluxing_Velocity, Vof_n, Volume_Flux_Tot, Volume_Flux_Sub)

          ! Compute the acceptor fluxes.
          call FLUX_ACCEPTOR (Volume_Flux_Sub)
   
          ! Compute BC (inflow) fluxes.
          call FLUX_BC (Fluxing_Velocity, Vof_n, Volume_Flux_Sub)

          ! Add the volume fluxes from this subcycle (Volume_Flux_Sub) to the
          ! total flux array (Volume_Flux_Tot), and update the volume fraction
          ! array (Vof).
          call VOLUME_ADVANCE (Volume_Flux_Sub, Volume_Flux_Tot, Vof)

       else

          ! Evaluate volume changes using a continuum advection scheme.
          call ADVECT_CONTINUUM (Fluxing_Velocity, Vof_n, Volume_Flux_Tot)

          ! Compute the acceptor fluxes.  Note that here we send in the total
          ! volume flux for the timestep, because there's no subcycling.
          call FLUX_ACCEPTOR (Volume_Flux_Tot)
   
          ! Compute BC (inflow) fluxes.  Note that here we send in the total
          ! volume flux for the timestep, because there's no subcycling.
          call FLUX_BC (Fluxing_Velocity, Vof_n, Volume_Flux_Tot)

          ! Update Vof.
          call CONTINUUM_ADVANCE (Volume_Flux_Tot, Vof)

       end if

       ! Make sure volume fractions of a particular material are within
       ! the allowed range (0 <= Vof <= 1) and that all materials sum to one.
       call VOF_BOUNDS (Vof, Volume_Flux_Tot)

    end do FLUXING_PASSES

    if (volume_track_interfaces) DEALLOCATE (Volume_Flux_Sub)

    ! Stop the volume advection timer.
    call stop_timer("Volume Tracking")

  END SUBROUTINE ADVECT_VOLUME

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><>

  SUBROUTINE ADVECT_CONTINUUM (Fluxing_Velocity, Vof_n, Volume_Flux_Tot)
    !=======================================================================
    ! Purpose(s):
    !
    !   Evaluate the volumetric material movement across every cell face
    !   using continuum advection methods (no geometric reconstruction)
    !
    !       Jim Sicilian (CCS-2)   October 2003
    !=======================================================================
    use fluid_data_module,      only: fluidVof, isImmobile
    use legacy_mesh_api,        only: ncells, nfc, Cell
    use parameter_module,       only: nmat
    use vof_data_module,        only: adv_dt

    ! Arguments
    real(r8), dimension(nfc,ncells),      intent(IN)    :: Fluxing_Velocity
    real(r8), dimension(nmat,ncells),     intent(IN)    :: Vof_n
    real(r8), dimension(nmat,nfc,ncells), intent(INOUT) :: Volume_Flux_Tot

    ! Local Variables
    integer :: n, f, m
    real(r8) :: Flux_Vol

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
   
    ! Calculate material volume fluxes out of cells, including through mesh
    ! boundary faces, in proportion to the fluid volumes in the donor cell.

    MATERIALS: do m = 1,nmat
       if (isImmobile(m)) cycle MATERIALS

       do n = 1,ncells
          if (fluidVof(n) > 0.0_r8) then
             do f = 1,nfc
                if (Fluxing_Velocity(f,n) > 0) then
                   Flux_Vol = adv_dt*Fluxing_Velocity(f,n)*Cell(n)%Face_Area(f)
                   Volume_Flux_Tot(m,f,n) = Flux_Vol*Vof_n(m,n)/fluidVof(n)
                end if
             end do
          end if
       end do

    end do MATERIALS

  END SUBROUTINE ADVECT_CONTINUUM

  SUBROUTINE CONTINUUM_ADVANCE (Volume_Flux_Tot, Vof)
    !=======================================================================
    ! Purpose(s):
    !
    !   Update the Vof array.
    !
    !=======================================================================
    use fluid_data_module,    only: isImmobile
    use legacy_mesh_api,      only: ncells, nfc, Cell
    use parameter_module,     only: nmat
    implicit none

    ! Arguments
    real(r8), dimension(nmat,nfc,ncells), intent(IN)    :: Volume_Flux_Tot
    real(r8), dimension(nmat,ncells),     intent(INOUT) :: Vof

    ! Local Variables
    integer :: f, m

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Loop over materials
    MATERIALS: do m = 1,nmat
       if(isImmobile(m)) cycle MATERIALS
       do f = 1,nfc
          Vof(m,:) = Vof(m,:) - Volume_Flux_Tot(m,f,:) / Cell%Volume
       end do
    end do MATERIALS

  END SUBROUTINE CONTINUUM_ADVANCE

  SUBROUTINE FLUX_ACCEPTOR (Volume_Flux_Sub)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute acceptor (negative) volume fluxes in this subcycle.
    !
    !=======================================================================
    use legacy_mesh_api,  only: ncells, nfc, EE_GATHER
    use parameter_module, only: nmat

    ! Arguments
    real(r8), dimension(nmat,nfc,ncells), intent(INOUT) :: Volume_Flux_Sub

    ! Local Variables
    integer :: m
    real(r8), dimension(nfc,ncells) :: acceptor_flux

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do m = 1,nmat
       call EE_GATHER (acceptor_Flux, Volume_Flux_Sub(m,:,:))
       where (acceptor_Flux > 0.0_r8) Volume_Flux_Sub(m,:,:) = - acceptor_Flux
    end do

  END SUBROUTINE FLUX_ACCEPTOR

  SUBROUTINE FLUX_BC (Fluxing_Velocity, Vof_n, Volume_Flux_Sub)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute inflow volume fluxes.
    !
    !=======================================================================
    use bc_module,         only: BC_Mat, IN_FLOW
    use input_utilities,   only: NULL_I
    use fluid_data_module, only: isImmobile
    use legacy_mesh_api,   only: ncells, nfc, Cell
    use parameter_module,  only: nmat
    use pgslib_module,     only: PGSLIB_GLOBAL_ANY
    use vof_data_module,   only: adv_dt

    ! Arguments
    real(r8), dimension(nfc,ncells),      intent(IN)    :: Fluxing_Velocity
    real(r8), dimension(nmat,ncells),     intent(IN)    :: Vof_n
    real(r8), dimension(nmat,nfc,ncells), intent(INOUT) :: Volume_Flux_Sub

    ! Local Variables
    real(r8) :: Sum_Vof_n
    real(r8), dimension(ncells) :: Flux_Vol
    integer :: f, n, m
    logical, dimension(ncells) :: Inflow_Mask

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    INFLOW_LOOP: do f = 1,nfc

       ! Find inflow faces.
       Inflow_Mask = IN_FLOW (f, Fluxing_Velocity)

       if (PGSLIB_GLOBAL_ANY(Inflow_Mask)) then

          ! Calculate the inflow flux volumes.
          where (Inflow_Mask) Flux_Vol = adv_dt*Fluxing_Velocity(f,:)*Cell%Face_Area(f)
   
          ! Zero out the subcycle volume fluxes at inflow faces.
          do m = 1,nmat
             where (Inflow_Mask) Volume_Flux_Sub(m,f,:) = 0.0_r8
          end do
   
          ! If inflow material specified as a BC, assign it.
          do n = 1,ncells
             if (Inflow_Mask(n) .and. BC_Mat(f,n) /= NULL_I) &
                Volume_Flux_Sub(BC_Mat(f,n),f,n) = Flux_Vol(n)
          end do
   
          ! If the user didn't specify an inflow material, assume that what's flowing in
          ! is more of what was in the cell at the beginning of the timestep.
          do n = 1,ncells
             if (Inflow_Mask(n) .and. BC_Mat(f,n) == NULL_I) then
                ! Sum the fluid Vof in the cell.
                Sum_Vof_n = 0.0_r8
                do m = 1,nmat
                   if (.not. isImmobile(m)) Sum_Vof_n = Sum_Vof_n + Vof_n(m,n)
                end do
                ! Set Volume_Flux_Sub in proportion to Vof_n/Sum_Vof_n.
                ! NNC, March 2013.  I think this is an error.  We only should be doing this for fluid materials.
                do m = 1,nmat
                   Volume_Flux_Sub(m,f,n) = Flux_Vol(n)*Vof_n(m,n)/Sum_Vof_n
                end do
             end if
          end do

       end if

    end do INFLOW_LOOP

  END SUBROUTINE FLUX_BC

  SUBROUTINE FLUX_RENORM (Fluxing_Velocity, Vof_n, Volume_Flux_Tot, Volume_Flux_Sub)
    !=======================================================================
    ! Purpose(s):
    !
    !   Scan all faces with an outward flux and determine if any
    !   material is over-exhausted from this cell.  If so lower
    !   the fluxes until the material is just exhausted.  Then
    !   loop over the faces and balance the individual material
    !   fluxes with the total face flux. The sum of the material
    !   volume fluxes (Volume_Flux_Sub) for each face should sum
    !   to the total volume flux for that face.
    !
    !=======================================================================
    use cutoffs_module,       only: cutvof
    use fluid_data_module,    only: isImmobile
    use legacy_mesh_api,      only: ncells, nfc, Cell
    use parameter_module,     only: nmat
    use vof_data_module,      only: adv_dt
 
    ! Arguments
    real(r8), dimension(nfc,ncells),      intent(IN)    :: Fluxing_Velocity
    real(r8), dimension(nmat,ncells),     intent(IN)    :: Vof_n
    real(r8), dimension(nmat,nfc,ncells), intent(INOUT) :: Volume_Flux_Tot
    real(r8), dimension(nmat,nfc,ncells), intent(INOUT) :: Volume_Flux_Sub

    ! Local Variables
    real(r8) :: Ratio, Sum, Cumul_Sum, Sum_not_maxed, Total_Face_Flux
    integer  :: norm_iter, f, m, n, number_not_maxed
    logical  :: Done_Renorm
    logical, dimension(nmat) :: Maxed

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    CELLS: do n = 1,ncells

       Maxed = .False.

       ! Loop over the renorm_loop a maximum of nmat - 1 times, to resolve all instances
       ! of fluxing more material than was in the cell at the beginning of the timestep

       RENORM_LOOP: do norm_iter = 1, nmat+1

          Done_Renorm = .True.

       ! We stay in this loop until the sum of material fluxes from each face of
       ! the cell equals the total face volume flux.  Where the cumulative sum of 
       ! individual material fluxes (from this and previous volume_track_subcycles) 
       ! exceeds the volume of a particular material originally within a 
       ! cell, we decrease those fluxes to equal the volume of material still available
       ! to be fluxed, and increase other fluxes appropriately.  If this increase
       ! leads to the over-exhaustion of other materials, we work our way through
       ! this loop again, and again, and again, and ... , until we're done.

       ! The first step is to determine if any material is being over-exhausted from 
       ! a cell.  If so mark it as MAXED and lower the Volume_Flux_Sub's so that the 
       ! material is just exhausted.

          MAT_LOOP: do m = 1,nmat

             ! Remember that at this stage of the code Face_Flux'es are only positive
             ! or zero.  Sum is the volume of material m attempting to leave the cell
             ! in this volume_track_subcycle; Cumul_Sum is Sum plus the material that
             ! has already left in previous volume_track_subcycles.
             Sum = 0.0_r8
             Cumul_Sum = 0.0_r8
             do f = 1,nfc
                Sum = Sum + Volume_Flux_Sub(m,f,n)
                Cumul_Sum = Cumul_Sum + MAX(Volume_Flux_Tot(m,f,n),0.0_r8)
             end do
             if (Sum == 0.0_r8) CYCLE MAT_LOOP
             Cumul_Sum = Cumul_Sum + Sum
 
             ! If the CUMULATIVE sum of outward fluxes across faces (Cumul_Sum)
             ! exceeds the amount of material ORIGINALLY in the cell (from Vof_n),
             ! calculate the 'Ratio' of fluid material volume still allowed to be
             ! fluxed to the flux volume, and note that we're not 'Done'
             Ratio = 0.0_r8  ! if none of this material was originally in the cell
             ! Update the Ratio for fluid materials; if the material isImmobile, 
             ! leave Ratio = 0.0_r8
             if (.not. isImmobile(m)) then
                Ratio = (Vof_n(m,n)*Cell(n)%Volume - (Cumul_Sum-Sum)) / Sum
             end if
             if (Ratio < 1.0_r8) then
                Done_Renorm = .False.
                Maxed(m) = .True.
             end if

             ! If Ratio < 1, lower the fluxes to match the material volume within
             ! the cell, and flag the cell and material number with 'Maxed'.
             if (Ratio < 1.0_r8) then
                do f = 1,nfc
                   Volume_Flux_Sub(m,f,n) = Ratio * Volume_Flux_Sub(m,f,n)
                end do
             end if
 
          end do MAT_LOOP

          if (Done_Renorm) exit RENORM_LOOP

       ! This cell had one/more fluxes reduced.  For each of the faces, if the sum
       ! of material fluxes is less than Total_Face_Flux, multiply all non-maxed 
       ! fluxes by another 'Ratio' (this time > 1) that restores the flux balance.  
       ! This may in turn over-exhaust one or more of these materials, and so from
       ! the bottom of this loop, we head back to the top.

          do f = 1,nfc

             ! Calculate the total flux volume through the cell face (is this already
             ! available elsewhere?), and if the flux volume is greater than zero, 
             ! then concern ourselves with adjusting individual material fluxes.

             Total_Face_Flux = adv_dt*Fluxing_Velocity(f,n)*Cell(n)%Face_Area(f)
             if (Total_Face_Flux > cutvof*Cell(n)%Volume) then
 
                ! Add up the sum of material fluxes at a face (Sum), and the sum of 
                ! un-maxed material fluxes (Sum_not_maxed).
                Sum = 0.0_r8
                Sum_not_maxed = 0.0_r8
                do m = 1,nmat
                   Sum = Sum + Volume_Flux_Sub(m,f,n)
                   if (.not. Maxed(m)) Sum_not_maxed = Sum_not_maxed + Volume_Flux_Sub(m,f,n)
                end do

                ! Ratio as defined below, when used to multiply the non-maxed fluxes at 
                ! a face, will restore the flux balance.
                if (Sum_not_maxed > 0.0_r8) then
                ! jms Note:  Ratio = (Total_Face_Flux - Maxed_Face_Flux) / Sum_not_maxed
                   Ratio = 1.0_r8 + (Total_Face_Flux - Sum) / Sum_not_maxed
                   do m = 1,nmat
                      if (.not. Maxed(m)) Volume_Flux_Sub(m,f,n) = Ratio * Volume_Flux_Sub(m,f,n)
                   end do
                else
                   number_not_maxed = 0
                   do m = 1,nmat
                      if (.not. Maxed(m) .and. .not.isImmobile(m)) number_not_maxed = number_not_maxed + 1
                   end do
                   if (number_not_maxed == 0) then
                       call TLS_panic ('FLUX_RENORM: cannot reassign face flux to any other material')
                   endif
                   Ratio = (Total_Face_Flux - Sum) / number_not_maxed
                   do m = 1,nmat
                      if (.not. Maxed(m).and. .not.isImmobile(m)) Volume_Flux_Sub(m,f,n) = Ratio
                   end do
                end if

             end if

          end do

       end do RENORM_LOOP

    end do CELLS

  END SUBROUTINE FLUX_RENORM
 

  SUBROUTINE VOLUME_ADVANCE (Volume_Flux_Sub, Volume_Flux_Tot, Vof)
    !=======================================================================
    ! Purpose(s):
    !
    !   Add the subcycle fluxes to the total, and update the Vof array.
    !
    !=======================================================================
    use fluid_data_module,    only: isImmobile
    use legacy_mesh_api,      only: ncells, nfc, Cell
    use parameter_module,     only: nmat

    ! Arguments
    real(r8), dimension(nmat,nfc,ncells), intent(IN)    :: Volume_Flux_Sub
    real(r8), dimension(nmat,nfc,ncells), intent(INOUT) :: Volume_Flux_Tot
    real(r8), dimension(nmat,ncells),     intent(INOUT) :: Vof

    ! Local Variables
    integer :: f, m

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Volume_Flux_Tot = Volume_Flux_Tot + Volume_Flux_Sub

    do m = 1,nmat
      if(isImmobile(m)) cycle
       do f = 1,nfc
          Vof(m,:) = Vof(m,:) - Volume_Flux_Sub(m,f,:) / Cell%Volume
       end do
    end do

  END SUBROUTINE VOLUME_ADVANCE

  SUBROUTINE VOF_BOUNDS (Vof, Volume_Flux_Tot)
    !=======================================================================
    ! Purpose(s):
    !
    !   Make sure volume fractions are within bounds:  0 <= Vof <= 1.
    !   If not, remove the overshoots (Vof > 1) and undershoots (Vof < 0).
    !
    !=======================================================================
    use cutoffs_module,     only: cutvof
    use fluid_data_module,  only: Void_Material_Exists, Void_Material_Index, &
                                  Void_Material_Count, isImmobile
    use legacy_mesh_api,    only: ncells, nfc
    use parameter_module,   only: nmat

    ! Arguments 
    real(r8), dimension(nmat,ncells),     intent(INOUT) :: Vof
    real(r8), dimension(nmat,nfc,ncells), intent(INOUT) :: Volume_Flux_Tot

    ! Local Variables
    integer  :: m, n, void_m, vmc
    real(r8) :: Ftot, Ftot_m1, void_volume
    real(r8), dimension(nmat) :: Delta_Vol
    logical  :: found

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    CELLS: do n = 1,ncells

       Ftot = 0.0_r8
       do m = 1,nmat

          if (.not.isImmobile(m)) then

             ! If volume fraction is > 1.0 - cutvof, round to one.
             if (Vof(m,n) > (1.0_r8-cutvof) .and. Vof(m,n) /= 1.0_r8) then
                found = .false.
                do vmc = 1,Void_Material_Count
                   if (m == Void_Material_Index(vmc)) then
                      Vof(m,n) = 1.0_r8
                      found = .true.
                      exit
                   end if
                end do
                if (.not. found) then
                   call ADJUST_FLUX_MATL (Volume_Flux_Tot, n, m, Vof(m,n), 1.0_r8)
                   Vof(m,n) = 1.0_r8
                end if
             end if

             ! If volume fraction is < cutvof; round to zero.
             if (Vof(m,n) < cutvof .and. Vof(m,n) /= 0.0_r8) then
                do vmc = 1,Void_Material_Count
                   if (m == Void_Material_Index(vmc)) then
                      Vof(m,n) = 0.0_r8
                      exit
                   end if
                end do
                if (Vof(m,n) /= 0.0_r8) then
                   call ADJUST_FLUX_MATL (Volume_Flux_Tot, n, m, Vof(m,n), 0.0_r8)
                   Vof(m,n) = 0.0_r8
                end if
             end if

          end if   ! end of the isImmobile test
   
          Ftot = Ftot + Vof(m,n)

       end do

       if (Ftot == 1.0_r8) cycle CELLS

       ! Renormalize the liquid volume fractions.
       if (Void_Material_Exists) then
          ! Check to see if void is already in the cell.
          void_m = 0
          void_volume = 0.0_r8
          do vmc = 1,Void_Material_Count
             if (Vof(Void_Material_Index(vmc),n) > 0.0_r8) then
                void_m = Void_Material_Index(vmc)
                void_volume = void_volume + Vof(void_m,n)
             end if
          end do
          if (void_m > 0) then
             if (Ftot > 1.0_r8) then
                ! Is there enough to balance the cell?
                if (void_volume > Ftot-1.0_r8) then
                   ! There is enough void ...
                   Ftot_m1 = Ftot - 1.0_r8
                   VMC_LOOP: do vmc = 1,Void_Material_Count
                      if (Vof(Void_Material_Index(vmc),n) > Ftot_m1) then
                         Vof(Void_Material_Index(vmc),n) = Vof(Void_Material_Index(vmc),n) - Ftot_m1
                         exit VMC_LOOP
                      else
                         Ftot_m1 = Ftot_m1 - Vof(Void_Material_Index(vmc),n)
                         Vof(Void_Material_Index(vmc),n) = 0.0_r8
                      end if
                   end do VMC_LOOP
                   cycle CELLS
                else
                   ! There isn't enough void ...
                   do vmc = 1,Void_Material_Count
                      void_m = Void_Material_Index(vmc)
                      Ftot = Ftot - Vof(void_m,n)
                      Vof(void_m,n) = 0.0_r8
                   end do
                   call ADJUST_FLUX_TOTAL (Volume_Flux_Tot, n, Ftot, Delta_Vol)
                   call ADJUST_VOFS (Vof, n, Delta_Vol)
                   cycle CELLS
                end if
             else
                ! Ftot < 1, and there's void already in the cel.
                Vof(void_m,n) = Vof(void_m,n) + 1.0_r8 - Ftot
                cycle CELLS
             end if
          end if
       end if

       ! Got to this line of code only if there's no void in this calculation,
       ! or at least there's no void in this cell.
       call ADJUST_FLUX_TOTAL (Volume_Flux_Tot, n, Ftot, Delta_Vol)
       call ADJUST_VOFS (Vof, n, Delta_Vol)
       cycle CELLS

    end do CELLS
 
  END SUBROUTINE VOF_BOUNDS

  SUBROUTINE ADJUST_VOFS (Vof, n, Delta_Vol)
    !======================================================================
    !  Purpose(s):
    !
    !     Adjust the Vof values in one cell to reflect the changes
    !     in volume transfer evaluated by ADJUST_FLUX_TOTAL
    !
    !    Jim Sicilian, CCS-2, December 2002
    !
    !======================================================================
    use legacy_mesh_api,  only: ncells, Cell
    use parameter_module, only: nmat

    ! Arguments
    real(r8), dimension(nmat,ncells), intent(INOUT) :: Vof
    integer, intent(in) :: n
    real(r8), dimension(nmat), intent(In) :: Delta_Vol
 
    ! Local Variables
    integer :: m

    do m = 1,nmat
       Vof(m,n) = Vof(m,n) + Delta_Vol(m)/Cell(n)%Volume
    end do

  END SUBROUTINE ADJUST_VOFS

  SUBROUTINE ADJUST_FLUX_MATL (Volume_Flux_Tot, n, MatID, Current_Material_Vof, Target_Material_Vof)
    !=======================================================================
    ! Purpose(s):
    !
    !  Adjust the material volume fluxes on the faces of a single cell to
    !  match the evaluated material volume to a target value.
    !
    !      Jim Sicilian,   CCS-2,   October 2002
    !
    !=======================================================================
    use legacy_mesh_api,        only: ncells, nfc, Cell
    use parameter_module,       only: nmat
    use projection_data_module, only: Boundary_Flag

    ! Arguments
    real(r8), dimension(nmat,nfc,ncells), intent(INOUT) :: Volume_Flux_Tot
    integer,  intent(IN) :: n, MatID
    real(r8), intent(IN) :: Current_Material_Vof
    real(r8), intent(IN) :: Target_Material_Vof

    ! Local Variables
    integer :: f
    real(r8) :: Inflow_Volume, Outflow_Volume, Volume_Change, Total_Flow, Change_Fraction
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Determine incoming and outgoing volumes changes of material MatID, for cell n.
    Inflow_Volume = 0.0_r8
    Outflow_Volume = 0.0_r8
    do f = 1, nfc
       ! Don't change dirichlet velocity BCs.
       if (Boundary_Flag(f,n)==2) cycle
       if (Volume_Flux_Tot(MatID,f,n) > 0.0_r8) then
          Outflow_Volume = Outflow_Volume + Volume_Flux_Tot(MatID,f,n)
       else
          Inflow_Volume = Inflow_Volume - Volume_Flux_Tot(MatID,f,n)
       end if
    end do

    ! Calculate the fractional change needed to adjust the current volume 
    ! to the target value.  The same fractional increase/decrease is applied
    ! to incoming and outgoing flows.
    Total_Flow = Inflow_Volume + Outflow_Volume
    Change_Fraction = Target_Material_Vof - Current_Material_Vof
    if(Total_Flow /= 0.0_r8) then
        Volume_Change = Change_Fraction*Cell(n)%Volume
        Change_Fraction = Volume_Change/Total_Flow
     else
        if(Change_Fraction < 0.0_r8) then
            ! jms Note:   If the material is to be removed from the cell
           ! look for a face that doesn't have incoming material, and 
           ! flux it out through that face
           do f = 1,nfc
              if (Boundary_Flag(f,n)==2 .or.        &
                  Volume_Flux_Tot(MatID,f,n) < 0.0_r8   ) cycle
                 Volume_Flux_Tot(MatID,f,n) = Volume_Flux_Tot(MatID,f,n) - Change_Fraction*Cell(n)%Volume
                 exit
           enddo
        else
           if (ABS(Change_Fraction) > WMaxMatU) then 
              WCountMatU = 0
              WMaxMatU = ABS(Change_Fraction)
              write(message,'(2(a,i0),a,es11.3)') 'unable to adjust material volume: cell=', n, &
                ', matid=', MatId, ', desired change fraction=', Change_Fraction
              call TLS_warn (message)
           elseif (WCountMatU < Wlimit) then
              WCountMatU = WCountMatU + 1
              write(message,'(2(a,i0),a,es11.3)') 'unable to adjust material volume: cell=', n, &
                ', matid=', MatId, ', desired change fraction=', Change_Fraction
              call TLS_warn (message)
           endif
        endif
        return
    endif

    ! Warning Message if the fractional change in volume is excessive
    if(ABS(Change_Fraction) > WMaxMat .and. Total_Flow > 1.0e-4*Cell(n)%Volume) then
       WCountMat = 0
       WMaxMat   = ABS(Change_Fraction)
       write(message,'(2(a,i0),a,es11.3)') 'excessive volume adjustment: cell=', n, &
             ', matid=', MatId, ', adjustment fraction=', Change_Fraction
       call TLS_warn (message)
    elseif(WCountMat < Wlimit .and. Total_Flow > 1.0e-4*Cell(n)%Volume) then
       WCountMat = WCountMat + 1
       write(message,'(2(a,i0),a,es11.3)') 'excessive volume adjustment: cell=', n, &
             ', matid=', MatId, ', adjustment fraction=', Change_Fraction
       call TLS_warn (message)
    endif

    ! Adjust the volume changes.
    do f = 1,nfc
       ! Don't change dirichlet velocity BCs.
       if (Boundary_Flag(f,n)==2) cycle
       if (Volume_Flux_Tot(MatID,f,n) > 0.0_r8) then
          Volume_Flux_Tot(MatID,f,n) = (1.0_r8-Change_Fraction)*Volume_Flux_Tot(MatID,f,n)
       else
          Volume_Flux_Tot(MatID,f,n) = (1.0_r8+Change_Fraction)*Volume_Flux_Tot(MatID,f,n)
       end if
    end do
    
  END SUBROUTINE ADJUST_FLUX_MATL

  SUBROUTINE ADJUST_FLUX_TOTAL (Volume_Flux_Tot, n, Current_Vof, Delta_Vol)
    !=======================================================================
    ! Purpose(s):
    !
    !  Adjust the material fluxes on the faces of a single cell to match
    !  the evaluated total material volume to a target value
    !
    !      Jim Sicilian,   CCS-2,   October 2002
    !
    !=======================================================================
    use fluid_data_module,      only: Void_material_Exists, Void_Material_Index, Void_Material_Count
    use legacy_mesh_api,        only: ncells, nfc, Cell
    use parameter_module,       only: nmat
    use projection_data_module, only: Boundary_Flag

    ! Arguments
    real(r8), dimension(nmat,nfc,ncells), intent(INOUT) :: Volume_Flux_Tot
    integer, intent(IN) :: n
    real(r8), intent(IN) :: Current_Vof
    real(r8), dimension(nmat), intent(OUT) :: Delta_Vol

    ! Local Variables
    integer :: f, m, v
    real(r8) :: Inflow_Volume, Outflow_Volume, Volume_Change, Total_Flow, Change_Fraction
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Determine incoming and outgoing volumes changes of material MatID, for cell n.
    Inflow_Volume = 0.0_r8
    Outflow_Volume = 0.0_r8
    Delta_Vol = 0.0_r8
    MAT_LOOP: do m = 1, nmat
       if(Void_material_Exists) then
          do v = 1, Void_Material_Count
             if(m == Void_Material_index(v)) cycle MAT_LOOP
          enddo
       endif
          do f = 1, nfc
            ! Don't change dirichlet velocity BCs.
            if (Boundary_Flag(f,n)==2 .or.    &
                Volume_Flux_Tot(m,f,n) == 0.0_r8  ) cycle
             if (Volume_Flux_Tot(m,f,n) > 0.0_r8) then
                Outflow_Volume = Outflow_Volume + Volume_Flux_Tot(m,f,n)
             else
                Inflow_Volume = Inflow_Volume - Volume_Flux_Tot(m,f,n)
             end if
          enddo
    end do MAT_LOOP

    ! Calculate the fractional change needed to adjust the current volume 
    ! to the target value.  The same fractional increase/decrease is applied
    ! to incoming and outgoing flows.
    Total_Flow = Inflow_Volume + Outflow_Volume
    Volume_Change = 1.0_r8 - Current_Vof
    if (Total_Flow /= 0.0_r8) then
       Volume_Change = Volume_Change*Cell(n)%Volume
       Change_Fraction = Volume_Change/Total_Flow
    else
       if (ABS(Volume_Change) > WMaxTotU) then 
          WCountTotU = 0
          WMaxTotU   = ABS(Volume_Change)
       elseif (WCountTotU < Wlimit) then
          WCountTotU = WCountTotU + 1
       endif
       write(message,'(a,i0,a,es11.3)') 'unable to adjust total volume: cell=', n, &
                                        ', desired change fraction=', Volume_Change
       call TLS_warn (message)
       return
    end if

    ! Warning Message if the fractional change in volume is excessive
    if(ABS(Change_Fraction) > WMaxTot .and. Total_Flow > 1.0e-4*Cell(n)%Volume) then
       WCountTot = 0
       WMaxTot = ABS(Change_Fraction)
       write(message,'(a,i0,a,es11.3)') 'excessive volume adjustment: cell=', n, &
                                        ', adjustment fraction=', Change_Fraction
       call TLS_warn (message)
    elseif (WCountTot < Wlimit .and. Total_Flow > 1.0e-4*Cell(n)%Volume) then
       WCountTot = WCountTot + 1
       write(message,'(a,i0,a,es11.3)') 'excessive volume adjustment: cell=', n, &
                                        ', adjustment fraction=', Change_Fraction
       call TLS_warn (message)
    endif

    ! Adjust the Volume Changes by face and accumulate them by cell.
    MAT_LOOP2: do m = 1, nmat
       if(Void_material_Exists) then
          do v = 1, Void_Material_Count
             if(m == Void_Material_index(v)) cycle MAT_LOOP2
          enddo
       endif
          do f = 1, nfc
             ! Don't change dirichlet velocity BCs.
             if (Boundary_Flag(f,n)==2) cycle
                if (Volume_Flux_Tot(m,f,n) == 0.0_r8) cycle
                if (Volume_Flux_Tot(m,f,n) > 0.0_r8) then
                   ! Adjust outgoing material transfers.
                   Delta_Vol(m) = Delta_Vol(m) + Change_Fraction*Volume_Flux_Tot(m,f,n)
                   Volume_Flux_Tot(m,f,n) = (1.0_r8-Change_Fraction)*Volume_Flux_Tot(m,f,n)
                else
                   ! Adjust incoming material transfers.
                   Delta_Vol(m) = Delta_Vol(m) - Change_Fraction*Volume_Flux_Tot(m,f,n)
                   Volume_Flux_Tot(m,f,n) = (1.0_r8+Change_Fraction)*Volume_Flux_Tot(m,f,n)
                end if
          end do
    end do MAT_LOOP2
    
  END SUBROUTINE ADJUST_FLUX_TOTAL

END MODULE ADVECT_VOLUME_MODULE
