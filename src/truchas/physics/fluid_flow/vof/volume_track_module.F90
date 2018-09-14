!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE VOLUME_TRACK_MODULE
  !=======================================================================
  ! Purpose:
  !
  !   Define the main volume-tracking procedures.
  !
  !   Public Interface:
  !
  !     * call VOLUME_TRACK (Vof, Fluxing_Velocity, Volume_Flux_Sub)
  !
  !       Control routine for a piecewise-linear volume tracking of
  !       material interfaces, in which material interfaces are
  !       reconstructed as planes from local volume fraction data.
  !
  ! Contains: VOLUME_TRACK
  !           INT_NORMAL
  !
  ! Author(s): Stewart J. Mosso (sjm@lanl.gov)
  !            Douglas B. Kothe (dbk@lanl.gov)
  !            Matthew Williams (mww@lanl.gov)
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private
 
  public :: VOLUME_TRACK

CONTAINS
 
  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE VOLUME_TRACK (Vof, Fluxing_Velocity, Volume_Flux_Sub)
    !=======================================================================
    ! Purpose(s):
    !
    !   Control routine for a volume-tracking interface treatment,
    !   in which material interfaces are reconstructed as planes
    !   from local volume fraction data.  The plane interface
    !   locations are then used to find material advection volumes
    !   at all relevant cell faces.
    !
    !=======================================================================
    use cutoffs_module,            only: cutvof, alittle
    use fluid_data_module,         only: FluidRho
    use flux_volume_module,        only: FLUX_VOL_QUANTITY, FLUX_VOL_VERTICES
    use interface_module,          only: Int_Geom, Int_Flux, SETUP_INTERFACE_GEOM, &
                                         SETUP_INTERFACE_FLUX
    use interface_output_module,   only: INTERFACE_OUTPUT, time_for_int_dump
    use interface_triangle_module, only: INTERFACE_TRIANGLES
    use locate_plane_module,       only: LOCATE_PLANE
    use legacy_mesh_api,           only: ncells, ndim, nfc, nvc, Cell, Vertex, EN_GATHER
    use legacy_matl_api,           only: nmat
    use parameter_module,          only: nicells
    use pgslib_module,             only: PGSLib_Global_MAXVAL
    use property_data_module,      only: Matpri
    use truncate_volume_module,    only: Trunc_Vol, TRUNCATE_VOLUME, FACE_PARAM
    use vof_data_module,           only: adv_dt, interface_area
    use truchas_timers
 
    ! Arguments
    real(r8), dimension(nmat,ncells),     intent(IN)    :: Vof
    real(r8), dimension(nfc,ncells),      intent(IN)    :: Fluxing_Velocity
    real(r8), dimension(nmat,nfc,ncells), intent(INOUT) :: Volume_Flux_Sub
 
    ! Local Variables
    integer :: status
    integer :: interfaces, materials, m, n, f, ff, ni, v, p
    integer,  dimension(ncells) :: Mat, last
    integer,  dimension(:,:), allocatable  :: Pri_Ptr
    real(r8), dimension(ncells) :: G1, G2, G3, Vp, Vofm, Vofint, MVF
    real(r8), dimension(nfc,ncells) :: Flux_Vol_Sum
    real(r8), dimension(ndim,nvc,ncells) :: Xv
    real(r8), dimension(nmat,ndim,ncells) :: Grad
    real(r8), dimension(:,:,:), allocatable :: Pack_vert
    logical,  dimension(ncells) :: Mask
    type(FLUX_VOL_QUANTITY), dimension(ncells) :: Flux_Vol

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the volume track timer
    call start_timer ("Reconstruct/Advect")

    ! Count the number of materials in each cell
    Mat = 0
    do m = 1,nmat
       where (Vof(m,:) > 0.0_r8) Mat = Mat + 1
    end do
 
    ! Number of interfaces to process, and the maximum number of materials in a cell
    materials = PGSLIB_Global_MAXVAL(Mat)
    interfaces = materials - 1

    ! Materials are advected in either ascending or descending priority.
    ! We choose ascending.  Form a pointer Pri_Ptr(materials,ncells) that
    ! contains a list of each cells' materials in ascending order.  This
    ! code presumes that for nmat materials, that priorities have been
    ! assigned between 1 and nmat.

    ALLOCATE (Pri_Ptr(materials,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('VOLUME_TRACK: Pri_Ptr(materials,ncells) allocation failed')
    Pri_Ptr = 0

    ! Use the last(ncells) array to keep track of the first index of Pri_Ptr
    last = 1

    ! Fill Pri_Ptr for each cell with a list of materials in priority order.  Recall
    ! that Matpri(1) contains the material of priority 1, Matpri(2) the material of 
    ! priority 2, ...
    do p = 1,nmat
       m = Matpri(p)
       do n = 1,ncells
          if (Vof(m,n) > 0.0_r8) then
             Pri_Ptr(last(n),n) = m
             last(n) = last(n) + 1
          end if
       end do
    end do

    ! Compute interface normal vectors for all the materials.
    call INT_NORMAL (Vof, Grad, MVF)

    ! Initialize the sum of volume fluxes on each donor face.
    Flux_Vol_Sum = 0.0_r8
 
    ! Loop over the interfaces in priority order; if interfaces <= 0,
    ! we don't ever go into this loop.
    INTERFACE_LOOP: do ni = 1, interfaces

       ! Accumulate the volume fraction of this material and
       ! the volume fraction of materials with lower priorities.
       Vofint = 0.0_r8
       do m = 1,ni
          do n = 1,ncells
             if (Pri_Ptr(m,n) > 0) Vofint(n) = Vofint(n) + Vof(Pri_Ptr(m,n),n)
          end do
       end do
 
       ! Force 0.0 <= Vofint <= 1.0
       Vofint = MIN(MAX(Vofint, 0.0_r8), 1.0_r8)
 
       ! Retrieve this material's interface normal and volume fraction from the Grad array.
       G1   = 0.0_r8
       G2   = 0.0_r8
       G3   = 0.0_r8
       Vofm = 0.0_r8
       do n = 1,ncells
          if (Pri_Ptr(ni,n) > 0) then
             G1(n) = Grad(Pri_Ptr(ni,n),1,n)
             G2(n) = Grad(Pri_Ptr(ni,n),2,n)
             G3(n) = Grad(Pri_Ptr(ni,n),3,n)
             Vofm(n) = Vof(Pri_Ptr(ni,n),n)
          end if
       end do
 
       ! Count the number of mixed donor cells
       Mask = Vofint > (0.0_r8 + cutvof) .and. Vofint < (1.0_r8 - cutvof)
       where(fluidRho < alittle) Mask = .false.
       nicells = COUNT(Mask)

       ! Allocate the Interface derived type
       ALLOCATE (Int_Geom(nicells), STAT = status)
       if (status /= 0) call TLS_panic ('VOLUME_TRACK: Int_Geom(nicells) allocation failed')
       ALLOCATE (Int_Flux(nicells), STAT = status)
       if (status /= 0) call TLS_panic ('VOLUME_TRACK: Int_Flux(nicells) allocation failed')
       ALLOCATE (Trunc_Vol(nicells,nfc), STAT = status)
       if (status /= 0) call TLS_panic ('VOLUME_TRACK: Trunc_Vol(nicells,nfc) allocation failed')

       if (interface_area) then

          ALLOCATE (Pack_vert(ndim,nvc,nicells), STAT = status)
          if (status /= 0) call TLS_panic ('VOLUME_TRACK: Pack_vert(ndim,nvc,nicells) allocation failed')

          ! Gather vertex coordinates into Xv then pack into Pack_vert
          do m = 1,ndim
             call EN_GATHER (Xv(m,:,:),Vertex%Coord(m))
          end do

          do m = 1,ndim
             do v = 1,nvc
                Pack_vert(m,v,:) = PACK(Xv(m,v,:),Mask)
             end do
          end do

       end if

       ! Pack interface data into the Interface derived type
       call SETUP_INTERFACE_GEOM (Mask, G1, G2, G3, Vofint)

       ! Locate each interface plane by computing rho, the plane constant.
       call LOCATE_PLANE ()

       ! Write out interface data if necessary.
       !-mf-jim Oct 12 2005 modified 
       if (time_for_int_dump) then
         call INTERFACE_OUTPUT ()
         if (n == interfaces) then
           time_for_int_dump=.false.
         endif
       endif

       ! Find Interface Areas
       if (interface_area) call INTERFACE_TRIANGLES (Pack_vert)
 
       ! Calculate delta advection volumes (Volume_Flux_Sub) for this material
       ! at each donor face and accumulate the sum.
       FACE_INTERFACE_LOOP: do f = 1,nfc

          ! Flux volumes
          Flux_Vol%Fc  = f
          Flux_Vol%Vol = adv_dt*Fluxing_Velocity(f,:)*Cell%Face_Area(f)
          where (Flux_Vol%Vol <= cutvof*Cell%Volume) 
             Flux_Vol%Fc = 0
             Flux_Vol%Vol = 0.0_r8
          endwhere

          call FLUX_VOL_VERTICES (f, Mask, Fluxing_Velocity, Flux_Vol) 

          call SETUP_INTERFACE_FLUX (Flux_Vol, Mask)

          ! Now compute the volume truncated by interface planes in each flux volumes.
          do ff = 1,nfc
             call FACE_PARAM ('flux_cell', ff)
          end do
          call TRUNCATE_VOLUME(Int_Flux%Advected_Volume)

          ! Use G1 as a temporary to hold the delta advection volumes.  Where
          ! the donor cell doesn't have this material, zero the advected volume
          ! for this material.
          G1 = 0.0_r8
 
          ! For clean donor cells, the entire Flux volume goes to the single donor material.
          ! The following line was replaced to deal with the possibility that the volume fraction
          ! of this material was increased to within cutvof of 1 in an earlier pass, converting the
          ! cell from an interface cell to a pure cell.  This is necessary because VOF_BOUNDS does
          ! not increase such volume fractions to 1.   Jim Sicilian   July 2002.

!!$       where ((.not.Mask) .and. Vofm == 1.0_r8 ) G1 = ABS(Flux_Vol%Vol)
          where ((.not.Mask) .and. Vofm >= (1.0_r8-cutvof) ) G1 = ABS(Flux_Vol%Vol)
 
          ! For clean donor cells, the face flux has been put in G1.
          ! For mixed donor cells, the face flux is in Int_Flux%Advection_Volume.
          Vp = 0.0_r8
          Vp = UNPACK(Int_Flux%Advected_Volume, Mask, G1)
 
          ! If Vp is close to 0 set it to 0.  If it is close
          ! to 1 set it to 1. This will avoid numerical round-off.
          where (Vp > (1.0_r8-cutvof)*ABS(Flux_Vol%Vol)) Vp = ABS(Flux_Vol%Vol)
 
          ! Make sure that the current material-integrated advection
          ! volume hasn't decreased from its previous value.  (This
          ! can happen if the interface significantly changed its
          ! orientation and now crosses previous interfaces.)  Also
          ! limit Volume_Flux_Sub to take no more than the material
          ! occupied volume in the donor cell.
          G2 = 0.0_r8
          where (ABS(Flux_Vol%Vol) > 0.0_r8) G2 = MIN(MAX(Vp - Flux_Vol_Sum(f,:), 0.0_r8), Vofm*Cell%Volume)
 
          ! Now gather advected volume information into Volume_Flux_Sub
          do n = 1,ncells
             if (Pri_Ptr(ni,n) > 0 .and. Flux_Vol(n)%Vol > 0.0_r8) then
                Volume_Flux_Sub(Pri_Ptr(ni,n),f,n) = G2(n)
                Flux_Vol_Sum(f,n) = Flux_Vol_Sum(f,n) + G2(n)
             end if
          end do

       end do FACE_INTERFACE_LOOP

       ! Deallocate the interface data
       DEALLOCATE (Int_Geom)
       DEALLOCATE (Int_Flux)
       DEALLOCATE (Trunc_Vol)
 
       ! Deallocate interface cell vertices if interface areas have been computed
       if (ALLOCATED(Pack_vert)) DEALLOCATE (Pack_vert)

    end do INTERFACE_LOOP

    ! If I'm right, then the INTERFACE_LOOP, which looped over the number of interfaces
    ! in each cell, did so one less times than the number of materials in each cell.  And
    ! so this last loop looks at the lowest-priority material, and so-to-speak, reconciles
    ! the books.  Yes?   - Markus

    ! Make sure the volume flux sum has been initialized if there are no interfaces.
    if (interfaces == 0) Flux_Vol_Sum = 0.0_r8

    ! Compute the advection volume for the last material.
    last = Mat
 
    ! Fetch the volume fraction of the donor cell's last material.
    do n = 1,ncells
       Vofm(n) = Vof(Pri_Ptr(last(n),n),n)
    end do

    FACE_LAST_LOOP: do f = 1,nfc

       ! Recalculate the total flux volume for this face.
       Flux_Vol%Vol = adv_dt*Fluxing_Velocity(f,:)*Cell%Face_Area(f)
       where (Flux_Vol%Vol <= cutvof*Cell%Volume) Flux_Vol%Vol = 0.0_r8

       ! The volume flux of the last material shouldn't be less than
       ! zero nor greater than the volume of this material in the donor cell.
       G2   = 0.0_r8
       where (ABS(Flux_Vol%Vol) > 0.0_r8) &
          G2 = MIN(MAX(ABS(Flux_Vol%Vol - Flux_Vol_Sum(f,:)), 0.0_r8), Vofm*Cell%Volume)
 
       ! Store the last material's volume flux.
       do n = 1,ncells
          if (G2(n) > cutvof*Cell(n)%Volume) then
             Volume_Flux_Sub(Pri_Ptr(last(n),n),f,n) = G2(n)
             Flux_Vol_Sum(f,n) = Flux_Vol_Sum(f,n) + G2(n)
          end if
       end do

       ! For donor cells containing only one material, assign the total flux.
       do n = 1,ncells
          if (Mat(n) == 1 .and. Flux_Vol(n)%Vol > 0.0_r8) then
             Volume_Flux_Sub(Pri_Ptr(last(n),n),f,n) = Flux_Vol(n)%Vol
             Flux_Vol_Sum(f,n) = Flux_Vol(n)%Vol
          end if
       end do

    end do FACE_LAST_LOOP

    DEALLOCATE (Pri_Ptr)

    ! Stop the volume track timer
    call stop_timer ("Reconstruct/Advect")
 
  END SUBROUTINE VOLUME_TRACK
 
  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE INT_NORMAL (Vof, Grad, MVF)
    !=======================================================================
    ! Purpose(s):
    !
    !   This routine calculates the interface normal for each
    !   material by computing the gradient of material volume fractions
    !
    !=======================================================================
    use cutoffs_module,       only: alittle
    use discrete_op_module,   only: GRADIENT
    use interface_module,     only: interface_topology_model
    use mollify,              only: MOLLIFY_CONV_SAVEMEM
    use legacy_matl_api,      only: nmat
    use legacy_mesh_api,      only: ncells, ndim
    use property_data_module, only: Matpri
    use vof_data_module,      only: interface_geometry
 
    ! Arguments
    real(r8), dimension(nmat,ncells),      intent(IN)  :: Vof
    real(r8), dimension(nmat,ndim,ncells), intent(OUT) :: Grad
    real(r8), dimension(ncells),           intent(OUT) :: MVF
 
    ! Local Variables
    integer :: status
    integer :: m, n, mp, mi, nc
    real(r8), dimension(ndim,nmat,ncells) :: Mat_Normal
    real(r8), dimension(ncells) :: Tmp
    real(r8), dimension(:,:), allocatable :: Nx, Ny, Nz
   
    integer, dimension(ncells) :: Mat    ! to count number of material present in a cell
    integer, dimension(nmat,ncells) :: mid    ! to identify material id in a given cell
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
    ! Initialize arrays
    Mat_Normal = 0.0_r8
    Grad  = 0.0_r8

    ! Use the Nx,Ny,Nz arrays for the convolution method.
    if (interface_topology_model == 'convolution model') then
       ALLOCATE(Nx(nmat,ncells),Ny(nmat,ncells),Nz(nmat,ncells), STAT=status)
       call TLS_fatal_if_any ((status /= 0), 'INT_NORMAL: error allocating Nx,Ny,Nz')
       Nx = 0.0_r8
       Ny = 0.0_r8
       Nz = 0.0_r8
    end if 

    ! Loop over materials
    do m = 1, nmat

       ! I also think that we
       ! could do much of this without Mat_Normal and without Nx, Ny, and Nz.  - Markus

       ! Calculate the volume fraction gradient for this material

       TOPOLOGY: if (interface_topology_model == 'least squares model') then

          call GRADIENT (Mat_Normal(1,m,:), Mat_Normal(2,m,:), Mat_Normal(3,m,:), &
               Vof(m,:), method = 'Green-Gauss')
 
          ! The interface normal is in the opposite sense of the gradient
          do n = 1,ndim
             Mat_Normal(n,m,:) = - Mat_Normal(n,m,:)
          end do

       else if (interface_topology_model == 'convolution model') then

          if (nmat == 2 .and. m == 2) then
             do n = 1,ndim
                Mat_Normal(n,m,:) = - Mat_Normal(n,1,:)
             end do
          else
             ! Why does this call need m if it is on a per material basis?
             call MOLLIFY_CONV_SAVEMEM (Tmp,MVF,Nx,Ny,Nz,m)

             Mat_Normal(1,m,:) = Nx(m,:)
             Mat_Normal(2,m,:) = Ny(m,:)
             Mat_Normal(3,m,:) = Nz(m,:)
          end if

       end if TOPOLOGY

       ! Kill all but one gradient component if reconstruction is piecewise constant

       if (TRIM(interface_geometry) == 'piecewise constant') then
          where (ABS(Mat_Normal(1,m,:)) >= ABS(Mat_Normal(2,m,:)) .and. &
                 ABS(Mat_Normal(1,m,:)) >= ABS(Mat_Normal(3,m,:)))
             Mat_Normal(2,m,:) = 0.0_r8 
             Mat_Normal(3,m,:) = 0.0_r8
          end where
          where (ABS(Mat_Normal(2,m,:)) >= ABS(Mat_Normal(1,m,:)) .and. &
                 ABS(Mat_Normal(2,m,:)) >= ABS(Mat_Normal(3,m,:)))
             Mat_Normal(1,m,:) = 0.0_r8 
             Mat_Normal(3,m,:) = 0.0_r8
          end where
          where (ABS(Mat_Normal(3,m,:)) >= ABS(Mat_Normal(1,m,:)) .and. &
                 ABS(Mat_Normal(3,m,:)) >= ABS(Mat_Normal(2,m,:)))
             Mat_Normal(1,m,:) = 0.0_r8 
             Mat_Normal(2,m,:) = 0.0_r8
          end where
       end if

    end do

    ! mmfran 07/22/11 --- begin changes
    ! bug fix for material order independent results
    ! consider special case the cell contains only two materials
    ! if cell has two materials only, their normal should be consistent
    Mat = 0
    mid = 0
    do m=1,nmat
      where(Vof(m,:) > 0.0_r8) Mat = Mat + 1
    enddo

    do nc = 1,ncells
      
      ! if there are more than 2 materials in the cell
      ! Sum the gradients in priority order.  This is equivalent to calculating the
      !  interface normal for a material composed of the first few materials
      if (Mat(nc) >= 3) then
        mp = Matpri(1)
        do mi = 2, nmat
           m = Matpri(mi)
           do n = 1,ndim
             if (Vof(m,nc) >  alittle)  Mat_Normal(n,m,nc) = Mat_Normal(n,m,nc) + Mat_Normal(n,mp,nc)
             if (Vof(m,nc) <  alittle)  Mat_Normal(n,m,nc) = Mat_Normal(n,mp,nc)
           end do
        end do
      endif

    ! if there are only two materials in the cell
      if (Mat(nc) == 2) then
        ! identify the material ID in the cells
        mi=1
        do m=1,nmat
           if (Vof(m,nc) > 0.0_r8) then
             mid(mi,nc)=m
             mi=mi+1
           endif
        end do
        Mat_Normal(:,mid(2,nc),nc) = - Mat_Normal(:,mid(1,nc),nc)
      endif

    end do !enddo nc loop

    ! note: some improvement could be considered by checking 
    ! which material has the highest priority
    ! mmfran 07/22/11 ---- end of changes

    ! Eliminate noise from the gradient
    do n = 1,ndim
       where (ABS(Mat_Normal(n,:,:)) < alittle) Mat_Normal(n,:,:) = 0.0_r8
    end do
 
    ! Save Mat_Normal into Grad
    do m = 1,nmat
 
       Grad(m,1,:) = Mat_Normal(1,m,:)
       Grad(m,2,:) = Mat_Normal(2,m,:)
       Grad(m,3,:) = Mat_Normal(3,m,:)
 
       ! Normalize the gradient
       Tmp = SQRT(Grad(m,1,:)**2 + Grad(m,2,:)**2 + Grad(m,3,:)**2)

       where (ABS(Tmp) > alittle)
          Tmp = 1.0_r8/Tmp
       elsewhere
          Tmp = 1.0_r8
       end where

       Grad(m,1,:) = Grad(m,1,:) * Tmp
       Grad(m,2,:) = Grad(m,2,:) * Tmp
       Grad(m,3,:) = Grad(m,3,:) * Tmp
 
       ! Set tiny components to zero, and set all components to 1.0 in pure cells.
       where (ABS(Grad(m,1,:)) < 1.0d-6) Grad(m,1,:) = 0.0_r8
       where (ABS(Grad(m,2,:)) < 1.0d-6) Grad(m,2,:) = 0.0_r8
       where (ABS(Grad(m,3,:)) < 1.0d-6) Grad(m,3,:) = 0.0_r8
 
       where (Grad(m,1,:) == 0.0_r8 .and. Grad(m,2,:) == 0.0_r8 .and. Grad(m,3,:) == 0.0_r8)
          Grad(m,1,:) = 1.0_r8
          Grad(m,2,:) = 1.0_r8
          Grad(m,3,:) = 1.0_r8
       end where
 
    end do

    ! Cleanup arrays allocated for the convolution method.
    if (ALLOCATED(Nx)) DEALLOCATE(Nx,Ny,Nz)
 
  END SUBROUTINE INT_NORMAL
 
END MODULE VOLUME_TRACK_MODULE
