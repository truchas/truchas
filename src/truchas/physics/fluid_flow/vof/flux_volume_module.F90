!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE FLUX_VOLUME_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define the variables and procedures associated with 
  !   computing the advection flux volume.
  !
  ! Public Interface(s):
  !
  !   * call FLUX_VOL_VERTICES (Fluxing_Velocity, Flux_Vol)
  !
  !     Compute the vertices that describe the flux volume.
  !
  ! Contains: FLUX_VOL_VERTICES
  !           SCREWED_VOLUME
  !
  ! Author(s): Douglas B. Kothe (LANL Group T-3, dbk@lanl.gov)
  !            S. Jay Mosso (LANL Group X-HM, sjm@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use legacy_mesh_api, only: ndim, nvc
  use truchas_logging_services
  implicit none
  private

  public :: FLUX_VOL_QUANTITY, FLUX_VOL_VERTICES

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Maximum allowable iterations for the flux volume vertex location
  integer, public, save :: flux_vol_iter_max

  ! Define FLUX structure
  type FLUX_VOL_QUANTITY
  
     ! Face number through which we this cell is fluxing
     integer :: Fc

     ! Volume of the flux
     real(r8) :: Vol

     ! Vertices of the flux volume
     real(r8), dimension(nvc,ndim) :: Xv
     
  end type FLUX_VOL_QUANTITY

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE FLUX_VOL_VERTICES (face, Mask, Fluxing_Velocity, Flux_Vol)
    !=======================================================================
    ! Purpose(s):
    !
    !   Given the value of Flux_Vol (the volume of material that moves
    !   through the current advection cell face), what are the vertices
    !   that describe this volume.  Four of these vertices will be the ones
    !   that describe the advection cell face.  The other four vertices will
    !   lie approximately "DIST" away from the advection cell face.  These
    !   are only approximately "DIST" away because the cell cross-sectional
    !   area may increase or decrease as one moves away from the advection
    !   cell face.  The value used is varied from "DIST" such that the 
    !   vertices describe a hexagonal volume that matches the value of
    !   Flux_Vol.
    !
    !=======================================================================
    use cutoffs_module,   only: alittle, cutvof
    use legacy_mesh_api,  only: ncells, nfc, nvf, Cell, Mesh, gather_vertex_coord
    use legacy_mesh_api,  only: CELL_TET, CELL_PYRAMID, CELL_PRISM, CELL_HEX
    use vof_data_module,  only: adv_dt

    ! Arguments
    integer, intent(IN) :: face
    logical,  dimension(ncells),     intent(IN) :: Mask
    real(r8), dimension(nfc,ncells), intent(IN) :: Fluxing_Velocity
    type(FLUX_VOL_QUANTITY), dimension(ncells), intent(INOUT) :: Flux_Vol

    ! Local Variables
    integer :: i, n, v, ia, ib, e, iter
    real(r8), dimension(nvf)         :: Percnt
    real(r8), dimension(nvc,ncells)  :: Vtx
    real(r8), dimension(nvf,ndim)    :: Uedge
    real(r8) :: Volume, Mult, Tmp, Dist
    integer, dimension(2,nvf,nfc) :: Edge_ends, &
        HEX_Edge_ends     = RESHAPE((/ 3,2,  4,1,  7,6,  8,5, &     ! face one edges
                                       1,4,  2,3,  5,8,  6,7, &     ! face two edges
                                       1,2,  4,3,  5,6,  8,7, &     ! face three edges
                                       2,1,  3,4,  6,5,  7,8, &     ! face four edges
                                       1,5,  2,6,  3,7,  4,8, &     ! face five edges
                                       5,1,  6,2,  7,3,  8,4 /), &  ! face six edges
                                       (/2,nvf,nfc/)), &
        PRISM_Edge_ends   = RESHAPE((/ 3,2,  4,1,  7,2,  8,1, &
                                       1,4,  2,3,  5,4,  6,3, &
                                       1,2,  4,3,  5,6,  8,7, &
                                       2,1,  3,4,  6,5,  7,8, &
                                       1,5,  2,6,  3,7,  4,8, &
                                       5,1,  6,2,  7,3,  8,4 /), &
                                       (/2,nvf,nfc/)), &
        PYRAMID_Edge_ends = RESHAPE((/ 3,2,  4,1,  7,2,  8,1, &
                                       1,4,  2,3,  5,4,  6,3, &
                                       1,2,  4,3,  5,2,  8,3, &
                                       2,1,  3,4,  6,1,  7,4, &
                                       1,5,  2,6,  3,7,  4,8, &
                                       5,1,  6,2,  7,3,  8,4 /), &
                                       (/2,nvf,nfc/)), &
        TET_Edge_ends     = RESHAPE((/ 3,2,  4,1,  7,2,  8,1, &
                                       1,4,  2,3,  5,4,  6,3, &
                                       1,3,  4,3,  5,3,  8,3, &
                                       2,4,  3,4,  6,4,  7,4, &
                                       1,5,  2,6,  3,7,  4,8, &
                                       5,1,  6,2,  7,3,  8,4 /), &
                                       (/2,nvf,nfc/))
    logical :: converged
    character(256) :: errmsg

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! First gather the vertex coordinates and store them in Vtx
    do n = 1,ndim
       call gather_vertex_coord (Vtx, dim=n)
       do v = 1,nvc
          Flux_Vol%Xv(v,n) = Vtx(v,:)
       end do
    end do

    do i = 1,ncells

       ! Need only go through this exercise for Mask = .true., and for outgoing fluxes
       if (.not.Mask(i)) cycle
       if (Fluxing_Velocity(face,i)*Cell(i)%Face_Area(face)*adv_dt &
              <= cutvof*Cell(i)%Volume) cycle

       ! An initial guess for the depth of the flux volume is 'Dist'
       Dist = Fluxing_Velocity(face,i) * adv_dt

       ! Define start and end vertices for the 'nvf' edges going back from 'face'
       CELL_TYPE:  select case (Mesh(i)%Cell_shape)
          case (CELL_TET)
             Edge_ends = TET_Edge_ends
          case (CELL_PYRAMID)
             Edge_ends = PYRAMID_Edge_ends
          case (CELL_PRISM)
             Edge_ends = PRISM_Edge_ends
          case (CELL_HEX)
             Edge_ends = HEX_Edge_ends
       end select CELL_TYPE

       ! Compute the edge unit vectors
       do e = 1, nvf
          ia  = Edge_ends(1,e,face)
          ib  = Edge_ends(2,e,face)

          Tmp = 0.0_r8
          do n = 1,ndim
             Uedge(e,n) = Flux_Vol(i)%Xv(ib,n) - Flux_Vol(i)%Xv(ia,n)
             Tmp = Tmp + Cell(i)%Face_Normal(n,face) * Uedge(e,n)
          end do
          Percnt(e) = -Dist/(Tmp+alittle)
       end do

       do e = 1, nvf
          if (Percnt(e) < 0.0_r8 .or. Percnt(e) > 1.0_r8) then
             write(errmsg,'(a,3es13.6)') 'FLUX_VOL_VERTICES: invalid flux volume or inverted element; cell centroid =', Cell(i)%Centroid
             call TLS_panic (errmsg)
          end if
       end do

       ! Initialize the scale factor to unity
       Mult = 1.0_r8

       converged = .false.

       ! Iterate to find the 4 other vertices that define the back end of the flux volume
       do iter = 1,flux_vol_iter_max

          ! Loop over edges to adjust the vertices
          do e = 1,nvf

             ia = Edge_ends(1,e,face)
             ib = Edge_ends(2,e,face)

             do n = 1,ndim
                Flux_Vol(i)%Xv(ib,n) = Flux_Vol(i)%Xv(ia,n) + Mult*Percnt(e)*Uedge(e,n)
             end do

          end do

          ! Compute the flux volume bounded by the computed vertices
          call SCREWED_VOLUME (Flux_Vol, i, Volume)

          ! Now compare this volume with the actual Flux_Vol
          if (ABS(Flux_Vol(i)%Vol - Volume) < cutvof*Cell(i)%Volume) then
             converged = .true.
             exit
          else
             Mult = Mult * Flux_Vol(i)%Vol/Volume
          endif

       end do

       ! Print out a warning message if we iterated up to the maximum
       if (.not.converged) then
          write(errmsg,'(a,i0,a,es12.5,a,i0)') &
              'Flux volume vertex iteration did not converge in ', &
              flux_vol_iter_max, '. Maximum flux volume difference is ', &
              ABS(Flux_Vol(i)%Vol-Volume), ' in cell ', i
          call TLS_warn (errmsg)
       end if

    end do

  END SUBROUTINE FLUX_VOL_VERTICES

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE SCREWED_VOLUME (Flux_Vol, icell, Volume)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute a hexahedral cell volume, according to the prescription
    !   given by J. Dukowicz, JCP 74: 493-496 (1988).
    !
    !=======================================================================
    use legacy_mesh_api, only: ncells

    ! Arguments
    type(FLUX_VOL_QUANTITY), dimension(ncells), intent(IN) :: Flux_Vol
    integer, intent(IN) :: icell
    real(r8), intent(OUT) :: Volume

    ! Local Variables
    real(r8), dimension(ndim) :: X1, X2, X3
    integer :: f, v1, v2, v3, v4, v5, v6, n
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Volume = 0.0_r8

    ! Loop over the six faces of the flux volume
    do f = 1,6

       select case (f)
          case (1)
             ! Side 1 (vertices 4-8-7-3)
             v1 = 8; v2 = 4; v3 = 7; v4 = 8; v5 = 3; v6 = 4
          case (2)
             ! Side 2 (vertices 1-2-6-5)
             v1 = 6; v2 = 2; v3 = 5; v4 = 6; v5 = 1; v6 = 2
          case (3)
             ! Side 3 (vertices 4-1-5-8)
             v1 = 5; v2 = 1; v3 = 8; v4 = 5; v5 = 4; v6 = 1
          case (4)
             ! Side 4 (vertices 3-7-6-2)
             v1 = 7; v2 = 3; v3 = 6; v4 = 7; v5 = 2; v6 = 3
          case (5)
             ! Side 5 (vertices 4-3-2-1)
             v1 = 3; v2 = 4; v3 = 2; v4 = 3; v5 = 1; v6 = 4
          case (6)
             ! Side 6 (vertices 8-5-6-7)
             v1 = 6; v2 = 5; v3 = 7; v4 = 6; v5 = 8; v6 = 5
       end select

       do n = 1,ndim
          X1(n) = Flux_Vol(icell)%Xv(v1,n) + Flux_Vol(icell)%Xv(v2,n)
          X2(n) = Flux_Vol(icell)%Xv(v3,n) + Flux_Vol(icell)%Xv(v4,n)
          X3(n) = Flux_Vol(icell)%Xv(v5,n) + Flux_Vol(icell)%Xv(v6,n)
       end do

       Volume = Volume + X1(1)*(X2(2)*X3(3) - X3(2)*X2(3)) &
                       + X1(2)*(X3(1)*X2(3) - X2(1)*X3(3)) &
                       + X1(3)*(X2(1)*X3(2) - X3(1)*X2(2))

    end do

    Volume = Volume / 12.0_r8

    ! Punt if Volume < 0
    if (Volume < 0.0_r8) then
       write(message,'(a,i0,a)') 'SCREWED_VOLUME: cell ', icell, ' contains a negative flux volume'
       call TLS_panic (message)
    end if

  END SUBROUTINE SCREWED_VOLUME

END MODULE FLUX_VOLUME_MODULE
