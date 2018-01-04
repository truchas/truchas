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

module flux_volume_modulez
  use kinds,  only: r8
  use truchas_logging_services
  implicit none
  private

  public :: flux_vol_quantity, flux_vol_vertices

  ! maximum allowable iterations for the flux volume vertex location
  integer,  parameter :: flux_vol_iter_max = 10

  type flux_vol_quantity
    integer  :: Fc           ! face number through which we this cell is fluxing
    real(r8) :: Vol          ! Volume of the flux
    real(r8) :: Xv(3,8) ! Vertices of the flux volume
  end type flux_vol_quantity

contains

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
  subroutine flux_vol_vertices (face, cell, is_mixed_donor_cell, dist, Flux_Vol, cutvof)
    use hex_types,     only: cell_data
    use cell_geometry, only: eval_hex_volumes

    integer,                 intent(in)    :: face
    type(cell_data),         intent(in)    :: cell
    logical,                 intent(in)    :: is_mixed_donor_cell
    real(r8),                intent(in)    :: dist, cutvof
    type(flux_vol_quantity), intent(inout) :: flux_vol

    integer :: n, v, ia, ib, e, iter, &
        Edge_ends(2,4,6) = reshape( [ &
        1,4,  2,3,  6,7,  5,8,   & ! face 1 edges
        2,1,  3,4,  7,8,  6,5,   & ! face 2 edges
        3,2,  4,1,  8,5,  7,6,   & ! face 3 edges
        1,2,  5,6,  8,7,  4,3,   & ! face 4 edges
        1,5,  4,8,  3,7,  2,6,   & ! face 5 edges
        5,1,  6,2,  7,3,  8,4] , & ! face 6 edges
        [2,4,6] )
    real(r8)       :: Percnt(4), tmp(8)
    real(r8)       :: Uedge(3,4)
    real(r8)       :: Volume, Mult, ndotuedge
    logical        :: converged
    character(256) :: errmsg

    ! first gather the vertex coordinates and store them in flux_vol%Xv
    Flux_Vol%Xv = cell%node

    ! need only go through this exercise for mixed donor cells and for outgoing fluxes
    if (.not.is_mixed_donor_cell .or. dist*cell%face_area(face) <= cutvof*cell%volume) return

    ! compute the edge unit vectors
    do e = 1, 4
      ia  = Edge_ends(1,e,face) ! front node
      ib  = Edge_ends(2,e,face) ! back node
      Uedge(:,e) = Flux_Vol%Xv(:,ib) - Flux_Vol%Xv(:,ia)

      ndotuedge = sum(cell%face_normal(:,face) * Uedge(:,e))
      Percnt(e) = -Dist/(ndotuedge+epsilon(1.0_r8))
    end do

    if (any(Percnt < 0.0_r8) .or. any(Percnt > 1.0_r8)) then
       do e = 1,8
         write(*,'(a,3es15.4)') 'flux nodes:',flux_vol%xv(:,e)
       end do
       write(*,*) dist
       do e = 1,4
         ndotuedge = sum(cell%face_normal(:,face) * Uedge(:,e))
         write(*,'(a,3es15.4)') 'ndotuedge: ',ndotuedge, -dist/(ndotuedge+epsilon(1.0_r8)),percnt(e)
       end do
      write(errmsg,'(a,3es13.6)') 'FLUX_VOL_VERTICES: invalid flux volume or inverted element; cell centroid =', &
           sum( cell%node, dim=2 ) / real(8,r8)
      call TLS_fatal (errmsg)
    end if

    ! iterate to find the 4 other vertices that define the back end of the flux volume
    mult = 1.0_r8; iter = 1; volume = 1e10_r8;
    do while (abs(flux_vol%vol - volume) > cutvof*cell%volume .and. iter<flux_vol_iter_max)
      ! loop over edges to adjust the vertices
      ! !$omp simd
      do e = 1,4
        ia = Edge_ends(1,e,face)
        ib = Edge_ends(2,e,face)
        Flux_Vol%Xv(:,ib) = Flux_Vol%Xv(:,ia) + Mult*Percnt(e)*Uedge(:,e)
      end do
      ! !$omp end simd

      ! compute the flux volume bounded by the computed vertices
      call eval_hex_volumes (flux_vol%xv, volume, tmp)

      ! increment multiplier for next iteration
      Mult = Mult * Flux_Vol%Vol/volume
      iter = iter+1
    end do

    ! print out a warning message if we iterated up to the maximum
    if (iter==flux_vol_iter_max) then
      write(errmsg,'(a,i0,a,es12.5,a,i0)') &
           'Flux volume vertex iteration did not converge in ', flux_vol_iter_max,&
           '. Maximum flux volume difference is ', abs(Flux_Vol%Vol-Volume) !,&
      !' in cell '!, cell_index
      call TLS_warn (errmsg)
    end if
  end subroutine flux_vol_vertices

end module flux_volume_modulez
