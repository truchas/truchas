!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE INTERFACE_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define interface derived types containing all interface
  !   and interface cell data necessary to perform a piecewise
  !   planar reconstruction.
  !
  ! Public Interface(s):
  !
  ! Contains: INTERFACE_DATA
  !           INTERFACE_FLUX_DATA
  !           INTERFACE_MODEL_DEFAULT
  !           SETUP_INTERFACE_GEOM
  !           SETUP_INTERFACE_FLUX
  !
  ! Author(s): S. Jay Mosso (LANL Group X-HM, sjm@lanl.gov)
  !            Matthew Williams (MST-8, mww@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use legacy_mesh_api, only: ncells, ndim, nvc
  implicit none
  private

  ! Public Procedures and Types
  public :: INTERFACE_DATA,           &
            INTERFACE_FLUX_DATA,      &
            SETUP_INTERFACE_GEOM,     &
            SETUP_INTERFACE_FLUX

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! character string in NUMERICS namelist
  character(LEN = 80), public, save :: interface_topology_model

  type INTERFACE_DATA
 
     ! The interface normal.
     real(r8) :: Normal(ndim)
 
     ! The interface plane constant.
     real(r8) :: Rho
 
     ! Volume fraction of the material(s) behind this interface.
     real(r8) :: Vof

     ! Area of the interface for material(s)
     real(r8) :: Area
 
     ! The donor cell volume.
     real(r8) :: Cell_Volume
 
     ! The physical coordinates of the cell vertices of the interface cells.
     real(r8) :: Cell_Coord(ndim,nvc)
 
  end type INTERFACE_DATA

  ! Declare an array of INTERFACE_DATA types
  type(INTERFACE_DATA), pointer, public, dimension(:) :: Int_Geom => null()

  type INTERFACE_FLUX_DATA

     ! The flux volume.
     real(r8) :: Flux_Vol

     ! The physical coordinates of the flux volume vertices of the interface cells.
     real(r8) :: Flux_Vol_Coord(ndim,nvc)
 
     ! The face identifier for this advection sweep.
     integer :: Face
 
     ! The advected volume of the materials behind this interface.
     real(r8) :: Advected_Volume

     ! The iteration count required to locate this plane.
     integer :: Iter

   end type INTERFACE_FLUX_DATA

  ! Declare an array of INTERFACE_FLUX_DATA types
  type(INTERFACE_FLUX_DATA), pointer, public, dimension(:) :: Int_Flux => null()


CONTAINS

  SUBROUTINE SETUP_INTERFACE_GEOM (Mask, G1, G2, G3, Vof)
    !=======================================================================
    !  Purpose(s):
    !
    !    Perform a masked pack of interface parameters into 
    !    the INTERFACE_DATA derived type.
    !
    !=======================================================================
    use legacy_mesh_api, only: Cell, gather_vertex_coord
 
    ! Arguments
    logical,  dimension(ncells), intent(IN) :: Mask
    real(r8), dimension(ncells), intent(IN) :: G1, G2, G3
    real(r8), dimension(ncells), intent(IN) :: Vof
 
    ! Local Variables
    integer :: n, v
    real(r8), dimension(ncells)     :: Tmp
    real(r8), dimension(nvc,ncells) :: Vtx1
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
    ! Pack the interface normal vector.
    Int_Geom%Normal(1) = PACK (G1, Mask)
    Int_Geom%Normal(2) = PACK (G2, Mask)
    Int_Geom%Normal(3) = PACK (G3, Mask)
 
    ! Pack the donor volume fraction.
    Int_Geom%Vof = PACK(Vof, Mask)
 
    ! The correct donor cell volume was put into Volume
    ! in the code above.  Just pack it into the arrays.
    Int_Geom%Cell_Volume = PACK (Cell%Volume, Mask)
 
    ! Gather the donor cell's vertices
    do n = 1,ndim
 
       ! Added the BOUNDARY argument.  Now off-pe data motion only done
       ! on first time through.
       call gather_vertex_coord (Vtx1, dim=n)
 
       do v = 1, nvc

          Tmp = Vtx1(v,:) 
          Int_Geom%Cell_Coord(n,v) = PACK (Tmp, Mask)
 
       end do
 
    end do
 
  END SUBROUTINE SETUP_INTERFACE_GEOM

  SUBROUTINE SETUP_INTERFACE_FLUX (Flux_vol, Mask)
    !=======================================================================
    !  Purpose(s):
    !
    !    Perform a masked pack of interface parameters into 
    !    the INTERFACE_DATA derived type.
    !
    !=======================================================================
    use flux_volume_module, only: FLUX_VOL_QUANTITY
 
    ! Arguments
    type(FLUX_VOL_QUANTITY), dimension(ncells), intent(IN) :: Flux_Vol
    logical, dimension(ncells), intent(IN) :: Mask
 
    ! Local Variables
    integer :: n, v
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Load the Flux_Vol
    Int_Flux%Flux_Vol = PACK (Flux_Vol%Vol, Mask)

    ! Gather the donor cell's flux volume vertices
    do v = 1, nvc
       do n = 1,ndim
          Int_Flux%Flux_Vol_Coord(n,v) = PACK(Flux_Vol%Xv(v,n),  Mask)
       end do
    end do

    ! This is needed for a directionally-swept algorithm.
!    Int_Flux%Face = PACK (Flux_Vol%Fc, Mask)
 
    ! Advected volume behind the interface (initialize to zero).
    Int_Flux%Advected_Volume = 0.0_r8
 
  END SUBROUTINE SETUP_INTERFACE_FLUX

END MODULE INTERFACE_MODULE
