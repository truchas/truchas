!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

MODULE OVERWRITE_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures that overwrite the standard initialization of the
  !   Matl, and Zone derived types.  These are null routines.  If you
  !   would like to see some examples, look in overwrite.F90.examples.
  !
  ! Public Interface:
  !
  !   * call OVERWRITE_MATL ()
  !
  !       Overwrite the standard initialization of the Matl derived
  !       type.
  !
  !   * call OVERWRITE_VEL ()
  !
  !       Overwrite the standard initialization of the velocity part of
  !       the Zone derived type (Zone%Vc).
  !
  !   * call OVERWRITE_ZONE ()
  !
  !       Overwrite the standard initialization of the Zone derived
  !       type.
  !
  ! Contains: OVERWRITE_MATL
  !           OVERWRITE_VEL
  !           OVERWRITE_ZONE
  !
  ! Author(s): Bryan Lally, (lally@lanl.gov)
  !
  !=======================================================================
  use truchas_logging_services
  implicit none
  private

  ! Public Procedures
  public :: OVERWRITE_MATL, OVERWRITE_VEL, OVERWRITE_ZONE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE OVERWRITE_MATL ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Overwrite the standard initialization of cell-centered
    !   material quantities computed in MATL_INIT
    !
    !=======================================================================
  END SUBROUTINE OVERWRITE_MATL

  SUBROUTINE OVERWRITE_VEL ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Overwrite the standard initialization of cell-centered
    !   velocities computed in VELOCITY_INIT
    !   velocity fields for testing the behavior of nonsoleniodal velocities
    !   and divergence free velocity fields have been added, also a velocity
    !   field for a 3D vortex can be used for rotating rigid bodies
    !   (with modifications to other modules)
    !   These velocity fields are currently commented out.
    !
    !=======================================================================
!    use constants_module, only: pi
!    use legacy_mesh_api,  only: ncells, nnodes, ndim
!    use zone_module,      only: Zone
!    use legacy_mesh_api,  only: Mesh, Cell, Vertex
!    use matl_module,      only: Matl, max_slots

    ! Arguments

    ! Local Variables
!    integer :: n, tes
!    real(r8), dimension(ncells)  :: D

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize XY-Plane Vortex .............................
    ! Streamfunction Given in BCG (JCP 85:257,1989)
!    do n = 1, ndim
!       select case (n)
!       case (1) ! x-component
!          Zone%Vc(n) =  -2.0*SIN(pi*Cell%Centroid(1))*SIN(pi*Cell%Centroid(1)) &
!                           *SIN(pi*Cell%Centroid(2))*COS(pi*Cell%Centroid(2))
!       case (2) ! y-component
!          Zone%Vc(n) =  2.0*SIN(pi*Cell%Centroid(2))*SIN(pi*Cell%Centroid(2)) &
!                           *SIN(pi*Cell%Centroid(1))*COS(pi*Cell%Centroid(1))
!       case (3) ! z-component
!          Zone%Vc(n) = 0.0
!       end select
!    end do

    ! Initialize YZ-Plane Vortex .............................
    ! Streamfunction Given in BCG (JCP 85:257,1989)
!    do n = 1, ndim
!       select case (n)
!       case (1) ! x-component
!          Zone%Vc(n) = 0.0
!       case (2) ! y-component
!          Zone%Vc(n) = -2.0 *SIN(pi*(Cell%Centroid(2) ))   &
!                           *SIN(pi*(Cell%Centroid(2) ))    &
!                           *SIN(pi*(Cell%Centroid(3) ))    &
!                           *COS(pi*(Cell%Centroid(3) ))
!       case (3) ! z-component
!          Zone%Vc(n) =  2.0*SIN(pi*(Cell%Centroid(3) ))    &
!                           *SIN(pi*(Cell%Centroid(3) ))    &
!                           *SIN(pi*(Cell%Centroid(2) ))    &
!                           *COS(pi*(Cell%Centroid(2) ))
!       end select
!    end do

    ! Initialize Tilted Divergence Free Rotational Field
    ! Streamfunction Given in Puckett/Saltzman (Physica D. Nov 1992)
!    do n = 1, ndim
!       select case (n)
!       case (1) ! x-component
!          Zone%Vc(n) = (2.0*pi*Dsqrt(3)/3.0)*((Cell%Centroid(3)-0.5) &
!                              - (Cell%Centroid(2)-0.5))

!       case (2) ! y-component
!          Zone%Vc(n) = (2.0*pi*Dsqrt(3)/3.0)*((Cell%Centroid(1)-0.5) &
!                              - (Cell%Centroid(3)-0.5))
!       case (3) ! z-component
!          Zone%Vc(n) = (2.0*pi*Dsqrt(3)/3.0)*((Cell%Centroid(2)-0.5) &
!                              - (Cell%Centroid(1)-0.5))
!       end select
!    end do


! ... Add nonsolenoidal noise (all velocities entering every other vertex)

!      do n=1,ncells
!       tes = Int((Cell(n)%Centroid(1)/0.125)+2.0)
!      Zone(n)%Vc(1) =  Zone(n)%Vc(1)  +  &
!             0.5*(-1)**tes
!       tes = Int((Cell(n)%Centroid(2)/0.125)+2.0)
!      Zone(n)%Vc(2) =  Zone(n)%Vc(2)  +  &
!             0.5*(-1)**tes
!       tes = Int((Cell(n)%Centroid(3)/0.125)+2.0)
!      Zone(n)%Vc(3) =  Zone(n)%Vc(3) +  &
!             0.5*(-1)**tes
!       enddo

! ... Add nonsolenoidal noise (diagonally coupled velocity field)

!      do n=1,ncells
!      Zone(n)%Vc(1) =  Zone(n)%Vc(1)  +  &
!             0.5*(-1)**(n)*Sin(8*pi*Cell(n)%Centroid(2))
!      Zone(n)%Vc(2) =  Zone(n)%Vc(2)  +  &
!             0.5*(-1)**(n+1)*Sin(8*pi*Cell(n)%Centroid(2))
!      Zone(n)%Vc(3) =  Zone(n)%Vc(3) +  &
!             0.0
!       enddo

! ... Add linear velocity field for debugging and translating

!       do n=1,ncells
!       Zone(n)%Vc(1) =  0.0
!       Zone(n)%Vc(2) =  -1.0
!       Zone(n)%Vc(3) =   1.0
!        enddo


    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  END SUBROUTINE OVERWRITE_VEL

  SUBROUTINE OVERWRITE_ZONE ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Overwrite the standard initialization of cell-centered
    !   quantities computed in ZONE_INIT.
    !
    !=======================================================================

!!$    use error_module,         only: ERROR_CHECK
!!$    use matl_module,          only: GATHER_VOF, max_slots, mat_slot, nmat
!!$    use legacy_mesh_api,      only: Mesh, Cell, Vertex
!!$    use legacy_mesh_api,      only: ncells, nnodes
!!$    use thermo,               only: Volume_Fraction
!!$    use zone_module,          only: Zone
!!$    use thermo_iterative,     only: h_of_t
!!$
!!$    ! Example of how to specify an initial temperature profile for all
!!$    ! mold materials based on a gradient from top to bottom.
!!$
!!$    ! local variables
!!$    real(r8), dimension(ncells)   :: metal_vof, slope
!!$    real(r8)                      :: temperature_bottom, temperature_top
!!$    real(r8)                      :: y_bottom, y_top
!!$    integer                    :: metal_material_id
!!$    real(r8), allocatable, dimension(:) :: Vof
!!$    integer                    :: n
!!$
!!$    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!!$    ALLOCATE (Volume_Fraction(0:nmat,ncells), STAT=n)
!!$    call ERROR_CHECK ((n /= 0), (/'FATAL: error allocating Volume_Fraction'/), 'OVERWRITE_ZONE')
!!$    ALLOCATE (Vof(ncells), STAT=n)
!!$    call ERROR_CHECK ((n /= 0), (/'FATAL: error allocating Vof'/), 'OVERWRITE_ZONE')
!!$
!!$    metal_material_id  = 1
!!$    y_bottom           = 0.000
!!$    y_top              = 0.146
!!$    temperature_bottom = 1000.
!!$    temperature_top    = 1200.
!!$
!!$    call GATHER_VOF (metal_material_id, metal_vof)
!!$
!!$    where (metal_vof < 0.5)
!!$       slope = (temperature_top - temperature_bottom) / (y_top - y_bottom)
!!$       Zone%Temp = temperature_bottom + slope * (Cell%Centroid(2) - y_bottom)
!!$    end where
!!$
!!$    ! Volume Fraction is need later inside of the h_of_t call
!!$    Volume_Fraction = 0.D0
!!$    do n = 1,nmat
!!$       call GATHER_VOF (n, Vof)
!!$       Volume_Fraction(n,:) = Vof
!!$    end do
!!$    !After temperature are overwritten we need to calculate new enthalpies
!!$    ! note: if OVERWRITE_MATL has been altered this won't work, you also need
!!$    ! to update the volume_fraction array defined in the thermo module
!!$    Zone(:)%Enthalpy = h_of_t(Volume_Fraction,Zone(:)%Temp)
!!$    DEALLOCATE (Volume_Fraction, STAT=n)
!!$    call ERROR_CHECK ((n /= 0), (/'FATAL: error deallocating Volume_Fraction'/), 'OVERWRITE_ZONE')
!!$    DEALLOCATE (Vof, STAT=n)
!!$    call ERROR_CHECK ((n /= 0), (/'FATAL: error deallocating Vof'/), 'OVERWRITE_ZONE')

  END SUBROUTINE OVERWRITE_ZONE

END MODULE OVERWRITE_MODULE
