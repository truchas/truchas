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
  !   BC, Matl, and Zone derived types.  These are null routines.  If you
  !   would like to see some examples, look in overwrite.F90.examples.
  !
  ! Public Interface:
  !
  !   * call OVERWRITE_BC ()
  !       Overwrite the standard initialization of the BC derived
  !       type.
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
  ! Contains: OVERWRITE_BC
  !           OVERWRITE_MATL
  !           OVERWRITE_VEL
  !           OVERWRITE_ZONE
  !
  ! Author(s): Bryan Lally, (lally@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  ! Public Procedures
  public :: OVERWRITE_BC,  OVERWRITE_MATL, &
            OVERWRITE_VEL, OVERWRITE_ZONE, &
            PRESCRIBE_VELOCITY, CREATE_PLUME

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE OVERWRITE_BC ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Overwrite the standard initialization of the BC derived type
    !   as computed in BC_INIT
    !
    !=======================================================================
  END SUBROUTINE OVERWRITE_BC

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
!    use parameter_module, only: max_slots
!    use legacy_mesh_api,  only: ncells, nnodes, ndim
!    use zone_module,      only: Zone
!    use legacy_mesh_api,  only: Mesh, Cell, Vertex
!    use matl_module,      only: Matl

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
!!$    use matl_module,          only: GATHER_VOF
!!$    use legacy_mesh_api,      only: Mesh, Cell, Vertex
!!$    use parameter_module,     only: max_slots, mat_slot, nmat
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

  !-----------------------------------------------------------------------------

  subroutine BoundaryBetweenMaterials (mask, material_1, material_2)
     !--------------------------------------------------------------------------
     ! Find all faces that have a cell made of material_1 on one side
     ! and material_2 on the other side.
     !
     ! mask should be declared logical, (nfc,ncells) in the caller. 
     ! It is set to true for the boundary faces.
     !
     ! THRESHOLD is in case a cell isn't completely made of
     ! material_1 or material_2.  Set as close to 1.0 as desired.
     !--------------------------------------------------------------------------

     use legacy_mesh_api,   only: ncells, nfc, EE_GATHER
     use matl_module,       only: gather_vof

     ! arguments
     logical, dimension(:,:) :: mask
     integer :: material_1
     integer :: material_2

     ! local variables
     real(r8), allocatable, dimension(:)   :: vof_1
     real(r8), allocatable, dimension(:)   :: vof_2
     real(r8), allocatable, dimension(:,:) :: vof_tmp
     real(r8), parameter :: THRESHOLD = 0.99
     integer :: c, f, status

     !--------------------------------------------------------------------------

     ! make sure mask is the right size
     ASSERT(size(mask,1) == nfc)
     ASSERT(size(mask,2) == ncells)

     ! get some working space
     allocate(vof_1(ncells),STAT=status)
     if (status /= 0) call TLS_panic ('BoundaryBetweenMaterials: allocate failed: vof_1')
     allocate(vof_2(ncells),STAT=status)
     if (status /= 0) call TLS_panic ('BoundaryBetweenMaterials: allocate failed: vof_2')
     allocate(vof_tmp(nfc,ncells),STAT=status)
     if (status /= 0) call TLS_panic ('BoundaryBetweenMaterials: allocate failed: vof_tmp')

     ! start with a null boundary
     mask = .false.

     ! get vof of material 1 in each cell
     call gather_vof (material_1,vof_1)

     ! get vof of material 2 in each cell
     call gather_vof (material_2,vof_2)

     ! gather vof of material 2 across cell faces
     call ee_gather (vof_tmp, vof_2)

     ! look for cells made of (mostly) material 1
     do c = 1, ncells
        if (vof_1(c) > THRESHOLD) then
           ! that touch cells made of (mostly) material 2 across faces
           do f = 1, nfc
              if (vof_tmp(f,c) > THRESHOLD) then
                 ! call this part of the boundary
                 mask(f,c) = .true.
              end if
           end do
        end if
     end do

     ! gather vof of material 1 across cell faces
     call ee_gather (vof_tmp, vof_1)

     ! look for cells made of (mostly) material 2
     do c = 1, ncells
        if (vof_2(c) > THRESHOLD) then
           ! that touch cells made of (mostly) material 1 across faces
           do f = 1, nfc
              if (vof_tmp(f,c) > THRESHOLD) then
                 ! call this part of the boundary
                 mask(f,c) = .true.
              end if
           end do
        end if
     end do

     ! release working space
     deallocate(vof_tmp)
     deallocate(vof_2)
     deallocate(vof_1)

  end subroutine BoundaryBetweenMaterials

  !-----------------------------------------------------------------------------

  subroutine BoundaryTExternalUncovered (mask)
     !--------------------------------------------------------------------------
     ! Find all external faces that don't have thermal boundary conditions.
     !
     ! mask should be declared logical, (nfc,ncells) in the caller. 
     ! It is set to true for the uncovered boundary faces.
     !--------------------------------------------------------------------------

     use legacy_mesh_api,   only: ncells, nfc, Mesh
     use bc_module,         only: Boundary, BC_T

     ! arguments
     logical, dimension(:,:) :: mask

     ! local variables
     integer :: f

     !--------------------------------------------------------------------------

     ! make sure mask is the right size
     ASSERT(size(mask,1) == nfc)
     ASSERT(size(mask,2) == ncells)

     do f = 1, nfc
        mask(f,:) = Mesh%Ngbr_cell(f) == 0 .and. .not. Boundary(BC_T,f)
     end do

  end subroutine BoundaryTExternalUncovered

  SUBROUTINE CREATE_PLUME(Phi,xc,yc,zc,r1,type)

    use legacy_mesh_api, only: ncells, ndim, Cell

    ! Arguments...
    real(r8), dimension(:), intent(INOUT) :: Phi
    real(r8), intent(IN) :: xc,yc,zc,r1
    character(*), intent(IN) :: type

    integer :: nc, n
    real(r8) :: sum, r, std, ro
    real(r8), dimension(ndim) :: xcent

    xcent(1) = xc
    xcent(2) = yc
    xcent(3) = zc

    std = 3.0*.02;

    ro = r1

    select case(type)

    case ('bounded')

       do nc = 1,ncells

          Phi(nc) = 0.0_r8
          sum     = 0.0_r8
          r       = 0.0_r8
          ! calculate distance...
          do n = 1, ndim
             r = r + (Cell(nc)%Centroid(n) - xcent(n))**2;
          end do
          r = sqrt(r)

          ! a bounded plume...
          if (r < ro) then
             Phi(nc) = 5.0 * (1.0 + cos(r * 3.141592654 / ro)) / 10.00;
          end if

       end do ! ncells loop...

    case('gaussian')

       do nc = 1,ncells
          Phi(nc) = 0.0_r8
          sum     = 0.0_r8
          r       = 0.0_r8
          ! calculate distance...
          do n = 1, ndim
             r = r + (Cell(nc)%Centroid(n) - xcent(n))**2;
          end do
          r = sqrt(r)

          ! a guassian plume...
          Phi(nc) = exp(-r * r / 2.0 / std / std);

       end do

    case ('box')

       do nc = 1,ncells

          Phi(nc) = 0.0_r8
          ! a step...
          if (Cell(nc)%Centroid(1) < 0.5 .AND. Cell(nc)%Centroid(1) > 0.1) then
             if (Cell(nc)%Centroid(3) < 0.8 .AND. Cell(nc)%Centroid(3) > 0.2) then
                Phi(nc) = 1.0
             end if
          end if



          ! a square...
          if (Cell(nc)%Centroid(1) < 0.4 .AND. Cell(nc)%Centroid(1) > 0.1) then

             if (Cell(nc)%Centroid(3) < 0.4 .AND. Cell(nc)%Centroid(3) > 0.1) then
                Phi(nc) = 1.0
             end if

          end if

       end do
    case ('phiequalsx')

       do nc = 1,ncells
          Phi(nc) = Cell(nc)%Centroid(1)
       end do

    case ('phiequalsy')

       do nc = 1,ncells
          Phi(nc) = Cell(nc)%Centroid(2)
       end do

    case ('phiequalsz')

       do nc = 1,ncells
          Phi(nc) = Cell(nc)%Centroid(3)
       end do

    end select

  END SUBROUTINE CREATE_PLUME

  SUBROUTINE PRESCRIBE_VELOCITY (type)
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
    use constants_module, only: pi
    use fluid_data_module, only: Fluxing_Velocity
    use legacy_mesh_api,  only: ncells, ndim, nfc, Cell
    use zone_module,      only: Zone

    ! Arguments

    character(*), intent(IN) :: type

    ! Local Variables
    integer  :: n, f,i,j,k
    real(r8) :: Vfx, Vfy, Vfz
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    select case(type)

    case('diagxy')

       do n=1,ncells
          ! cell centers...
          Zone(n)%Vc(1) =  1.0d0
          Zone(n)%Vc(2) =  1.0d0
          Zone(n)%Vc(3) =  0.0d0
          ! face velocities used for fluxing...
          do f=1,nfc
             Vfx = 1.0d0
             Vfy = 1.0d0
             Vfz = 0.0d0
             Fluxing_Velocity(f,n) = Vfx*Cell(n)%Face_Normal(1,f) + &
                                     Vfy*Cell(n)%Face_Normal(2,f) + &
                                     Vfz*Cell(n)%Face_Normal(3,f) 
          end do
       enddo

    case('diagxz')

       do n=1,ncells
          ! cell centers...
          Zone(n)%Vc(1) =  1.0d0
          Zone(n)%Vc(2) =  0.0d0
          Zone(n)%Vc(3) =  1.0d0
          ! face velocities used for fluxing...
          do f=1,nfc
             Vfx = 1.0d0
             Vfy = 0.0d0
             Vfz = 1.0d0
             Fluxing_Velocity(f,n) = Vfx*Cell(n)%Face_Normal(1,f) + &
                                     Vfy*Cell(n)%Face_Normal(2,f) + &
                                     Vfz*Cell(n)%Face_Normal(3,f) 
          end do
       enddo

    case('diagyz')

       do n=1,ncells
          ! cell centers...
          Zone(n)%Vc(1) =  0.0d0
          Zone(n)%Vc(2) =  1.0d0
          Zone(n)%Vc(3) =  1.0d0
          ! face velocities used for fluxing...
          do f=1,nfc
             Vfx = 0.0d0
             Vfy = 1.0d0
             Vfz = 1.0d0
             Fluxing_Velocity(f,n) = Vfx*Cell(n)%Face_Normal(1,f) + &
                                     Vfy*Cell(n)%Face_Normal(2,f) + &
                                     Vfz*Cell(n)%Face_Normal(3,f) 
          end do
       enddo

    case('singlevortexxy')

       ! Initialize XY-Plane Vortex .............................
       ! Streamfunction Given in BCG (JCP 85:257,1989)
       do n = 1, ndim
          select case (n)
          case (1) ! x-component
             Zone%Vc(n) = -2.0 *SIN(pi*(Cell%Centroid(1) ))   &
                  *SIN(pi*(Cell%Centroid(1) ))    &
                  *SIN(pi*(Cell%Centroid(2) ))    &
                  *COS(pi*(Cell%Centroid(2) ))
          case (2) ! y-component
             Zone%Vc(n) =  2.0*SIN(pi*(Cell%Centroid(2) ))    &
                  *SIN(pi*(Cell%Centroid(2) ))    &
                  *SIN(pi*(Cell%Centroid(1) ))    &
                  *COS(pi*(Cell%Centroid(1) ))
          case (3) ! z-component
             Zone%Vc(n) = 0.0d0
          end select
       end do

       ! now for faces...
       do i=1,ncells
          do j=1,nfc
             do k=1,ndim
                select case(k)
                case (1)
                   Vfx =  -2.0 *SIN(pi*(Cell(i)%Face_Centroid(1,j) ))   &
                               *SIN(pi*(Cell(i)%Face_Centroid(1,j) ))   &
                               *SIN(pi*(Cell(i)%Face_Centroid(2,j) ))   &
                               *COS(pi*(Cell(i)%Face_Centroid(2,j) ))
                case (2)
                   Vfy =   2.0 *SIN(pi*(Cell(i)%Face_Centroid(2,j) ))   &
                               *SIN(pi*(Cell(i)%Face_Centroid(2,j) ))   &
                               *SIN(pi*(Cell(i)%Face_Centroid(1,j) ))   &
                               *COS(pi*(Cell(i)%Face_Centroid(1,j) ))
                case (3)
                   Vfz = 0.0d0
                end select
             end do
             Fluxing_Velocity(j,i) = Vfx*Cell(i)%Face_Normal(1,j) + &
                                     Vfy*Cell(i)%Face_Normal(2,j) + &
                                     Vfz*Cell(i)%Face_Normal(3,j) 
          end do
       end do

    case('singlevortexxz')

       ! Initialize XZ-Plane Vortex .............................
       ! Streamfunction Given in BCG (JCP 85:257,1989)
       do n = 1, ndim
          select case (n)
          case (1) ! x-component
             Zone%Vc(n) = -2.0 *SIN(pi*(Cell%Centroid(1) ))   &
                  *SIN(pi*(Cell%Centroid(1) ))    &
                  *SIN(pi*(Cell%Centroid(3) ))    &
                  *COS(pi*(Cell%Centroid(3) ))
          case (2) ! y-component
             Zone%Vc(n) = 0.0d0
          case (3) ! z-component
             Zone%Vc(n) =  2.0*SIN(pi*(Cell%Centroid(3) ))    &
                  *SIN(pi*(Cell%Centroid(3) ))    &
                  *SIN(pi*(Cell%Centroid(1) ))    &
                  *COS(pi*(Cell%Centroid(1) ))
          end select
       end do

       ! now for faces...
       do i=1,ncells
          do j=1,nfc
             do k=1,ndim
                select case(k)
                case (1)
                   Vfx =  -2.0 *SIN(pi*(Cell(i)%Face_Centroid(1,j) ))   &
                               *SIN(pi*(Cell(i)%Face_Centroid(1,j) ))   &
                               *SIN(pi*(Cell(i)%Face_Centroid(3,j) ))   &
                               *COS(pi*(Cell(i)%Face_Centroid(3,j) ))
                case (2)
                   Vfy = 0.0d0
                case (3)
                   Vfz =   2.0 *SIN(pi*(Cell(i)%Face_Centroid(3,j) ))   &
                               *SIN(pi*(Cell(i)%Face_Centroid(3,j) ))   &
                               *SIN(pi*(Cell(i)%Face_Centroid(1,j) ))   &
                               *COS(pi*(Cell(i)%Face_Centroid(1,j) ))
                end select
             end do
             Fluxing_Velocity(j,i) = Vfx*Cell(i)%Face_Normal(1,j) + &
                                     Vfy*Cell(i)%Face_Normal(2,j) + &
                                     Vfz*Cell(i)%Face_Normal(3,j) 
          end do
       end do

    case('singlevortexyz')

       ! Initialize XY-Plane Vortex .............................
       ! Streamfunction Given in BCG (JCP 85:257,1989)
       do n = 1, ndim
          select case (n)
          case (1) ! x-component
             Zone%Vc(n) = 0.0d0
          case (2) ! y-component
             Zone%Vc(n) =  2.0*SIN(pi*(Cell%Centroid(2) ))    &
                  *SIN(pi*(Cell%Centroid(2) ))    &
                  *SIN(pi*(Cell%Centroid(3) ))    &
                  *COS(pi*(Cell%Centroid(3) ))
          case (3) ! z-component
             Zone%Vc(n) = -2.0 *SIN(pi*(Cell%Centroid(3) ))   &
                  *SIN(pi*(Cell%Centroid(3) ))    &
                  *SIN(pi*(Cell%Centroid(2) ))    &
                  *COS(pi*(Cell%Centroid(2) ))
          end select
       end do

       ! now for faces...
       do i=1,ncells
          do j=1,nfc
             do k=1,ndim
                select case(k)
                case (1)
                   Vfx = 0.0d0
                case (2)
                   Vfy =   2.0 *SIN(pi*(Cell(i)%Face_Centroid(2,j) ))   &
                               *SIN(pi*(Cell(i)%Face_Centroid(2,j) ))   &
                               *SIN(pi*(Cell(i)%Face_Centroid(3,j) ))   &
                               *COS(pi*(Cell(i)%Face_Centroid(3,j) ))
                case (3)
                   Vfz =  -2.0 *SIN(pi*(Cell(i)%Face_Centroid(3,j) ))   &
                               *SIN(pi*(Cell(i)%Face_Centroid(3,j) ))   &
                               *SIN(pi*(Cell(i)%Face_Centroid(2,j) ))   &
                               *COS(pi*(Cell(i)%Face_Centroid(2,j) ))
                end select 
             end do
             Fluxing_Velocity(j,i) = Vfx*Cell(i)%Face_Normal(1,j) + &
                                     Vfy*Cell(i)%Face_Normal(2,j) + &
                                     Vfz*Cell(i)%Face_Normal(3,j) 
          end do
       end do

    case('rotationxy')

       do n = 1, ndim
          select case (n)
          case (1) ! x-component
             Zone%Vc(n) = -2.0 * (Cell%Centroid(2) - 0.5d0)
          case (2) ! y-component
             Zone%Vc(n) =  2.0 * (Cell%Centroid(1) - 0.5d0)
          case (3) ! z-component
             Zone%Vc(n) = 0.0d0
          end select
       end do

       ! now for faces...
       do i=1,ncells
          do j=1,nfc
             do k=1,ndim
                select case(k)
                case (1)
                   Vfx = -2.0d0*(Cell(i)%Face_Centroid(2,j) - 0.5d0)
                case (2)
                   Vfy =  2.0d0*(Cell(i)%Face_Centroid(1,j) - 0.5d0)
                case (3)
                   Vfz = 0.0d0
                end select
             end do
             Fluxing_Velocity(j,i) = Vfx*Cell(i)%Face_Normal(1,j) + &
                                     Vfy*Cell(i)%Face_Normal(2,j) + &
                                     Vfz*Cell(i)%Face_Normal(3,j) 
          end do
       end do

    case('rotationxz')

       do n = 1, ndim
          select case (n)
          case (1) ! x-component
             Zone%Vc(n) = -2.0 * (Cell%Centroid(3) - 0.5d0)
          case (2) ! y-component
             Zone%Vc(n) = 0.0d0
          case (3) ! z-component
             Zone%Vc(n) =  2.0 * (Cell%Centroid(1) - 0.5d0)
          end select
       end do

       ! now for faces...
       do i=1,ncells
          do j=1,nfc
             do k=1,ndim
                select case(k)
                case (1)
                   Vfx = -2.0d0*(Cell(i)%Face_Centroid(3,j) - 0.5d0)
                case (2)
                   Vfy = 0.0d0
                case (3)
                   Vfz =  2.0d0*(Cell(i)%Face_Centroid(1,j) - 0.5d0)
                end select
             end do
             Fluxing_Velocity(j,i) = Vfx*Cell(i)%Face_Normal(1,j) + &
                                     Vfy*Cell(i)%Face_Normal(2,j) + &
                                     Vfz*Cell(i)%Face_Normal(3,j) 
          end do
       end do

    case('rotationyz')

       do n = 1, ndim
          select case (n)
          case (1) ! x-component
             Zone%Vc(n) = 0.0d0
          case (2) ! y-component
             Zone%Vc(n) =  2.0 * (Cell%Centroid(3) - 0.5d0)
          case (3) ! z-component
             Zone%Vc(n) = -2.0 * (Cell%Centroid(2) - 0.5d0)
          end select
       end do

       ! now for faces...
       do i=1,ncells
          do j=1,nfc
             do k=1,ndim
                select case(k)
                case (1)
                   Vfx = 0.0d0
                case (2)
                   Vfy =  2.0d0*(Cell(i)%Face_Centroid(3,j) - 0.5d0)
                case (3)
                   Vfz = -2.0d0*(Cell(i)%Face_Centroid(2,j) - 0.5d0)
                end select
             end do
             Fluxing_Velocity(j,i) = Vfx*Cell(i)%Face_Normal(1,j) + &
                                     Vfy*Cell(i)%Face_Normal(2,j) + &
                                     Vfz*Cell(i)%Face_Normal(3,j) 
          end do
       end do

    end select

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  END SUBROUTINE PRESCRIBE_VELOCITY

END MODULE OVERWRITE_MODULE
