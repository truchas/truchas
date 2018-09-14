!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module surface_tension_module
  !=======================================================================
  ! Purpose:
  !
  !   Surface tension modeling variables and procedures.
  !
  !   Public Interface(s):
  !
  !     * call CSF (Cell, Mesh, Vertex, dt, Mom_Delta)
  !         Compute the cell-centered CSF force due to the 
  !         tangential component of surface tension
  !         and increment the momentum delta array.
  ! 
  !     * call CSF_FACE (Fcsf)
  !         Compute the face_centered CSF force due to 
  !         the normal component of surface tension (Fcsf).
  !        
  !
  ! Contains: CSF
  !           DELTA_FUNCTION
  !           CSF_FACE 
  !           DIVERGENCE_INTNORMAL
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use scalar_func_class
  use scalar_func_containers
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  public :: read_surface_tension_namelist
  Public :: CSF, CSF_FACE

  integer, parameter :: MAX_NAME_LEN = 31, MAX_FACE_SET_IDS = 32

  !! This flag is enabled by the PHYSICS namelist.
  logical, public, save :: surface_tension = .false.
  
  !! The remaining variables are defined by READ_SURFACE_TENSION_NAMELIST.
  logical, public, save :: csf_normal = .false.
  logical, public, save :: csf_tangential = .false.
  logical, public, save :: csf_boundary = .false.
  real,    public, save :: dsig_dT
  integer, allocatable, public :: face_set_ids(:)
  !! Volume mesh cell size in z direction for cells adjacent to side set on
  !! which "csf_boundary"-type surface tension is to be applied. Computed by
  !! FLUID_INIT().
  real(r8), allocatable, public :: csf_z(:)

  integer, public, save :: surfmat1, surfmat2
  class(scalar_func), allocatable, public :: sigma_func ! surface tension coefficient
  integer :: kernel = 0  ! initialized to one of the following values
  integer, parameter :: WILLIAMS_KERNEL = 1
  integer, parameter :: RUDMAN_KERNEL   = 2

contains

  SUBROUTINE CSF (dt, Mom_Delta)
    !=======================================================================
    ! Purpose(s):
    !   Compute cell-centered forces with the CSF (Continuum Surface
    !   Force) model for surface tension, as detailed in the Journal
    !   of Computational Physics, 100: 335-354 (1993). Increment the
    !   momentum delta array by dt*(Fx,Fy,Fz), where (Fx,Fy,Fz)
    !   is the CSF force.
    !=======================================================================
    use cutoffs_module,       only: alittle, cutvof
    use discrete_op_module,   only: GRADIENT
    use legacy_matl_api,      only: GATHER_VOF, nmat
    use legacy_mesh_api,      only: ncells, ndim
    use property_module,      only: density_material
    use zone_module,          only: Zone
    use legacy_mesh_api,      only: ncells, mesh_face_set

    real(r8), intent(in) :: dt
    real(r8), intent(inout) :: Mom_Delta(:,:)

    integer :: n, m, j, fminz, fmaxz, nssets
    real(r8) :: state(1), cz
    real(r8), dimension(ncells) :: Color, Delta
    real(r8), dimension(ncells) :: Fx, Fy, Fz
    real(r8), dimension(ncells) :: Kappa, Sigma
    real(r8), dimension(ncells) :: dC_dx, dC_dy, dC_dz
    real(r8), dimension(ncells) :: dS_dx, dS_dy, dS_dz

    ASSERT(csf_tangential .or. csf_boundary)
    ASSERT(size(Mom_Delta,1) == ndim)
    ASSERT(size(Mom_Delta,2) == ncells)

    ! start the surface tension timer
    call start_timer('Predictor Surface Tension')

    ! Initialize.
    Fx = 0.0_r8; Fy = 0.0_r8; Fz = 0.0_r8

    ! NNC, Dec 2012.  Formerly, surface tension was triggered when the first
    ! material having the surf10 property was encountered.  This property also
    ! pointed to the other material.  The COLOR field below was the volume
    ! fraction sum of all materials except this other material. I have preserved
    ! this behavior treating surfmat1 as the first material and surfmat2 as
    ! the other, though I do not understand why COLOR was not just taken to be
    ! the volume fraction of the first (or other) material.
    Color = 0.0_r8

    !- tangential surface tension applied on a boundary surface given by set id
    if (csf_boundary) then

      Fx=0.0_r8
      Fy=0.0_r8
      Fz=0.0_r8

      ! use variable dC_dx,dC_dy,dC_dz temporarly as the temperature gradient
      call GRADIENT(dC_dx,dC_dy, dC_dz,Zone%Temp,method='Green-Gauss')

      ! compute tangential surface tension force only in cells adjacent to bc
      do j = 1, ncells
        if (csf_z(j) > 0.0) then
          Fx(j) = dsig_dT * dC_dx(j) / csf_z(j)
          Fy(j) = dsig_dT * dC_dy(j) / csf_z(j)
          Fz(j) = 0.0_r8
        end if
      end do

    else

      do m = 1, nmat
        if (m == surfmat2) cycle
        call gather_vof (m, Kappa)
        Color = Color + Kappa
      end do

      ! Compute the cell-centered color gradient.
      call GRADIENT (dC_dx, dC_dy, dC_dz, Color, method='Green-Gauss')

      ! Compute the surface delta function.
      call DELTA_FUNCTION (Color, dC_dx, dC_dy, dC_dz, Delta)

      ! Normalize the color gradient (use Sigma as a temporary).
      Sigma = SQRT(dC_dx*dC_dx + dC_dy*dC_dy + dC_dz*dC_dz)
      Sigma = MERGE(0.0_r8, 1.0_r8/Sigma, Sigma <= alittle)
      dC_dx = Sigma*dC_dx; dC_dy = Sigma*dC_dy; dC_dz = Sigma*dC_dz

      ! Get the cell-centered surface tension coefficient.
      !call PROPERTY (Zone%Temp_Old, m, 'surface_tension', Value = Sigma)
      do j = 1, ncells
        state(1) = Zone(j)%Temp_Old
        sigma(j) = sigma_func%eval(state)
      end do

      ! Density weighting of the delta function
      Delta=Delta*2.*Zone%Rho/(density_material(surfmat1)+density_material(surfmat2))

      ! Compute the cell-centered sigma gradient
      call GRADIENT (dS_dx, dS_dy, dS_dz, Sigma, method = 'Green-Gauss')

      ! Construct the CSF tangential component. 
      where (abs(dC_dx).gt.cutvof.or.abs(dC_dy).gt.cutvof & 
             .or.abs(dC_dz).gt.cutvof)
         Fx = Fx + Delta*((1.0_r8 - dC_dx*dC_dx)*dS_dx - &
                                 dC_dy*dC_dx *dS_dy - &
                                 dC_dz*dC_dx *dS_dz)
         Fy = Fy + Delta*((1.0_r8 - dC_dy*dC_dy)*dS_dy - &
                                 dC_dz*dC_dy *dS_dz - &
                                 dC_dx*dC_dy *dS_dx)
         Fz = Fz + Delta*((1.0_r8 - dC_dz*dC_dz)*dS_dz - &
                                 dC_dx*dC_dz *dS_dx - &
                                 dC_dy*dC_dz *dS_dy)
      end where

    endif

    ! Increment Momentum Delta
    Mom_Delta(1,:) = Mom_Delta(1,:) + dt*Fx
    Mom_Delta(2,:) = Mom_Delta(2,:) + dt*Fy
    Mom_Delta(3,:) = Mom_Delta(3,:) + dt*Fz

    !stop the surface tension timer
    call stop_timer ('Predictor Surface Tension')

  END SUBROUTINE CSF

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

  SUBROUTINE DELTA_FUNCTION (Color, Nx, Ny, Nz, Delta)
    !=======================================================================
    ! Purpose(s):
    !   Compute a scalar cell-centered surface delta function for an
    !   interface which is defined by a characteristic (color) function.
    !   The color function is assumed to be unity on one side of the
    !   interface and zero on the other side.  (Nx, Ny, Nz) is provided
    !   as the cell-centered color gradient. 
    !=======================================================================
    use legacy_mesh_api, only: ncells

    ! Argument List
    real(r8), dimension(ncells), intent(IN)  :: Color, Nx, Ny, Nz 
    real(r8), dimension(ncells), intent(OUT) :: Delta

    ! Currently the delta function is just the color gradient magnitude.
    Delta = SQRT(Nx*Nx + Ny*Ny + Nz*Nz)

  END SUBROUTINE DELTA_FUNCTION

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE CSF_FACE (Fcsf_face)
    !=======================================================================
    ! Purpose(s):
    !   Compute face-centered forces with the CSF (Continuum Surface
    !   Force) model for the normal component of surface tension
    !  
    ! Note: The surface tension forces are computed at FACE-centered 
    !   to have a force balanced consistent formulation 
    !   with the pressure gradient                
    ! 
    ! References: 
    !   J. Comput. Physics, 213: 141-173 (2006) Francois et al.
    !   LA-UR-03-2128 report (ASME paper FEDS 2003) Francois et al.
    !   J. Comput. Physics, 100: 335-354 (1992) Brackbill et al.
    !        
    ! Author: Marianne M. Francois, LANL CCS-2 (mmfran@lanl.gov)
    !=======================================================================
    use cutoffs_module,              only: alittle, cutvof
    use discrete_derivatives,        only: GRADIENT_FACE
    use discrete_op_module,          only: GRADIENT,VERTEX_AVG 
    use discrete_ops_data,           only: use_ortho_face_gradient
    use legacy_matl_api,             only: GATHER_VOF, nmat
    use mollify,                     only: MOLLIFY_CONV_SAVEMEM,     &
                                           interface_smoothing_length
    use legacy_mesh_api,             only: ncells, ndim, nfc, nvc, nnodes
    use legacy_mesh_api,             only: EN_GATHER, LINEAR_PROP 
    use zone_module,                 only: Zone 
    use kernel_interpolation_module, only: KERN_CONVOLUTION_CENTER
 
    ! Argument List
    real(r8), dimension(ndim,nfc,ncells), intent(INOUT) :: Fcsf_face
 
   ! Local Variables
    integer :: j, m, n, f, status
    real(r8) :: state(1)
    !integer :: curvmodel
    real(r8)   :: d    
    real(r8), dimension(:),     allocatable :: Weight
    real(r8), dimension(:,:,:), allocatable :: dC
    real(r8), dimension(:,:),   allocatable :: Kappa_face
    real(r8), dimension(:,:),   allocatable :: Kappa_v
    real(r8), dimension(:),     allocatable :: Kappa_vtx
    real(r8), dimension(:,:),   allocatable :: Sigma_face
    real(r8), dimension(:,:),   allocatable :: Sigma_v
    real(r8), dimension(:),     allocatable :: Sigma_vtx 
    real(r8), dimension(:),     allocatable :: Kappa,Sigma,Tmp
    real(r8), dimension(:),     allocatable :: Color,mo_col
    real(r8), dimension(:,:),   allocatable :: mo_color
    real(r8), dimension(:,:),   allocatable :: dC_dx,dC_dy,dC_dz
    real(r8), dimension(:,:),   allocatable :: Fx,Fy,Fz
    real(r8), dimension(:,:),   allocatable :: dC1,dC2,dC3
    real(r8), dimension(:),     allocatable :: dCn1,dCn2,dCn3
    real(r8), dimension(:,:),   allocatable :: dCcc     
    logical,  dimension(:,:),   allocatable :: face_flag 
   
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
 
    ALLOCATE (Color(ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: Color(ncells) allocation failed')

    ALLOCATE (Kappa(ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: Kappa(ncells) allocation failed')
    
    ALLOCATE (dC(ndim,nfc,ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: dC(ndim,nfc,ncells) allocation failed')

    ALLOCATE (dC_dx(nfc,ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: dC_dx(ndim,nfc,ncells) allocation failed')
     
    ALLOCATE (dC_dy(nfc,ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: dC_dy(ndim,nfc,ncells) allocation failed')
 
    ALLOCATE (dC_dz(nfc,ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: dC_dz(ndim,nfc,ncells) allocation failed')
     
    ALLOCATE (Weight(ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: Weight(ncells) allocation failed')

    ALLOCATE (Sigma(ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: Sigma(ncells) allocation failed')

    ALLOCATE (Kappa_face(nfc,ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: Kappa_face(nfc,ncells) allocation failed')

    ALLOCATE (Sigma_face(nfc,ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: Sigma_face(nfc,ncells) allocation failed')

    ALLOCATE (Fx(nfc,ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: Fx(nfc,ncells) allocation failed')

    ALLOCATE (Fy(nfc,ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: Fy(nfc,ncells) allocation failed')

    ALLOCATE (Fz(nfc,ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: Fz(nfc,ncells) allocation failed')

    !-----------------------------------------------
    ! Start the surface tension timer
    call  start_timer ("Projection Surface Tension")

    ! Initialize.
    Fx=0.0_r8; Fy=0.0_r8; Fz=0.0_r8

    ! NNC, Dec 2012.  Formerly, surface tension was triggered when the first
    ! material having the surf10 property was encountered.  This property also
    ! pointed to the other material.  The COLOR field below was the volume
    ! fraction sum of all materials except this other material. I have preserved
    ! this behavior treating surfmat1 as the first material and surfmat2 as
    ! the other, though I don't understand why COLOR wasn't just taken to be
    ! the volume fraction of the first (or other) material.
    Color = 0.0_r8
    do m = 1, nmat
      if (m == surfmat2) cycle
      call gather_vof (m, Kappa)
      Color = Color + Kappa
    end do

    !-----------------------------------------------------
    !Compute the face-centered color gradient.          
    dC=0.0_r8
    dC_dx=0.0_r8; dC_dy=0.0_r8; dC_dz=0.0_r8
    Weight=1.0_r8
    call GRADIENT_FACE (PHI=Color,GRAD=dC,WEIGHT=Weight, &
                        USE_ORTHO=use_ortho_face_gradient)
    dC_dx=dC(1,:,:)
    dC_dy=dC(2,:,:)
    dC_dz=dC(3,:,:)

    dC_dx=MERGE(0.0_r8,dC_dx,abs(dC_dx) <= cutvof)
    dC_dy=MERGE(0.0_r8,dC_dy,abs(dC_dy) <= cutvof)
    dC_dz=MERGE(0.0_r8,dC_dz,abs(dC_dz) <= cutvof)

    !------------------------------------------------------
    !Compute curvature Kappa at cell-centers
    !------------------------------------------------------
    select case (kernel)
    case (WILLIAMS_KERNEL)

      ALLOCATE (dCcc(ndim,ncells),STAT=status)
      if (status/=0) call TLS_panic ('CSF_FACE: dCcc(ndim,ncells) allocation failed')

      ALLOCATE (dC1(nmat,ncells),STAT=status)
      if (status/=0) call TLS_panic ('CSF_FACE: dC1(nmat,ncells) allocation failed')

      ALLOCATE (dC2(nmat,ncells),STAT=status)
      if (status/=0) call TLS_panic ('CSF_FACE: dC2(nmat,ncells) allocation failed')

      ALLOCATE (dC3(nmat,ncells),STAT=status)
      if (status/=0) call TLS_panic ('CSF_FACE: dC3(nmat,ncells) allocation failed')

      ALLOCATE (mo_color(nmat,ncells),STAT=status)
      if (status/=0) call TLS_panic ('CSF_FACE: mo_color(nmat,ncells) allocation failed')

      call MOLLIFY_CONV_SAVEMEM(Color,mo_color,dC1,dC2,dC3,surfmat1)

      dCcc(1,:)=dC1(surfmat1,:)
      dCcc(2,:)=dC2(surfmat1,:)
      dCcc(3,:)=dC3(surfmat1,:)

      call DIVERGENCE_INTNORMAL(Kappa,dCcc)

      DEALLOCATE(dCcc)
      DEALLOCATE(dC1,dC2,dC3)
      DEALLOCATE(mo_color)

    case (RUDMAN_KERNEL)

      ALLOCATE (mo_col(ncells),STAT=status)
      if (status /=0) call TLS_panic ('CSF_FACE: mo_col(ncells) allocation failed')

      ALLOCATE(dCn1(ncells),STAT=status)
      if (status /=0) call TLS_panic ('CSF_FACE: dCn1(ncells) allocation failed')
      dCn1=0.0_r8

      ALLOCATE(dCn2(ncells),STAT=status)
      if (status /=0) call TLS_panic ('CSF_FACE: dCn2(ncells) allocation failed')
      dCn2=0.0_r8

      ALLOCATE(dCn3(ncells),STAT=status)
      if (status /=0) call TLS_panic ('CSF_FACE: dCn3(ncells) allocation failed')
      dCn3=0.0_r8

      ALLOCATE(dCcc(ndim,ncells),STAT=status)
      if (status /=0) call TLS_panic ('CSF_FACE: dCcc(ndim,ncells) allocation failed')
      dCcc=0.0_r8

      ALLOCATE(Tmp(ncells),STAT=status)
      if (status /=0) call TLS_panic ('CSF_FACE: Tmp(ncells) allocation failed')
      Tmp=0.0_r8

      d=interface_smoothing_length 

      call KERN_CONVOLUTION_CENTER(Color,mo_col,d)

      ! compute the gradient of the mollified Color function
      ! to obtain the normal vectors

      call GRADIENT(dCn1,dCn2,dCn3,mo_col)

      ! normalize the normal vectors

      Tmp=sqrt(dCn1**2+dCn2**2+dCn3**2)
      Tmp=MERGE(0.0_r8,1.0_r8/Tmp, Tmp<=alittle)
      where (mo_col > cutvof .and. mo_col < 1.0_r8 - cutvof)
        dCcc(1,:)=dCn1*Tmp
        dCcc(2,:)=dCn2*Tmp
        dCcc(3,:)=dCn3*Tmp
      end where

      ! compute the divergence of the normal vector
      ! to obtain curvature Kappa
      call DIVERGENCE_INTNORMAL(Kappa,dCcc)

      Kappa=-Kappa

      DEALLOCATE(dCn1,dCn2,dCn3,Tmp,dCcc,mo_col)

    case default
      INSIST(.false.)
    end select

    !---------------------------------------------
    ! Interpolate cell-centered curvature (Kappa) 
    ! at face-centers Kappa_face

    ! Allocate variables

    ALLOCATE (Kappa_vtx(nnodes),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: Kappa_vtx(nnodes) allocation failed')

    ALLOCATE (Kappa_v(nvc,ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: Kappa_v(nvc,ncells) allocation failed')

    Kappa_face=0.0_r8

    !-averaged Kappa at faces
    ! first averaged at vertices
    ! then at faces (Kappa is a scalar)

     call VERTEX_AVG(Kappa,Kappa_vtx) 
     call EN_GATHER(Kappa_v,Kappa_vtx)
     do f=1,nfc
       call LINEAR_PROP(f,Kappa_v,Kappa_face(f,:))
     enddo

    DEALLOCATE(Kappa_vtx,Kappa_v)

    !-Flag out faces that do not contribute to surface tension force 
    ALLOCATE (face_flag(nfc,ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: face_flag(nfc,ncells) allocation failed')
    face_flag=.false.

    do n=1,ncells
      do f=1,nfc
        if (abs(dC_dx(f,n)).gt.cutvof.or.abs(dC_dy(f,n)).gt.cutvof &
            .or.abs(dC_dz(f,n)).gt.cutvof) then
          face_flag(f,n)=.true.
        endif 
        if (.not.face_flag(f,n)) then
          Kappa_face(f,n)=0.0_r8
        endif
      enddo
    enddo

    DEALLOCATE (face_flag)

    !---------------------------------------------
    ! Estimate surface tension coefficient

    !call PROPERTY (Zone%Temp_Old, surfmat1, 'surface_tension', Value = Sigma)
    do j = 1, ncells
      state(1) = Zone(j)%Temp_Old
      Sigma(j) = sigma_func%eval(state)
    end do

    !-Interpolate surface tension coefficient at faces

    ! Allocate variables

    ALLOCATE (Sigma_vtx(nnodes),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: Sigma_vtx(nnodes) allocation failed')

    ALLOCATE (Sigma_v(nvc,ncells),STAT=status)
    if (status /=0) call TLS_panic ('CSF_FACE: Sigma_v(nvc,ncells) allocation failed')

    Sigma_face=0.0_r8

    !-averaged Sigma at faces
    ! first averaged at vertices
    ! then at faces (Sigma is a scalar)

     call VERTEX_AVG(Sigma,Sigma_vtx) 
     call EN_GATHER(Sigma_v,Sigma_vtx)
     do f=1,nfc
       call LINEAR_PROP(f,Sigma_v,Sigma_face(f,:))
     enddo

    DEALLOCATE(Sigma_vtx,Sigma_v)

    !-------------------------
    ! CSF normal component.

    NORMAL: if (csf_normal) then
       ! Construct the CSF normal component, nonzero only
       do f=1,nfc 
         Fx(f,:)=Fx(f,:)+Sigma_face(f,:)*Kappa_face(f,:)*dC_dx(f,:)
         Fy(f,:)=Fy(f,:)+Sigma_face(f,:)*Kappa_face(f,:)*dC_dy(f,:)
         Fz(f,:)=Fz(f,:)+Sigma_face(f,:)*Kappa_face(f,:)*dC_dz(f,:)
       enddo 
    end if NORMAL

  !--------------------------

    Fcsf_face(1,:,:)=Fx(:,:)
    Fcsf_face(2,:,:)=Fy(:,:)
    Fcsf_face(3,:,:)=Fz(:,:)
   
    ! Stop the surface tension timer.
    call stop_timer ("Projection Surface Tension")
 
    DEALLOCATE (Kappa)     
    DEALLOCATE (Sigma)
    DEALLOCATE (Color)
    DEALLOCATE (Kappa_face,Sigma_face)
    DEALLOCATE (dC_dx,dC_dy,dC_dz)
    DEALLOCATE (Weight) 
    DEALLOCATE (dC)
    DEALLOCATE (Fx,Fy,Fz) 
 
  END SUBROUTINE CSF_FACE

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE DIVERGENCE_INTNORMAL(D,NV)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the cell-centered divergence (D) of a cell-centered
    !   vector NV with a discrete approximation
    !   to Gauss theorem, whereby the integral of Div*NV over the
    !   cell volume is converted to an integral of NV over the cell
    !   surface area vector.  The area integral is approximated as a
    !   discrete sum over cell faces. This control volume formulation
    !   is discretely conservative, i.e., adjacent face contributions
    !   telescope, leaving only boundary contributions.
    !
    ! Note: this subroutine is based on the subroutine DIVERGENCE, 
    !       and has different boundary condition
    ! 
    !=======================================================================
    use cutoffs_module,    only: alittle
    use legacy_mesh_api,   only: Cell, Mesh, Vrtx_Face, Vertex
    use legacy_mesh_api,   only: ncells, ndim, nfc, nfv, nnodes, nvc
    use legacy_mesh_api,   only: EN_GATHER, EN_SUM_SCATTER, LINEAR_PROP

    ! Arguments
    real(r8), dimension(ndim,ncells),    intent(IN)  :: NV
    real(r8), dimension(ncells),         intent(OUT) :: D

    ! Local Variables
    logical, dimension(ncells)           :: Mask
    integer                              :: f, face, n, v, i
    real(r8), dimension(ncells)          :: Tmp, Normal
    real(r8), dimension(ndim,nfc,ncells) :: NV_Face
    real(r8), dimension(nvc,ncells)      :: NV_v
    real(r8), dimension(nnodes)          :: NV_vtx
     
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Loop over each velocity component, computing a face-centered
    ! velocity component from cell-centered velocities
    do n = 1,ndim
       ! Loop over each vertex, computing the portion of the
       ! cell-centered normal component apportioned to each vertex
       NV_v = 0.0_r8
       do v = 1,nvc
          Normal = NV(n,:)

          ! Loop over those faces associated with this vertex,
          ! and subtract out the face-normal part of the normal
          ! component for those faces that are boundary faces. This
          ! operation mocks up the presence of a ghost cell across
          ! from each boundary face that has an equal and opposite
          ! cell-centered normal. Note: nfv - number of faces/vertex,
          ! Vrtx_Face(v,f) = face number associated with vertex f
          ! This assume symmetric BC 
  
           do f = 1,nfv
             face = Vrtx_Face(v,f)

             ! Set the Mask to true on boundary faces, 
             ! except for those faces with a specified pressure,
             ! where flow can go in/out
             Mask = Mesh%Ngbr_cell(face) == 0

             ! Normal part of the velocity
             Tmp = 0.0_r8
             do i = 1,ndim
               Tmp = Tmp + NV(i,:)*Cell%Face_Normal(i,face)
             end do

             ! Subtract out the normal part
             where (Mask)
                Normal = Normal - Tmp*Cell%Face_Normal(n,face)
             end where
          end do

          ! Weight the contribution of the normal to this
          ! vertex with the inverse cell volume
          NV_v(v,:) = Normal/Cell%Volume
       end do

       ! Sum-Scatter the velocity vector to each vertex
       call EN_SUM_SCATTER (NV_vtx, NV_v)

       ! Multiply by the inverse volume sum
       NV_Vtx = NV_Vtx * Vertex%Rsum_Rvol

       ! Gather V_Vtx into cell vector NV_v
       call EN_GATHER (NV_v, NV_Vtx)

       ! Loop over faces, computing a face-centered
       ! vector from the vertex-averaged velocities
       do f = 1,nfc
          ! Interpolate vertex values to this face
          call LINEAR_PROP (f, NV_v, NV_face(n,f,:))
       end do
    end do

    ! normalize NV_face
    do f=1,nfc 
      Tmp=sqrt(NV_face(1,f,:)**2+NV_face(2,f,:)**2+NV_face(3,f,:)**2)        
      Tmp=MERGE(0.0_r8,1/Tmp,Tmp<=alittle)
      NV_face(1,f,:)=Tmp*NV_face(1,f,:)
      NV_face(2,f,:)=Tmp*NV_face(2,f,:) 
      NV_face(3,f,:)=Tmp*NV_face(3,f,:)
    enddo
  
    ! Loop over faces and accumulate the divergence
    ! from these face-centered velocities
    D = 0.0_r8
    do f = 1,nfc
       Tmp = 0.0_r8
       ! Apply BC if wall adhesion
 
       ! Dot product of face interface normal and unit normals
       do n = 1,ndim
         Tmp = Tmp + NV_Face(n,f,:)*Cell%Face_Normal(n,f)
       end do

       ! Multiply by area
       D = D + Tmp*Cell%Face_Area(f)
    end do

    ! Normalize by cell volume
    D = D/Cell%Volume

    ! Eliminate noise
    D = MERGE(0.0_r8, D, ABS(D) <= alittle)

  END SUBROUTINE DIVERGENCE_INTNORMAL

  subroutine read_surface_tension_namelist (lun)

    use input_utilities
    use string_utilities, only: i_to_c, raise_case
    use parallel_communication, only: is_IOP, broadcast
    use scalar_func_factories, only: alloc_const_scalar_func
    use scalar_func_table, only: lookup_func
    use property_module, only: get_truchas_material_id
    use fluid_data_module, only: IsImmobile

    integer, intent(in) :: lun

    integer :: ios, n
    logical :: found

    integer :: bndry_face_set_ids(MAX_FACE_SET_IDS)
    integer  :: interface_materials(2)
    real(r8) :: sigma_constant
    character(32) :: smoothing_kernel, sigma_function
    character(len=8+MAX_NAME_LEN) :: label
    namelist /surface_tension/ csf_normal, csf_tangential, smoothing_kernel, &
        interface_materials, sigma_constant, sigma_function, dsig_dT, &
        csf_boundary, bndry_face_set_ids
    
    call TLS_info ('')
    call TLS_info ('Reading SURFACE_TENSION namelist ...')

    !! Locate the SURFACE_TENSION namelist (first occurence).
    if (is_IOP) then
      rewind lun
      call seek_to_namelist (lun, 'SURFACE_TENSION', found, iostat=ios)
    end if

    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('Error reading input file: iostat=' // i_to_c(ios))

    call broadcast (found)
    if (.not.found) call TLS_fatal ('SURFACE_TENSION namelist not found')

    !! Read the namelist.
    if (is_IOP) then
      csf_normal = .false.
      csf_tangential = .false.
      csf_boundary = .false.
      smoothing_kernel = NULL_C
      interface_materials = NULL_I
      sigma_constant = NULL_R
      sigma_function = NULL_C
      dsig_dT = 0.0
      bndry_face_set_ids = NULL_I
      read(lun,nml=surface_tension,iostat=ios)
    end if

    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('Error reading SURFACE_TENSION namelist.')

    !! Broadcast the namelist variables.
    call broadcast (csf_normal)
    call broadcast (csf_tangential)
    call broadcast (csf_boundary)
    call broadcast (smoothing_kernel)
    call broadcast (interface_materials)
    call broadcast (sigma_constant)
    call broadcast (sigma_function)
    call broadcast (dsig_dT)
    call broadcast (bndry_face_set_ids)

    if (.not.(csf_normal .or. csf_tangential .or. csf_boundary)) then
      call TLS_fatal ('At least one of csf_normal, csf_tangential, or csf_boundary  must be enabled.')
    endif

    if (csf_boundary .and. csf_normal) call TLS_fatal ('csf_boundary and csf_normal are mutually exclusive options')
if (csf_boundary .and. csf_tangential) call TLS_fatal ('csf_boundary and csf_tangential are mutually exclusive options')

    if (csf_normal) then
      smoothing_kernel = raise_case(smoothing_kernel)
      select case (smoothing_kernel)
      case ('WILLIAMS')
        kernel = WILLIAMS_KERNEL
      case ('RUDMAN')
        kernel = RUDMAN_KERNEL
      case (NULL_C)
        !call TLS_fatal ('SMOOTHING_KERNEL must be assigned a value when CSF_NORMAL is enabled.')
        kernel = RUDMAN_KERNEL
        call TLS_info ('  using default SMOOTHING_KERNEL value: "RUDMAN"')
      case default
        call TLS_fatal ('Unknown value for SMOOTHING_KERNEL: "' // trim(smoothing_kernel) // '"')
      end select
    end if
    
    if (csf_boundary) then

      !! Check for a non-empty FACE_SET_IDS.
      if (count(bndry_face_set_ids /= NULL_I) == 0) then
        call TLS_fatal('no values assigned to FACE_SET_IDS')
      endif

      allocate(face_set_ids(count(bndry_face_set_ids/=NULL_I)))
      face_set_ids = pack(bndry_face_set_ids, mask=(bndry_face_set_ids/=NULL_I))

    else

      !! Verify that the interface material numbers refer to distinct defined fluids
      if (any(interface_materials == NULL_I)) then
        call TLS_fatal ('INTERFACE_MATERIALS must be assigned two values.')
      else
        surfmat1 = get_truchas_material_id(interface_materials(1))
        if (surfmat1 < 1) call TLS_fatal ('Unknown material number INTERFACE_MATERIALS(1): ' &
                                          // i_to_c(interface_materials(1)))
        surfmat2 = get_truchas_material_id(interface_materials(2))
        if (surfmat2 < 1) call TLS_fatal ('Unknown material number INTERFACE_MATERIALS(2): ' &
                                          // i_to_c(interface_materials(2)))
        if (surfmat1 == surfmat2) call TLS_fatal ('INTERFACE_MATERIALS must be assigned distinct material numbers')
        if (IsImmobile(surfmat1)) call TLS_fatal ('INTERFACE_MATERIALS(1) is not a fluid')
        if (IsImmobile(surfmat2)) call TLS_fatal ('INTERFACE_MATERIALS(2) is not a fluid')
      end if

    end if
    
    !! Verify that only one of SIGMA_CONSTANT and SIGMA_FUNCTION was specified.
    if (sigma_constant == NULL_R .eqv. sigma_function == NULL_C) then
      call TLS_fatal ('Exactly one of SIGMA_CONSTANT and SIGMA_FUNCTION must be assigned a value.')
    end if
    
    !! Get or create the specified surface tension coefficient function.
    if (sigma_constant /= NULL_R) then
      if (sigma_constant < 0.0_r8) then
        call TLS_fatal ('SIGMA_CONSTANT must be >= 0')
      end if
      call alloc_const_scalar_func (sigma_func, sigma_constant)
    else
      call lookup_func (sigma_function, sigma_func)
      if (.not.allocated(sigma_func)) then
        call TLS_fatal ('Unknown SIGMA_FUNCTION name: "' // trim(sigma_function) // '"')
      end if
    end if

  end subroutine read_surface_tension_namelist

end module surface_tension_module
