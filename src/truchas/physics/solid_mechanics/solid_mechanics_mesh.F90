!!
!! Solid Mechanics Mesh Module
!!      
!! This module holds a copy of the mesh pointer and defines subroutines
!! and functions related to the mesh required throughout the solid mechanics physics
!! package.
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PUBLIC Data 
!!
!!    SM_MESH
!!    A pointer to solid mechanics dist_mesh 
!!
!!    SM_MESH_TO_OLD_CELL, SM_OLD_TO_MESH_CELL
!!    SM_MESH_TO_OLD_NODE, SM_OLD_TO_MESH_NODE
!!    Cell and node based data structures that permute the old mesh to dist_mesh
!!
!!    SM_GAP_ELEMENT
!!    Integer pointer. I have no idea what this is. Indices map to gap elements
!!    old mesh? I was mimicking what I saw in the heat transfer module.
!!
!! PUBLIC Parameters
!!
!!    These are duplicate parameters that originally came from the
!!    mesh module. Preserved here to maintain the package capability to
!!    handle meshes that are not hex meshes. 
!!    SM_CELL_TET
!!    SM_CELL_PYRAMID
!!    SM_CELL_PRISM
!!    SM_CELL_HEX
!!    SM_GAP_ELEMENT_1      
!!    SM_GAP_ELEMENT_3      
!!    SM_GAP_ELEMENT_5
!!
!!    These are duplicate parameters from the parameter module. Solid mechanics
!!    modules should these parameters instead of the PARAMETER_MODULE parameters
!!    These should be removed one the change to DIST_MESH is complete.
!!    NDIM
!!    NFC
!!    NVC
!!    NVF
!!    NEC
!!    NROT
!!    NCOMPS
!!
!! PUBLIC Subroutine/Functions
!!
!!    SM_MESH_ENABLE
!!    A wrapper call to DIST_MESH's enable_mesh. Only used when initializing
!! 
!!    SM_MESH_INIT
!!    Initialize (copy) SM_MESH pointer and other public data 
!!
!!    SM_CELL_SHAPE
!!    A function that return cell shape parameter. Only checks the number of
!!    nodes of the first cell element. THIS WILL FAIL with mixed element meshes
!!
!!    SM_FACE_CENTROID
!!    Compute the physical face centroid. Based on FACE_CENTROID_PHYSICAL in the cell geometry
!!    module.
!!
!!    SM_CELL_CENTROID
!!    Compute the physical cell centroid. Based on CELL_CENTROID_PHYSICAL in the cell geometry
!!    module.
!!
!!    SM_CELL_CENTROID_LOGICAL
!!    Based on CELL_CENTROID_LOGICAL in the cell geometry module.

#include "f90_assert.fpp"
      

module solid_mechanics_mesh

  use unstr_mesh_type
  use parallel_permutations, only: par_perm

  implicit none
  private

  !! Solid mechanics DIST_MESH Pointer
  type(unstr_mesh), pointer, public :: sm_mesh => null()

  !! Permutation structures to transfer data to and from old mesh
  type(par_perm), target, save, public :: sm_mesh_to_old_cell, sm_old_to_mesh_cell
  type(par_perm), target, save, public :: sm_mesh_to_old_node, sm_old_to_mesh_node

  !! I don't know if this should be a module variable or local in sm_mesh_init
  integer, pointer, public, save :: sm_gap_elements(:) => null()

  !!
  public sm_mesh_init, sm_mesh_enable!, sm_cell_shape

!NNC: not yet used -- disabling
!NNC  !! Imported Cell geometry subroutines
!NNC  public :: sm_face_centroid, sm_cell_centroid, sm_face_centroid_logical

  !! From the old mesh module a data block defining the edge ordering of a HEX
  !! element.
  integer, public :: Cell_Edge(2,12)
  data Cell_Edge/1,2, 2,3, 3,4, 4,1, 2,6, 3,7, 4,8, 1,5, 5,6, 6,7, 7,8, 8,5/

  
  !! Parameter definitions originally from the old mesh module
  integer, parameter, public :: SM_CELL_TET     = 1
  integer, parameter, public :: SM_CELL_PYRAMID = 2
  integer, parameter, public :: SM_CELL_PRISM   = 3
  integer, parameter, public :: SM_CELL_HEX     = 4
  !! define another cell shape for "zero volume" gap elements
  !! The gap element number corresponds to the first non-degenerate cell face
  integer, parameter, public :: SM_GAP_ELEMENT_1  = 5
  integer, parameter, public :: SM_GAP_ELEMENT_3  = 6
  integer, parameter, public :: SM_GAP_ELEMENT_5  = 7

  !! These are parameters orignally from the parameter module
  !! Eventually these will be removed and replaced with 
  !! Number of physical dimensions
  integer, parameter, public :: ndim = 3

  !! Number of faces per cell
  integer, parameter, public :: nfc = 2*ndim

  !! Number of vertices per cell
  integer, parameter, public :: nvc = 2**ndim

  !! Number of vertices per face
  integer, parameter, public :: nvf = 2**(ndim - 1)

  !! Number of edges per cell
  integer, parameter, public :: nec = ndim*2**(ndim - 1)

  !! Number of rotation axes
  integer, parameter, public :: nrot = 2*ndim - 3

  !! Number of stress/strain components
  integer, parameter, public :: ncomps = (ndim*(ndim + 1))/2
  
  !! Namelist mesh name, used in sm_mesh_init and sm_mesh_enable
  character(4) :: sm_mesh_name = 'MAIN'


contains

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 

  ! Purpose: Initialize the solid mechanics mesh pointer and the permutation
  !          structures.
  subroutine sm_mesh_init
    use mesh_manager, only: unstr_mesh_ptr
    use mesh_module, only: unpermute_mesh_vector, unpermute_vertex_vector
    use parallel_permutations, only: create_par_perm

    integer, pointer :: dummy(:) => null()
    logical :: found

    ! Initialize the mesh pointer
    sm_mesh => unstr_mesh_ptr(sm_mesh_name)
    INSIST(associated(sm_mesh))

    ! Initialize the permutation structures -- I've copied
    ! the generate_mesh_mappings routine from  mesh_interop.F90
    call create_par_perm(unpermute_mesh_vector, &
                         sm_mesh%xcell(:sm_mesh%ncell_onp), &
                         sm_old_to_mesh_cell, sm_gap_elements, &
                         sm_mesh_to_old_cell, dummy)
    INSIST(size(dummy) == 0)
    deallocate(dummy)
    INSIST(are_gap_elements(sm_gap_elements))

    call create_par_perm(unpermute_vertex_vector, &
                         sm_mesh%xnode(:sm_mesh%nnode_onp), &
                         sm_old_to_mesh_node, sm_gap_elements, &
                         sm_mesh_to_old_node,dummy)
    INSIST(size(dummy) == 0)
    deallocate(dummy)
    INSIST(are_gap_elements(sm_gap_elements))

    contains
      logical function are_gap_elements (list)
      use mesh_module, only: mesh, GAP_ELEMENT_1
      integer, intent(in) :: list(:)
      integer :: j
      are_gap_elements = .false.
      do j = 1, size(list)
        if (mesh(list(j))%cell_shape < GAP_ELEMENT_1) return
      end do
      are_gap_elements = .true.
    end function are_gap_elements
  end subroutine sm_mesh_init

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Purpose: A wrapper to the DIST_MESH call to enable the mesh pointer
  subroutine sm_mesh_enable
    use mesh_manager, only : enable_mesh
    logical :: found
    call enable_mesh(sm_mesh_name,found)
    INSIST(found)
  end subroutine sm_mesh_enable
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  
!NNC: not yet used -- disabling
!NNC  function sm_cell_shape (mesh) result (cell_shape)
!NNC    type(unstr_mesh), intent(in) :: mesh
!NNC    integer :: cell_shape
!NNC
!NNC    select case ( size(mesh%cnode,dim=1) )
!NNC      case (4)
!NNC        cell_shape=SM_CELL_TET
!NNC      case (8)
!NNC        cell_shape=SM_CELL_HEX
!NNC      case default
!NNC        cell_shape=SM_CELL_HEX
!NNC    end select
!NNC
!NNC  end function sm_cell_shape
!NNC  
!NNC  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!NNC
!NNC  subroutine sm_face_centroid(mesh,centroids)
!NNC    use kinds,         only: r8
!NNC    use linear_module, only: linear_prop
!NNC
!NNC    type(dist_mesh), intent(in) :: mesh
!NNC    real(r8), intent(out) :: centroids(ndim,nfc,mesh%ncell)
!NNC
!NNC    integer :: j, k, inode 
!NNC    integer :: f, n
!NNC    real(r8), allocatable :: coord(:,:,:)
!NNC
!NNC    allocate(coord(ndim,size(mesh%cnode,dim=1),mesh%ncell))
!NNC
!NNC    !! Replace the EN_GATHER call in the NDIM_LOOP
!NNC    do j = 1, mesh%ncell
!NNC      do k = 1, size(mesh%cnode,dim=1)
!NNC        inode = mesh%cnode(k,j)
!NNC        coord(:,k,j) = mesh%x(:,inode)
!NNC      end do  
!NNC    end do
!NNC
!NNC    ! Now loop through the dimensions and faces
!NNC    do n = 1, ndim
!NNC      do f = 1,nfc
!NNC        call linear_prop(f,coord(n,:,:),centroids(n,f,:))
!NNC      end do 
!NNC    end do 
!NNC
!NNC    ! Deallocate
!NNC    deallocate(coord)
!NNC
!NNC   
!NNC  end subroutine sm_face_centroid  
!NNC  
!NNC  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!NNC
!NNC  subroutine sm_face_centroid_logical(mesh,centroids)
!NNC    !=======================================================================
!NNC    ! Purpose(s):
!NNC    !
!NNC    !   Replacement for Cell%Face_Centroid_L. Lifted from face_centroid_logical
!NNC    !   in cell_geometry_module.F90
!NNC    !
!NNC    !   Compute the face centroid (Face_centroid), which is part
!NNC    !   of the Cell structure. A face centroid is assigned for
!NNC    !   each face of the cell, and is a LOGICAL coordinate.
!NNC    !     Input  - Vertex%Coord(ndim)
!NNC    !     Output - Cell%Face_centroid_L(ndim,nfc)
!NNC    !
!NNC    !=======================================================================
!NNC    use kinds,                only: r8
!NNC    use ArrayAllocate_Module, only: ARRAYCREATE, ARRAYDESTROY
!NNC    use cutoffs_module,       only: alittle
!NNC
!NNC    ! Input/Output
!NNC    type(dist_mesh), intent(in)  :: mesh
!NNC    real(r8),        intent(out) :: centroids(ndim,nfc,mesh%ncell)
!NNC
!NNC    ! Local Variables
!NNC    integer :: j, k, inode, iface
!NNC    integer :: f, i, i1, i2, n, v11, v12, v13, v14,    &
!NNC                                v21, v22, v23, v24, v31, v32, v33, v34, &
!NNC                                nn
!NNC    integer :: n1 = 1, n2 = 2, n3 = 3
!NNC    real(r8), dimension(ndim)   :: Coef
!NNC    !real(r8), dimension(ncells) :: Face_Area
!NNC    real(r8), allocatable :: Face_Area(:,:)
!NNC    real(r8), allocatable :: Face_Normal(:,:,:)
!NNC    real(r8), pointer, dimension(:,:,:) :: Xv
!NNC
!NNC    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!NNC
!NNC    ! Only do this if it is 3-D.
!NNC    if (ndim == 3) then
!NNC       ! Allocate temporaries
!NNC       !call ARRAYCREATE (Xv, 1, ndim, 1, nvc, 1, ncells, 'Array Xv(ndim,nvc,ncells)')
!NNC       call ARRAYCREATE (Xv, 1, ndim, 1, nvc, 1, mesh%ncell, 'Array Xv(ndim,nvc,ncells)')
!NNC
!NNC       ! Gather vertex coordinates
!NNC       !do i = 1,ndim
!NNC       !   call EN_GATHER (Xv(i,:,:), Vertex%Coord(i), BOUNDARY=Vrtx_Bdy(i)%Data)
!NNC       !end do
!NNC       do j = 1, mesh%ncell
!NNC         do k = 1, size(mesh%cnode,dim=1)
!NNC           inode = mesh%cnode(k,j)
!NNC           Xv(:,k,j) = mesh%x(:,inode)
!NNC         end do  
!NNC       end do  
!NNC    end if
!NNC
!NNC    ! Compute the Face_Area and the Face_Normal arrays
!NNC    allocate(Face_Area(nfc,mesh%ncell))
!NNC    allocate(Face_Normal(ndim,nfc,mesh%ncell))
!NNC    Face_Area = 0.0_r8
!NNC    do j = 1, mesh%ncell
!NNC      do f = 1, size(mesh%cface,dim=1)
!NNC        iface = mesh%cface(f,j)
!NNC        Face_Area(f,j) = mesh%area(iface)
!NNC        Face_Normal(:,f,j) = mesh%normal(:,iface)
!NNC      end do 
!NNC    end do    
!NNC
!NNC    ! Loop over faces
!NNC    FACE_LOOP: do f = 1,nfc
!NNC
!NNC       ! Preinitialize the face centroid
!NNC       do n = 1,ndim
!NNC          !Cell%Face_Centroid_L(n,f) = 0.5_r8
!NNC          CENTROIDS(n,f,:)          = 0.5_r8
!NNC          Coef(n)                   = 1
!NNC       end do
!NNC
!NNC       ! We are done if this is 2-D.
!NNC       if (ndim == 2) cycle FACE_LOOP
!NNC
!NNC       ! Defined outside the loop  
!NNC       !Face_Area = 1.0 / (12*(Cell%Face_Area(f) + alittle))
!NNC
!NNC       select case (f)
!NNC
!NNC          case (1)
!NNC
!NNC          ! Side 1 (vertices 4-8-7-3)
!NNC          v11 = 1; v12 = 1; v13 = 1; v14 = 1
!NNC          v21 = 7; v22 = 8; v23 = 3; v24 = 4
!NNC          v31 = 8; v32 = 4; v33 = 7; v34 = 3
!NNC          nn  = 1
!NNC
!NNC          !Cell%Face_centroid_L(nn,f) = 0; Coef(nn) = 0
!NNC          CENTROIDS(nn,f,:) = 0; Coef(nn) = 0
!NNC
!NNC          case (2)
!NNC
!NNC          ! Side 2 (vertices 1-2-6-5)
!NNC          v11 = 1; v12 = 1; v13 = 1; v14 = 1
!NNC          v21 = 2; v22 = 1; v23 = 6; v24 = 5
!NNC          v31 = 6; v32 = 2; v33 = 5; v34 = 1
!NNC          nn  = 1
!NNC
!NNC          !Cell%Face_centroid_L(nn,f) = 1; Coef(nn) = 0
!NNC          CENTROIDS(nn,f,:) = 1; Coef(nn) = 0
!NNC
!NNC          case (3)
!NNC
!NNC          ! Side 3 (vertices 4-1-5-8)
!NNC          v11 = 1; v12 = 4; v13 = 5; v14 = 8
!NNC          v21 = 1; v22 = 1; v23 = 1; v24 = 1
!NNC          v31 = 5; v32 = 1; v33 = 8; v34 = 4
!NNC          nn  = 2
!NNC
!NNC          !Cell%Face_centroid_L(nn,f) = 0; Coef(nn) = 0
!NNC          CENTROIDS(nn,f,:) = 0; Coef(nn) = 0
!NNC
!NNC          case (4)
!NNC
!NNC          ! Side 4 (vertices 3-7-6-2)
!NNC          v11 = 6; v12 = 7; v13 = 2; v14 = 3
!NNC          v21 = 1; v22 = 1; v23 = 1; v24 = 1
!NNC          v31 = 7; v32 = 3; v33 = 6; v34 = 2
!NNC          nn  = 2
!NNC
!NNC          !Cell%Face_centroid_L(nn,f) = 1; Coef(nn) = 0
!NNC          CENTROIDS(nn,f,:) = 1; Coef(nn) = 0
!NNC
!NNC          case (5)
!NNC
!NNC          ! Side 5 (vertices 4-3-2-1)
!NNC          v11 = 2; v12 = 3; v13 = 1; v14 = 4
!NNC          v21 = 3; v22 = 4; v23 = 2; v24 = 1
!NNC          v31 = 1; v32 = 1; v33 = 1; v34 = 1
!NNC          nn  = 3
!NNC
!NNC          !Cell%Face_centroid_L(nn,f) = 0; Coef(nn) = 0
!NNC          CENTROIDS(nn,f,:) = 0; Coef(nn) = 0
!NNC
!NNC          case (6)
!NNC
!NNC          ! Side 6 (vertices 8-5-6-7)
!NNC          v11 = 5; v12 = 8; v13 = 6; v14 = 7
!NNC          v21 = 6; v22 = 5; v23 = 7; v24 = 8
!NNC          v31 = 1; v32 = 1; v33 = 1; v34 = 1
!NNC          nn  = 3
!NNC
!NNC          !Cell%Face_centroid_L(nn,f) = 1; Coef(nn) = 0
!NNC          CENTROIDS(nn,f,:) = 1; Coef(nn) = 0
!NNC
!NNC       end select
!NNC
!NNC       NDIM_LOOP: do i = 1,ndim
!NNC
!NNC          select case (i)
!NNC             case (1)
!NNC                i1 = 2; i2 = 3
!NNC             case (2)
!NNC                i1 = 3; i2 = 1
!NNC             case (3)
!NNC                i1 = 1; i2 = 2
!NNC          end select
!NNC
!NNC          !Cell%Face_centroid_L(n1,f) = Cell%Face_centroid_L(n1,f) + &
!NNC          !      Coef(n1)*Face_Area*Cell%Face_Normal(i,f)* &
!NNC          !      ((Xv(i1,v11,:) - Xv(i1,v12,:))*(Xv(i2,v13,:) - Xv(i2,v14,:)) &
!NNC          !     - (Xv(i2,v11,:) - Xv(i2,v12,:))*(Xv(i1,v13,:) - Xv(i1,v14,:)))
!NNC          CENTROIDS(n1,f,:) = CENTROIDS(n1,f,:) + &
!NNC                Coef(n1)*Face_Area(f,:)*Face_Normal(i,f,:)* &
!NNC                ((Xv(i1,v11,:) - Xv(i1,v12,:))*(Xv(i2,v13,:) - Xv(i2,v14,:)) &
!NNC               - (Xv(i2,v11,:) - Xv(i2,v12,:))*(Xv(i1,v13,:) - Xv(i1,v14,:)))
!NNC
!NNC
!NNC          !Cell%Face_centroid_L(n2,f) = Cell%Face_centroid_L(n2,f) + &
!NNC          !     Coef(n2)*Face_Area*Cell%Face_Normal(i,f)* &
!NNC          !     ((Xv(i1,v21,:) - Xv(i1,v22,:))*(Xv(i2,v23,:) - Xv(i2,v24,:)) &
!NNC          !    - (Xv(i2,v21,:) - Xv(i2,v22,:))*(Xv(i1,v23,:) - Xv(i1,v24,:)))
!NNC          CENTROIDS(n2,f,:) = CENTROIDS(n2,f,:) + &
!NNC               Coef(n2)*Face_Area(f,:)*Face_Normal(i,f,:)* &
!NNC               ((Xv(i1,v21,:) - Xv(i1,v22,:))*(Xv(i2,v23,:) - Xv(i2,v24,:)) &
!NNC              - (Xv(i2,v21,:) - Xv(i2,v22,:))*(Xv(i1,v23,:) - Xv(i1,v24,:)))
!NNC
!NNC          !Cell%Face_centroid_L(n3,f) = Cell%Face_centroid_L(n3,f) + &
!NNC          !     Coef(n3)*Face_Area*Cell%Face_Normal(i,f)* &
!NNC          !     ((Xv(i1,v31,:) - Xv(i1,v32,:))*(Xv(i2,v33,:) - Xv(i2,v34,:)) &
!NNC          !    - (Xv(i2,v31,:) - Xv(i2,v32,:))*(Xv(i1,v33,:) - Xv(i1,v34,:)))
!NNC          CENTROIDS(n3,f,:) = CENTROIDS(n3,f,:) + &
!NNC               Coef(n3)*Face_Area(f,:)*Face_Normal(i,f,:)* &
!NNC               ((Xv(i1,v31,:) - Xv(i1,v32,:))*(Xv(i2,v33,:) - Xv(i2,v34,:)) &
!NNC              - (Xv(i2,v31,:) - Xv(i2,v32,:))*(Xv(i1,v33,:) - Xv(i1,v34,:)))
!NNC
!NNC       end do NDIM_LOOP
!NNC
!NNC    end do FACE_LOOP
!NNC
!NNC    ! Deallocate temporaries
!NNC    if (ASSOCIATED(Xv)) call ARRAYDESTROY (Xv, 'Array Xv(ndim,nvc,ncells)')
!NNC
!NNC    ! Deallocate the face and normal temporary arrays
!NNC    deallocate(Face_Area)
!NNC    deallocate(Face_Normal)
!NNC
!NNC
!NNC
!NNC  end subroutine sm_face_centroid_logical
!NNC
!NNC  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!NNC
!NNC  subroutine sm_cell_centroid(mesh,centroids)
!NNC    use kinds,                only : r8
!NNC    use ArrayAllocate_Module, only : ARRAYCREATE, ARRAYDESTROY
!NNC    !=======================================================================
!NNC    ! Purpose(s):
!NNC    !
!NNC    !   Compute the cell centroids, given a distributed mesh type
!NNC    !     Originally defined in cell_geometry_module.F90 in subroutine
!NNC    !     cell_centroid
!NNC    !
!NNC    !=======================================================================
!NNC
!NNC    ! Input/Output
!NNC    type(dist_mesh), intent(in)  :: mesh
!NNC    real(r8),        intent(out) :: centroids(ndim,mesh%ncell)
!NNC
!NNC
!NNC
!NNC    ! Local Variables
!NNC    integer :: j, k, inode, ncells
!NNC    integer :: i, i1, i2, v1, v2, v3, v4, v5, v6, v7, v8
!NNC    real(r8), pointer, dimension(:,:) :: L, M, N, LxD3, MxD2, NxD1
!NNC    real(r8), pointer, dimension(:,:) :: Tmp, D1, D2, D3, Dv, D1xDv, D2xDv, D3xDv
!NNC    real(r8), pointer, dimension(:,:,:) :: Xv
!NNC
!NNC    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!NNC
!NNC    ! Define the mesh parameters
!NNC    ncells = mesh%ncell
!NNC
!NNC    ! Explicitly allocate temporaries
!NNC    call ARRAYCREATE (L, 1, ndim, 1, ncells, 'Array L(ndim,ncells)')
!NNC    call ARRAYCREATE (M, 1, ndim, 1, ncells, 'Array M(ndim,ncells)')
!NNC    call ARRAYCREATE (N, 1, ndim, 1, ncells, 'Array N(ndim,ncells)')
!NNC    call ARRAYCREATE (Xv, 1, ndim, 1, nvc, 1, ncells, 'Array Xv(ndim,nvc,ncells)')
!NNC    call ARRAYCREATE (LxD3, 1, nrot, 1, ncells, 'Array LxD3(nrot,ncells)')
!NNC    call ARRAYCREATE (MxD2, 1, nrot, 1, ncells, 'Array MxD2(nrot,ncells)')
!NNC    call ARRAYCREATE (NxD1, 1, nrot, 1, ncells, 'Array NxD1(nrot,ncells)')
!NNC
!NNC    ! Only need these for 3-D.
!NNC    if (ndim == 3) then
!NNC       call ARRAYCREATE (Tmp, 1, ndim, 1, ncells, 'Array Tmp(ndim,ncells)')
!NNC       call ARRAYCREATE (D1, 1, ndim, 1, ncells, 'Array D1(ndim,ncells)')
!NNC       call ARRAYCREATE (D2, 1, ndim, 1, ncells, 'Array D2(ndim,ncells)')
!NNC       call ARRAYCREATE (D3, 1, ndim, 1, ncells, 'Array D3(ndim,ncells)')
!NNC       call ARRAYCREATE (Dv, 1, ndim, 1, ncells, 'Array Dv(ndim,ncells)')
!NNC       call ARRAYCREATE (D1xDv, 1, nrot, 1, ncells, 'Array D1xDv(nrot,ncells)')
!NNC       call ARRAYCREATE (D2xDv, 1, nrot, 1, ncells, 'Array D2xDv(nrot,ncells)')
!NNC       call ARRAYCREATE (D3xDv, 1, nrot, 1, ncells, 'Array D3xDv(nrot,ncells)')
!NNC    end if
!NNC
!NNC    ! Set vertex numbers
!NNC    select case(ndim)
!NNC       case (2)
!NNC          v1 = 1; v2 = 2; v3 = 3; v4 = 4
!NNC       case (3)
!NNC          v1 = 1; v2 = 2; v3 = 3; v4 = 4
!NNC          v5 = 5; v6 = 6; v7 = 7; v8 = 8
!NNC    end select
!NNC
!NNC    ! Gather vertex coordinates
!NNC    !do i = 1,ndim
!NNC    !   call EN_GATHER (Xv(i,:,:), Vertex%Coord(i), BOUNDARY=Vrtx_Bdy(i)%Data)
!NNC    !end do
!NNC    do j = 1, ncells
!NNC      do k = 1, nvc
!NNC        inode = SM_MESH%CNODE(k,j)
!NNC        Xv(:,k,j) = SM_MESH%X(:,inode)
!NNC      end do  
!NNC    end do  
!NNC
!NNC    ! Compute quantities needed for the centroid coordinates
!NNC    do i = 1,ndim
!NNC
!NNC       select case (ndim)
!NNC
!NNC          case (2)
!NNC
!NNC             L(i,:) = Xv(i,v1,:) - Xv(i,v4,:)
!NNC             M(i,:) = Xv(i,v3,:) - Xv(i,v4,:)
!NNC             N(i,:) = - Xv(i,v1,:) + Xv(i,v2,:) - Xv(i,v3,:) + Xv(i,v4,:)
!NNC
!NNC          case (3)
!NNC
!NNC             L(i,:) = Xv(i,v1,:) + Xv(i,v2,:) + Xv(i,v6,:) + Xv(i,v5,:) &
!NNC                    - Xv(i,v3,:) - Xv(i,v4,:) - Xv(i,v8,:) - Xv(i,v7,:)
!NNC             M(i,:) = Xv(i,v2,:) + Xv(i,v3,:) + Xv(i,v7,:) + Xv(i,v6,:) &
!NNC                    - Xv(i,v4,:) - Xv(i,v1,:) - Xv(i,v5,:) - Xv(i,v8,:)
!NNC             N(i,:) = Xv(i,v8,:) + Xv(i,v5,:) + Xv(i,v6,:) + Xv(i,v7,:) &
!NNC                    - Xv(i,v3,:) - Xv(i,v2,:) - Xv(i,v1,:) - Xv(i,v4,:)
!NNC             D1(i,:) = Xv(i,v2,:) + Xv(i,v6,:) + Xv(i,v4,:) + Xv(i,v8,:) &
!NNC                     - Xv(i,v1,:) - Xv(i,v5,:) - Xv(i,v3,:) - Xv(i,v7,:)
!NNC             D2(i,:) = Xv(i,v3,:) + Xv(i,v4,:) + Xv(i,v5,:) + Xv(i,v6,:) &
!NNC                     - Xv(i,v1,:) - Xv(i,v2,:) - Xv(i,v7,:) - Xv(i,v8,:)
!NNC             D3(i,:) = Xv(i,v1,:) + Xv(i,v4,:) + Xv(i,v6,:) + Xv(i,v7,:) &
!NNC                     - Xv(i,v2,:) - Xv(i,v3,:) - Xv(i,v5,:) - Xv(i,v8,:)
!NNC             Dv(i,:) = Xv(i,v1,:) + Xv(i,v3,:) + Xv(i,v6,:) + Xv(i,v8,:) &
!NNC                     - Xv(i,v2,:) - Xv(i,v4,:) - Xv(i,v5,:) - Xv(i,v7,:)
!NNC
!NNC       end select
!NNC
!NNC    end do
!NNC
!NNC    select case (ndim)
!NNC
!NNC       case (2)
!NNC
!NNC          i1 = 1; i2 = 2 
!NNC          do i = 1,nrot
!NNC             LxD3(i,:) = L(i1,:)*M(i2,:) - L(i2,:)*M(i1,:)
!NNC             MxD2(i,:) = L(i1,:)*N(i2,:) - L(i2,:)*N(i1,:)
!NNC             NxD1(i,:) = N(i1,:)*M(i2,:) - N(i2,:)*M(i1,:)
!NNC          end do
!NNC
!NNC       case (3)
!NNC
!NNC          L = 0.25_r8*L
!NNC          M = 0.25_r8*M
!NNC          N = 0.25_r8*N
!NNC          Tmp = 0
!NNC
!NNC          do i = 1,nrot
!NNC
!NNC             select case (i)
!NNC             case (1)
!NNC                i1 = 2; i2 = 3
!NNC             case (2)
!NNC                i1 = 3; i2 = 1
!NNC             case (3)
!NNC                i1 = 1; i2 = 2
!NNC             end select
!NNC
!NNC             LxD3(i,:) = L(i1,:)*D3(i2,:) - L(i2,:)*D3(i1,:)
!NNC             MxD2(i,:) = M(i1,:)*D2(i2,:) - M(i2,:)*D2(i1,:)
!NNC             NxD1(i,:) = N(i1,:)*D1(i2,:) - N(i2,:)*D1(i1,:)
!NNC
!NNC             D1xDv(i,:) = D1(i1,:)*Dv(i2,:) - D1(i2,:)*Dv(i1,:)
!NNC             D2xDv(i,:) = D2(i1,:)*Dv(i2,:) - D2(i2,:)*Dv(i1,:)
!NNC             D3xDv(i,:) = D3(i1,:)*Dv(i2,:) - D3(i2,:)*Dv(i1,:)
!NNC
!NNC          end do
!NNC
!NNC          do i = 1,ndim
!NNC
!NNC             Tmp(v1,:) = Tmp(v1,:) + L(i,:)*(MxD2(i,:) - NxD1(i,:)) &
!NNC                        + (N(i,:)*D2xDv(i,:) - M(i,:)*D1xDv(i,:))/12
!NNC             Tmp(v2,:) = Tmp(v2,:) + M(i,:)*(NxD1(i,:) - LxD3(i,:)) &
!NNC                        + (L(i,:)*D1xDv(i,:) - N(i,:)*D3xDv(i,:))/12
!NNC             Tmp(v3,:) = Tmp(v3,:) + N(i,:)*(LxD3(i,:) - MxD2(i,:)) &
!NNC                        + (M(i,:)*D3xDv(i,:) - L(i,:)*D2xDv(i,:))/12
!NNC
!NNC          end do
!NNC
!NNC          do i = 1,ndim
!NNC             !Tmp(i,:) = 0.5_r8 + Tmp(i,:)/(24*Cell%Volume)
!NNC             Tmp(i,:) = 0.5_r8 + Tmp(i,:)/(24*mesh%volume)
!NNC          end do
!NNC
!NNC    end select
!NNC
!NNC    ! Compute the centroid
!NNC    do i = 1,ndim
!NNC
!NNC       select case (ndim)
!NNC
!NNC          case (2)
!NNC
!NNC             !Cell%Centroid(i) = LxD3(1,:)*(12*Xv(i,v4,:) + 6*(L(i,:) + M(i,:)) + 3*N(i,:)) + &
!NNC             !                   MxD2(1,:)*(6*Xv(i,v4,:) + 4*L(i,:) + 3*M(i,:) + 2*N(i,:)) + &
!NNC             !                   NxD1(1,:)*(6*Xv(i,v4,:) + 3*L(i,:) + 4*M(i,:) + 2*N(i,:))
!NNC             centroids(i,:) = LxD3(1,:)*(12*Xv(i,v4,:) + 6*(L(i,:) + M(i,:)) + 3*N(i,:)) + &
!NNC                                MxD2(1,:)*(6*Xv(i,v4,:) + 4*L(i,:) + 3*M(i,:) + 2*N(i,:)) + &
!NNC                                NxD1(1,:)*(6*Xv(i,v4,:) + 3*L(i,:) + 4*M(i,:) + 2*N(i,:))
!NNC             !Cell%Centroid(i) = Cell%Centroid(i)/(12*Cell%Volume)
!NNC             centroids(i,:) = centroids(i,:)/(12*mesh%volume)
!NNC
!NNC          case (3)
!NNC
!NNC             !Cell%Centroid(i) = Xv(i,v4,:) &
!NNC             !                 + Tmp(v1,:)*(Xv(i,v1,:) - Xv(i,v4,:)) &
!NNC             !                 + Tmp(v2,:)*(Xv(i,v3,:) - Xv(i,v4,:)) &
!NNC             !                 + Tmp(v1,:)*Tmp(v2,:)*(Xv(i,v2,:) + Xv(i,v4,:) - Xv(i,v1,:) - Xv(i,v3,:)) &
!NNC             !                 + Tmp(v3,:)*(Xv(i,v8,:) - Xv(i,v4,:) &
!NNC             !                 + Tmp(v1,:)*(Xv(i,v4,:) + Xv(i,v5,:) - Xv(i,v1,:) - Xv(i,v8,:)) &
!NNC             !                 + Tmp(v2,:)*(Xv(i,v4,:) + Xv(i,v7,:) - Xv(i,v3,:) - Xv(i,v8,:)) &
!NNC             !                 + Tmp(v1,:)*Tmp(v2,:)*Dv(i,:))
!NNC             centroids(i,:) = Xv(i,v4,:) &
!NNC                              + Tmp(v1,:)*(Xv(i,v1,:) - Xv(i,v4,:)) &
!NNC                              + Tmp(v2,:)*(Xv(i,v3,:) - Xv(i,v4,:)) &
!NNC                              + Tmp(v1,:)*Tmp(v2,:)*(Xv(i,v2,:) + Xv(i,v4,:) - Xv(i,v1,:) - Xv(i,v3,:)) &
!NNC                              + Tmp(v3,:)*(Xv(i,v8,:) - Xv(i,v4,:) &
!NNC                              + Tmp(v1,:)*(Xv(i,v4,:) + Xv(i,v5,:) - Xv(i,v1,:) - Xv(i,v8,:)) &
!NNC                              + Tmp(v2,:)*(Xv(i,v4,:) + Xv(i,v7,:) - Xv(i,v3,:) - Xv(i,v8,:)) &
!NNC                              + Tmp(v1,:)*Tmp(v2,:)*Dv(i,:))
!NNC
!NNC
!NNC       end select
!NNC
!NNC    end do
!NNC
!NNC    ! Explicitly deallocate temporaries
!NNC    call ARRAYDESTROY (L, 'Array L(ndim,ncells)')
!NNC    call ARRAYDESTROY (M, 'Array M(ndim,ncells)')
!NNC    call ARRAYDESTROY (N, 'Array N(ndim,ncells)')
!NNC    call ARRAYDESTROY (Xv, 'Array Xv(ndim,nvc,ncells)')
!NNC    call ARRAYDESTROY (LxD3, 'Array LxD3(nrot,ncells)')
!NNC    call ARRAYDESTROY (MxD2, 'Array MxD2(nrot,ncells)')
!NNC    call ARRAYDESTROY (NxD1, 'Array NxD1(nrot,ncells)')
!NNC    if (ndim == 3) then
!NNC       call ARRAYDESTROY (Tmp, 'Array Tmp(nvc,ncells)')
!NNC       call ARRAYDESTROY (D1, 'Array D1(ndim,ncells)')
!NNC       call ARRAYDESTROY (D2, 'Array D2(ndim,ncells)')
!NNC       call ARRAYDESTROY (D3, 'Array D3(ndim,ncells)')
!NNC       call ARRAYDESTROY (Dv, 'Array Dv(ndim,ncells)')
!NNC       call ARRAYDESTROY (D1xDv, 'Array D1xDv(nrot,ncells)')
!NNC       call ARRAYDESTROY (D2xDv, 'Array D2xDv(nrot,ncells)')
!NNC       call ARRAYDESTROY (D3xDv, 'Array D3xDv(nrot,ncells)')
!NNC    end if
!NNC
!NNC
!NNC  end subroutine sm_cell_centroid

End module solid_mechanics_mesh 
