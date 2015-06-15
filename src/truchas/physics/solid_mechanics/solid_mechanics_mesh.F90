!!
!!
!! Solid Mechanics Mesh Module
!!      
!! This module holds a copy of the mesh pointer and defines subroutines
!! and functions related to the mesh required throughout the solid mechanics physics
!! package.
!! 
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

  use dist_mesh_type
  use parallel_permutations, only: par_perm

  implicit none
  private

  !! Solid mechanics DIST_MESH Pointer
  type(dist_mesh), pointer, public :: sm_mesh => null()

  !! Permutation structures to transfer data to and from old mesh
  type(par_perm), target, save, public :: sm_mesh_to_old_cell, sm_old_to_mesh_cell
  type(par_perm), target, save, public :: sm_mesh_to_old_node, sm_old_to_mesh_node

  !! I don't know if this should be a module variable or local in sm_mesh_init
  integer, pointer, public, save :: sm_gap_elements(:) => null()

  !!
  public sm_mesh_init, sm_mesh_enable, sm_cell_shape

  !! Imported Cell geometry subroutines
  public :: sm_face_centroid, sm_cell_centroid, sm_face_centroid_logical

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
    use mesh_broker, only: dist_mesh_ptr
    use mesh_module, only: unpermute_mesh_vector, unpermute_vertex_vector
    use parallel_permutations, only: create_par_perm

    integer, pointer :: dummy(:) => null()
    logical :: found

    ! Initialize the mesh pointer
    sm_mesh => dist_mesh_ptr(sm_mesh_name)
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
    use mesh_broker, only : enable_mesh
    logical :: found
    call enable_mesh(sm_mesh_name,found)
    INSIST(found)
  end subroutine sm_mesh_enable
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  
  function sm_cell_shape (mesh) result (cell_shape)
    type(dist_mesh), intent(in) :: mesh
    integer :: cell_shape

    select case ( size(mesh%cnode,dim=1) )
      case (4)
        cell_shape=SM_CELL_TET
      case (8)
        cell_shape=SM_CELL_HEX
      case default
        cell_shape=SM_CELL_HEX
    end select

  end function sm_cell_shape
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  subroutine sm_face_centroid(mesh,centroids)
    use kinds,         only: r8
    use linear_module, only: linear_prop

    type(dist_mesh), intent(in) :: mesh
    real(r8), intent(out) :: centroids(ndim,nfc,mesh%ncell)

    integer :: j, k, inode 
    integer :: f, n
    real(r8), allocatable :: coord(:,:,:)

    allocate(coord(ndim,size(mesh%cnode,dim=1),mesh%ncell))

    !! Replace the EN_GATHER call in the NDIM_LOOP
    do j = 1, mesh%ncell
      do k = 1, size(mesh%cnode,dim=1)
        inode = mesh%cnode(k,j)
        coord(:,k,j) = mesh%x(:,inode)
      end do  
    end do

    ! Now loop through the dimensions and faces
    do n = 1, ndim
      do f = 1,nfc
        call linear_prop(f,coord(n,:,:),centroids(n,f,:))
      end do 
    end do 

    ! Deallocate
    deallocate(coord)

   
  end subroutine sm_face_centroid  
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  subroutine sm_face_centroid_logical(mesh,centroids)
    !=======================================================================
    ! Purpose(s):
    !
    !   Replacement for Cell%Face_Centroid_L. Lifted from face_centroid_logical
    !   in cell_geometry_module.F90
    !
    !   Compute the face centroid (Face_centroid), which is part
    !   of the Cell structure. A face centroid is assigned for
    !   each face of the cell, and is a LOGICAL coordinate.
    !     Input  - Vertex%Coord(ndim)
    !     Output - Cell%Face_centroid_L(ndim,nfc)
    !
    !=======================================================================
    use kinds,                only: r8
    use ArrayAllocate_Module, only: ARRAYCREATE, ARRAYDESTROY
    use cutoffs_module,       only: alittle

    ! Input/Output
    type(dist_mesh), intent(in)  :: mesh
    real(r8),        intent(out) :: centroids(ndim,nfc,mesh%ncell)

    ! Local Variables
    integer :: j, k, inode, iface
    integer :: f, i, i1, i2, n, v11, v12, v13, v14,    &
                                v21, v22, v23, v24, v31, v32, v33, v34, &
                                nn
    integer :: n1 = 1, n2 = 2, n3 = 3
    real(r8), dimension(ndim)   :: Coef
    !real(r8), dimension(ncells) :: Face_Area
    real(r8), allocatable :: Face_Area(:,:)
    real(r8), allocatable :: Face_Normal(:,:,:)
    real(r8), pointer, dimension(:,:,:) :: Xv

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Only do this if it is 3-D.
    if (ndim == 3) then
       ! Allocate temporaries
       !call ARRAYCREATE (Xv, 1, ndim, 1, nvc, 1, ncells, 'Array Xv(ndim,nvc,ncells)')
       call ARRAYCREATE (Xv, 1, ndim, 1, nvc, 1, mesh%ncell, 'Array Xv(ndim,nvc,ncells)')

       ! Gather vertex coordinates
       !do i = 1,ndim
       !   call EN_GATHER (Xv(i,:,:), Vertex%Coord(i), BOUNDARY=Vrtx_Bdy(i)%Data)
       !end do
       do j = 1, mesh%ncell
         do k = 1, size(mesh%cnode,dim=1)
           inode = mesh%cnode(k,j)
           Xv(:,k,j) = mesh%x(:,inode)
         end do  
       end do  
    end if

    ! Compute the Face_Area and the Face_Normal arrays
    allocate(Face_Area(nfc,mesh%ncell))
    allocate(Face_Normal(ndim,nfc,mesh%ncell))
    Face_Area = 0.0_r8
    do j = 1, mesh%ncell
      do f = 1, size(mesh%cface,dim=1)
        iface = mesh%cface(f,j)
        Face_Area(f,j) = mesh%area(iface)
        Face_Normal(:,f,j) = mesh%normal(:,iface)
      end do 
    end do    

    ! Loop over faces
    FACE_LOOP: do f = 1,nfc

       ! Preinitialize the face centroid
       do n = 1,ndim
          !Cell%Face_Centroid_L(n,f) = 0.5_r8
          CENTROIDS(n,f,:)          = 0.5_r8
          Coef(n)                   = 1
       end do

       ! We are done if this is 2-D.
       if (ndim == 2) cycle FACE_LOOP

       ! Defined outside the loop  
       !Face_Area = 1.0 / (12*(Cell%Face_Area(f) + alittle))

       select case (f)

          case (1)

          ! Side 1 (vertices 4-8-7-3)
          v11 = 1; v12 = 1; v13 = 1; v14 = 1
          v21 = 7; v22 = 8; v23 = 3; v24 = 4
          v31 = 8; v32 = 4; v33 = 7; v34 = 3
          nn  = 1

          !Cell%Face_centroid_L(nn,f) = 0; Coef(nn) = 0
          CENTROIDS(nn,f,:) = 0; Coef(nn) = 0

          case (2)

          ! Side 2 (vertices 1-2-6-5)
          v11 = 1; v12 = 1; v13 = 1; v14 = 1
          v21 = 2; v22 = 1; v23 = 6; v24 = 5
          v31 = 6; v32 = 2; v33 = 5; v34 = 1
          nn  = 1

          !Cell%Face_centroid_L(nn,f) = 1; Coef(nn) = 0
          CENTROIDS(nn,f,:) = 1; Coef(nn) = 0

          case (3)

          ! Side 3 (vertices 4-1-5-8)
          v11 = 1; v12 = 4; v13 = 5; v14 = 8
          v21 = 1; v22 = 1; v23 = 1; v24 = 1
          v31 = 5; v32 = 1; v33 = 8; v34 = 4
          nn  = 2

          !Cell%Face_centroid_L(nn,f) = 0; Coef(nn) = 0
          CENTROIDS(nn,f,:) = 0; Coef(nn) = 0

          case (4)

          ! Side 4 (vertices 3-7-6-2)
          v11 = 6; v12 = 7; v13 = 2; v14 = 3
          v21 = 1; v22 = 1; v23 = 1; v24 = 1
          v31 = 7; v32 = 3; v33 = 6; v34 = 2
          nn  = 2

          !Cell%Face_centroid_L(nn,f) = 1; Coef(nn) = 0
          CENTROIDS(nn,f,:) = 1; Coef(nn) = 0

          case (5)

          ! Side 5 (vertices 4-3-2-1)
          v11 = 2; v12 = 3; v13 = 1; v14 = 4
          v21 = 3; v22 = 4; v23 = 2; v24 = 1
          v31 = 1; v32 = 1; v33 = 1; v34 = 1
          nn  = 3

          !Cell%Face_centroid_L(nn,f) = 0; Coef(nn) = 0
          CENTROIDS(nn,f,:) = 0; Coef(nn) = 0

          case (6)

          ! Side 6 (vertices 8-5-6-7)
          v11 = 5; v12 = 8; v13 = 6; v14 = 7
          v21 = 6; v22 = 5; v23 = 7; v24 = 8
          v31 = 1; v32 = 1; v33 = 1; v34 = 1
          nn  = 3

          !Cell%Face_centroid_L(nn,f) = 1; Coef(nn) = 0
          CENTROIDS(nn,f,:) = 1; Coef(nn) = 0

       end select

       NDIM_LOOP: do i = 1,ndim

          select case (i)
             case (1)
                i1 = 2; i2 = 3
             case (2)
                i1 = 3; i2 = 1
             case (3)
                i1 = 1; i2 = 2
          end select

          !Cell%Face_centroid_L(n1,f) = Cell%Face_centroid_L(n1,f) + &
          !      Coef(n1)*Face_Area*Cell%Face_Normal(i,f)* &
          !      ((Xv(i1,v11,:) - Xv(i1,v12,:))*(Xv(i2,v13,:) - Xv(i2,v14,:)) &
          !     - (Xv(i2,v11,:) - Xv(i2,v12,:))*(Xv(i1,v13,:) - Xv(i1,v14,:)))
          CENTROIDS(n1,f,:) = CENTROIDS(n1,f,:) + &
                Coef(n1)*Face_Area(f,:)*Face_Normal(i,f,:)* &
                ((Xv(i1,v11,:) - Xv(i1,v12,:))*(Xv(i2,v13,:) - Xv(i2,v14,:)) &
               - (Xv(i2,v11,:) - Xv(i2,v12,:))*(Xv(i1,v13,:) - Xv(i1,v14,:)))


          !Cell%Face_centroid_L(n2,f) = Cell%Face_centroid_L(n2,f) + &
          !     Coef(n2)*Face_Area*Cell%Face_Normal(i,f)* &
          !     ((Xv(i1,v21,:) - Xv(i1,v22,:))*(Xv(i2,v23,:) - Xv(i2,v24,:)) &
          !    - (Xv(i2,v21,:) - Xv(i2,v22,:))*(Xv(i1,v23,:) - Xv(i1,v24,:)))
          CENTROIDS(n2,f,:) = CENTROIDS(n2,f,:) + &
               Coef(n2)*Face_Area(f,:)*Face_Normal(i,f,:)* &
               ((Xv(i1,v21,:) - Xv(i1,v22,:))*(Xv(i2,v23,:) - Xv(i2,v24,:)) &
              - (Xv(i2,v21,:) - Xv(i2,v22,:))*(Xv(i1,v23,:) - Xv(i1,v24,:)))

          !Cell%Face_centroid_L(n3,f) = Cell%Face_centroid_L(n3,f) + &
          !     Coef(n3)*Face_Area*Cell%Face_Normal(i,f)* &
          !     ((Xv(i1,v31,:) - Xv(i1,v32,:))*(Xv(i2,v33,:) - Xv(i2,v34,:)) &
          !    - (Xv(i2,v31,:) - Xv(i2,v32,:))*(Xv(i1,v33,:) - Xv(i1,v34,:)))
          CENTROIDS(n3,f,:) = CENTROIDS(n3,f,:) + &
               Coef(n3)*Face_Area(f,:)*Face_Normal(i,f,:)* &
               ((Xv(i1,v31,:) - Xv(i1,v32,:))*(Xv(i2,v33,:) - Xv(i2,v34,:)) &
              - (Xv(i2,v31,:) - Xv(i2,v32,:))*(Xv(i1,v33,:) - Xv(i1,v34,:)))

       end do NDIM_LOOP

    end do FACE_LOOP

    ! Deallocate temporaries
    if (ASSOCIATED(Xv)) call ARRAYDESTROY (Xv, 'Array Xv(ndim,nvc,ncells)')

    ! Deallocate the face and normal temporary arrays
    deallocate(Face_Area)
    deallocate(Face_Normal)



  end subroutine sm_face_centroid_logical

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  subroutine sm_cell_centroid(mesh,centroids)
    use kinds,                only : r8
    use ArrayAllocate_Module, only : ARRAYCREATE, ARRAYDESTROY
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the cell centroids, given a distributed mesh type
    !     Originally defined in cell_geometry_module.F90 in subroutine
    !     cell_centroid
    !
    !=======================================================================

    ! Input/Output
    type(dist_mesh), intent(in)  :: mesh
    real(r8),        intent(out) :: centroids(ndim,mesh%ncell)



    ! Local Variables
    integer :: j, k, inode, ncells
    integer :: i, i1, i2, v1, v2, v3, v4, v5, v6, v7, v8
    real(r8), pointer, dimension(:,:) :: L, M, N, LxD3, MxD2, NxD1
    real(r8), pointer, dimension(:,:) :: Tmp, D1, D2, D3, Dv, D1xDv, D2xDv, D3xDv
    real(r8), pointer, dimension(:,:,:) :: Xv

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Define the mesh parameters
    ncells = mesh%ncell

    ! Explicitly allocate temporaries
    call ARRAYCREATE (L, 1, ndim, 1, ncells, 'Array L(ndim,ncells)')
    call ARRAYCREATE (M, 1, ndim, 1, ncells, 'Array M(ndim,ncells)')
    call ARRAYCREATE (N, 1, ndim, 1, ncells, 'Array N(ndim,ncells)')
    call ARRAYCREATE (Xv, 1, ndim, 1, nvc, 1, ncells, 'Array Xv(ndim,nvc,ncells)')
    call ARRAYCREATE (LxD3, 1, nrot, 1, ncells, 'Array LxD3(nrot,ncells)')
    call ARRAYCREATE (MxD2, 1, nrot, 1, ncells, 'Array MxD2(nrot,ncells)')
    call ARRAYCREATE (NxD1, 1, nrot, 1, ncells, 'Array NxD1(nrot,ncells)')

    ! Only need these for 3-D.
    if (ndim == 3) then
       call ARRAYCREATE (Tmp, 1, ndim, 1, ncells, 'Array Tmp(ndim,ncells)')
       call ARRAYCREATE (D1, 1, ndim, 1, ncells, 'Array D1(ndim,ncells)')
       call ARRAYCREATE (D2, 1, ndim, 1, ncells, 'Array D2(ndim,ncells)')
       call ARRAYCREATE (D3, 1, ndim, 1, ncells, 'Array D3(ndim,ncells)')
       call ARRAYCREATE (Dv, 1, ndim, 1, ncells, 'Array Dv(ndim,ncells)')
       call ARRAYCREATE (D1xDv, 1, nrot, 1, ncells, 'Array D1xDv(nrot,ncells)')
       call ARRAYCREATE (D2xDv, 1, nrot, 1, ncells, 'Array D2xDv(nrot,ncells)')
       call ARRAYCREATE (D3xDv, 1, nrot, 1, ncells, 'Array D3xDv(nrot,ncells)')
    end if

    ! Set vertex numbers
    select case(ndim)
       case (2)
          v1 = 1; v2 = 2; v3 = 3; v4 = 4
       case (3)
          v1 = 1; v2 = 2; v3 = 3; v4 = 4
          v5 = 5; v6 = 6; v7 = 7; v8 = 8
    end select

    ! Gather vertex coordinates
    !do i = 1,ndim
    !   call EN_GATHER (Xv(i,:,:), Vertex%Coord(i), BOUNDARY=Vrtx_Bdy(i)%Data)
    !end do
    do j = 1, ncells
      do k = 1, nvc
        inode = SM_MESH%CNODE(k,j)
        Xv(:,k,j) = SM_MESH%X(:,inode)
      end do  
    end do  

    ! Compute quantities needed for the centroid coordinates
    do i = 1,ndim

       select case (ndim)

          case (2)

             L(i,:) = Xv(i,v1,:) - Xv(i,v4,:)
             M(i,:) = Xv(i,v3,:) - Xv(i,v4,:)
             N(i,:) = - Xv(i,v1,:) + Xv(i,v2,:) - Xv(i,v3,:) + Xv(i,v4,:)

          case (3)

             L(i,:) = Xv(i,v1,:) + Xv(i,v2,:) + Xv(i,v6,:) + Xv(i,v5,:) &
                    - Xv(i,v3,:) - Xv(i,v4,:) - Xv(i,v8,:) - Xv(i,v7,:)
             M(i,:) = Xv(i,v2,:) + Xv(i,v3,:) + Xv(i,v7,:) + Xv(i,v6,:) &
                    - Xv(i,v4,:) - Xv(i,v1,:) - Xv(i,v5,:) - Xv(i,v8,:)
             N(i,:) = Xv(i,v8,:) + Xv(i,v5,:) + Xv(i,v6,:) + Xv(i,v7,:) &
                    - Xv(i,v3,:) - Xv(i,v2,:) - Xv(i,v1,:) - Xv(i,v4,:)
             D1(i,:) = Xv(i,v2,:) + Xv(i,v6,:) + Xv(i,v4,:) + Xv(i,v8,:) &
                     - Xv(i,v1,:) - Xv(i,v5,:) - Xv(i,v3,:) - Xv(i,v7,:)
             D2(i,:) = Xv(i,v3,:) + Xv(i,v4,:) + Xv(i,v5,:) + Xv(i,v6,:) &
                     - Xv(i,v1,:) - Xv(i,v2,:) - Xv(i,v7,:) - Xv(i,v8,:)
             D3(i,:) = Xv(i,v1,:) + Xv(i,v4,:) + Xv(i,v6,:) + Xv(i,v7,:) &
                     - Xv(i,v2,:) - Xv(i,v3,:) - Xv(i,v5,:) - Xv(i,v8,:)
             Dv(i,:) = Xv(i,v1,:) + Xv(i,v3,:) + Xv(i,v6,:) + Xv(i,v8,:) &
                     - Xv(i,v2,:) - Xv(i,v4,:) - Xv(i,v5,:) - Xv(i,v7,:)

       end select

    end do

    select case (ndim)

       case (2)

          i1 = 1; i2 = 2 
          do i = 1,nrot
             LxD3(i,:) = L(i1,:)*M(i2,:) - L(i2,:)*M(i1,:)
             MxD2(i,:) = L(i1,:)*N(i2,:) - L(i2,:)*N(i1,:)
             NxD1(i,:) = N(i1,:)*M(i2,:) - N(i2,:)*M(i1,:)
          end do

       case (3)

          L = 0.25_r8*L
          M = 0.25_r8*M
          N = 0.25_r8*N
          Tmp = 0

          do i = 1,nrot

             select case (i)
             case (1)
                i1 = 2; i2 = 3
             case (2)
                i1 = 3; i2 = 1
             case (3)
                i1 = 1; i2 = 2
             end select

             LxD3(i,:) = L(i1,:)*D3(i2,:) - L(i2,:)*D3(i1,:)
             MxD2(i,:) = M(i1,:)*D2(i2,:) - M(i2,:)*D2(i1,:)
             NxD1(i,:) = N(i1,:)*D1(i2,:) - N(i2,:)*D1(i1,:)

             D1xDv(i,:) = D1(i1,:)*Dv(i2,:) - D1(i2,:)*Dv(i1,:)
             D2xDv(i,:) = D2(i1,:)*Dv(i2,:) - D2(i2,:)*Dv(i1,:)
             D3xDv(i,:) = D3(i1,:)*Dv(i2,:) - D3(i2,:)*Dv(i1,:)

          end do

          do i = 1,ndim

             Tmp(v1,:) = Tmp(v1,:) + L(i,:)*(MxD2(i,:) - NxD1(i,:)) &
                        + (N(i,:)*D2xDv(i,:) - M(i,:)*D1xDv(i,:))/12
             Tmp(v2,:) = Tmp(v2,:) + M(i,:)*(NxD1(i,:) - LxD3(i,:)) &
                        + (L(i,:)*D1xDv(i,:) - N(i,:)*D3xDv(i,:))/12
             Tmp(v3,:) = Tmp(v3,:) + N(i,:)*(LxD3(i,:) - MxD2(i,:)) &
                        + (M(i,:)*D3xDv(i,:) - L(i,:)*D2xDv(i,:))/12

          end do

          do i = 1,ndim
             !Tmp(i,:) = 0.5_r8 + Tmp(i,:)/(24*Cell%Volume)
             Tmp(i,:) = 0.5_r8 + Tmp(i,:)/(24*mesh%volume)
          end do

    end select

    ! Compute the centroid
    do i = 1,ndim

       select case (ndim)

          case (2)

             !Cell%Centroid(i) = LxD3(1,:)*(12*Xv(i,v4,:) + 6*(L(i,:) + M(i,:)) + 3*N(i,:)) + &
             !                   MxD2(1,:)*(6*Xv(i,v4,:) + 4*L(i,:) + 3*M(i,:) + 2*N(i,:)) + &
             !                   NxD1(1,:)*(6*Xv(i,v4,:) + 3*L(i,:) + 4*M(i,:) + 2*N(i,:))
             centroids(i,:) = LxD3(1,:)*(12*Xv(i,v4,:) + 6*(L(i,:) + M(i,:)) + 3*N(i,:)) + &
                                MxD2(1,:)*(6*Xv(i,v4,:) + 4*L(i,:) + 3*M(i,:) + 2*N(i,:)) + &
                                NxD1(1,:)*(6*Xv(i,v4,:) + 3*L(i,:) + 4*M(i,:) + 2*N(i,:))
             !Cell%Centroid(i) = Cell%Centroid(i)/(12*Cell%Volume)
             centroids(i,:) = centroids(i,:)/(12*mesh%volume)

          case (3)

             !Cell%Centroid(i) = Xv(i,v4,:) &
             !                 + Tmp(v1,:)*(Xv(i,v1,:) - Xv(i,v4,:)) &
             !                 + Tmp(v2,:)*(Xv(i,v3,:) - Xv(i,v4,:)) &
             !                 + Tmp(v1,:)*Tmp(v2,:)*(Xv(i,v2,:) + Xv(i,v4,:) - Xv(i,v1,:) - Xv(i,v3,:)) &
             !                 + Tmp(v3,:)*(Xv(i,v8,:) - Xv(i,v4,:) &
             !                 + Tmp(v1,:)*(Xv(i,v4,:) + Xv(i,v5,:) - Xv(i,v1,:) - Xv(i,v8,:)) &
             !                 + Tmp(v2,:)*(Xv(i,v4,:) + Xv(i,v7,:) - Xv(i,v3,:) - Xv(i,v8,:)) &
             !                 + Tmp(v1,:)*Tmp(v2,:)*Dv(i,:))
             centroids(i,:) = Xv(i,v4,:) &
                              + Tmp(v1,:)*(Xv(i,v1,:) - Xv(i,v4,:)) &
                              + Tmp(v2,:)*(Xv(i,v3,:) - Xv(i,v4,:)) &
                              + Tmp(v1,:)*Tmp(v2,:)*(Xv(i,v2,:) + Xv(i,v4,:) - Xv(i,v1,:) - Xv(i,v3,:)) &
                              + Tmp(v3,:)*(Xv(i,v8,:) - Xv(i,v4,:) &
                              + Tmp(v1,:)*(Xv(i,v4,:) + Xv(i,v5,:) - Xv(i,v1,:) - Xv(i,v8,:)) &
                              + Tmp(v2,:)*(Xv(i,v4,:) + Xv(i,v7,:) - Xv(i,v3,:) - Xv(i,v8,:)) &
                              + Tmp(v1,:)*Tmp(v2,:)*Dv(i,:))


       end select

    end do

    ! Explicitly deallocate temporaries
    call ARRAYDESTROY (L, 'Array L(ndim,ncells)')
    call ARRAYDESTROY (M, 'Array M(ndim,ncells)')
    call ARRAYDESTROY (N, 'Array N(ndim,ncells)')
    call ARRAYDESTROY (Xv, 'Array Xv(ndim,nvc,ncells)')
    call ARRAYDESTROY (LxD3, 'Array LxD3(nrot,ncells)')
    call ARRAYDESTROY (MxD2, 'Array MxD2(nrot,ncells)')
    call ARRAYDESTROY (NxD1, 'Array NxD1(nrot,ncells)')
    if (ndim == 3) then
       call ARRAYDESTROY (Tmp, 'Array Tmp(nvc,ncells)')
       call ARRAYDESTROY (D1, 'Array D1(ndim,ncells)')
       call ARRAYDESTROY (D2, 'Array D2(ndim,ncells)')
       call ARRAYDESTROY (D3, 'Array D3(ndim,ncells)')
       call ARRAYDESTROY (Dv, 'Array Dv(ndim,ncells)')
       call ARRAYDESTROY (D1xDv, 'Array D1xDv(nrot,ncells)')
       call ARRAYDESTROY (D2xDv, 'Array D2xDv(nrot,ncells)')
       call ARRAYDESTROY (D3xDv, 'Array D3xDv(nrot,ncells)')
    end if


  end subroutine sm_cell_centroid

End module solid_mechanics_mesh 
