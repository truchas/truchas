!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  unsplit_geometry_volume_tracker_type
!! 
!!  Author: Robert Chiodi (robertchiodi@lanl.gov)
!!  June 2019
!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module unsplit_geometric_volume_tracker_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64, i8 => int64
  use volume_tracker_class
  use truchas_logging_services
  use truchas_timers
  use unstr_mesh_type
  use index_partitioning
  use irl_fortran_interface
  implicit none
  private

  integer, parameter, private :: advect_band = 5 ! Size of band to do advection

  type, extends(volume_tracker), public :: unsplit_geometric_volume_tracker
    private
    type(unstr_mesh), pointer :: mesh ! unowned reference
    integer :: location_iter_max ! maximum number of iterations to use in fitting interface
    integer :: subcycles
    logical :: nested_dissection
    real(r8), allocatable :: normal(:,:,:)
    real(r8), allocatable :: face_flux(:,:)
    ! node/face/cell workspace
    real(r8), allocatable :: w_node(:,:)
    real(r8), allocatable :: velocity_weightings(:)
    real(r8), allocatable :: projected_nodes(:,:)
    integer, allocatable :: priority(:), bc_index(:), local_face(:), inflow_mat(:)
    integer :: nrealfluid, nfluid, nmat ! # of non-void fluids, # of fluids incl. void, # of materials
    integer, allocatable :: boundary_recon_to_cell(:) ! Mapping of boundary reconstruction ID to neighboring inside cell ID
    ! Start IRL objects
    type(ObjServer_PlanarLoc_type) :: object_server_planar_localizer
    type(ObjServer_PlanarSep_type) :: object_server_planar_separator
    type(ObjServer_LocSepLink_type) :: object_server_localized_separator_link
    type(PlanarLoc_type), allocatable :: planar_localizer(:)
    type(PlanarSep_type), allocatable :: planar_separator(:)
    type(LocSepLink_type), allocatable :: localized_separator_link(:)
    type(Tet_type) :: IRL_tet
    type(Pyrmd_type) :: IRL_pyramid
    type(TriPrism_type) :: IRL_wedge
    type(Hex_type) :: IRL_hex
    type(CapDod_LLLL_type) :: IRL_CapDod_LLLL
    type(CapDod_LLLT_type) :: IRL_CapDod_LLLT
    type(CapDod_LTLT_type) :: IRL_CapDod_LTLT
    type(CapDod_LLTT_type) :: IRL_CapDod_LLTT
    type(CapDod_LTTT_type) :: IRL_CapDod_LTTT
    type(CapDod_TTTT_type) :: IRL_CapDod_TTTT
    type(CapOcta_LLL_type) :: IRL_CapOcta_LLL
    type(CapOcta_LLT_type) :: IRL_CapOcta_LLT
    type(CapOcta_LTT_type) :: IRL_CapOcta_LTT
    type(CapOcta_TTT_type) :: IRL_CapOcta_TTT            
    ! End IRL objects
    ! Used for IRL advection
    integer, allocatable :: flux_geometry_class(:)
    integer, allocatable :: flux_node(:)
  contains
    procedure, public  :: init
    procedure, public  :: flux_volumes
    procedure, public  :: set_inflow_material
    procedure, private :: init_irl_mesh
    procedure, private :: generate_flux_classes
    procedure, private :: compute_velocity_weightings
    procedure, private :: set_irl_interfaces
    procedure, private :: adjust_planes_match_VOF
    procedure, private :: reset_volume_moments
    procedure, private :: compute_node_velocities
    procedure, private :: compute_projected_nodes
    procedure, private :: compute_fluxes
    procedure, private :: update_vof
    procedure, private :: normals
    procedure, private :: getBCMaterialFractions
  end type unsplit_geometric_volume_tracker

contains

  subroutine init(this, mesh, nrealfluid, nfluid, nmat, liq_matid, params)

    use parameter_list_type
    use property_module, only: get_truchas_material_id
    use f08_intrinsics, only: findloc

    class(unsplit_geometric_volume_tracker), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    integer, intent(in) :: nrealfluid, nfluid, nmat, liq_matid(:)
    type(parameter_list), intent(inout) :: params
    integer :: i, j, k
    
    this%mesh => mesh
    call this%mesh%init_cell_centroid
    call this%mesh%init_face_centroid
    call this%mesh%init_face_normal_dist
    this%nrealfluid = nrealfluid
    this%nfluid = nfluid
    this%nmat = nmat

    call params%get('location_iter_max', this%location_iter_max, default=40)
    call params%get('cutoff', this%cutoff, default=1.0e-8_r8)
    call params%get('subcycles', this%subcycles, default=2)
    print*,'Subcycling is currently disabled for unsplit transport'
    call params%get('nested_dissection', this%nested_dissection, default=.true.)

    ! convert user material ids to array index
    if (params%is_vector('material_priority')) then
      call params%get('material_priority', this%priority)
      do i = 1,size(this%priority)
        if (this%priority(i) < 1) cycle ! solid (-1) is handled later
        this%priority(i) = findloc(liq_matid, get_truchas_material_id(this%priority(i)))

        ! make sure we found a liquid material
        ! TODO: need better error handling here
        INSIST(this%priority(i) > 0)
      end do

      ! The current expectation is that a user will
      ! use a material number of -1 to indicate solid.
      ! Internally, if solid is present it is the last material
      if (this%nmat > this%nfluid) then
        where (this%priority == -1) this%priority = this%nmat
      end if
    else
      this%priority = [(i, i=1,this%nmat)]
    end if
    ASSERT(size(this%priority) == this%nmat)
    ASSERT(all(this%priority > 0) .and. all(this%priority <= this%nmat))
    ! TODO: assert that each material appears exactly once
    
    ! Allocate face fluxes we'll store
    allocate(this%face_flux(this%nmat, this%mesh%nface))

    allocate(this%normal(3,this%nmat,mesh%ncell))
    allocate(this%w_node(3,mesh%nnode))
    allocate(this%velocity_weightings(mesh%xcnode(mesh%ncell+1)-1))
    allocate(this%projected_nodes(3, mesh%nnode))

    ! list of boundary face ids
    j = count(this%mesh%fcell(2,:this%mesh%nface) == 0)
    allocate(this%bc_index(j), this%local_face(j), this%inflow_mat(j))
    j = 1
    do i = 1, this%mesh%nface
      if (this%mesh%fcell(2,i) == 0) then
        this%bc_index(j) = i
        k = this%mesh%fcell(1,i)
        this%local_face(j) = &
            findloc(this%mesh%cface(this%mesh%xcface(k):this%mesh%xcface(k+1)-1), i)
        j = j + 1
      end if
    end do
    this%inflow_mat = 0
    
    ! Initialize IRL arrays needed for advection
    call this%init_irl_mesh()
    
    ! Initialize IRL cell types we will use
    call new(this%IRL_tet)
    call new(this%IRL_pyramid)
    call new(this%IRL_wedge)
    call new(this%IRL_hex)
    call new(this%IRL_CapDod_LLLL)
    call new(this%IRL_CapDod_LLLT)
    call new(this%IRL_CapDod_LTLT)
    call new(this%IRL_CapDod_LLTT)
    call new(this%IRL_CapDod_LTTT)
    call new(this%IRL_CapDod_TTTT)
    call new(this%IRL_CapOcta_LLL)
    call new(this%IRL_CapOcta_LLT)
    call new(this%IRL_CapOcta_LTT)
    call new(this%IRL_CapOcta_TTT)            

    ! Reorganize and generate needed face-flux information
    call this%generate_flux_classes

    ! Calculate weightings for interpolation of cell center velocity to node
    call this%compute_velocity_weightings

  end subroutine init

  ! flux volumes routine assuming flux_vol is a cface-like array
  ! flux volumes routine assuming vel is stored on faces, length 1:mesh%nfaces
  subroutine flux_volumes(this, vel, vel_cc, vof_n, vof, flux_vol, fluids, void, dt, a_vof_band)
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: vel(:), vel_cc(:,:), vof_n(:,:), dt
    real(r8), intent(out) :: flux_vol(:,:), vof(:,:)
    integer, intent(in) :: fluids, void
    integer, intent(in) :: a_vof_band(:,:)

    integer :: i,j

    flux_vol = 0.0_r8
    vof = vof_n

    call start_timer('reconstruction')
    call this%normals(vof)    
    call this%set_irl_interfaces(vof)    
    call this%reset_volume_moments(vof)
    call stop_timer('reconstruction')

    call start_timer('advection')
    call this%compute_node_velocities(vel_cc)

    call this%compute_projected_nodes(dt)
    
    call this%compute_fluxes(vel, dt, vof, a_vof_band)
    
    ! Need to communicate fluxes
    call gather_boundary(this%mesh%face_ip, this%face_flux)
          
    call this%update_vof(a_vof_band, flux_vol, vof)
   
    call stop_timer('advection')
    
  end subroutine flux_volumes

  !! Set the inflow material for the given boundary faces. A material index
  !! of 0 will result in materials being fluxed in proportion to the material
  !! fractions present in the cell. This is the preset default. The BC_INDEX
  !! array component is ordered, so we use a binary search to locate the faces
  !! in the the array, so we can set the corresponding inflow material index.
  !! TODO: If FACES is also ordered (likely) the search can be improved further.

  subroutine set_inflow_material(this, mat, faces)
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    integer, intent(in) :: mat  ! material index
    integer, intent(in) :: faces(:) ! face indices
    integer :: i
    ASSERT(mat >= 0)
    do i = 1, size(faces)
      block
        integer :: k, k1, k2
        k1 = 1; k2 = size(this%bc_index)
        do while (k1 < k2)
          k = (k1 + k2)/2
          if (this%bc_index(k) < faces(i)) then
            k1 = k + 1
          else
            k2 = k
          end if
        end do
        ASSERT(k1 == k2)
        ASSERT(this%bc_index(k1) == faces(i))
        this%inflow_mat(k1) = mat
      end block
    end do
  end subroutine set_inflow_material
    
  subroutine init_irl_mesh(this)
  
    use cell_topology, only : get_face_nodes
    use irl_interface_helper
    
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    
    integer :: j, f, face_index, cn
    integer :: number_of_cell_faces
    integer :: total_recon_needed
    real(r8) :: tmp_plane(4)
    integer :: face_to_boundary_mapping(this%mesh%nface)
    integer :: found_internal, found_boundary, total_internal

    ! QUESTION : face_to_... would be much better as a
    ! hash table, as most indices are empty/not used.
    
    call setVFBounds(this%cutoff)                                         
    call setVFTolerance_IterativeDistanceFinding(1.0e-12_r8)
    
    ! Allocate reconstruction for each cell and one for each boundary face.
    ! We're going to kind-of fake a ghost cell
    total_recon_needed = this%mesh%ncell+size(this%bc_index)

    ! Mapping of boundary_id to the corresponding cell that is
    ! actually inside the domain. Index for the boundary_id
    ! in the array is boundary_id - this%mesh%ncell
    allocate(this%boundary_recon_to_cell(size(this%bc_index)))

    do j = 1, size(this%bc_index)
       this%boundary_recon_to_cell(j) = this%mesh%fcell(1,this%bc_index(j))
       face_to_boundary_mapping(this%bc_index(j)) = this%mesh%ncell + j
    end do
    
    ! Allocate block storage to increase locality of objects we will be creating
    call new(this%object_server_planar_localizer, int(total_recon_needed, i8))
    call new(this%object_server_planar_separator, int(total_recon_needed, i8))
    call new(this%object_server_localized_separator_link, int(total_recon_needed, i8))   
        
    ! Allocate the Fortran memory
    allocate(this%planar_localizer(total_recon_needed))
    allocate(this%planar_separator(total_recon_needed))
    allocate(this%localized_separator_link(total_recon_needed))
    
    ! Now allocate on the IRL side the different planar reconstructions
    do j = 1, total_recon_needed
      call new(this%planar_localizer(j), this%object_server_planar_localizer)
      call new(this%planar_separator(j), this%object_server_planar_separator)
      call new(this%localized_separator_link(j), &
               this%object_server_localized_separator_link, &
               this%planar_localizer(j), &
               this%planar_separator(j) )
    end do

    ! Now setup PlanarLocalizers (one per cell) and link together.
    do j = 1, this%mesh%ncell    
      number_of_cell_faces = this%mesh%xcface(j+1)-this%mesh%xcface(j)
                                                                    
      call setNumberOfPlanes(this%planar_localizer(j), number_of_cell_faces)
      call setId(this%localized_separator_link(j), j)
            
      associate(face_id => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1), &
                cn => this%mesh%cnhbr(this%mesh%xcnhbr(j):this%mesh%xcnhbr(j+1)-1))
        
        total_internal = count(cn /= 0)
        found_internal = 0
        found_boundary = 0        
        do f = 1, number_of_cell_faces
          ! Set PlanarLocalizer Plane for this face
          face_index = face_id(f)
          tmp_plane(1:3) = this%mesh%normal(:,face_index) / this%mesh%area(face_index)
          tmp_plane(4) = dot_product(tmp_plane(1:3), this%mesh%x(:,this%mesh%fnode(this%mesh%xfnode(face_index))))
          if(btest(this%mesh%cfpar(j), f)) then 
            ! Want outward oriented normal
            tmp_plane = -tmp_plane
          end if
         
          ! Create and Link this plane in the LocalizedSeparatorLink
          if(cn(f) /= 0) then
            found_internal = found_internal + 1 ! IRL is 0-based indexed
            call setPlane(this%planar_localizer(j), found_internal-1, tmp_plane(1:3), tmp_plane(4))         
            call setEdgeConnectivity(this%localized_separator_link(j), found_internal-1, &
                                     this%localized_separator_link(cn(f))) 
          else
            ! Connect to correct "outside" marking localized_separator_link
            ! Making these last in reconstruction should keep all inside domain volume inside the domain
            found_boundary = found_boundary + 1 ! IRL is 0-based index
            call setPlane(this%planar_localizer(j), total_internal + found_boundary-1, tmp_plane(1:3), tmp_plane(4))
            call setEdgeConnectivity(this%localized_separator_link(j), total_internal + found_boundary-1, &
                                     this%localized_separator_link(face_to_boundary_mapping(face_index)))
          end if
                
        end do
      end associate
      
    end do
      
    ! Make "outside" reconstructions consume all volume given to them.
    do j = this%mesh%ncell+1,total_recon_needed
      call setNumberOfPlanes(this%planar_separator(j), 1)
      call setPlane(this%planar_separator(j), 0, [0.0_r8, 0.0_r8, 0.0_r8], 1.0_r8)
      call setNumberOfPlanes(this%planar_localizer(j), 1)
      call setPlane(this%planar_localizer(j), 0, [0.0_r8, 0.0_r8, 0.0_r8], 1.0_r8)
      call setEdgeConnectivityNull(this%localized_separator_link(j), 0) ! No link back to domain
      call setId(this%localized_separator_link(j), j)
    end do

  end subroutine init_irl_mesh

  subroutine compute_velocity_weightings(this)

    class(unsplit_geometric_volume_tracker), intent(inout) :: this

    integer :: j, n, nodes_in_cell
    real(r8) :: squared_distance, magnitude

    ! Compute average length scale for each node
    this%w_node(1:2,:) = 0.0_r8
    do j = 1, this%mesh%ncell
      nodes_in_cell = this%mesh%xcnode(j+1) - this%mesh%xcnode(j)
      associate( node_id => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
        do n = 1, nodes_in_cell
          squared_distance = dot_product(this%mesh%x(:,node_id(n))-this%mesh%cell_centroid(:,j), &
                                         this%mesh%x(:,node_id(n))-this%mesh%cell_centroid(:,j))
          this%w_node(1,node_id(n)) = this%w_node(1,node_id(n)) + squared_distance
          this%w_node(2,node_id(n)) = this%w_node(2,node_id(n)) + 1.0_r8
          this%velocity_weightings(this%mesh%xcnode(j)+n-1) = squared_distance
        end do

      end associate

    end do

    ! Now normalize for  average squared distance
    this%w_node(1,:) = this%w_node(1,:) / (this%w_node(2,:) + tiny(1.0_r8))

    ! Compute weightings using previously stored squared distance and average length scale
    this%w_node(2,:) = 0.0_r8
    do j = 1, this%mesh%ncell
      nodes_in_cell = this%mesh%xcnode(j+1) - this%mesh%xcnode(j)
      associate( node_id => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
        do n = 1, nodes_in_cell
          this%velocity_weightings(this%mesh%xcnode(j)+n-1) = &
               exp(-this%velocity_weightings(this%mesh%xcnode(j)+n-1)/this%w_node(1,node_id(n)))
          this%w_node(2,node_id(n)) = this%w_node(2,node_id(n)) + this%velocity_weightings(this%mesh%xcnode(j)+n-1)
        end do
      end associate
    end do

    ! Normalize weightings to have an L2 norm of 1.
    do j = 1, this%mesh%ncell
      nodes_in_cell = this%mesh%xcnode(j+1) - this%mesh%xcnode(j)
      associate( node_id => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
        do n = 1, nodes_in_cell
          this%velocity_weightings(this%mesh%xcnode(j)+n-1) = &
               this%velocity_weightings(this%mesh%xcnode(j)+n-1) / (this%w_node(2,node_id(n)) + tiny(1.0_r8))
        end do
      end associate
    end do
    
  end subroutine compute_velocity_weightings
  
  subroutine normals(this, vof)

    use flow_operators, only: gradient_cc
    use f08_intrinsics, only: findloc
    intrinsic :: norm2

    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in)  :: vof(:,:)

    real(r8) :: mag
    integer :: i,j,k,c
    logical :: hasvof(size(vof,dim=1))

    do i = 1, this%nmat
      call gradient_cc(this%normal(1,i,:), this%normal(2,i,:), this%normal(3,i,:), &
          vof(i,:), this%w_node(1,:), this%w_node(2,:))
    end do

    do i = 1, this%mesh%ncell_onP
      hasvof = vof(:,i) > 0.0_r8
      c = count(hasvof)

      this%normal(:,:,i) = -this%normal(:,:,i)
      ! enforce consistency for two materials
      if (c == 2) then
        j = findloc(hasvof,.true.)
        k = findloc(hasvof,.true.,back=.true.)
        this%normal(:,k,i) = -this%normal(:,j,i)
      endif
      
      ! normalize and remove smallish components due to robustness issues in nested disection
      do j = 1 , this%nmat
        if (vof(j,i) < this%cutoff .or. vof(j,i) > 1.0_r8 - this%cutoff) then
          this%normal(:,j,i) = 0.0_r8
          cycle
        end if
        ! remove small values
        do k = 1, 3
          if (abs(this%normal(k,j,i)) < epsilon(1.0_r8)) then
            this%normal(k,j,i) = 0.0_r8
          end if
        end do
        ! normalize if possible
        mag = norm2(this%normal(:,j,i))
        if (mag > epsilon(1.0_r8)) then
          this%normal(:,j,i) = this%normal(:,j,i)/mag
        else
        ! QUESTION : What is the correct thing to do here? Needs to be valid normal.
          this%normal(:,j,i) = 1.0_r8 / sqrt(3.0_r8)
        end if
      end do
    end do
    ! will need normals for vof reconstruction in ghost cells
    call gather_boundary(this%mesh%cell_ip, this%normal)

  end subroutine normals
  
  subroutine set_irl_interfaces(this, vof)
    
    use cell_geom_type
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: vof(:,:)
    
    
    integer :: j
    real(r8) :: distance_guess
    
    ASSERT(this%nmat == 2) ! Will aleviate after code working with two phases
    do j = 1, this%mesh%ncell
    
      call setNumberOfPlanes(this%planar_separator(j), 1)
      if(vof(1,j) < this%cutoff .or. vof(1,j) > 1.0_r8 - this%cutoff) then
        ! Full phase cell
        distance_guess = sign(1.0_r8, vof(1,j)-0.5_r8)
        call setPlane(this%planar_separator(j),0, [0.0_r8, 0.0_r8, 0.0_r8], distance_guess)
      else
        ! Mixed cell
        distance_guess = dot_product(this%normal(:,1,j), this%mesh%cell_centroid(:,j))
        call setPlane(this%planar_separator(j), 0, this%normal(:,1,j), distance_guess)
        
        associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
          
        call this%adjust_planes_match_VOF(this%mesh%x(:,cn), vof(:,j), &
                                          this%planar_separator(j) )        
        end associate               
      
      end if
          
    end do
    
  end subroutine set_irl_interfaces
  
  subroutine adjust_planes_match_VOF(this, a_cell, a_vof, a_planar_separator)
  
    use irl_interface_helper
    use cell_geometry
  
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: a_cell(:,:)
    real(r8), intent(in) :: a_vof(:)
    type(PlanarSep_type), intent(inout) :: a_planar_separator

    real(r8), parameter :: VOF_tolerance = 1.0e-14_r8
    
    ! PlanarSeparator should already be setup with a valid normal and
    ! guess for the distance (atleast somewhere inside cell so we can calculate
    ! a gradient in IRL for the optimization)    
    select case(size(a_cell,dim=2))
      case (4) ! tet      
        call truchas_poly_to_irl(a_cell, this%IRL_tet)
        call matchVolumeFraction(this%IRL_tet, a_vof(1), a_planar_separator, VOF_tolerance)
      
      case (5) ! pyramid
        call truchas_poly_to_irl(a_cell, this%IRL_pyramid)
        call matchVolumeFraction(this%IRL_pyramid, a_vof(1), a_planar_separator, VOF_tolerance)
      
      case (6) ! Wedge
        call truchas_poly_to_irl(a_cell, this%IRL_wedge)
        call matchVolumeFraction(this%IRL_wedge, a_vof(1), a_planar_separator, VOF_tolerance)
    
      case (8) ! Hex
        call truchas_poly_to_irl(a_cell, this%IRL_hex)
        call matchVolumeFraction(this%IRL_hex, a_vof(1), a_planar_separator, VOF_tolerance)
    
      case default
        call TLS_fatal('Unknown Truchas cell type during plane-distance setting')
      end select  
    
  end subroutine adjust_planes_match_VOF

 subroutine reset_volume_moments(this, a_vof)
  
    use irl_interface_helper
    use cell_geometry
  
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(out) :: a_vof(:,:)

    integer :: j
    real(r8) :: new_volume
    
    do j = 1, this%mesh%ncell_onP
      associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
          
        select case(size(cn))
        case (4) ! tet      
          call truchas_poly_to_irl(this%mesh%x(:,cn), this%IRL_tet)
          call getNormMoments(this%IRL_tet, this%planar_separator(j), new_volume)
          a_vof(1,j) = new_volume / calculateVolume(this%IRL_tet)
          a_vof(2,j) = 1.0_r8 - a_vof(1,j)

        case (5) ! pyramid
          call truchas_poly_to_irl(this%mesh%x(:,cn), this%IRL_pyramid)
          call getNormMoments(this%IRL_pyramid, this%planar_separator(j), new_volume)
          a_vof(1,j) = new_volume / calculateVolume(this%IRL_pyramid)
          a_vof(2,j) = 1.0_r8 - a_vof(1,j)        

        case (6) ! Wedge
          call truchas_poly_to_irl(this%mesh%x(:,cn), this%IRL_wedge)
          call getNormMoments(this%IRL_wedge, this%planar_separator(j), new_volume)
          a_vof(1,j) = new_volume / calculateVolume(this%IRL_wedge)
          a_vof(2,j) = 1.0_r8 - a_vof(1,j)        

        case (8) ! Hex
          call truchas_poly_to_irl(this%mesh%x(:,cn), this%IRL_hex)
          call getNormMoments(this%IRL_hex, this%planar_separator(j), new_volume)
          a_vof(1,j) = new_volume / calculateVolume(this%IRL_hex)
          a_vof(2,j) = 1.0_r8 - a_vof(1,j)        

        case default
          call TLS_fatal('Unknown Truchas cell type during plane-distance setting')
        end select

      end associate

    end do
    
  call gather_boundary(this%mesh%cell_ip, a_vof)
    
  end subroutine reset_volume_moments
    
  subroutine compute_node_velocities(this, a_cell_centered_vel)
  
      class(unsplit_geometric_volume_tracker), intent(inout) :: this
      real(r8), intent(in) :: a_cell_centered_vel(:,:)
      
      integer :: j, nodes_in_cell, n
      
      ! Compute velocity at each node as the surface-area weighted average of 
      ! all connected face velocities.
      this%w_node = 0.0_r8
      do j = 1, this%mesh%ncell      
        nodes_in_cell = this%mesh%xcnode(j+1) - this%mesh%xcnode(j)
        associate( node_id => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))

          ! Weights already normalized, so don't need to do it again
          do n = 1, nodes_in_cell
            this%w_node(1:3,node_id(n)) = this%w_node(1:3,node_id(n)) + &
                 a_cell_centered_vel(:,j) * this%velocity_weightings(this%mesh%xcnode(j)+n-1)
          end do
        end associate
      
      end do
      
      call gather_boundary(this%mesh%node_ip, this%w_node)
      
  end subroutine compute_node_velocities

  subroutine compute_projected_nodes(this, a_dt)

    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: a_dt

    integer :: n

    do n = 1, this%mesh%nnode
      this%projected_nodes(:,n) = this%mesh%x(:,n) + this%w_node(1:3,n)*(-a_dt)
    end do

  end subroutine compute_projected_nodes
  
  subroutine compute_fluxes(this, a_face_vel, a_dt, a_vof, a_vof_band)
  
    use irl_interface_helper
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: a_face_vel(:)
    real(r8), intent(in) :: a_dt
    real(r8), intent(in) :: a_vof(:,:)
    integer, intent(in)  :: a_vof_band(:,:)
        
    integer :: f, v, i, t, k
    integer ::  number_of_nodes, neighbor_cell_index
    real(r8) :: cell_nodes(3,9), phase_volume(this%nmat), cell_volume
    type(TagAccVM_SepVol_type) :: tagged_sepvol
    integer :: current_tag, min_band
    
    ASSERT(this%nmat == 2) ! Will alleviate later

    call new(tagged_sepvol)
    call getMoments_setMethod(2)
      
    this%face_flux = 0.0_r8
    do f = 1, this%mesh%nface_onP

      neighbor_cell_index = this%mesh%fcell(1, f)
      if(this%mesh%fcell(2,f) /= 0) then
        cell_volume = min(this%mesh%volume(neighbor_cell_index), &
                          this%mesh%volume(this%mesh%fcell(2,f)))
        min_band = min(minval(abs(a_vof_band(:,neighbor_cell_index))),minval(abs(a_vof_band(:,this%mesh%fcell(2,f)))))       
      else
        cell_volume = this%mesh%volume(neighbor_cell_index)
        min_band = minval(abs(a_vof_band(:,neighbor_cell_index)))        
      end if       

      if(min_band > advect_band) then
         cycle
       end if

      call setMinimumVolToTrack(cell_volume*1.0e-15_r8)
    
      number_of_nodes = this%mesh%xfnode(f+1)-this%mesh%xfnode(f)

      ! Grab face nodes we will extrude from
      ! Will place face nodes in latter half of cell_nodes so that
      ! if (dot_product(flux_vel, face_norm) > 0), the volume is > 0
      associate( node_index => this%flux_node( &
           this%mesh%xfnode(f):this%mesh%xfnode(f+1)-1))
        ! Face nodes
        cell_nodes(:,1:number_of_nodes) = &
             this%mesh%x(:,node_index(1:number_of_nodes))

        ! Projected nodes
        cell_nodes(:,number_of_nodes+1:2*number_of_nodes) = this%projected_nodes(:, node_index(1:number_of_nodes))

        ! Initial guess for volume conservative cap vertex
        cell_nodes(:,2*number_of_nodes+1) = sum(cell_nodes(:,number_of_nodes+1:2*number_of_nodes),2)/real(number_of_nodes,r8)

        select case(this%flux_geometry_class(f))
          case(1) ! Triangular Face -> Octahedron Volume
          
            call construct(this%IRL_CapOcta_LLL, cell_nodes)
            call adjustCapToMatchVolume(this%IRL_CapOcta_LLL, a_dt*a_face_vel(f)*this%mesh%area(f))
            if(f == 7659) then
              print*,'Volume for face', f
              print*,calculateVolume(this%IRL_CapOcta_LLL)
            end if            
            call getMoments(this%IRL_CapOcta_LLL, &
                            this%localized_separator_link(neighbor_cell_index), &
                            tagged_sepvol)

          case(2) ! Triangular Face -> Octahedron Volume
          
            call construct(this%IRL_CapOcta_LLT, cell_nodes)            
            call adjustCapToMatchVolume(this%IRL_CapOcta_LLT, a_dt*a_face_vel(f)*this%mesh%area(f))
            call getMoments(this%IRL_CapOcta_LLT, &
                            this%localized_separator_link(neighbor_cell_index), &
                            tagged_sepvol)
            
          case(3) ! Triangular Face -> Octahedron Volume
          
            call construct(this%IRL_CapOcta_LTT, cell_nodes)
            call adjustCapToMatchVolume(this%IRL_CapOcta_LTT, a_dt*a_face_vel(f)*this%mesh%area(f))       
            call getMoments(this%IRL_CapOcta_LTT, &
                            this%localized_separator_link(neighbor_cell_index), &
                            tagged_sepvol)

          case(4) ! Triangular Face -> Octahedron Volume
          
            call construct(this%IRL_CapOcta_TTT, cell_nodes)            
            call adjustCapToMatchVolume(this%IRL_CapOcta_TTT, a_dt*a_face_vel(f)*this%mesh%area(f))
            call getMoments(this%IRL_CapOcta_TTT, &
                            this%localized_separator_link(neighbor_cell_index), &
                            tagged_sepvol)      
            
          case(5) ! Quad Face -> Dodecahedron Volume
            
            call construct(this%IRL_CapDod_LLLL, cell_nodes)
            call adjustCapToMatchVolume(this%IRL_CapDod_LLLL, a_dt*a_face_vel(f)*this%mesh%area(f))
            call getMoments(this%IRL_CapDod_LLLL, &
                            this%localized_separator_link(neighbor_cell_index), &
                            tagged_sepvol)

          case(6) ! Quad Face -> Dodecahedron Volume
            
            call construct(this%IRL_CapDod_LLLT, cell_nodes)            
            call adjustCapToMatchVolume(this%IRL_CapDod_LLLT, a_dt*a_face_vel(f)*this%mesh%area(f))            
            call getMoments(this%IRL_CapDod_LLLT, &
                            this%localized_separator_link(neighbor_cell_index), &
                            tagged_sepvol)

          case(7) ! Quad Face -> Dodecahedron Volume
            
            call construct(this%IRL_CapDod_LTLT, cell_nodes)            
            call adjustCapToMatchVolume(this%IRL_CapDod_LTLT, a_dt*a_face_vel(f)*this%mesh%area(f))            
            call getMoments(this%IRL_CapDod_LTLT, &
                            this%localized_separator_link(neighbor_cell_index), &
                            tagged_sepvol)

          case(8) ! Quad Face -> Dodecahedron Volume
            
            call construct(this%IRL_CapDod_LLTT, cell_nodes)            
            call adjustCapToMatchVolume(this%IRL_CapDod_LLTT, a_dt*a_face_vel(f)*this%mesh%area(f))

            call getMoments(this%IRL_CapDod_LLTT, &
                            this%localized_separator_link(neighbor_cell_index), &
                            tagged_sepvol)

          case(9) ! Quad Face -> Dodecahedron Volume
            
            call construct(this%IRL_CapDod_LTTT, cell_nodes)            
            call adjustCapToMatchVolume(this%IRL_CapDod_LTTT, a_dt*a_face_vel(f)*this%mesh%area(f))            
            call getMoments(this%IRL_CapDod_LTTT, &
                            this%localized_separator_link(neighbor_cell_index), &
                            tagged_sepvol)

          case(10) ! Quad Face -> Dodecahedron Volume
            
            call construct(this%IRL_CapDod_TTTT, cell_nodes)            
            call adjustCapToMatchVolume(this%IRL_CapDod_TTTT, a_dt*a_face_vel(f)*this%mesh%area(f))            
            call getMoments(this%IRL_CapDod_TTTT, &
                            this%localized_separator_link(neighbor_cell_index), &
                            tagged_sepvol)

          case default
            call TLS_fatal('Unknown face type. How was this not found in this%generate_flux_classes()?')          
        end select
          
      end associate

      phase_volume = 0.0_r8
      do t = 0, getSize(tagged_sepvol)-1
         current_tag = getTagForIndex(tagged_sepvol, t)
         if(current_tag > this%mesh%ncell) then
            ! Is a boundary condition, all volume in phase 0
            phase_volume(:) = &
                 phase_volume(:) + this%getBCMaterialFractions(current_tag, a_vof) * getVolumeAtIndex(tagged_sepvol, t, 0)
          else
            ! Is inside domain, trust actual volumes
            phase_volume(1) = phase_volume(1) + getVolumeAtIndex(tagged_sepvol, t, 0)
            phase_volume(2) = phase_volume(2) + getVolumeAtIndex(tagged_sepvol, t, 1)
         end if
      end do
      this%face_flux(:,f) = phase_volume

      if(abs(a_dt*a_face_vel(f)*this%mesh%area(f) - sum(phase_volume)) > 1.0e-13 .and. &
         abs(a_dt*a_face_vel(f)*this%mesh%area(f)) > epsilon(1.0_r8)  ) then
        print*,'Face failed', f, this%mesh%fcell(:,f)
        print*,'Case: ', this%flux_geometry_class(f)
        print*,a_dt*a_face_vel(f)*this%mesh%area(f), sum(phase_volume)
        print*,phase_volume
      end if
      
    end do
    
  end subroutine compute_fluxes

  pure function getBCMaterialFractions(this, a_tag, a_vof_n) result(a_fractions)

    class(unsplit_geometric_volume_tracker), intent(in) :: this
    integer, intent(in) :: a_tag
    real(r8), intent(in) :: a_vof_n(:,:)
    real(r8) :: a_fractions(this%nmat)

    a_fractions = 0.0_r8
    if(this%inflow_mat(a_tag - this%mesh%ncell) > 0) then
       ! Known material coming in
       a_fractions(this%inflow_mat(a_tag - this%mesh%ncell)) = 1.0_r8
    else
       ! Take as Neumann on neighbor cell inside domain
       a_fractions(:) = a_vof_n(:,this%boundary_recon_to_cell(a_tag - this%mesh%ncell))
    end if

  end function getBCMaterialFractions
  
  subroutine update_vof(this, a_vof_band, a_flux_vol, a_vof)
  
    class(unsplit_geometric_volume_tracker), intent(in) :: this
    integer, intent(in) :: a_vof_band(:,:)
    real(r8), intent(inout) :: a_flux_vol(:,:)
    real(r8), intent(inout) :: a_vof(:,:)
    
    integer :: j, f, k
    integer :: number_of_faces
    real(r8) :: cell_volume_sum
    logical :: vof_modified

    real(r8) :: tmp_sum(this%nmat), tmp_sum_prev(this%nmat), summed_removed(this%nmat)
    
    ! Fill the ragged a_flux_vol array
    a_flux_vol = 0.0_r8
    do j = 1, this%mesh%ncell  
      if(minval(abs(a_vof_band(:,j))) > advect_band) then
        cycle
      end if
      associate( fn => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1) )
      do f = 1, size(fn)
        if(btest(this%mesh%cfpar(j), f)) then
          a_flux_vol(:,this%mesh%xcface(j)+f-1) = -this%face_flux(:,fn(f))
        else
          a_flux_vol(:,this%mesh%xcface(j)+f-1) =  this%face_flux(:,fn(f))
        end if
      end do
            
      end associate
      
    end do

    tmp_sum = 0.0_r8
    do j = 1, this%mesh%ncell_onP
      tmp_sum = tmp_sum + a_vof(:,j)*this%mesh%volume(j)
    end do
    print*,'Before update', tmp_sum
    tmp_sum_prev = tmp_sum
    
    ! Use a_flux_vol to update a_vof now
    summed_removed = 0.0_r8
    do j = 1, this%mesh%ncell_onP
      if(minval(abs(a_vof_band(:,j))) > advect_band) then
        cycle
      end if      
      a_vof(:,j) = a_vof(:,j) * this%mesh%volume(j)
      number_of_faces = this%mesh%xcface(j+1) - this%mesh%xcface(j)
      cell_volume_sum = this%mesh%volume(j)
      do f = 1, number_of_faces
        a_vof(:,j) = a_vof(:,j) - a_flux_vol(:, this%mesh%xcface(j)+f-1)
        cell_volume_sum = cell_volume_sum - sum(a_flux_vol(:, this%mesh%xcface(j)+f-1))
      end do

      ! DEBUG
      if(abs(this%mesh%volume(j) - cell_volume_sum) .gt. 1.0e-14*this%mesh%volume(j)) then
        print*,'Unconservative fluxing in volume!!', j, this%mesh%volume(j), cell_volume_sum
        print*,a_vof(:,j)
        print*,sum(abs(a_vof(:,j)))
      end if
      
      a_vof(:,j) = a_vof(:,j) / cell_volume_sum
      vof_modified = .false.
      do k = 1, this%nmat
        if(a_vof(k,j) < this%cutoff .and. abs(a_vof(k,j)) > EPSILON(1.0_r8)) then
          vof_modified = .true.
          summed_removed(k) = summed_removed(k) - a_vof(k,j)
          a_vof(k,j) = 0.0_r8
          ! TODO Also update IRL planar_separator plane to reflect this
          if(k==1) then
            call setPlane(this%planar_separator(j),0,[0.0_r8, 0.0_r8, 0.0_r8],-1.0_r8)
          end if
          if(k==2) then
            call setPlane(this%planar_separator(j),0,[0.0_r8, 0.0_r8, 0.0_r8], 1.0_r8)
          end if          
        else if(a_vof(k,j) > 1.0_r8 - this%cutoff .and. abs(a_vof(k,j)) < 1.0_r8 - EPSILON(1.0_r8)) then
          vof_modified = .true.
          summed_removed(k) = summed_removed(k) + (1.0_r8 - a_vof(k,j))
          a_vof(k,j) = 1.0_r8
          ! TODO Also update IRL planar_separator plane to reflect this
          if(k==1) then
            call setPlane(this%planar_separator(j),0,[0.0_r8, 0.0_r8, 0.0_r8], 1.0_r8)
          end if
          if(k==2) then
            call setPlane(this%planar_separator(j),0,[0.0_r8, 0.0_r8, 0.0_r8],-1.0_r8)
          end if                    
        end if
      end do
    
      if(vof_modified) then
        ! Rescale VOF to be valid (sum to 1)
        a_vof(:,j) = a_vof(:,j) / sum(a_vof(:,j))
      end if
    end do

    call gather_boundary(this%mesh%cell_ip, a_vof)       

    tmp_sum = 0.0_r8
    do j = 1, this%mesh%ncell_onP
      tmp_sum = tmp_sum + a_vof(:,j)*this%mesh%volume(j)
    end do
    print*,'After update', tmp_sum
    print*,'Change ', tmp_sum - tmp_sum_prev
    print*,'Clipped ', summed_removed

  end subroutine update_vof    

  ! Routine to generate consistently approximated flux volumes
  !
  ! This routine makes heavy use of lookup tables, which is why it
  ! is down here hidden from the above prettier code. These
  ! lookup tables will be used to reorder face node indices and
  ! classify the necessary flux volume geometry class to use
  ! from IRL. By doing this, a completely conformal mesh should
  ! be formed during the back-advection phase, with no gaps
  ! or overlapping regions.
  !
  ! The integer Ids for this%flux_geometry_class correspond to:
  ! Octa_LLL ( 1)
  ! Octa_LLT ( 2)
  ! Octa_LTT ( 3)
  ! Octa_TTT ( 4)
  ! Dod_LLLL ( 5)
  ! Dod_LLLT ( 6)
  ! Dod_LTLT ( 7)
  ! Dod_LLTT ( 8)
  ! Dod_LTTT ( 9)
  ! Dod_TTTT (10)
  !
  ! The lookup tables will encode this information, as well as
  ! the number of shifts of vertices from the mesh%fnode order
  ! to the consistent flux_node order, using lookup tables
  ! based on an ID by accounting for the edge having a leading
  ! diagonalization (L) as a 0 (unset) bit, and a trailing
  ! diagonalization(T) as a 1 (set) bit. For example,
  ! the case of LTLT is 0+2+0+8 = 10.
  !
  ! The diagonalization of each edge will be set from the faces
  ! on a first-come, first-served basis. Each face in
  ! this%mesh is looped over, and if the diagonalization of the
  ! edge has not yet been set, will be set to L (arbitrary choice).
  ! Note, setting L for this edge will make it a T for the opposite
  ! direction traversal of the edge.
  
  subroutine generate_flux_classes(this)

    class(unsplit_geometric_volume_tracker), intent(inout) :: this

    integer, parameter :: &
         Octa_geometry_class(8) = &
         [1, 2, 2, 3, &
          2, 3, 3, 4]

    integer, parameter :: &
         Octa_shift(8) = &
         [0, 2, 1, 1, &
          0, 2, 0, 0]
    
    integer, parameter :: &
         Dod_geometry_class(16) = &
         [5, 6, 6, 8, &
          6, 7, 8, 9, &
          6, 8, 7, 9, &
          8, 9, 9, 10]

    integer, parameter :: &
         Dod_shift(16) = &
         [0, 1, 2, 2, &
          1, 1, 3, 1, &
          0, 3, 0, 2,&
          0, 3, 0, 0]
    integer, allocatable :: edges(:,:) ! Temporary working array. Would be better as hash_map
    integer :: f
    integer :: n, number_of_nodes    
    integer :: edge_end
    integer :: lookup_case
    integer :: shift_amount, r, tmp_node, reordered_vertices(4)
    integer :: irl_ordering(4)
    
    allocate(this%flux_geometry_class(this%mesh%nface))
    allocate(this%flux_node(this%mesh%xfnode(this%mesh%nface+1)-1))
    allocate(edges(this%mesh%nnode,this%mesh%nnode))


    ! Now loop through domain and set the edge orientations
    ! on first-come, first-served basis
    edges = -1 ! Sentinel value to mark as not-yet assigned.
    do f = 1, this%mesh%nface
      associate(node_id => this%mesh%fnode(this%mesh%xfnode(f):this%mesh%xfnode(f+1)-1))
        irl_ordering = reorder_node_id(node_id)
        number_of_nodes = size(node_id)
        do n = 1, number_of_nodes
          edge_end = irl_ordering(mod(n,number_of_nodes)+1)
          if(edges(irl_ordering(n),edge_end) == -1) then
            ! Not yet set, so then set it.
            edges(irl_ordering(n), edge_end) = 0 ! This direction leading
            edges(edge_end, irl_ordering(n)) = 1 ! Opposite direction trailing
          end if
          
        end do        
        
      end associate
    end do

    ! For each face, now store correct flux class and
    ! node ordering to be consistent
    do f = 1, this%mesh%nface
      associate(node_id => this%mesh%fnode(this%mesh%xfnode(f):this%mesh%xfnode(f+1)-1))
        irl_ordering = reorder_node_id(node_id)        
        number_of_nodes = size(node_id)
        ASSERT(number_of_nodes == 3 .or. number_of_nodes == 4)
        lookup_case = 1
        do n = 1, number_of_nodes          
          edge_end = irl_ordering(mod(n,number_of_nodes)+1)
          ASSERT(edges(irl_ordering(n),edge_end) == 0 .or. &
                 edges(irl_ordering(n),edge_end) == 1)
          lookup_case = lookup_case + edges(irl_ordering(n), edge_end) * 2**(n-1)
        end do

        select case(number_of_nodes)
          case(3) ! Triangle -> Octahedron
            ASSERT(lookup_case > 0 .and. lookup_case <= 8)          
            this%flux_geometry_class(f) = Octa_geometry_class(lookup_case)
            shift_amount = Octa_shift(lookup_case)

          case(4) ! Quad -> Dodecahedron
            ASSERT(lookup_case > 0 .and. lookup_case <= 16)
            this%flux_geometry_class(f) = Dod_geometry_class(lookup_case)
            shift_amount = Dod_shift(lookup_case)
            
          case default
            call TLS_fatal('Cannot create flux geometry face with nodes!=[3,4]')
        end select

        reordered_vertices(1:number_of_nodes) = irl_ordering(1:number_of_nodes)
        do r = 1, shift_amount
          tmp_node = reordered_vertices(number_of_nodes)
          do n = number_of_nodes, 2, -1
            reordered_vertices(n) = reordered_vertices(n-1)
          end do
          reordered_vertices(1) = tmp_node
        end do

        this%flux_node(this%mesh%xfnode(f):this%mesh%xfnode(f+1)-1) = reordered_vertices(1:number_of_nodes)
        
      end associate

    end do
    
    deallocate(edges)

  contains

    function reorder_node_id(a_node_list) result(a_irl_ordered_list)
      integer, intent(in) :: a_node_list(:)
      integer :: a_irl_ordered_list(4)

      ASSERT(size(a_node_list) == 3 .or. size(a_node_list) == 4)
      
      if(size(a_node_list) == 3) then
        a_irl_ordered_list = [a_node_list(2), a_node_list(3), a_node_list(1), -1]
      end if

      if(size(a_node_list) == 4) then
        a_irl_ordered_list = [a_node_list(2), a_node_list(3), a_node_list(4), a_node_list(1)]        
      end if
      
      return
    end function reorder_node_id

  end subroutine generate_flux_classes

end module unsplit_geometric_volume_tracker_type
