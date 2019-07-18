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
  use parameter_module, only : string_len
  implicit none
  private

  integer, parameter, private :: advect_band = 5 ! Size of band to do advection

  type, extends(volume_tracker), public :: unsplit_geometric_volume_tracker
    private
    type(unstr_mesh), pointer :: mesh ! unowned reference
    integer :: location_iter_max ! maximum number of iterations to use in fitting interface
    integer :: subcycles
    logical :: nested_dissection
    character(string_len) :: interface_reconstruction_name
    real(r8), allocatable :: normal(:,:,:)
    real(r8), allocatable :: cell_flux(:,:)
    real(r8), allocatable :: phase_centroids(:,:,:)
    ! node/face/cell workspace
    real(r8), allocatable :: w_node(:,:)
    real(r8), allocatable :: node_velocity(:,:)
    real(r8), allocatable :: velocity_weightings(:)
    real(r8), allocatable :: lsq_weightings(:,:)
    real(r8), allocatable :: projected_nodes(:,:)
    integer, allocatable :: priority(:), bc_index(:), local_face(:), inflow_mat(:)
    integer :: nrealfluid, nfluid, nmat ! # of non-void fluids, # of fluids incl. void, # of materials
    integer, allocatable :: boundary_recon_to_cell(:) ! Mapping of boundary reconstruction ID to neighboring inside cell ID
    real(r8), allocatable :: cell_bounding_box(:,:,:)
    real(r8), allocatable :: correction_vertex(:,:)
    ! Start IRL objects
    type(ObjServer_PlanarLoc_type) :: object_server_planar_localizer
    type(ObjServer_PlanarSep_type) :: object_server_planar_separator
    type(ObjServer_LocSepLink_type) :: object_server_localized_separator_link
    type(PlanarLoc_type), allocatable :: planar_localizer(:)
    type(PlanarSep_type), allocatable :: planar_separator(:)
    type(LocSepLink_type), allocatable :: localized_separator_link(:)
    type(Poly_type), allocatable :: interface_polygons(:) 
    type(Tet_type) :: IRL_tet
    type(Pyrmd_type) :: IRL_pyramid
    type(TriPrism_type) :: IRL_wedge
    type(Hex_type) :: IRL_hex
    type(SymTet_type) :: IRL_sym_tet
    type(SymPyrmd_type) :: IRL_sym_pyramid
    type(SymTriPrism_type) :: IRL_sym_wedge
    type(SymHex_type) :: IRL_sym_hex    
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
    procedure, public :: write_interface
    procedure, private :: init_irl_mesh
    procedure, private :: generate_flux_classes
    procedure, private :: generate_cell_bounding_boxes
    procedure, private :: compute_velocity_weightings
    procedure, private :: compute_lsq_weightings
    procedure, private :: set_irl_interfaces
    procedure, private :: adjust_planes_match_VOF
    procedure, private :: construct_interface_polygons
    procedure, private :: construct_nonzero_polygon
    procedure, private :: compute_projected_nodes
    procedure, private :: get_velocity_at_point
    procedure, private :: compute_effective_cfl
    procedure, private :: compute_correction_vertices
    procedure, private :: compute_fluxes
    procedure, private :: set_irl_volume_geometry
    procedure, private :: refined_advection_check
    procedure, private :: moments_from_geometric_cutting
    procedure, private :: getBCMaterialFractions    
    procedure, private :: interface_reconstruction    
    procedure, private :: normals_youngs
    procedure, private :: normals_newgrad_youngs    
    procedure, private :: normals_swartz
    procedure, private :: normals_mof    
    procedure, private :: normals_lvira
    procedure, private :: normals_lvira_hex_execute
    procedure, private :: normals_lvira_tet_execute    
    procedure, private :: identify_reconstruction_neighborhood
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
    character(:), allocatable :: interface_recon
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
    call TLS_warn('Subcycling is currently disabled for unsplit transport')
    call params%get('interface_reconstruction', interface_recon, default='Youngs')
    this%interface_reconstruction_name = interface_recon
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
    ASSERT(this%nmat == 2)
    allocate(this%cell_flux(this%nmat, this%mesh%ncell))
    allocate(this%phase_centroids(3, this%nmat, this%mesh%ncell))
    do j = 1, this%mesh%ncell
      do i = 1, this%nmat
        this%phase_centroids(:,i,j) = this%mesh%cell_centroid(:,j)
      end do
    end do

    allocate(this%normal(3,this%nmat,mesh%ncell))
    allocate(this%w_node(3,mesh%nnode))
    allocate(this%node_velocity(3,mesh%nnode))    
    allocate(this%velocity_weightings(mesh%xndcell(mesh%nnode_onP+1)-1))
    allocate(this%lsq_weightings(3,mesh%xcnc(mesh%ncell_onP+1)-1))
    allocate(this%projected_nodes(3, mesh%nnode))
    allocate(this%cell_bounding_box(3,2,this%mesh%ncell))
    allocate(this%correction_vertex(3,this%mesh%nface))

    ! list of boundary face ids
    j = count(this%mesh%fcell(:,:this%mesh%nface) == 0)
    allocate(this%bc_index(j), this%local_face(j), this%inflow_mat(j))
    j = 1
    do i = 1, this%mesh%nface
      if (this%mesh%fcell(2,i) == 0) then
        this%bc_index(j) = i
        k = this%mesh%fcell(1,i)
        this%local_face(j) = &
            findloc(this%mesh%cface(this%mesh%xcface(k):this%mesh%xcface(k+1)-1), i)
        j = j + 1
     else if(this%mesh%fcell(1,i) == 0) then
        this%bc_index(j) = i
        k = this%mesh%fcell(2,i)
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
    call new(this%IRL_sym_tet)
    call new(this%IRL_sym_pyramid)
    call new(this%IRL_sym_wedge)
    call new(this%IRL_sym_hex)    
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

    ! Initialize and allocate polygons   
    allocate(this%interface_polygons(this%mesh%ncell))
    do j = 1, this%mesh%ncell
      call new(this%interface_polygons(j))
    end do

    ! Reorganize and generate needed face-flux information
    call this%generate_flux_classes

    ! Generate cell bounding boxes
    call this%generate_cell_bounding_boxes

    ! Calculate weightings for interpolation of cell center velocity to node
    call this%compute_velocity_weightings

    ! Calculate weights for a LSQ gradient based on face-neighbors
    call this%compute_lsq_weightings

  end subroutine init

  ! flux volumes routine assuming flux_vol is a cface-like array
  ! flux volumes routine assuming vel is stored on faces, length 1:mesh%nfaces
  subroutine flux_volumes(this, vel, vel_cc, vof_n, vof, flux_vol, fluids, void, dt, a_interface_band)
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: vel(:), vel_cc(:,:), vof_n(:,:), dt
    real(r8), intent(out) :: flux_vol(:,:), vof(:,:)
    integer, intent(in) :: fluids, void
    integer, intent(in) :: a_interface_band(:)

    integer :: i,j    

    ! DEBUGGING
    integer :: f
    real(r8) :: tmp

    flux_vol = 0.0_r8
    vof = vof_n

    call start_timer('reconstruction')
    call this%interface_reconstruction(vof)    
    call this%construct_interface_polygons(a_interface_band)
    call stop_timer('reconstruction')

    call start_timer('advection')
    call this%compute_projected_nodes(vel_cc, dt)    
    call this%compute_effective_cfl(dt)

    call this%compute_correction_vertices(vel, dt, a_interface_band)

    call this%compute_fluxes(vel, dt, vof_n, vof, a_interface_band)

    call stop_timer('advection')

  end subroutine flux_volumes

  subroutine write_interface(this, t, dt, cycle_number)

    use truchas_phase_interface_output, only : TPIO_write_mesh
    
    class(unsplit_geometric_volume_tracker), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: dt
    integer, intent(in) :: cycle_number
    
    call TPIO_write_mesh(this%interface_polygons, t, dt, cycle_number)

    return
  end subroutine write_interface

  
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

    face_to_boundary_mapping = -10000
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

    integer :: j, n
    real(r8) :: squared_distance

    ! Node weightings
    this%w_node(1,:) = 0.0_r8
    do n = 1, this%mesh%nnode_onP
      associate(cn => this%mesh%ndcell(this%mesh%xndcell(n):this%mesh%xndcell(n+1)-1))
        do j = 1, size(cn)
          squared_distance = sum((this%mesh%x(:,n)-this%mesh%cell_centroid(:,cn(j)))**2)
          this%velocity_weightings(this%mesh%xndcell(n)+j-1) = 1.0_r8 / sqrt(squared_distance)
        end do
      end associate
      ! Normalize weighting
      this%velocity_weightings(this%mesh%xndcell(n):this%mesh%xndcell(n+1)-1) = &
           this%velocity_weightings(this%mesh%xndcell(n):this%mesh%xndcell(n+1)-1) / &
           sum(this%velocity_weightings(this%mesh%xndcell(n):this%mesh%xndcell(n+1)-1))      
    end do
    
  end subroutine compute_velocity_weightings

  subroutine compute_lsq_weightings(this)

    class(unsplit_geometric_volume_tracker), intent(inout) :: this

    integer :: j, n, neighbors, q, r
    real(r8) :: squared_distance, weight, dx(3)
    real(r8) :: A(3,3), Ainv(3,3), determinant
    integer :: max_stencil_size
    real(r8), allocatable :: B(:,:)

    max_stencil_size = -1
    do j = 1, this%mesh%ncell_onP
      max_stencil_size = max(max_stencil_size, this%mesh%xcnc(j+1)-this%mesh%xcnc(j))
    end do
    allocate(B(3,max_stencil_size))
    
    ! Compute average length scale for each node
    do j = 1, this%mesh%ncell_onP
      neighbors = this%mesh%xcnc(j+1) - this%mesh%xcnc(j)
      associate( cn => this%mesh%cnc(this%mesh%xcnc(j):this%mesh%xcnc(j+1)-1))

        A = 0.0_r8
        B = 0.0_r8
        do n = 1, neighbors
          if(cn(n) /= 0) then
             squared_distance = dot_product(this%mesh%cell_centroid(:,cn(n))-this%mesh%cell_centroid(:,j), &
                                            this%mesh%cell_centroid(:,cn(n))-this%mesh%cell_centroid(:,j))
             weight = 1.0_r8 / sqrt(squared_distance)
             dx = this%mesh%cell_centroid(:,cn(n)) - this%mesh%cell_centroid(:,j)
             do q = 1, 3
               do r = 1,3
                 A(q,r) = A(q,r) + weight * dx(q)*dx(r)
               end do
             end do
             do q = 1, 3
               B(q,n) = weight*dx(q)
             end do
          end if
        end do     
        
        determinant = A(1,1) * (A(3,3) * A(2,2) - A(2,3)**2) &
             - A(2,1) * (A(3,3) * A(1,2) - A(2,3) * A(1,3)) &
             + A(1,3) * (A(2,3) * A(1,2) - A(2,2) * A(1,3))
        ASSERT(determinant > 0.0_r8)
        Ainv(1,1) = A(3,3) * A(2,2) - A(2,3)**2  
        Ainv(1,2)  = A(1,3) * A(2,3) - A(3,3) * A(1,2)  
        Ainv(1,3)  = A(1,2) * A(2,3) - A(1,3) * A(2,2)  
        
        Ainv(2,2) = A(3,3) * A(1,1) - A(1,3)**2  
        Ainv(2,3) =  A(1,2) * A(1,3) - A(1,1) * A(2,3)
        
        Ainv(3,3) = A(1,1) * A(2,2) - A(1,2)**2
      
        Ainv(2,1) = Ainv(1,2)
        Ainv(3,1) = Ainv(1,3)
        Ainv(3,2) = Ainv(2,3)
        
        determinant = (A(1,1) * Ainv(1,1)) + (A(1,2) * Ainv(1,2)) + (A(1,3) * Ainv(1,3))
        ASSERT(determinant > 0.0_r8)
        Ainv = Ainv / determinant
        
        do n = 1, neighbors
          do q = 1, 3            
            this%lsq_weightings(q,this%mesh%xcnc(j)+n-1) = dot_product(Ainv(q,:),B(:,n))
          end do
        end do

      end associate

    end do

    deallocate(B)
    
  end subroutine compute_lsq_weightings
  
  subroutine interface_reconstruction(this, vof)

    use parallel_communication, only : global_sum
    use parameter_module, only : string_len

    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in)  :: vof(:,:)

    integer :: j,nrecon
    character(string_len) :: message        

    nrecon = 0
    do j = 1, this%mesh%ncell_onP
      if(vof(1,j) < this%cutoff .or. vof(1,j) > 1.0_r8 - this%cutoff) then
         cycle
      end if
      nrecon = nrecon + 1
    end do
    write(message, '(a,i10)') 'Reconstructions being performed', global_sum(nrecon)
    call TLS_info(message)
    
    select case(trim(this%interface_reconstruction_name))

      case('Youngs')
         call this%normals_youngs(vof)

      case('NewGrad Youngs')
         call this%normals_newgrad_youngs(vof)
         
      case('Swartz')
        call this%normals_youngs(vof)
        call this%normals_swartz(vof)

     case('MOF')
        ! Don't have valid centroids on first iteration
        call TLS_Fatal('Need to think of what to do on first iteration when no centroids.')
        call this%normals_mof(vof)

      case('LVIRA')
        call this%normals_youngs(vof)
        call this%set_irl_interfaces(vof)
        call this%normals_lvira(vof)
        
      case default
        call TLS_fatal('Unknown reconstruction type requested for interface_reconstruction')
      end select

    call this%set_irl_interfaces(vof)            
    
  end subroutine interface_reconstruction
  
  subroutine normals_youngs(this, vof)    

    use flow_operators, only: gradient_cc
    use f08_intrinsics, only: findloc
    intrinsic :: norm2

    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in)  :: vof(:,:)

    real(r8) :: mag
    integer :: i,k,c
    logical :: hasvof(size(vof,dim=1))


   call gradient_cc(this%normal(1,1,:), this%normal(2,1,:), this%normal(3,1,:), &
       vof(1,:), this%w_node(1,:), this%w_node(2,:))    

    this%normal(:,1,1:this%mesh%ncell_onP) = -this%normal(:,1,1:this%mesh%ncell_onP)
    do i = 1, this%mesh%ncell_onP      
      
      ! normalize and remove smallish components due to robustness issues in nested disection
      if (vof(1,i) < this%cutoff .or. vof(1,i) > 1.0_r8 - this%cutoff) then
         this%normal(:,1,i) = 0.0_r8
         cycle
      end if
      ! remove small values
      do k = 1, 3
        if (abs(this%normal(k,1,i)) < epsilon(1.0_r8)) then
           this%normal(k,1,i) = 0.0_r8
        end if
      end do
      ! normalize if possible
      mag = norm2(this%normal(:,1,i))
      if (mag > epsilon(1.0_r8)) then
         this%normal(:,1,i) = this%normal(:,1,i)/mag
      else
         ! QUESTION : What is the correct thing to do here? Needs to be valid normal.
         this%normal(:,1,i) = 1.0_r8 / sqrt(3.0_r8)
      end if
    end do

    this%normal(:,2,1:this%mesh%ncell_onP) = -this%normal(:,1,1:this%mesh%ncell_onP)
    call gather_boundary(this%mesh%cell_ip, this%normal)
    
  end subroutine normals_youngs

  subroutine normals_newgrad_youngs(this, vof)    

    use flow_operators, only: gradient_cc
    use f08_intrinsics, only: findloc
    intrinsic :: norm2

    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in)  :: vof(:,:)

    real(r8) :: mag
    integer :: i,k,c
    logical :: hasvof(size(vof,dim=1))

    this%cell_flux = vof
    !call filter_vof(this%mesh,this%cell_flux(:,:), 1)
    do i = 1, this%mesh%ncell_onP      
      this%normal(:,1,i) = 0.0_r8
      if (vof(1,i) >= this%cutoff .and. vof(1,i) <= 1.0_r8 - this%cutoff) then
         ! Need a normal for reconstruction
         associate( cn => this%mesh%cnc(this%mesh%xcnc(i):this%mesh%xcnc(i+1)-1))      
           do c = 1, size(cn)
             if(cn(c) /= 0) then
                this%normal(:,1,i) = this%normal(:,1,i) + this%lsq_weightings(:,this%mesh%xcnc(i)+c-1)&
                     *(this%cell_flux(1,cn(c))-this%cell_flux(1,i))
             end if
           end do
         end associate
      end if
    end do    

    this%normal(:,1,1:this%mesh%ncell_onP) = -this%normal(:,1,1:this%mesh%ncell_onP)      
    do i = 1, this%mesh%ncell_onP      
      
      ! Already addressed, set to zero, just cycle
      if (vof(1,i) < this%cutoff .or. vof(1,i) > 1.0_r8 - this%cutoff) then
         cycle
      end if
      ! remove small values
      do k = 1, 3
        if (abs(this%normal(k,1,i)) < epsilon(1.0_r8)) then
           this%normal(k,1,i) = 0.0_r8
        end if
      end do
      ! normalize if possible
      mag = norm2(this%normal(:,1,i))
      if (mag > epsilon(1.0_r8)) then
         this%normal(:,1,i) = this%normal(:,1,i)/mag
      else
         ! QUESTION : What is the correct thing to do here? Needs to be valid normal.
         this%normal(:,1,i) = 1.0_r8 / sqrt(3.0_r8)
      end if
    end do

    this%normal(:,2,1:this%mesh%ncell_onP) = -this%normal(:,1,1:this%mesh%ncell_onP)
    call gather_boundary(this%mesh%cell_ip, this%normal)

  contains

    subroutine filter_vof(a_mesh, a_vof, a_iter)

      type(unstr_mesh), intent(in) :: a_mesh
      real(r8), intent(inout) :: a_vof(:,:)
      integer, intent(in) :: a_iter

      integer :: j, n, next, last, active
      real(r8) :: vol_sum      

      next = 2
      last = 1
      do i = 1, a_iter
        do j = 1, a_mesh%ncell_onP
          associate( cn => this%mesh%cnhbr(this%mesh%xcnhbr(j):this%mesh%xcnhbr(j+1)-1))
            a_vof(next,j) = 0.0_r8
            vol_sum = 0.0_r8
            do n = 1, size(cn)
              if(cn(n) /= 0) then
                 vol_sum = vol_sum + a_mesh%volume(cn(n))
                 a_vof(next,j) = a_vof(next,j) + a_vof(last, cn(n))*a_mesh%volume(cn(n))
              end if
            end do
            a_vof(next,j) = a_vof(next,j) / (vol_sum + tiny(1.0_r8))
          end associate
        end do
        call gather_boundary(a_mesh%cell_ip, a_vof(next,:))          
        next = 2-(next-1)
        last = 2-(last-1)
      end do      

      if(last == 2) then
         a_vof(1,:) = a_vof(2,:)
      end if
      
    end subroutine filter_vof
    
  end subroutine normals_newgrad_youngs
  
  subroutine normals_swartz(this, a_vof)

    use irl_interface_helper
    use constants_module, only : pi
    use cell_geometry, only: cross_product
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in)  :: a_vof(:,:)    

    integer, parameter :: max_iter = 4
    
    logical :: done
    integer :: j, n, outer_iter
    real(r8), allocatable :: polygon_center(:,:), initial_normal(:,:), old_normal(:,:)
    real(r8) :: distance_guess
    real(r8) :: angle_difference, rotation_axis(3), normal_center_line(3)
    real(r8) :: rotation_quaternion(4), inv_rotation_quaternion(4), rotated_normal(4)
    real(r8) :: rotation_axis_mag
    integer :: number_of_active_neighbors

    ASSERT(this%nmat == 2)
    
    allocate(polygon_center(3,this%mesh%ncell))
    allocate(initial_normal(3,this%mesh%ncell))
    allocate(old_normal(3,this%mesh%ncell))    
    initial_normal = this%normal(:,1,:)
    done = .false.
    outer_iter = 0
    do while(.not. done)
      outer_iter = outer_iter + 1

      ! Set current interfaces in IRL (with satisfting volume fractions) and polygons
      polygon_center = 0.0_r8      
      do j = 1, this%mesh%ncell_onP
        if(a_vof(1,j) > this%cutoff .and. a_vof(1,j) < 1.0_r8 - this%cutoff) then
          distance_guess = dot_product(this%normal(:,1,j), this%mesh%cell_centroid(:,j))
          call setNumberOfPlanes(this%planar_separator(j), 1)
          call setPlane(this%planar_separator(j), 0, this%normal(:,1,j), distance_guess)
          
          associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
            
            call this%adjust_planes_match_VOF(this%mesh%x(:,cn), a_vof(:,j), &
                 this%planar_separator(j) )

            call this%construct_nonzero_polygon(this%mesh%x(:,cn), this%planar_separator(j), &
                 0, this%interface_polygons(j))
          end associate

          do n = 1, getNumberOfVertices(this%interface_polygons(j))
            polygon_center(:,j) = polygon_center(:,j) + getPt(this%interface_polygons(j), n-1)
          end do
          polygon_center(:,j) = polygon_center(:,j) / real(getNumberOfVertices(this%interface_polygons(j)), r8)         
        end if
      end do

      call gather_boundary(this%mesh%cell_ip, polygon_center)
      call gather_boundary(this%mesh%cell_ip, this%normal(:,1,:))

      ! Now do Jacobi type update on the normal orientation guesses.
      old_normal = this%normal(:,1,:)
      this%normal(:,1,:) = 0.0_r8
      do j = 1, this%mesh%ncell_onP
        if(a_vof(1,j) > this%cutoff .and. a_vof(1,j) < 1.0_r8 - this%cutoff) then

          associate(cn => this%mesh%cnhbr(this%mesh%xcnhbr(j):this%mesh%xcnhbr(j+1)-1))
            number_of_active_neighbors = 0
            do n = 1, size(cn)
              if(cn(n) /= 0) then
                if(a_vof(1,cn(n)) > this%cutoff .and. a_vof(1,cn(n)) < 1.0_r8 - this%cutoff) then
                  angle_difference = acos(max(-1.0_r8,min(1.0_r8, &
                       dot_product(old_normal(:,j), old_normal(:,cn(n))))))
                  if(abs(angle_difference) < pi/6.0_r8) then
                    ! Use this neighbor
                    normal_center_line = polygon_center(:,cn(n))-polygon_center(:,j)
                    normal_center_line = normal_center_line / sqrt(sum(normal_center_line**2))
                    if(abs(dot_product(normal_center_line, old_normal(:,j))) < 1.0_r8 - 1.0e-4_r8) then 
                      rotation_axis = cross_product(normal_center_line, old_normal(:,j))
                      rotation_axis = rotation_axis / sqrt(sum(rotation_axis**2))
                    else                      
                      cycle
                    end if
                    number_of_active_neighbors = number_of_active_neighbors + 1                    
                    rotation_quaternion(1) = cos(0.5_r8*0.5_r8*pi)
                    rotation_quaternion(2:4) = rotation_axis * sin(0.5_r8*0.5_r8*pi)
                    inv_rotation_quaternion(1) = rotation_quaternion(1)
                    inv_rotation_quaternion(2:4) = -rotation_quaternion(2:4)
                    ! Product of rotation_quaternion and normal_center_line
                    rotated_normal(1) = -dot_product(rotation_quaternion(2:4),normal_center_line)
                    rotated_normal(2:4) = rotation_quaternion(1)*normal_center_line + &
                         cross_product(rotation_quaternion(2:4),normal_center_line)
                    ! Product of (rotated_quaternion*normal_center_line) and inv_rotation_quaternion
                    rotation_quaternion = rotated_normal
                    rotated_normal(1) = rotation_quaternion(1)*inv_rotation_quaternion(1) &
                         -dot_product(rotation_quaternion(2:4),inv_rotation_quaternion(2:4))
                    rotated_normal(2:4) = rotation_quaternion(1)*inv_rotation_quaternion(2:4) &
                         + inv_rotation_quaternion(1)*rotation_quaternion(2:4) &
                         + cross_product(rotation_quaternion(2:4),inv_rotation_quaternion(2:4))
                    rotated_normal(2:4) = rotated_normal(2:4) / sqrt(sum(rotated_normal(2:4)**2))
                    this%normal(:,1,j) = this%normal(:,1,j) + rotated_normal(2:4)
                  end if                  
                end if
              end if

            end do

          end associate

          if(number_of_active_neighbors > 0) then
            this%normal(:,1,j) = this%normal(:,1,j) / sqrt(sum(this%normal(:,1,j)**2))
          else
            this%normal(:,1,j) = old_normal(:,j)
          end if

        end if

      end do
      
      ! Exit criteria
      if(outer_iter == max_iter) then
        done = .true.
      end if
      
    end do

    this%normal(:,2,:) = -this%normal(:,1,:)

    call gather_boundary(this%mesh%cell_ip, this%normal)
    
    deallocate(polygon_center)
    deallocate(initial_normal)
    deallocate(old_normal)
    
  end subroutine normals_swartz

  subroutine normals_mof(this, a_vof)

    use irl_interface_helper
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this

    real(r8), intent(in) :: a_vof(:,:)

    integer :: j
    type(SepVM_type) :: sepvm
    real(r8) :: tmp_plane(4)

    call new(sepvm)

    do j = 1, this%mesh%ncell_onP
      if(a_vof(1,j) < this%cutoff .or. a_vof(1,j) > 1.0_r8 - this%cutoff) then
        cycle
      end if

      call construct(sepvm, [this%cell_flux(1,j), this%phase_centroids(1:3,1,j), &
                             this%cell_flux(2,j), this%phase_centroids(1:3,2,j)])
      
      associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))          
        select case(size(cn))
        case (4) ! tet      
          call truchas_poly_to_irl(this%mesh%x(:,cn), this%IRL_tet)
          call reconstructMOF3D(this%IRL_tet, sepvm, this%planar_separator(j))

        case (5) ! pyramid
          call truchas_poly_to_irl(this%mesh%x(:,cn), this%IRL_pyramid)
          call TLS_FATAL('Not yet implemented') 
          !call reconstructMOF3D(this%IRL_pyramid, sepvm, this%planar_separator(j))

        case (6) ! Wedge
          call truchas_poly_to_irl(this%mesh%x(:,cn), this%IRL_wedge)
          call TLS_FATAL('Not yet implemented') 
          !call reconstructMOF3D(this%IRL_wedge, sepvm, this%planar_separator(j))
          
        case (8) ! Hex
          call truchas_poly_to_irl(this%mesh%x(:,cn), this%IRL_hex)
          call reconstructMOF3D(this%IRL_hex, sepvm, this%planar_separator(j))
          
        case default
          call TLS_fatal('Unknown Truchas cell type during MOF reconstruction.')
        end select

      end associate

      tmp_plane = getPlane(this%planar_separator(j), 0)
      this%normal(:,1,j) = tmp_plane(1:3)
      
    end do    

    this%normal(:,2,:) = -this%normal(:,1,:)
    call gather_boundary(this%mesh%cell_ip, this%normal)
    
  end subroutine normals_mof
    
  subroutine normals_lvira(this, a_vof)

    use traversal_tracker_type, only : traversal_tracker
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: a_vof(:,:)

    integer :: j
    integer :: cell_types(4), stencil_type_case
    type(traversal_tracker) :: stencil

    ASSERT(this%nmat == 2)

    do j = 1, this%mesh%ncell_onP
      if(a_vof(1,j) < this%cutoff .or. a_vof(1,j) > 1.0_r8 - this%cutoff) then
        cycle
      end if
      call this%identify_reconstruction_neighborhood(j, 1.5_r8*(this%mesh%volume(j)**(1.0_r8/3.0_r8)), &
           cell_types, stencil)
      ASSERT(sum(cell_types) == stencil%get_number_of_visited_cells())
      
      if(sum(cell_types) == maxval(cell_types)) then
        ! Homogeneous stencil
        stencil_type_case = maxloc(cell_types,1)
      end if
      select case(stencil_type_case)
        case(1) ! Homogeneous tet mesh
          call this%normals_lvira_tet_execute(j, stencil, a_vof)          
          
        case(2) ! Homogeneous pyramid mesh
          call TLS_fatal('Not yet implemented')
          
        case(3) ! Homogeneous wedge mesh
          call TLS_fatal('Not yet implemented')
          
        case(4) ! Homogeneous hex  mesh
          call this%normals_lvira_hex_execute(j, stencil, a_vof)
          
        case  default ! Heterogeneous mesh
          call TLS_fatal('Not yet implemented')
      end select
      
    end do

    this%normal(:,2,:) = -this%normal(:,1,:)
    call gather_boundary(this%mesh%cell_ip, this%normal)

  end subroutine normals_lvira

  subroutine normals_lvira_tet_execute(this, a_center_index, a_stencil, a_vof)

    use irl_interface_helper
    use traversal_tracker_type, only : traversal_tracker    
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    integer, intent(in) :: a_center_index
    type(traversal_tracker), intent(in) :: a_stencil
    real(r8), intent(in) :: a_vof(:,:)

    integer :: n, cell_index
    type(LVIRANeigh_Tet_type) :: neighborhood
    type(Tet_type), allocatable :: IRL_tet(:)
    real(r8) :: tmp_plane(4)
    
    ! Create IRL cell geometries
    allocate(IRL_tet(a_stencil%get_number_of_visited_cells()))   

    ! Build IRL stencil
    call new(neighborhood)
    call setSize(neighborhood, a_stencil%get_number_of_visited_cells())
    do n = 1, a_stencil%get_number_of_visited_cells()      
      call new(IRL_tet(n))
      cell_index = a_stencil%get_visited_cell_index(n)
      associate (cn => this%mesh%cnode(this%mesh%xcnode(cell_index):this%mesh%xcnode(cell_index+1)-1))      
        call truchas_poly_to_irl(this%mesh%x(:,cn), IRL_tet(n))
      end associate
      call setMember(neighborhood, n-1, IRL_tet(n), a_vof(1,cell_index))
    end do
    call setCenterOfStencil(neighborhood, 0)

    ! Perform LVIRA
    call reconstructLVIRA3D(neighborhood, this%planar_separator(a_center_index))   

    ! Move normal to this%normal
    tmp_plane = getPlane(this%planar_separator(a_center_index),0)
    this%normal(:,1,a_center_index) = tmp_plane(1:3)

    deallocate(IRL_tet)

  end subroutine normals_lvira_tet_execute

  
  subroutine normals_lvira_hex_execute(this, a_center_index, a_stencil, a_vof)

    use irl_interface_helper
    use traversal_tracker_type, only : traversal_tracker    
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    integer, intent(in) :: a_center_index
    type(traversal_tracker), intent(in) :: a_stencil
    real(r8), intent(in) :: a_vof(:,:)

    integer :: n, cell_index
    type(LVIRANeigh_Hex_type) :: neighborhood
    type(Hex_type), allocatable :: IRL_hex(:)
    real(r8) :: tmp_plane(4)
    
    ! Create IRL cell geometries
    allocate(IRL_hex(a_stencil%get_number_of_visited_cells()))   

    ! Build IRL stencil
    call new(neighborhood)
    call setSize(neighborhood, a_stencil%get_number_of_visited_cells())
    do n = 1, a_stencil%get_number_of_visited_cells()      
      call new(IRL_hex(n))
      cell_index = a_stencil%get_visited_cell_index(n)
      associate (cn => this%mesh%cnode(this%mesh%xcnode(cell_index):this%mesh%xcnode(cell_index+1)-1))      
        call truchas_poly_to_irl(this%mesh%x(:,cn), IRL_hex(n))
      end associate
      call setMember(neighborhood, n-1, IRL_hex(n), a_vof(1,cell_index))
    end do
    call setCenterOfStencil(neighborhood, 0)

    ! Perform LVIRA
    call reconstructLVIRA3D(neighborhood, this%planar_separator(a_center_index))   

    ! Move normal to this%normal
    tmp_plane = getPlane(this%planar_separator(a_center_index),0)
    this%normal(:,1,a_center_index) = tmp_plane(1:3)

    deallocate(IRL_hex)

  end subroutine normals_lvira_hex_execute
    

  subroutine identify_reconstruction_neighborhood(this, a_starting_index, a_stencil_length, &
       a_cell_types, a_stencil)

    use traversal_tracker_type, only : traversal_tracker    

    class(unsplit_geometric_volume_tracker), intent(in) :: this
    integer, intent(in) :: a_starting_index
    real(r8), intent(in) :: a_stencil_length
    integer, intent(out) :: a_cell_types(:)
    type(traversal_tracker), intent(out) :: a_stencil

    integer :: current_cell, n
    real(r8) :: distance_to_centroid

    a_cell_types = 0
    call a_stencil%init([a_starting_index])

    ! Add starting_cell to cell_types
    select case(this%mesh%xcnode(a_starting_index+1)-this%mesh%xcnode(a_starting_index))
    case (4) ! tet      
      a_cell_types(1) = a_cell_types(1) + 1                
    case (5) ! pyramid
      a_cell_types(2) = a_cell_types(2) + 1                
    case (6) ! Wedge
      a_cell_types(3) = a_cell_types(3) + 1                
    case (8) ! Hex
      a_cell_types(4) = a_cell_types(4) + 1
    case default
      call TLS_fatal('Unknown Truchas cell type during stencil generation')
    end select

    do while(a_stencil%still_cells_to_visit())
      current_cell = a_stencil%get_next_cell()
      associate( cn => this%mesh%cnhbr(this%mesh%xcnhbr(current_cell):this%mesh%xcnhbr(current_cell+1)-1))
        do n = 1, size(cn)
          if(cn(n) /= 0) then
             if(a_stencil%cell_not_encountered(cn(n))) then
                distance_to_centroid = sqrt(sum((&
                     this%mesh%cell_centroid(:,cn(n))-this%mesh%cell_centroid(:,a_starting_index))**2))
                if(distance_to_centroid < a_stencil_length) then
                   call a_stencil%add_cell(cn(n))
                   select case(this%mesh%xcnode(cn(n)+1)-this%mesh%xcnode(cn(n)))
                   case (4) ! tet      
                      a_cell_types(1) = a_cell_types(1) + 1                
                   case (5) ! pyramid
                      a_cell_types(2) = a_cell_types(2) + 1                
                   case (6) ! Wedge
                      a_cell_types(3) = a_cell_types(3) + 1                
                   case (8) ! Hex
                      a_cell_types(4) = a_cell_types(4) + 1
                   case default
                      call TLS_fatal('Unknown Truchas cell type during stencil generation')
                   end select
                end if
             end if
          end if
        end do
      end associate
    end do

  end subroutine identify_reconstruction_neighborhood
  
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
    
  subroutine construct_interface_polygons(this, a_interface_band)

    use irl_interface_helper    
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    integer, intent(in) :: a_interface_band(:)

    integer :: j

    do j = 1, this%mesh%ncell
      if(a_interface_band(j) /= 0) then
        call zeroPolygon(this%interface_polygons(j))
        cycle
      end if
      associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))

        call this%construct_nonzero_polygon(this%mesh%x(:,cn), this%planar_separator(j),&
             0, this%interface_polygons(j))
        
      end associate
    end do


  end subroutine construct_interface_polygons

  subroutine construct_nonzero_polygon(this, a_cell, a_planar_separator, a_plane_index, a_polygon)

    use irl_interface_helper    

    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: a_cell(:,:)
    type(PlanarSep_type), intent(in) :: a_planar_separator
    integer, intent(in) :: a_plane_index
    type(Poly_type), intent(inout) :: a_polygon
    
    select case(size(a_cell,2))
    case (4) ! tet
      call truchas_poly_to_irl(a_cell, this%IRL_tet)
      call getPoly(this%IRL_tet, a_planar_separator, a_plane_index, a_polygon)          
    case (5) ! pyramid
      call truchas_poly_to_irl(a_cell, this%IRL_pyramid)          
      call TLS_fatal('Not yet implemented')
    case (6) ! Wedge
      call truchas_poly_to_irl(a_cell, this%IRL_wedge)          
      call TLS_fatal('Not yet implemented')
    case (8) ! Hex
      call truchas_poly_to_irl(a_cell, this%IRL_hex)          
      call getPoly(this%IRL_Hex, a_planar_separator, a_plane_index, a_polygon)
    case default
      call TLS_fatal('Unknown Truchas cell type during polygon construction setting')
    end select
    
  end subroutine construct_nonzero_polygon
      
  
  subroutine compute_projected_nodes(this, a_cell_centered_vel, a_dt)
  
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: a_cell_centered_vel(:,:)
    real(r8), intent(in) :: a_dt
    
    integer :: j, n, iter

    integer :: cell_id
    integer, parameter :: max_iter = 5
    real(r8) :: old_dist, current_dist, vel(3), projected_vel(3), v(3,4)
    real(r8), parameter :: tolerance = 1.0e-4_r8 ! 1% change in distance
    
    ! Compute velocity at each node as the surface-area weighted average of 
    ! all connected face velocities.
    this%node_velocity = 0.0_r8
    do n = 1, this%mesh%nnode_onP
      associate( cn => this%mesh%ndcell(this%mesh%xndcell(n):this%mesh%xndcell(n+1)-1))        
        ! Weights already normalized, so don't need to do it again
        do j = 1, size(cn)
          this%node_velocity(:,n) = this%node_velocity(:,n) + &
               a_cell_centered_vel(:,cn(j)) * this%velocity_weightings(this%mesh%xndcell(n)+j-1)
        end do
      end associate
      
    end do
    call gather_boundary(this%mesh%node_ip, this%node_velocity)
    
    
    ! First order explicit Euler
     this%projected_nodes(:,1:this%mesh%nnode_onP) = this%mesh%x(:,1:this%mesh%nnode_onP) + &
          this%node_velocity(:,1:this%mesh%nnode_onP)*(-a_dt)
     this%w_node(:,1:this%mesh%nnode_onP) = (this%mesh%x(:,1:this%mesh%nnode_onP) - &
          this%projected_nodes(:,1:this%mesh%nnode_onP)) &
          / a_dt
     call gather_boundary(this%mesh%node_ip, this%w_node )         
     call gather_boundary(this%mesh%node_ip, this%projected_nodes )
    
    ! ! Second Order (Spatially) Simplectic Integration
    ! node_loop : do n = 1, this%mesh%nnode_onP
    !   this%projected_nodes(:,n) = this%mesh%x(:,n)
    !   cell_id = this%mesh%ndcell(this%mesh%xndcell(n))
    !   integrate : do iter = 1, max_iter
    !     vel = this%get_velocity_at_point(0.5_r8*(this%projected_nodes(:,n)+this%mesh%x(:,n)), cell_id)
    !     if(cell_id > this%mesh%ncell) then ! Outside domain, use Forward Euler
    !        this%projected_nodes(:,n) = this%mesh%x(:,n) - a_dt*this%node_velocity(:,n)
    !        cycle node_loop
    !     end if
    !     this%projected_nodes(:,n) = this%mesh%x(:,n) - a_dt*vel
    !     current_dist = sum((this%projected_nodes(:,n) - this%mesh%x(:,n))**2)
    !     if(iter .ne. 1) then
    !        if( (current_dist - old_dist)/(old_dist+tiny(1.0_r8)) < tolerance) then
    !           cycle node_loop
    !        end if
    !     end if
    !     old_dist = current_dist
    !   end do integrate
    ! end do node_loop

    ! ! Store in w_node the linear velocity that would create the projected node in 1 step.
    ! this%w_node(:,1:this%mesh%nnode_onP) = (this%mesh%x(:,1:this%mesh%nnode_onP) - &
    !      this%projected_nodes(:,1:this%mesh%nnode_onP)) &
    !      / a_dt    
    ! call gather_boundary(this%mesh%node_ip, this%w_node )
    ! call gather_boundary(this%mesh%node_ip, this%projected_nodes )

    ! Explicit RK4
    ! node_loop : do n = 1, this%mesh%nnode_onP
    !   cell_id = this%mesh%ndcell(this%mesh%xndcell(n))
    !   v(:,1)=this%get_velocity_at_point(this%mesh%x(:,n), cell_id)
    !   if(cell_id > this%mesh%ncell) then ! Outside domain, use Forward Euler
    !      this%projected_nodes(:,n) = this%mesh%x(:,n) - a_dt*this%node_velocity(:,n)
    !      cycle node_loop
    !   end if
    !   v(:,2)=this%get_velocity_at_point(this%mesh%x(:,n)+0.5_r8*-a_dt*v(:,1),cell_id)
    !   if(cell_id > this%mesh%ncell) then ! Outside domain, use Forward Euler
    !      this%projected_nodes(:,n) = this%mesh%x(:,n) - a_dt*this%node_velocity(:,n)
    !      cycle node_loop
    !   end if
    !   v(:,3)=this%get_velocity_at_point(this%mesh%x(:,n)+0.5_r8*-a_dt*v(:,2),cell_id)
    !   if(cell_id > this%mesh%ncell) then ! Outside domain, use Forward Euler
    !      this%projected_nodes(:,n) = this%mesh%x(:,n) - a_dt*this%node_velocity(:,n)
    !      cycle node_loop
    !   end if
    !   v(:,4)=this%get_velocity_at_point(this%mesh%x(:,n)+       -a_dt*v(:,3),cell_id)
    !   if(cell_id > this%mesh%ncell) then ! Outside domain, use Forward Euler
    !      this%projected_nodes(:,n) = this%mesh%x(:,n) - a_dt*this%node_velocity(:,n)
    !      cycle node_loop
    !   end if
    !   this%projected_nodes(:,n) = this%mesh%x(:,n)-a_dt/6.0_r8*(v(:,1)+2.0_r8*v(:,2)+2.0_r8*v(:,3)+v(:,4))    
    ! end do node_loop

    ! ! Store in w_node the linear velocity that would create the projected node in 1 step.
    ! this%w_node(:,1:this%mesh%nnode_onP) = (this%mesh%x(:,1:this%mesh%nnode_onP) - &
    !      this%projected_nodes(:,1:this%mesh%nnode_onP)) &
    !      / a_dt    
    ! call gather_boundary(this%mesh%node_ip, this%w_node )
    ! call gather_boundary(this%mesh%node_ip, this%projected_nodes )   
    
      
  end subroutine compute_projected_nodes

  function get_velocity_at_point(this, a_pt, a_identified_cell) result(a_vel)

    class(unsplit_geometric_volume_tracker), intent(in) :: this
    real(r8), intent(in) :: a_pt(3)
    integer, intent(inout) :: a_identified_cell
    real(r8) :: a_vel(3)

    integer :: n, nodes_in_cell
    real(r8) :: characteristic_length, weights(8)

    ! Note, when passed in, a_identified_cell should be guess for which cell the point is in.
    ! When function returns, it will be the actual cell the point is located in.
    a_identified_cell =  locatePt(a_pt, this%localized_separator_link(a_identified_cell))
    if(a_identified_cell > this%mesh%ncell) then ! Projected outside domain
       a_vel = 0.0_r8
       return
    end if
    
    nodes_in_cell = this%mesh%xcnode(a_identified_cell+1) - this%mesh%xcnode(a_identified_cell)
    associate( node_id => this%mesh%cnode(this%mesh%xcnode(a_identified_cell):this%mesh%xcnode(a_identified_cell+1)-1))
      characteristic_length = 0.0_r8
      do n = 1, nodes_in_cell
        weights(n) = sum((a_pt-this%mesh%x(:,node_id(n)))**2)
      end do
      characteristic_length = sum(weights(1:nodes_in_cell)) / real(nodes_in_cell, r8)

      do n = 1, nodes_in_cell
        weights(n) =  sqrt(1.0 / (weights(n)/characteristic_length + tiny(1.0_r8)))
      end do
      weights(1:nodes_in_cell) = weights(1:nodes_in_cell) / sum(weights(1:nodes_in_cell))

      a_vel = 0.0_r8
      do n = 1, nodes_in_cell
        a_vel = a_vel + weights(n)*this%node_velocity(:,node_id(n))
      end do

    end associate

    return
  end function get_velocity_at_point

  subroutine compute_effective_cfl(this, a_dt)

    use parallel_communication, only : global_sum, global_maxval
    use parameter_module, only : string_len

    class(unsplit_geometric_volume_tracker), intent(in) :: this
    real(r8), intent(in) :: a_dt

    integer :: f, n
    integer :: number_of_nodes, edge_end
    real(r8) :: limiting_cfl
    real(r8) :: edge_length, vel1, vel2
    integer :: number_of_crossed_faces
    real(r8) :: projected_sum
    real(r8) :: max_vel,  ivey_cfl
    character(string_len) :: message    
    
    
    limiting_cfl = -huge(1.0_r8)
    ivey_cfl = -huge(1.0_r8)
    number_of_crossed_faces = 0
    do f = 1, this%mesh%nface_onP

      associate(node_id => this%mesh%fnode(this%mesh%xfnode(f):this%mesh%xfnode(f+1)-1))
        number_of_nodes = size(node_id)
        do n = 1, number_of_nodes
          edge_end = node_id(mod(n,number_of_nodes)+1)
          edge_length = sqrt(dot_product(this%mesh%x(:,node_id(n))-this%mesh%x(:,edge_end), &
               this%mesh%x(:,node_id(n))-this%mesh%x(:,edge_end)))
          ! Want velocity magnitude in edge direction          
          vel1 = dot_product(this%w_node(1:3,node_id(n)), (this%mesh%x(:,edge_end)-this%mesh%x(:,node_id(n)))/edge_length)
          vel2 = dot_product(this%w_node(1:3,edge_end), (this%mesh%x(:,node_id(n))-this%mesh%x(:,edge_end))/edge_length)
          max_vel = max(sqrt(dot_product(this%w_node(1:3,node_id(n)),this%w_node(1:3,node_id(n)))), &
               sqrt(dot_product(this%w_node(1:3,edge_end),this%w_node(1:3,edge_end))))
          ivey_cfl = max(ivey_cfl, max_vel*a_dt/edge_length)
          ! This is a necessary but not sufficient criteria for  detecting
          ! "hour-glass modes" that will be created during projection          
          projected_sum =  (max(vel1,0.0_r8)+max(vel2,0.0_r8))*a_dt/edge_length
          if(projected_sum > 1.0_r8) then
           number_of_crossed_faces = number_of_crossed_faces + 1
          end if

          ! This is some measure of CFL based on edge-parallel vertex travel
          limiting_cfl = max(limiting_cfl, max(abs(vel1), abs(vel2))*a_dt/edge_length)
        end do
      end associate
      
    end do

    write(message, '(a,es12.4)') 'Effective Unsplit CFL', global_maxval(limiting_cfl)
    call TLS_info(message)
    write(message, '(a,es12.4)') 'Ivey CFL', global_maxval(ivey_cfl)
    call TLS_info(message)
    number_of_crossed_faces =  global_sum(number_of_crossed_faces)
    if(number_of_crossed_faces > 0) then
       write(message, '(i4,a)') number_of_crossed_faces, 'face crossings might exist, could lose discrete conservation.'
       call TLS_warn(message)
    end if
      
  end subroutine compute_effective_cfl

  subroutine compute_correction_vertices(this, a_face_vel, a_dt, a_interface_band)

    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: a_face_vel(:)
    real(r8), intent(in) :: a_dt
    integer, intent(in)  :: a_interface_band(:)
    
    integer :: f, min_band,n
    integer :: number_of_nodes
    real(r8) :: correct_volume, cell_nodes(3,9)

    do f = 1, this%mesh%nface_onP

      if(this%mesh%fcell(2,f) /= 0) then
        min_band = min(abs(a_interface_band(this%mesh%fcell(1,f))),abs(a_interface_band(this%mesh%fcell(2,f))))
      else
        min_band = abs(a_interface_band(this%mesh%fcell(1,f)))
      end if       

      if(min_band > advect_band) then
         cycle
      end if

      correct_volume = a_dt*a_face_vel(f)*this%mesh%area(f)

      number_of_nodes = this%mesh%xfnode(f+1)-this%mesh%xfnode(f)
      
      associate( node_index => this%flux_node( &
           this%mesh%xfnode(f):this%mesh%xfnode(f+1)-1))
        ! Face nodes
        cell_nodes(:,1:number_of_nodes) = this%mesh%x(:,node_index)
        
        ! Projected nodes
        cell_nodes(:,number_of_nodes+1:2*number_of_nodes) = this%projected_nodes(:, node_index)
        
        ! Initial guess for volume conservative cap vertex
        cell_nodes(:,2*number_of_nodes+1) = sum(cell_nodes(:,number_of_nodes+1:2*number_of_nodes),2)/real(number_of_nodes,r8)
        
        select case(this%flux_geometry_class(f))
        case(1) ! Triangular Face -> Octahedron Volume
           
           call construct(this%IRL_CapOcta_LLL, cell_nodes)
           call adjustCapToMatchVolume(this%IRL_CapOcta_LLL, correct_volume)
           this%correction_vertex(:,f) = getPt(this%IRL_CapOcta_LLL, 6)
           
        case(2) ! Triangular Face -> Octahedron Volume
           
           call construct(this%IRL_CapOcta_LLT, cell_nodes)
           call adjustCapToMatchVolume(this%IRL_CapOcta_LLT, correct_volume)
           this%correction_vertex(:,f) = getPt(this%IRL_CapOcta_LLT, 6)
           
        case(3) ! Triangular Face -> Octahedron Volume
           
           call construct(this%IRL_CapOcta_LTT, cell_nodes)
           call adjustCapToMatchVolume(this%IRL_CapOcta_LTT, correct_volume)
           this%correction_vertex(:,f) = getPt(this%IRL_CapOcta_LTT, 6)
           
        case(4) ! Triangular Face -> Octahedron Volume
           
           call construct(this%IRL_CapOcta_TTT, cell_nodes)            
           call adjustCapToMatchVolume(this%IRL_CapOcta_TTT, correct_volume)
           this%correction_vertex(:,f) = getPt(this%IRL_CapOcta_TTT, 6)         
           
        case(5) ! Quad Face -> Dodecahedron Volume
           
           call construct(this%IRL_CapDod_LLLL, cell_nodes)
           call adjustCapToMatchVolume(this%IRL_CapDod_LLLL, correct_volume)
           this%correction_vertex(:,f) = getPt(this%IRL_CapDod_LLLL, 8)
           
        case(6) ! Quad Face -> Dodecahedron Volume
           
           call construct(this%IRL_CapDod_LLLT, cell_nodes)            
           call adjustCapToMatchVolume(this%IRL_CapDod_LLLT, correct_volume)
           this%correction_vertex(:,f) = getPt(this%IRL_CapDod_LLLT, 8)         
           
        case(7) ! Quad Face -> Dodecahedron Volume
           
           call construct(this%IRL_CapDod_LTLT, cell_nodes)            
           call adjustCapToMatchVolume(this%IRL_CapDod_LTLT, correct_volume)
           this%correction_vertex(:,f) = getPt(this%IRL_CapDod_LTLT, 8)         
           
        case(8) ! Quad Face -> Dodecahedron Volume
           
           call construct(this%IRL_CapDod_LLTT, cell_nodes)            
           call adjustCapToMatchVolume(this%IRL_CapDod_LLTT, correct_volume)
           this%correction_vertex(:,f) = getPt(this%IRL_CapDod_LLTT, 8)      
           
        case(9) ! Quad Face -> Dodecahedron Volume
           
           call construct(this%IRL_CapDod_LTTT, cell_nodes)            
           call adjustCapToMatchVolume(this%IRL_CapDod_LTTT, correct_volume)
           this%correction_vertex(:,f) = getPt(this%IRL_CapDod_LTTT, 8)         
           
        case(10) ! Quad Face -> Dodecahedron Volume
           
           call construct(this%IRL_CapDod_TTTT, cell_nodes)            
           call adjustCapToMatchVolume(this%IRL_CapDod_TTTT, correct_volume)
           this%correction_vertex(:,f) = getPt(this%IRL_CapDod_TTTT, 8)         
           
        case default
           call TLS_fatal('Unknown face type. How was this not found in this%generate_flux_classes()?')          
        end select
        
      end associate
      
    end do

    call gather_boundary(this%mesh%face_ip, this%correction_vertex)
    
  end subroutine compute_correction_vertices

  subroutine compute_fluxes(this, a_face_vel, a_dt, a_old_vof, a_new_vof, a_interface_band)
  
    use irl_interface_helper
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: a_face_vel(:)
    real(r8), intent(in) :: a_dt
    real(r8), intent(in) :: a_old_vof(:,:)
    real(r8), intent(inout) :: a_new_vof(:,:)    
    integer, intent(in)  :: a_interface_band(:)
        
    integer :: j, k, f, cell_ind
    integer ::  number_of_nodes, neighbor_cell_index
    real(r8) :: cell_nodes(3,14)
    real(r8) :: bounding_box(3,2)
    real(r8) :: moment_flux(8)    
    logical :: geometric_cutting_needed
    
    ASSERT(this%nmat == 2) ! Will alleviate later

    call getMoments_setMethod(1)
      
    this%cell_flux = 0.0_r8
    do j = 1, this%mesh%ncell_onP     

      if(abs(a_interface_band(j)) > advect_band) then
         cycle
      end if
    
      number_of_nodes = this%mesh%xcnode(j+1)-this%mesh%xcnode(j)

      associate( node_index => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
        
        ! Projected cell nodes
        cell_nodes(:,1:number_of_nodes) = this%projected_nodes(:,node_index)
        
        associate(face_index => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
          do f = 1, size(face_index)
            cell_nodes(:,number_of_nodes+f) = this%correction_vertex(:,face_index(f))
          end do
        end associate

        call this%set_irl_volume_geometry(cell_nodes, number_of_nodes, bounding_box)
        geometric_cutting_needed = this%refined_advection_check(bounding_box, &
             j, a_interface_band(j), a_interface_band)
      end associate

      if(geometric_cutting_needed) then
         ! Perform cutting
         call setMinimumVolToTrack(this%mesh%volume(j)*1.0e-15_r8)         
         call this%moments_from_geometric_cutting(number_of_nodes, &
              j, a_old_vof, moment_flux)
         this%cell_flux(1:2,j) = moment_flux(1:2)
         this%phase_centroids(:,1,j) = moment_flux(3:5)
         this%phase_centroids(:,2,j) = moment_flux(6:8)
         ! Update VOF
         ! Only time vof should change is when geometric cutting is needed.
         a_new_vof(:,j) = this%cell_flux(:,j) / sum(this%cell_flux(:,j))         
         if(a_new_vof(1,j) < this%cutoff) then
            this%cell_flux(1,j) = 0.0_r8
            this%cell_flux(2,j) = this%mesh%volume(j)
            a_new_vof(1,j) = 0.0_r8
            a_new_vof(2,j) = 1.0_r8
            this%phase_centroids(:,1,j) = this%mesh%cell_centroid(:,j)
            this%phase_centroids(:,2,j) = this%mesh%cell_centroid(:,j)            
         else if(a_new_vof(1,j) > 1.0_r8 - this%cutoff) then
            this%cell_flux(1,j) = this%mesh%volume(j)
            this%cell_flux(2,j) = 0.0_r8            
            a_new_vof(1,j) = 1.0_r8
            a_new_vof(2,j) = 0.0_r8
            this%phase_centroids(:,1,j) = this%mesh%cell_centroid(:,j)
            this%phase_centroids(:,2,j) = this%mesh%cell_centroid(:,j)            
         else
            if(trim(this%interface_reconstruction_name) == 'MOF') then
               cell_ind = j
               this%phase_centroids(:,1,j) = this%phase_centroids(:,1,j) + &
                    this%get_velocity_at_point(this%phase_centroids(:,1,j), cell_ind)*a_dt
               this%phase_centroids(:,2,j) = this%phase_centroids(:,2,j) + &
                    this%get_velocity_at_point(this%phase_centroids(:,2,j), cell_ind)*a_dt               
            end if
         end if
      end if
    end do
    
    ! Need to communicate fluxes
    call gather_boundary(this%mesh%cell_ip, a_new_vof)

    ! Don't need to communicate for MoF, but might use for other things later
    call gather_boundary(this%mesh%cell_ip, this%phase_centroids(:,1,:))
    call gather_boundary(this%mesh%cell_ip, this%phase_centroids(:,2,:))    
    
  end subroutine compute_fluxes

  subroutine set_irl_volume_geometry(this, a_cell, a_base_number_of_nodes, a_bounding_box)

    use irl_interface_helper
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this    
    real(r8), intent(in) :: a_cell(:,:)
    integer, intent(in) :: a_base_number_of_nodes
    real(r8), intent(out) :: a_bounding_box(3,2)

    integer :: n

    select case(a_base_number_of_nodes)
    case(4) ! Projected Tet
       call truchas_poly_to_irl(a_cell, this%IRL_sym_tet)
       call getBoundingPts(this%IRL_sym_tet, a_bounding_box(:,1), a_bounding_box(:,2))
       
    case(5) ! Projected Pyramid
       call truchas_poly_to_irl(a_cell, this%IRL_sym_pyramid)
       call getBoundingPts(this%IRL_sym_pyramid, a_bounding_box(:,1), a_bounding_box(:,2))       
       
    case(6) ! Projected Wedge
       call truchas_poly_to_irl(a_cell, this%IRL_sym_wedge)
       call getBoundingPts(this%IRL_sym_wedge, a_bounding_box(:,1), a_bounding_box(:,2))       
       
    case(8) ! Projected Hex
       call truchas_poly_to_irl(a_cell, this%IRL_sym_hex)
       call getBoundingPts(this%IRL_sym_hex, a_bounding_box(:,1), a_bounding_box(:,2))       
       
    case default
       call TLS_fatal('Unknown projected cell type in set_irl_volume_geometry.')          
    end select

  end subroutine set_irl_volume_geometry

  function refined_advection_check(this, a_bounding_box, a_starting_cell, a_starting_band, a_interface_band) &
       result(a_geometric_cutting_needed)

    use traversal_tracker_type, only : traversal_tracker        

    class(unsplit_geometric_volume_tracker), intent(in) :: this
    real(r8), intent(in) :: a_bounding_box(:,:)
    integer, intent(in) :: a_starting_cell
    integer, intent(in) :: a_starting_band
    integer, intent(in) :: a_interface_band(:)
    logical :: a_geometric_cutting_needed

    type(traversal_tracker) :: coverage
    integer :: current_cell, n, dim

    if(a_starting_band == 0) then
       a_geometric_cutting_needed = .true.
       return
    end if

    call coverage%init([a_starting_cell])

    do while(coverage%still_cells_to_visit())
      current_cell = coverage%get_next_cell()

      associate( cn => this%mesh%cnhbr(this%mesh%xcnhbr(current_cell):this%mesh%xcnhbr(current_cell+1)-1))
        neigh_loop : do n = 1, size(cn)
          if(cn(n) /= 0) then
             if(coverage%cell_not_encountered(cn(n))) then
                ! Now check if intersection of axis aligned bounding boxes
                ! X direction
                do dim = 1,3
                  if(this%cell_bounding_box(dim,1,cn(n)) > a_bounding_box(dim,2) .or. &
                       this%cell_bounding_box(dim,2,cn(n)) < a_bounding_box(dim,1)) then
                     cycle neigh_loop
                  end if
                end do

                ! Made it this far, bounding_box's intersect. See if would require cutting
                if(a_interface_band(cn(n))*a_starting_band <= 0) then
                   a_geometric_cutting_needed = .true.
                   return
                end if
                call coverage%add_cell(cn(n))
             end if
          end if
        end do neigh_loop
      end associate
    end do
    
    a_geometric_cutting_needed = .false.
    return
    
  end function refined_advection_check

  subroutine moments_from_geometric_cutting(this, a_base_number_of_nodes, a_starting_index, a_vof, a_cell_flux)

    class(unsplit_geometric_volume_tracker), intent(in) :: this
    integer, intent(in) :: a_base_number_of_nodes
    integer, intent(in) :: a_starting_index
    real(r8), intent(in) :: a_vof(:,:)
    real(r8), intent(inout) :: a_cell_flux(:)

    integer :: t
    type(TagAccVM_SepVM_type) :: tagged_sepvm
    type(SepVM_type) :: sepvm
    integer :: current_tag
    real(r8) :: boundary_vof(this%nmat)

    ! Note, the correct object type for a_geom_class should already
    ! be correctly construction in this%

    call new(tagged_sepvm)
      
    select case(a_base_number_of_nodes)
    case(4) ! Projected Tet

       call getMoments(this%IRL_sym_tet, &
            this%localized_separator_link(a_starting_index), &
            tagged_sepvm)
       
    case(5) ! Projected Pyramid

       call getMoments(this%IRL_sym_pyramid, &
            this%localized_separator_link(a_starting_index), &
            tagged_sepvm)       
       
    case(6) ! Projected Wedge

       call getMoments(this%IRL_sym_wedge, &
            this%localized_separator_link(a_starting_index), &
            tagged_sepvm)       
       
    case(8) ! Projected Hex
       
       call getMoments(this%IRL_sym_hex, &
            this%localized_separator_link(a_starting_index), &
            tagged_sepvm)
       
    case default
       call TLS_fatal('Unknown projected cell type in moments_from_geometric_cutting.')          
    end select
       
    a_cell_flux = 0.0_r8
    do t = 0, getSize(tagged_sepvm)-1
      current_tag = getTagForIndex(tagged_sepvm, t)
      call getSepVMAtIndex(tagged_sepvm, t, sepvm)      
      if(current_tag > this%mesh%ncell) then
         ! Is a boundary condition, all volume in phase 0
          boundary_vof =  this%getBCMaterialFractions(current_tag, a_vof)
          a_cell_flux(1:2) = &
               a_cell_flux(1:2) + boundary_vof * getVolume(sepvm, 0)
          a_cell_flux(3:5) = a_cell_flux(3:5) + boundary_vof(1) * getCentroid(sepvm,0)
          a_cell_flux(6:8) = a_cell_flux(6:8) + boundary_vof(2) * getCentroid(sepvm,0)           
      else
         ! Is inside domain, trust actual volumes
         a_cell_flux(1) = a_cell_flux(1) + getVolume(sepvm, 0)
         a_cell_flux(2) = a_cell_flux(2) + getVolume(sepvm, 1)
         a_cell_flux(3:5) = a_cell_flux(3:5) + getCentroid(sepvm, 0)
         a_cell_flux(6:8) = a_cell_flux(6:8) + getCentroid(sepvm, 1)         
      end if
    end do
    ! Normalize to get centroid location
    if(abs(a_cell_flux(1)) > tiny(1.0_r8)) then
       a_cell_flux(3:5) = a_cell_flux(3:5) / a_cell_flux(1)
    end if

    if(abs(a_cell_flux(2)) > tiny(1.0_r8)) then
       a_cell_flux(6:8) = a_cell_flux(6:8) / a_cell_flux(2)
    end if

    return
  end subroutine moments_from_geometric_cutting
  
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

    return
  end function getBCMaterialFractions
  
  ! subroutine update_vof(this, a_interface_band, a_flux_vol, a_vof)
  
  !   class(unsplit_geometric_volume_tracker), intent(in) :: this
  !   integer, intent(in) :: a_interface_band(:)
  !   real(r8), intent(inout) :: a_flux_vol(:,:)
  !   real(r8), intent(inout) :: a_vof(:,:)
    
  !   integer :: j, f, k
  !   integer :: number_of_faces
  !   real(r8) :: cell_volume_sum
  !   logical :: vof_modified
    
  !   ! Fill the ragged a_flux_vol array
  !   a_flux_vol = 0.0_r8
  !   do j = 1, this%mesh%ncell_onP
  !     if(abs(a_interface_band(j)) > advect_band) then
  !       cycle
  !     end if
  !     associate( fn => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1) )      
  !      do f = 1, size(fn)
  !       if(btest(this%mesh%cfpar(j), f)) then
  !         a_flux_vol(:,this%mesh%xcface(j)+f-1) = -this%face_flux(:,fn(f))
  !       else
  !         a_flux_vol(:,this%mesh%xcface(j)+f-1) =  this%face_flux(:,fn(f))
  !       end if
  !     end do
            
  !     end associate
      
  !   end do

  !   ! Use a_flux_vol to update a_vof now
  !   do j = 1, this%mesh%ncell_onP
  !     if(abs(a_interface_band(j)) > advect_band) then
  !       cycle
  !     end if
  !     a_vof(:,j) = a_vof(:,j) * this%mesh%volume(j)
  !     cell_volume_sum = this%mesh%volume(j)      
  !     number_of_faces = this%mesh%xcface(j+1) - this%mesh%xcface(j)
  !     do f = 1, number_of_faces
  !       a_vof(:,j) = a_vof(:,j) - a_flux_vol(:, this%mesh%xcface(j)+f-1)
  !       cell_volume_sum = cell_volume_sum - sum(a_flux_vol(:, this%mesh%xcface(j)+f-1))
  !     end do

  !      ! DEBUG
  !      if(abs(this%mesh%volume(j) - cell_volume_sum) .gt. 1.0e-14) then
  !        print*,'Unconservative fluxing in volume!!', j, this%mesh%volume(j), cell_volume_sum
  !        print*,a_vof(:,j)
  !        print*,sum(abs(a_vof(:,j)))
  !      end if
  !     a_vof(:,j) = a_vof(:,j) / cell_volume_sum
  !     vof_modified = .false.
  !     ! This will lead to inconsistency with face fluxes
  !     do k = 1, this%nmat
  !       if(a_vof(k,j) < this%cutoff .and. abs(a_vof(k,j)) > EPSILON(1.0_r8)) then
  !         vof_modified = .true.
  !         a_vof(k,j) = 0.0_r8
  !       else if(a_vof(k,j) > 1.0_r8 - this%cutoff .and. abs(1.0_r8 - a_vof(k,j)) >  EPSILON(1.0_r8)) then
  !         vof_modified = .true.
  !         a_vof(k,j) = 1.0_r8
  !       end if
  !     end do
    
  !     if(vof_modified) then
  !       ! Rescale VOF to be valid (sum to 1)
  !       a_vof(:,j) = a_vof(:,j) / sum(a_vof(:,j))
  !     end if
  !   end do

  !   call gather_boundary(this%mesh%cell_ip, a_vof)

  ! end subroutine update_vof    

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
  ! on a first-come, first-served basis based on the global Node ID
  ! in order to be globally consistent. Each face in
  ! this%mesh is looped over, and edges are set to have Leading
  ! diagonalizations (L) from the node with the lower global
  ! ID to the higher one.    
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
    integer :: edge_direction
    integer :: f
    integer :: n, number_of_nodes    
    integer :: edge_end
    integer :: lookup_case
    integer :: shift_amount, r, tmp_node, reordered_vertices(4)
    integer :: irl_ordering(4)
    
    allocate(this%flux_geometry_class(this%mesh%nface))
    allocate(this%flux_node(this%mesh%xfnode(this%mesh%nface+1)-1))

    ! Treat edges as being diagonalized from lower global node id to higher
    ! For each face, now store correct flux class and
    ! node ordering to be consistent
    do f = 1, this%mesh%nface
      associate(node_id => this%mesh%fnode(this%mesh%xfnode(f):this%mesh%xfnode(f+1)-1))
        number_of_nodes = size(node_id)
        ASSERT(number_of_nodes == 3 .or. number_of_nodes == 4)        
        irl_ordering = reorder_node_id(node_id)
        irl_ordering = make_smallest_global_ind_first(this%mesh, irl_ordering(1:number_of_nodes))
        lookup_case = 1
        do n = 1, number_of_nodes          
          edge_end = irl_ordering(mod(n,number_of_nodes)+1)
          if(this%mesh%xnode(irl_ordering(n)) < this%mesh%xnode(edge_end)) then
             edge_direction = 0
          else
             edge_direction = 1
          end if
          lookup_case = lookup_case + edge_direction * 2**(n-1)
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

  contains

    ! Reordering face to that used by IRL
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

    ! This is needed because it seems subdomains do not agree on face ordering
    ! of nodes? 
    function make_smallest_global_ind_first(a_mesh, a_node_list) result(a_new_list)
      type(unstr_mesh), intent(in) :: a_mesh
      integer, intent(in) :: a_node_list(:)      
      integer :: a_new_list(4)

      integer :: n, smallest_gid, sgid_ind, ind
      
      ASSERT(size(a_node_list) == 3 .or. size(a_node_list) == 4)

      smallest_gid = maxval(a_mesh%xnode(a_node_list))+1
      do n = 1, size(a_node_list)
        if(a_mesh%xnode(a_node_list(n)) < smallest_gid) then
           smallest_gid = a_mesh%xnode(a_node_list(n))
           sgid_ind = n
        end if
      end do

      do n = 1, size(a_node_list)
        ind = mod(sgid_ind+n-1,size(a_node_list))+1
        a_new_list(n) = a_node_list(ind)
      end do
      
      
    end function make_smallest_global_ind_first
    
  end subroutine generate_flux_classes

  subroutine generate_cell_bounding_boxes(this)

    class(unsplit_geometric_volume_tracker), intent(inout) :: this

    integer :: n, d

    do n = 1, this%mesh%ncell
      associate(cn => this%mesh%cnode(this%mesh%xcnode(n):this%mesh%xcnode(n+1)-1))
        do d = 1,3
          this%cell_bounding_box(d, 1, n) = minval(this%mesh%x(d,cn))
          this%cell_bounding_box(d, 2, n) = maxval(this%mesh%x(d,cn))        
        end do
      end associate
    end do

  end subroutine generate_cell_bounding_boxes

end module unsplit_geometric_volume_tracker_type
