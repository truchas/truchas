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
  use unsplit_volume_tracker_class
  use truchas_logging_services
  use truchas_timers
  use unstr_mesh_type
  use index_partitioning
  use irl_fortran_interface
  use parameter_module, only : string_len
  use cell_tagged_mm_volumes_type
  implicit none
  private

  integer, parameter, private :: advect_band = 5 ! Size of band to do advection

  type, extends(unsplit_volume_tracker), public :: unsplit_geometric_volume_tracker
    private
    type(unstr_mesh), pointer :: mesh ! unowned reference
    integer :: location_iter_max ! maximum number of iterations to use in fitting interface
    integer :: subcycles
    logical :: nested_dissection
    character(string_len) :: interface_reconstruction_name
    real(r8), allocatable :: normal(:,:,:)
    type(cell_tagged_mm_volumes), allocatable :: face_flux(:)
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
    type(PlanarLoc_type), allocatable :: planar_localizer(:)
    type(PlanarSep_type), allocatable :: planar_separator(:,:)
    type(PlanarSepPath_type) :: planar_separator_path
    type(PlanarSepPathGroup_type), allocatable :: planar_separator_path_group(:)
    type(LocSepGroupLink_type), allocatable :: localized_separator_group_link(:)
    type(Poly_type), allocatable, private :: interface_polygons(:) 
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
    procedure, public :: write_interface
    procedure, private :: init_irl_mesh
    procedure, private :: generate_flux_classes
    procedure, private :: compute_velocity_weightings
    procedure, private :: set_irl_priority_order
    procedure, private :: interface_reconstruction    
    procedure, private :: normals_youngs
    procedure, private :: normals_swartz
    procedure, private :: normals_lvira
    procedure, private :: normals_lvira_hex_execute
    procedure, private :: normals_lvira_tet_execute    
    procedure, private :: identify_reconstruction_neighborhood    
    procedure, private :: set_irl_interfaces
    procedure, private :: set_plane_distances
    procedure, private :: adjust_planes_match_VOF
    procedure, private :: reset_volume_moments
    procedure, private :: construct_interface_polygons
    procedure, private :: construct_nonzero_polygon
    procedure, private :: compute_node_velocities
    procedure, private :: compute_effective_cfl
    procedure, private :: compute_projected_nodes
    procedure, private :: compute_fluxes
    procedure, private :: update_vof
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
    character(:), allocatable :: interface_recon
    integer :: i, j, k
    integer :: material_present(nmat)
    
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

      ! QUESTION: Come back and visit solid. Will need to be first.
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
    material_present = 0
    do i = 1, this%nmat
      material_present(this%priority(i)) = material_present(this%priority(i)) + 1
    end do
    ! Check each material shows up exactly once
    ASSERT(maxval(material_present) == 1)
    ASSERT(minval(material_present) == 1)
    
    ! Allocate face fluxes we'll store
    allocate(this%face_flux(this%mesh%nface))

    allocate(this%normal(3,this%nmat,mesh%ncell))
    allocate(this%w_node(3,mesh%nnode))
    allocate(this%velocity_weightings(mesh%xcnode(mesh%ncell+1)-1))
    allocate(this%projected_nodes(3, mesh%nnode))

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

    ! Calculate weightings for interpolation of cell center velocity to node
    call this%compute_velocity_weightings

  end subroutine init

  ! flux volumes routine assuming flux_vol is a cface-like array
  ! flux volumes routine assuming vel is stored on faces, length 1:mesh%nfaces
  subroutine flux_volumes(this, vel, vel_cc, vof_n, vof, flux_vol, fluids, void, dt, a_interface_band)

    use integer_real8_tuple_vector_type
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: vel(:), vel_cc(:,:), vof_n(:,:), dt
    type(cell_tagged_mm_volumes), intent(out) :: flux_vol(:)
    real(r8), intent(out) :: vof(:,:)
    integer, intent(in) :: fluids, void
    integer, intent(in) :: a_interface_band(:)

    integer :: i,j,k

    ! DEBUGGING
    integer :: f
    real(r8) :: tmp
    type(integer_real8_tuple_vector) :: cell_local_volumes
    type(integer_real8_tuple_vector), pointer :: cell_local_volumes_ptr    



    ! ! Doing a simple test on these classes...
    ! call cell_local_volumes%resize(3)
    ! call cell_local_volumes%set(1, 3, 6.0_r8)
    ! call cell_local_volumes%set(2, 5, 10.0_r8)
    ! call cell_local_volumes%set(3, -2, -4.0_r8)
    ! call this%face_flux(1)%reserve(2)
    ! call this%face_flux(1)%add_cell_fluxes(15, cell_local_volumes)
    ! call cell_local_volumes%resize(0)
    ! call cell_local_volumes%push_back(15, 30.0_r8)
    ! call cell_local_volumes%push_back(22, 44.0_r8)
    ! call this%face_flux(1)%add_cell_fluxes(4096, cell_local_volumes)


    ! do f = 1, this%face_flux(1)%get_number_of_cells()
    !   cell_local_volumes_ptr => this%face_flux(1)%get_cell_fluxes(f)
    !   print*,'For cell', this%face_flux(1)%get_cell_id(f)
    !   do k = 1, cell_local_volumes_ptr%size()
    !     print*, k, cell_local_volumes_ptr%at_int(k), cell_local_volumes_ptr%at_r8(k)
    !   end do
    ! end do



    ! call TLS_fatal('Debug')

    do f = 1, this%mesh%nface
      call this%face_flux(f)%clear()
    end do

    vof = vof_n
    call start_timer('reconstruction')
    call this%set_irl_priority_order(vof)
    call this%interface_reconstruction(vof)
    call this%reset_volume_moments(vof)    
    call this%construct_interface_polygons(a_interface_band)
    call stop_timer('reconstruction')

    call start_timer('advection')
    call this%compute_node_velocities(vel_cc)
    call this%compute_effective_cfl(dt)

    call this%compute_projected_nodes(dt)

    call this%compute_fluxes(vel, dt, vof, a_interface_band)

    call this%update_vof(a_interface_band, flux_vol, vof)

    ! ! DEBUG
    ! do j = 1, this%mesh%ncell_onP
    !    tmp = 0.0_r8
    !    associate(fid => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
    !      do f = 1, size(fid)
    !         if(btest(this%mesh%cfpar(j), f)) then
    !           tmp = tmp + vel(fid(f))*dt*this%mesh%area(fid(f))
    !         else
    !           tmp = tmp - vel(fid(f))*dt*this%mesh%area(fid(f))
    !         end if
    !      end do
    !    end associate

    !    if(abs(tmp > 1.0e-12)) then
    !      print*,'Velocity field has divergence of ', tmp, ' for cell', j
    !      print*,'This is ',100.0_r8*tmp/this%mesh%volume(j),' % of the cell volume.'
    !    end if       
    ! end do
    ! ! END DEBUG

    call stop_timer('advection')

    do f = 1, this%mesh%nface
      call this%face_flux(f)%clear()
    end do    

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
    
    integer :: j, f, face_index, cn, k
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
    call new(this%object_server_planar_separator, int(this%nmat*total_recon_needed, i8))
        
    ! Allocate the Fortran memory
    allocate(this%planar_localizer(total_recon_needed))
    allocate(this%planar_separator(this%nmat, total_recon_needed))
    allocate(this%planar_separator_path_group(total_recon_needed))
    allocate(this%localized_separator_group_link(total_recon_needed))    
    call new(this%planar_separator_path)
    
    ! Now allocate on the IRL side the different planar reconstructions
    do j = 1, total_recon_needed
      call new(this%planar_localizer(j), this%object_server_planar_localizer)
      call new(this%planar_separator_path_group(j))      
      do k = 1, this%nmat
        call new(this%planar_separator(k, j), this%object_server_planar_separator)
        call construct(this%planar_separator_path, this%planar_separator(k,j))
        call addPlanarSeparatorPath(this%planar_separator_path_group(j), this%planar_separator_path, k)
      end do
      call new(this%localized_separator_group_link(j), &
               this%planar_localizer(j), &
               this%planar_separator_path_group(j) )           
    end do

    ! Now setup PlanarLocalizers (one per cell) and link together.
    do j = 1, this%mesh%ncell    
      number_of_cell_faces = this%mesh%xcface(j+1)-this%mesh%xcface(j)
                                                                    
      call setNumberOfPlanes(this%planar_localizer(j), number_of_cell_faces)
      call setId(this%localized_separator_group_link(j), j)      
            
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
            call setEdgeConnectivity(this%localized_separator_group_link(j), found_internal-1, &
                                     this%localized_separator_group_link(cn(f)))             
          else
            ! Connect to correct "outside" marking localized_separator_link
            ! Making these last in reconstruction should keep all inside domain volume inside the domain
            found_boundary = found_boundary + 1 ! IRL is 0-based index
            call setPlane(this%planar_localizer(j), total_internal + found_boundary-1, tmp_plane(1:3), tmp_plane(4))
            call setEdgeConnectivity(this%localized_separator_group_link(j), total_internal + found_boundary-1, &
                                     this%localized_separator_group_link(face_to_boundary_mapping(face_index)))            
          end if
                
        end do
      end associate
      
    end do
      
    ! Make "outside" reconstructions consume all volume given to them.
    do j = this%mesh%ncell+1,total_recon_needed
      call setNumberOfPlanes(this%planar_separator(1,j), 1)
      call setPlane(this%planar_separator(1,j), 0, [0.0_r8, 0.0_r8, 0.0_r8], 1.0_r8)
      call setNumberOfPlanes(this%planar_localizer(j), 1)
      call setPlane(this%planar_localizer(j), 0, [0.0_r8, 0.0_r8, 0.0_r8], 1.0_r8)
      call setEdgeConnectivityNull(this%localized_separator_group_link(j), 0) ! No link back to domain      
      call setId(this%localized_separator_group_link(j), j)
      call setPriorityOrder(this%planar_separator_path_group(j), 1, [1])
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

    call gather_boundary(this%mesh%node_ip, this%w_node(1,:))    

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

    ! Normalize weightings.
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

  subroutine set_irl_priority_order(this, a_vof)

    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in)  :: a_vof(:,:)

    integer :: j, k, priority_index
    integer :: irl_priority(this%nmat), valid_size

    do j = 1, this%mesh%ncell
      valid_size = 0
      do k = 1, this%nmat
        priority_index = this%priority(k)
        if(a_vof(priority_index, j) > this%cutoff) then
          valid_size = valid_size + 1
          irl_priority(valid_size) = priority_index
        end if
      end do
      ASSERT(valid_size > 0)
      call setPriorityOrder(this%planar_separator_path_group(j), valid_size, irl_priority)      
    end do

  end subroutine set_irl_priority_order

  subroutine interface_reconstruction(this, vof)


    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in)  :: vof(:,:)


    select case(trim(this%interface_reconstruction_name))

      case('Youngs')
        call this%normals_youngs(vof)

      case('Swartz')
        call this%normals_youngs(vof)
        call this%normals_swartz(vof)

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

    call gather_boundary(this%mesh%cell_ip, this%normal)
    
  end subroutine normals_youngs

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
          call setNumberOfPlanes(this%planar_separator(1,j), 1)
          call setPlane(this%planar_separator(1,j), 0, this%normal(:,1,j), distance_guess)
          
          associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
            
            call this%adjust_planes_match_VOF(this%mesh%x(:,cn), a_vof(1,j), &
                 this%planar_separator(1,j) )

            call this%construct_nonzero_polygon(this%mesh%x(:,cn), this%planar_separator(1,j), &
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
    call reconstructLVIRA3D(neighborhood, this%planar_separator(1,a_center_index))   

    ! Move normal to this%normal
    tmp_plane = getPlane(this%planar_separator(1,a_center_index),0)
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
    call reconstructLVIRA3D(neighborhood, this%planar_separator(1,a_center_index))   

    ! Move normal to this%normal
    tmp_plane = getPlane(this%planar_separator(1,a_center_index),0)
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
            distance_to_centroid = sqrt(sum((&
                 this%mesh%cell_centroid(:,cn(n))-this%mesh%cell_centroid(:,a_starting_index))**2))
            if(distance_to_centroid < a_stencil_length .and. a_stencil%cell_not_encountered(cn(n))) then
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
        end do
      end associate
    end do

  end subroutine identify_reconstruction_neighborhood
  
  subroutine set_irl_interfaces(this, a_vof)
    
    use cell_geom_type
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: a_vof(:,:)
    
    
    integer :: j, k, priority_index, priority_size
    real(r8) :: distance_guess
    
    do j = 1, this%mesh%ncell
      priority_size = getPriorityOrderSize(this%planar_separator_path_group(j))
      associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))        
        call this%set_plane_distances(this%mesh%x(:,cn), this%mesh%cell_centroid(:,j), &
             this%normal(:,:,j), a_vof(:,j), &
             priority_size, this%planar_separator(:,j), this%planar_separator_path_group(j))
      end associate
    end do
    
  end subroutine set_irl_interfaces

  subroutine set_plane_distances(this, a_cell, a_cell_centroid, a_normal, a_vof, a_priority_size, &
       a_planar_separator, a_path_group)

    use irl_interface_helper

    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: a_cell(:,:)
    real(r8), intent(in) :: a_cell_centroid(3)
    real(r8), intent(in) :: a_normal(:,:)
    real(r8), intent(in) :: a_vof(:)
    integer, intent(in) :: a_priority_size
    type(PlanarSep_type), intent(inout) :: a_planar_separator(:)
    type(PlanarSepPathGroup_type), intent(inout) :: a_path_group

    integer :: phase, priority_index
    real(r8) :: tmp_plane(4)
    real(r8) :: ordered_vof(a_priority_size)

    real(r8), parameter :: VOF_tolerance = 1.0e-14_r8    

    do phase = 0, a_priority_size-1
      priority_index = getPriorityOrderTag(a_path_group, phase)
      tmp_plane(1:3) = a_normal(:,priority_index)
      tmp_plane(4) = dot_product(tmp_plane(1:3), a_cell_centroid)
      call setNumberOfPlanes(a_planar_separator(priority_index), 1)
      call setPlane(a_planar_separator(priority_index),0, tmp_plane(1:3), tmp_plane(4))
      ordered_vof(phase+1) = a_vof(priority_index)
    end do

    ! PlanarSeparator should already be setup with a valid normal and
    ! guess for the distance (atleast somewhere inside cell so we can calculate
    ! a gradient in IRL for the optimization)    
    select case(size(a_cell,dim=2))
      case (4) ! tet      
        call truchas_poly_to_irl(a_cell, this%IRL_tet)
        call matchGroupVolumeFraction(this%IRL_tet, a_priority_size, ordered_vof, a_path_group, VOF_tolerance)
      
      case (5) ! pyramid
        call truchas_poly_to_irl(a_cell, this%IRL_pyramid)
        call matchGroupVolumeFraction(this%IRL_pyramid, a_priority_size, ordered_vof, a_path_group, VOF_tolerance)
      
      case (6) ! Wedge
        call truchas_poly_to_irl(a_cell, this%IRL_wedge)
        call matchGroupVolumeFraction(this%IRL_wedge, a_priority_size, ordered_vof, a_path_group, VOF_tolerance)
    
      case (8) ! Hex
        call truchas_poly_to_irl(a_cell, this%IRL_hex)
        call matchGroupVolumeFraction(this%IRL_hex, a_priority_size, ordered_vof, a_path_group, VOF_tolerance)
    
      case default
        call TLS_fatal('Unknown Truchas cell type during plane-distance setting')
      end select      

    
  end subroutine set_plane_distances
   
  
  subroutine adjust_planes_match_VOF(this, a_cell, a_vof, a_planar_separator)
  
    use irl_interface_helper
    use cell_geometry
  
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: a_cell(:,:)
    real(r8), intent(in) :: a_vof
    type(PlanarSep_type), intent(inout) :: a_planar_separator

    real(r8), parameter :: VOF_tolerance = 1.0e-14_r8
    
    ! PlanarSeparator should already be setup with a valid normal and
    ! guess for the distance (atleast somewhere inside cell so we can calculate
    ! a gradient in IRL for the optimization)    
    select case(size(a_cell,dim=2))
      case (4) ! tet      
        call truchas_poly_to_irl(a_cell, this%IRL_tet)
        call matchVolumeFraction(this%IRL_tet, a_vof, a_planar_separator, VOF_tolerance)
      
      case (5) ! pyramid
        call truchas_poly_to_irl(a_cell, this%IRL_pyramid)
        call matchVolumeFraction(this%IRL_pyramid, a_vof, a_planar_separator, VOF_tolerance)
      
      case (6) ! Wedge
        call truchas_poly_to_irl(a_cell, this%IRL_wedge)
        call matchVolumeFraction(this%IRL_wedge, a_vof, a_planar_separator, VOF_tolerance)
    
      case (8) ! Hex
        call truchas_poly_to_irl(a_cell, this%IRL_hex)
        call matchVolumeFraction(this%IRL_hex, a_vof, a_planar_separator, VOF_tolerance)
    
      case default
        call TLS_fatal('Unknown Truchas cell type during plane-distance setting')
      end select  
    
  end subroutine adjust_planes_match_VOF

 subroutine reset_volume_moments(this, a_vof)
  
    use irl_interface_helper
    use parallel_communication, only : global_sum
    use parameter_module, only : string_len    
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(inout) :: a_vof(:,:)

    integer :: j, phase, phase_id
    type(TagAccVM_Vol_type) :: new_volumes

    ! DEBUG/Diagnostics
    real(r8) :: volume_change(this%nmat)
    real(r8) :: tmp_vof    
    character(string_len) :: message, myformat

    volume_change = 0.0_r8
    call new(new_volumes)
    do j = 1, this%mesh%ncell_onP
      call setMinimumVolToTrack(this%mesh%volume(j)*1.0e-15_r8)
      associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))          
        select case(size(cn))
        case (4) ! tet      
          call truchas_poly_to_irl(this%mesh%x(:,cn), this%IRL_tet)
          call getNormMoments(this%IRL_tet, this%planar_separator_path_group(j), new_volumes)

        case (5) ! pyramid
          call truchas_poly_to_irl(this%mesh%x(:,cn), this%IRL_pyramid)
          call getNormMoments(this%IRL_pyramid, this%planar_separator_path_group(j), new_volumes)
          
        case (6) ! Wedge
          call truchas_poly_to_irl(this%mesh%x(:,cn), this%IRL_wedge)
          call getNormMoments(this%IRL_wedge, this%planar_separator_path_group(j), new_volumes)
          
        case (8) ! Hex
          call truchas_poly_to_irl(this%mesh%x(:,cn), this%IRL_hex)
          call getNormMoments(this%IRL_hex, this%planar_separator_path_group(j), new_volumes)
          
        case default
          call TLS_fatal('Unknown Truchas cell type during plane-distance setting')
        end select

      end associate

      do phase = 0, getSize(new_volumes)-1
        phase_id = getTagForIndex(new_volumes, phase)
        tmp_vof = a_vof(phase_id, j)
        a_vof(phase_id, j) = getVolumeAtIndex(new_volumes, phase) / this%mesh%volume(j)
        volume_change(phase_id) = volume_change(phase_id) + (a_vof(phase_id, j) - tmp_vof)*this%mesh%volume(j)
      end do      
    end do

    call gather_boundary(this%mesh%cell_ip, a_vof)

    do phase = 1, this%nmat
      volume_change(phase) = global_sum(volume_change(phase))
    end do

    write(myformat, '(a,i1,a)') '(a,',this%nmat,'es12.4)'
    write(message, trim(myformat)) 'Phase volumes lost during reconstruction', volume_change
    call TLS_info(message)
    
  end subroutine reset_volume_moments

  subroutine construct_interface_polygons(this, a_interface_band)

    use irl_interface_helper    
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    integer, intent(in) :: a_interface_band(:)

    integer :: j, priority_index

    ! QUESTION: How to produce polygons for multiple nested disection interfaces?
    
    do j = 1, this%mesh%ncell
      if(a_interface_band(j) /= 0) then
        call zeroPolygon(this%interface_polygons(j))
        cycle
      end if
      priority_index = getPriorityOrderTag(this%planar_separator_path_group(j), 0)
      associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
        call this%construct_nonzero_polygon(this%mesh%x(:,cn), this%planar_separator(priority_index,j),&
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
      
      call gather_boundary(this%mesh%node_ip, this%w_node(1:3,:))
      
  end subroutine compute_node_velocities

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
       write(message, '(i4,a)') number_of_crossed_faces, 'face crossings might exist, could lost discrete conservation.'
       call TLS_warn(message)
    end if
      
  end subroutine compute_effective_cfl

  subroutine compute_projected_nodes(this, a_dt)

    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: a_dt

    integer :: n

    do n = 1, this%mesh%nnode_onP
      this%projected_nodes(:,n) = this%mesh%x(:,n) + this%w_node(1:3,n)*(-a_dt)
    end do

    call gather_boundary(this%mesh%node_ip, this%projected_nodes)

  end subroutine compute_projected_nodes
  
  subroutine compute_fluxes(this, a_face_vel, a_dt, a_vof, a_interface_band)
  
    use irl_interface_helper
    use integer_real8_tuple_vector_type
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: a_face_vel(:)
    real(r8), intent(in) :: a_dt
    real(r8), intent(in) :: a_vof(:,:)
    integer, intent(in)  :: a_interface_band(:)
        
    integer :: f, v, i, t, k
    integer ::  number_of_nodes, neighbor_cell_index, phase_id, phase, phase_size
    real(r8) :: cell_nodes(3,9), phase_volume(this%nmat), cell_volume, correct_flux_vol
    type(TagAccVM2_Vol_type) :: irl_tagged_volumes
    type(TagAccVM_Vol_type) :: irl_cell_local_volumes
    type(integer_real8_tuple_vector) :: cell_local_volumes
    type(integer_real8_tuple_vector), pointer :: cell_local_volumes_ptr    
    integer :: current_tag, min_band
    
    call new(irl_tagged_volumes)
    call getMoments_setMethod(1)
    call cell_local_volumes%reserve(this%nmat)
      
    do f = 1, this%mesh%nface_onP
      call this%face_flux(f)%clear()

      neighbor_cell_index = this%mesh%fcell(1, f)
      if(this%mesh%fcell(2,f) /= 0) then
        cell_volume = min(this%mesh%volume(neighbor_cell_index), &
                          this%mesh%volume(this%mesh%fcell(2,f)))
        min_band = min(abs(a_interface_band(neighbor_cell_index)),abs(a_interface_band(this%mesh%fcell(2,f))))
      else
        cell_volume = this%mesh%volume(neighbor_cell_index)
        min_band = abs(a_interface_band(neighbor_cell_index))
      end if       

      ! Handle case if far from interface, means
      ! flux is just a single material. This is the
      ! simplification/model that we can use neighbor-face fluxing
      ! to calculate the flux values. This is true for
      ! material VOF advection, but an approximation for
      ! other quantities based on this (momentum, etc.)
      if(min_band > advect_band) then
        ! Single phase flux
        phase_id = maxloc(a_vof(:,neighbor_cell_index),1)
        correct_flux_vol = a_dt*a_face_vel(f)*this%mesh%area(f)
        call this%face_flux(f)%reserve(1)
        call cell_local_volumes%resize(1)
        ! Assuming neumann condition here if on boundary (and there's no interface nearby)
        ! Only prevents influxing of material into domain with that material not with
        ! advect_band cells
        call cell_local_volumes%set(1, phase_id, correct_flux_vol)
        ! Determine which cell from upwinding
        if(a_face_vel(f) > 0.0_r8) then
          call this%face_flux(f)%add_cell_fluxes(this%mesh%fcell(1,f), cell_local_volumes)
        else if(a_face_vel(f) < 0.0_r8) then
          if(this%mesh%fcell(2,f) /= 0) then
            call this%face_flux(f)%add_cell_fluxes(this%mesh%fcell(2,f), cell_local_volumes)            
          else ! Boundary condition
            ! Handle boundary somehow?
            ! -1 will be used as a sentinel value to signify it came from the boundary.
            call this%face_flux(f)%add_cell_fluxes(-1, cell_local_volumes)
          end if
        end if
        
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
            call getMoments(this%IRL_CapOcta_LLL, &
                            this%localized_separator_group_link(neighbor_cell_index), &
                            irl_tagged_volumes)

          case(2) ! Triangular Face -> Octahedron Volume
          
            call construct(this%IRL_CapOcta_LLT, cell_nodes)
            call adjustCapToMatchVolume(this%IRL_CapOcta_LLT, a_dt*a_face_vel(f)*this%mesh%area(f))
            call getMoments(this%IRL_CapOcta_LLT, &
                            this%localized_separator_group_link(neighbor_cell_index), &
                            irl_tagged_volumes)
            
          case(3) ! Triangular Face -> Octahedron Volume
          
            call construct(this%IRL_CapOcta_LTT, cell_nodes)
            call adjustCapToMatchVolume(this%IRL_CapOcta_LTT, a_dt*a_face_vel(f)*this%mesh%area(f))       
            call getMoments(this%IRL_CapOcta_LTT, &
                            this%localized_separator_group_link(neighbor_cell_index), &
                            irl_tagged_volumes)

          case(4) ! Triangular Face -> Octahedron Volume
          
            call construct(this%IRL_CapOcta_TTT, cell_nodes)            
            call adjustCapToMatchVolume(this%IRL_CapOcta_TTT, a_dt*a_face_vel(f)*this%mesh%area(f))
            call getMoments(this%IRL_CapOcta_TTT, &
                            this%localized_separator_group_link(neighbor_cell_index), &
                            irl_tagged_volumes) 
            
          case(5) ! Quad Face -> Dodecahedron Volume
            
            call construct(this%IRL_CapDod_LLLL, cell_nodes)
            call adjustCapToMatchVolume(this%IRL_CapDod_LLLL, a_dt*a_face_vel(f)*this%mesh%area(f))
            call getMoments(this%IRL_CapDod_LLLL, &
                            this%localized_separator_group_link(neighbor_cell_index), &
                            irl_tagged_volumes)

          case(6) ! Quad Face -> Dodecahedron Volume
            
            call construct(this%IRL_CapDod_LLLT, cell_nodes)            
            call adjustCapToMatchVolume(this%IRL_CapDod_LLLT, a_dt*a_face_vel(f)*this%mesh%area(f))            
            call getMoments(this%IRL_CapDod_LLLT, &
                            this%localized_separator_group_link(neighbor_cell_index), &
                            irl_tagged_volumes)

          case(7) ! Quad Face -> Dodecahedron Volume
            
            call construct(this%IRL_CapDod_LTLT, cell_nodes)            
            call adjustCapToMatchVolume(this%IRL_CapDod_LTLT, a_dt*a_face_vel(f)*this%mesh%area(f))            
            call getMoments(this%IRL_CapDod_LTLT, &
                            this%localized_separator_group_link(neighbor_cell_index), &
                            irl_tagged_volumes)

          case(8) ! Quad Face -> Dodecahedron Volume
            
            call construct(this%IRL_CapDod_LLTT, cell_nodes)            
            call adjustCapToMatchVolume(this%IRL_CapDod_LLTT, a_dt*a_face_vel(f)*this%mesh%area(f))
            call getMoments(this%IRL_CapDod_LLTT, &
                            this%localized_separator_group_link(neighbor_cell_index), &
                            irl_tagged_volumes)

          case(9) ! Quad Face -> Dodecahedron Volume
            
            call construct(this%IRL_CapDod_LTTT, cell_nodes)            
            call adjustCapToMatchVolume(this%IRL_CapDod_LTTT, a_dt*a_face_vel(f)*this%mesh%area(f))            
            call getMoments(this%IRL_CapDod_LTTT, &
                            this%localized_separator_group_link(neighbor_cell_index), &
                            irl_tagged_volumes)

          case(10) ! Quad Face -> Dodecahedron Volume
            
            call construct(this%IRL_CapDod_TTTT, cell_nodes)            
            call adjustCapToMatchVolume(this%IRL_CapDod_TTTT, a_dt*a_face_vel(f)*this%mesh%area(f))            
            call getMoments(this%IRL_CapDod_TTTT, &
                            this%localized_separator_group_link(neighbor_cell_index), &
                            irl_tagged_volumes)

          case default
            call TLS_fatal('Unknown face type. How was this not found in this%generate_flux_classes()?')          
        end select
          
      end associate

      call this%face_flux(f)%reserve(getSize(irl_tagged_volumes))
      do t = 0, getSize(irl_tagged_volumes)-1
        current_tag = getTagForIndex(irl_tagged_volumes, t)
        call getAtIndex(irl_tagged_volumes, t, irl_cell_local_volumes)
        call cell_local_volumes%resize(0)        
        if(current_tag > this%mesh%ncell) then
          ! Is a boundary condition, all volume in phase 0          
          phase_volume = this%getBCMaterialFractions(current_tag, a_vof) * getVolumeAtIndex(irl_cell_local_volumes, 0)
          do k = 1, this%nmat
            if(abs(phase_volume(k)) > tiny(1.0_r8)) then              
              call cell_local_volumes%push_back(k, phase_volume(k))           
            end if
          end do
        else
          ! Is inside domain, trust actual volumes
          phase_size = getSize(irl_cell_local_volumes)
          do phase = 0, phase_size - 1 
            phase_id = getTagForIndex(irl_cell_local_volumes, phase)
            if(abs(getVolumeAtIndex(irl_cell_local_volumes, phase)) > tiny(1.0_r8)) then
              call cell_local_volumes%push_back(phase_id, getVolumeAtIndex(irl_cell_local_volumes, phase))
            end if
          end do          
        end if
        if(cell_local_volumes%size() > 0) then
           call this%face_flux(f)%add_cell_fluxes(current_tag, cell_local_volumes)
        end if
      end do      
    end do

    ! Need to communicate fluxes
    call communicate_face_fluxes(this%mesh, this%face_flux)
    
  end subroutine compute_fluxes

  ! This is going to be a brutal and slow communication.
  ! Could probably be sped up in future with some
  ! specialization of the data structures.
  subroutine communicate_face_fluxes(a_mesh, a_face_flux)

    use integer_real8_tuple_vector_type
    use parallel_communication    
    
    type(unstr_mesh), pointer, intent(in) :: a_mesh
    type(cell_tagged_mm_volumes), intent(inout) :: a_face_flux(:)

    integer :: f, c, p, ci
    integer :: current_buffer_size, max_buffer_size
    integer ::nfcells, nphases, cell_id
    real(r8), allocatable :: face_flux_buffer(:,:)
    type(integer_real8_tuple_vector), pointer :: tagged_phase_array
    type(integer_real8_tuple_vector) :: tmp_phase_array


    ! Need to communicate cell ID using global cell ID and map back to local cell ID.
    ! The cell ID's stored in the cell_tagged moments is for THAT processor,
    ! not globally.
    call TLS_FATAL('Debug')
    
    ! Calculate needed buffer size
    max_buffer_size = -1    
    do f = 1, a_mesh%nface_onP      
      current_buffer_size = 1 ! Number of cells
      current_buffer_size = current_buffer_size +  a_face_flux(f)%get_number_of_cells() ! Cell id's
      do c = 1, a_face_flux(f)%get_number_of_cells()
        tagged_phase_array => a_face_flux(f)%get_cell_fluxes(c)
        current_buffer_size = current_buffer_size + 1 ! Number of phase Ids
        current_buffer_size = current_buffer_size + 2*tagged_phase_array%size() ! Phase data
      end do
      max_buffer_size = max(max_buffer_size, current_buffer_size)
    end do

    ! Allocate the buffer
    max_buffer_size = global_maxval(max_buffer_size)
    allocate(face_flux_buffer(max_buffer_size, a_mesh%nface))
    
    ! Fill the buffer. Note the ordering, needed for unpacking
    do f = 1, a_mesh%nface_onP
      ci = 0
      ci = ci + 1      
      face_flux_buffer(ci,f) = real(a_face_flux(f)%get_number_of_cells(), r8) ! Number of cells
      do c = 1, a_face_flux(f)%get_number_of_cells()
        ci = ci + 1
        face_flux_buffer(ci,f) = real(a_face_flux(f)%get_cell_id(c), r8) ! Cell id
        tagged_phase_array => a_face_flux(f)%get_cell_fluxes(c)        
        ci = ci + 1        
        face_flux_buffer(ci,f) = real(tagged_phase_array%size(), r8) ! Number of phases
        do p = 1, tagged_phase_array%size()
          ci = ci + 1
          face_flux_buffer(ci,f) = real(tagged_phase_array%at_int(p), r8) ! phase tag
          ci = ci + 1          
          face_flux_buffer(ci,f) = tagged_phase_array%at_r8(p) ! phase volume
        end do
      end do
    end do

    ! Communicate buffer
    call gather_boundary(a_mesh%face_ip, face_flux_buffer)

    ! Unpack the buffer
    do f = a_mesh%nface_onP+1, a_mesh%nface
      call a_face_flux(f)%clear()
      ci = 0
      ci = ci + 1      
      nfcells = NINT(face_flux_buffer(ci,f))
      print*,'Number of cells ', f, nfcells
      call a_face_flux(f)%reserve(nfcells)
      do c = 1, nfcells
        ci = ci + 1
        cell_id = NINT(face_flux_buffer(ci,f))
        print*,'Cell id', cell_id
        ci = ci + 1
        nphases = NINT(face_flux_buffer(ci,f))
        print*,'Number of phases', nphases
        call tmp_phase_array%resize(nphases)
        do p = 1, nphases
          ci = ci + 1
          print*,'phase_id and volume',NINT(face_flux_buffer(ci, f)), face_flux_buffer(ci+1,f)
          call tmp_phase_array%set(p, NINT(face_flux_buffer(ci, f)), face_flux_buffer(ci+1,f))
          ci = ci + 1
        end do
        call a_face_flux(f)%add_cell_fluxes(cell_id, tmp_phase_array)
      end do      
    end do

  end subroutine communicate_face_fluxes

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
  
  subroutine update_vof(this, a_interface_band, a_flux_vol, a_vof)

    use integer_real8_tuple_vector_type
  
    class(unsplit_geometric_volume_tracker), intent(in) :: this
    integer, intent(in) :: a_interface_band(:)
    type(cell_tagged_mm_volumes), intent(out) :: a_flux_vol(:)
    real(r8), intent(inout) :: a_vof(:,:)
    
    integer :: j, f, c, k
    integer :: number_of_faces
    real(r8) :: cell_volume_sum
    real(r8) :: vof_flux(this%nmat)
    type(integer_real8_tuple_vector), pointer :: fluxed_phases
    logical :: vof_modified
    
    ! Fill the ragged a_flux_vol array
    do j = 1, this%mesh%ncell
      associate( fn => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1) )        
        do f = 1, size(fn)

          
       ! do c = 1, this%face_flux(fn(f))%get_number_of_cells()
       !   fluxed_phases => this%face_flux(fn(f))%get_cell_fluxes(c)        
       !   do k = 1, fluxed_phases%size()
       !     if(fluxed_phases%at_int(k) == 0) then
       !       print*,'HOWWWW', j, fn(f),fluxed_phases%at_int(k)
       !     end if
       !   end do
       ! end do

          if(btest(this%mesh%cfpar(j), f)) then
            a_flux_vol(this%mesh%xcface(j)+f-1) = this%face_flux(fn(f))
            call a_flux_vol(this%mesh%xcface(j)+f-1)%negate_in_place()

            
            ! do c = 1, a_flux_vol(this%mesh%xcface(j)+f-1)%get_number_of_cells()
            !   fluxed_phases =>  a_flux_vol(this%mesh%xcface(j)+f-1)%get_cell_fluxes(c)
            !   do k = 1, fluxed_phases%size()
            !     if(fluxed_phases%at_int(k) == 0) then
            !       print*,'Top part',j,fn(f),fluxed_phases%at_int(k)
            !     end if
            !   end do
            ! end do

          else
            a_flux_vol(this%mesh%xcface(j)+f-1) = this%face_flux(fn(f))

            
            ! do c = 1, a_flux_vol(this%mesh%xcface(j)+f-1)%get_number_of_cells()
            !   fluxed_phases =>  a_flux_vol(this%mesh%xcface(j)+f-1)%get_cell_fluxes(c)
            !   do k = 1, fluxed_phases%size()
            !     if(fluxed_phases%at_int(k) == 0) then
            !       print*,'Bottom part',j,fn(f),fluxed_phases%at_int(k)
            !     end if
            !   end do
            ! end do

            
          end if
        end do
      end associate
    end do

    ! Use a_flux_vol to update a_vof now
    do j = 1, this%mesh%ncell_onP
      if(abs(a_interface_band(j)) > advect_band) then
        cycle
      end if
      a_vof(:,j) = a_vof(:,j) * this%mesh%volume(j)
      cell_volume_sum = this%mesh%volume(j)      
      number_of_faces = this%mesh%xcface(j+1) - this%mesh%xcface(j)
      do f = 1, number_of_faces
        vof_flux = 0.0_r8
        do c = 1, a_flux_vol(this%mesh%xcface(j)+f-1)%get_number_of_cells()
          fluxed_phases =>  a_flux_vol(this%mesh%xcface(j)+f-1)%get_cell_fluxes(c)          
          do k = 1, fluxed_phases%size()
            vof_flux(fluxed_phases%at_int(k)) = vof_flux(fluxed_phases%at_int(k)) + fluxed_phases%at_r8(k)
          end do
        end do
        a_vof(:,j) = a_vof(:,j) - vof_flux
        cell_volume_sum = cell_volume_sum - sum(vof_flux)
      end do
      a_vof(:,j) = a_vof(:,j) / cell_volume_sum
      vof_modified = .false.
      ! NOTE: This will lead to slight inconsistency with face fluxes
      do k = 1, this%nmat
        if(a_vof(k,j) < this%cutoff .and. abs(a_vof(k,j)) > EPSILON(1.0_r8)) then
          vof_modified = .true.
          a_vof(k,j) = 0.0_r8
        else if(a_vof(k,j) > 1.0_r8 - this%cutoff .and. abs(1.0_r8 - a_vof(k,j)) >  EPSILON(1.0_r8)) then
          vof_modified = .true.
          a_vof(k,j) = 1.0_r8
        end if
      end do
    

!      if(vof_modified) then
        ! Rescale VOF to be valid (sum to 1)
        a_vof(:,j) = a_vof(:,j) / sum(a_vof(:,j))
!      end if
    end do

    call gather_boundary(this%mesh%cell_ip, a_vof)

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
    do f = 1, this%mesh%nface_onP
      associate(node_id => this%mesh%fnode(this%mesh%xfnode(f):this%mesh%xfnode(f+1)-1))
        irl_ordering = reorder_node_id(node_id)        
        number_of_nodes = size(node_id)
        ASSERT(number_of_nodes == 3 .or. number_of_nodes == 4)
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
