!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

  type, extends(volume_tracker), public :: unsplit_geometric_volume_tracker
    private
    type(unstr_mesh), pointer :: mesh ! unowned reference
    integer :: location_iter_max ! maximum number of iterations to use in fitting interface
    integer :: subcycles
    logical :: nested_dissection
    real(r8) :: cutoff ! allow volume fraction {0,(cutoff,1]}
    real(r8), allocatable :: normal(:,:,:)
    real(r8), allocatable :: face_flux(:,:)
    ! node/face/cell workspace
    real(r8), allocatable :: w_node(:,:)
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
    type(Octa_type) :: IRL_octahedron
    type(Dod_type) :: IRL_dodecahedron
    ! End IRL objects
  contains
    procedure, public  :: init
    procedure, public  :: flux_volumes
    procedure, public  :: set_inflow_material
    procedure, private :: init_irl_mesh
    procedure, private :: set_irl_interfaces
    procedure, private :: adjust_planes_match_VOF
    procedure, private :: compute_node_velocities
    procedure, private :: compute_fluxes
    procedure, private :: project_vertex
    procedure, private :: update_vof
    procedure, private :: clean_VOF
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
    allocate(this%w_node(4,mesh%nnode))

    ! list of boundary face ids
    j = count(this%mesh%fcell(2,:this%mesh%nface_onP) == 0)
    allocate(this%bc_index(j), this%local_face(j), this%inflow_mat(j))
    j = 1
    do i = 1, this%mesh%nface_onP
      if (this%mesh%fcell(2,i) == 0) then
        this%bc_index(j) = i
        k = this%mesh%fcell(1,i)
        this%local_face(j) = &
            findloc(this%mesh%cface(this%mesh%xcface(k):this%mesh%xcface(k+1)-1), i)
        j = j + 1
      end if
    end do
    this%inflow_mat = 0
    
    ! Initialize IRL array needed for advection
    call this%init_irl_mesh()
    
    ! Initialize IRL cell types we will use
    call new(this%IRL_tet)
    call new(this%IRL_pyramid)
    call new(this%IRL_wedge)
    call new(this%IRL_hex)
    call new(this%IRL_octahedron)
    call new(this%IRL_dodecahedron)

  end subroutine init

  ! flux volumes routine assuming flux_vol is a cface-like array
  ! flux volumes routine assuming vel is stored on faces, length 1:mesh%nfaces
  subroutine flux_volumes(this, vel, vel_cc, vof_n, vof, flux_vol, fluids, void, dt)
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: vel(:), vel_cc(:,:), vof_n(:,:), dt
    real(r8), intent(out) :: flux_vol(:,:), vof(:,:)
    integer, intent(in) :: fluids, void

    integer :: i,j

    flux_vol = 0.0_r8
    vof = vof_n
    call this%normals(vof)
    
    call this%set_irl_interfaces(vof)
    
    call this%compute_node_velocities(vel_cc)
    
    call this%compute_fluxes(vel, dt, vof_n)
    
    ! Need to communicate fluxes
    call gather_boundary(this%mesh%face_ip, this%face_flux)
          
    call this%update_vof(flux_vol, vof)
    
    call this%clean_vof(vof)

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
      do f = 1, number_of_cell_faces
        ! Set PlanarLocalizer Plane for this face
        face_index = face_id(f)
        tmp_plane(1:3) = this%mesh%normal(:,face_index) / this%mesh%area(face_index)
        tmp_plane(4) = dot_product(tmp_plane(1:3), this%mesh%x(:,this%mesh%fnode(this%mesh%xfnode(face_index))))
        if(btest(this%mesh%cfpar(j), f)) then 
          ! Want outward oriented normal
          tmp_plane = -tmp_plane
        end if        
        call setPlane(this%planar_localizer(j), f-1, tmp_plane(1:3), tmp_plane(4))        
        
        ! Link this plane in the LocalizedSeparatorLink
        if(cn(f) /= 0) then
          call setEdgeConnectivity(this%localized_separator_link(j), f-1, &
                                   this%localized_separator_link(cn(f))) 
        else
          ! Connect to correct "outside" marking localized_separator_link
          call setEdgeConnectivity(this%localized_separator_link(j), f-1, &
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

!    call TLS_fatal('IMPLEMENT ORGANIZATION OF BOUNDARY CELL PLANES BEING LAST IN RECONSTRUCTION')

   end subroutine init_irl_mesh
  
  subroutine normals(this, vof)

    use flow_operators, only: gradient_cc
    use f08_intrinsics, only: findloc
    intrinsic :: norm2

    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in)  :: vof(:,:)

    real(r8) :: mag
    integer :: i,j,k,c
    logical :: hasvof(size(vof,dim=1))

    call start_timer('normals')
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

    call stop_timer('normals')

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
    
    ! PlanarSeparator should already be setup with a valid normal and
    ! guess for the distance (atleast somewhere inside cell so we can calculate
    ! a gradient in IRL for the optimization)    
    select case(size(a_cell,dim=2))
      case (4) ! tet      
        call truchas_tet_to_irl(a_cell, this%IRL_tet)
        call matchVolumeFraction(this%IRL_tet, a_vof(1), a_planar_separator)
      
      case (5) ! pyramid
        call truchas_pyramid_to_irl(a_cell, this%IRL_pyramid)
        call matchVolumeFraction(this%IRL_pyramid, a_vof(1), a_planar_separator)
      
      case (6) ! Wedge
        call truchas_wedge_to_irl(a_cell, this%IRL_wedge)
        call matchVolumeFraction(this%IRL_wedge, a_vof(1), a_planar_separator)
    
      case (8) ! Hex
        call truchas_hex_to_irl(a_cell, this%IRL_hex)
        call matchVolumeFraction(this%IRL_hex, a_vof(1), a_planar_separator)
    
      case default
        call TLS_fatal('Unknown Truchas cell type during plane-distance setting')
      end select  
    
  end subroutine adjust_planes_match_VOF
  
  subroutine compute_node_velocities(this, a_cell_centered_vel)
  
      class(unsplit_geometric_volume_tracker), intent(inout) :: this
      real(r8), intent(in) :: a_cell_centered_vel(:,:)
      
      integer :: j, nodes_in_cell, n
      real(r8) :: displacement(3), weight
      
      ! Compute velocity at each node as the surface-area weighted average of 
      ! all connected face velocities.
      this%w_node = 0.0_r8
      do j = 1, this%mesh%ncell      
        nodes_in_cell = this%mesh%xcnode(j+1) - this%mesh%xcnode(j)
        associate( node_id => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))

        do n = 1, nodes_in_cell
          displacement = this%mesh%x(:,node_id(n)) - this%mesh%cell_centroid(:,j)
          ! QUESTION : What would be the correct weighting here? 
          ! Something based on inverse distance to centroid?
          weight = this%mesh%volume(j) ! Weighting
          this%w_node(1:3,node_id(n)) = this%w_node(1:3,node_id(n)) + a_cell_centered_vel(:,j) * weight
          this%w_node(4  ,node_id(n)) = this%w_node(4  ,node_id(n)) + weight
        end do
        end associate
      
      end do
      
      ! Now normalize node velocities to get average value
      ! Question : Will need to sum up globally the this%w_node before averaging
      do n = 1, this%mesh%nnode_onP
        this%w_node(1:3,n) = this%w_node(1:3,n) / (this%w_node(4, n)+ tiny(1.0_r8))
      end do
      
      call gather_boundary(this%mesh%node_ip, this%w_node)
      
  end subroutine compute_node_velocities
  
  subroutine compute_fluxes(this, a_face_vel, a_dt, a_vof_n)
  
    use irl_interface_helper
    
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: a_face_vel(:)
    real(r8), intent(in) :: a_dt
    real(r8), intent(in) :: a_vof_n(:,:)
        
    integer :: f, v, i, t, k
    integer ::  number_of_nodes, neighbor_cell_index
    real(r8) :: cell_nodes(3,9), phase_volume(this%nmat), cell_volume
    type(TagAccVM_SepVol_type) :: tagged_sepvol
    integer :: current_tag
    
    ! FOR NOW, ADVECT EVERYWHERE.
    ! IN FUTURE, MAKE BAND AND JUST ADVECT IN THERE

    ASSERT(this%nmat == 2)

    call new(tagged_sepvol)
    call getMoments_setMethod(1)
    
    do f = 1, this%mesh%nface_onP
    
      number_of_nodes = this%mesh%xfnode(f+1)-this%mesh%xfnode(f)
      ! Grab face nodes we will extrude from
      ! Will place face nodes in latter half of cell_nodes so that
      ! if (dot_product(flux_vel, face_norm) > 0), the volume is > 0
      associate( node_index => this%mesh%fnode( &
                                 this%mesh%xfnode(f):this%mesh%xfnode(f+1)-1))
        cell_nodes(:,number_of_nodes+1:2*number_of_nodes) = &
        this%mesh%x(:,node_index(1:number_of_nodes))
        
        do v = 1, number_of_nodes
          cell_nodes(:,v) = this%project_vertex(cell_nodes(:,number_of_nodes+v), node_index(v), -a_dt)
        end do
                                                                 
        neighbor_cell_index = this%mesh%fcell(1, f)
        if(this%mesh%fcell(2,f) == 0) then
          cell_volume = this%mesh%volume(neighbor_cell_index)
        else
          cell_volume = 0.5_r8*(this%mesh%volume(neighbor_cell_index) + &
                                this%mesh%volume(this%mesh%fcell(2,f)))
        end if

        call setMinimumVolToTrack(cell_volume*1.0e-15_r8)    
        
        select case(number_of_nodes)
          case(3) ! Triangular Face -> Octahedron Volume
          
            call truchas_octa_to_irl(cell_nodes, this%IRL_octahedron)
            call getMoments(this%IRL_octahedron, &
                            this%localized_separator_link(neighbor_cell_index), &
                            tagged_sepvol)
            
            this%face_flux(:,f) = phase_volume          
          case(4) ! Quad Face -> Dodecahedron Volume
            
            call truchas_dod_to_irl(cell_nodes, this%IRL_dodecahedron)
            call getMoments(this%IRL_dodecahedron, &
                            this%localized_separator_link(neighbor_cell_index), &
                            tagged_sepvol)
            
            this%face_flux(:,f) = phase_volume
            
          case default
            call TLS_fatal('Face with unhandled amount of nodes found.')          
        end select
          
      end associate


      phase_volume = 0.0_r8
      do t = 0, getSize(tagged_sepvol)-1
         current_tag = getTagForIndex(tagged_sepvol, t)
         if(current_tag > this%mesh%ncell) then            
            ! Is a boundary condition, all volume in phase 0
            phase_volume(:) = &
                 phase_volume(:) + this%getBCMaterialFractions(current_tag, a_vof_n) * getVolumeAtIndex(tagged_sepvol, t, 0)
         else
            ! Is inside domain, trust actual volumes
            phase_volume(1) = phase_volume(1) + getVolumeAtIndex(tagged_sepvol, t, 0)
            phase_volume(2) = phase_volume(2) + getVolumeAtIndex(tagged_sepvol, t, 1)
         end if
      end do    
                     
      this%face_flux(:,f) = phase_volume
          
    end do
    
  end subroutine compute_fluxes

  pure function project_vertex(this, a_pt, a_node_index, a_dt) result(a_proj_pt)

    class(unsplit_geometric_volume_tracker), intent(in) :: this
    real(r8), intent(in) :: a_pt(3)
    integer, intent(in) :: a_node_index
    real(r8), intent(in) :: a_dt
    real(r8) :: a_proj_pt(3)
    
    a_proj_pt = a_pt + a_dt * this%w_node(1:3,a_node_index)
    return
  end function project_vertex

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
  
  subroutine update_vof(this, a_flux_vol, a_vof)
  
    class(unsplit_geometric_volume_tracker), intent(in) :: this
    real(r8), intent(inout) :: a_flux_vol(:,:)
    real(r8), intent(inout) :: a_vof(:,:)
    
    integer :: j, f
    integer :: number_of_faces
    real(r8) :: cell_volume_sum
    
    ! Fill the ragged a_flux_vol array
    do j = 1, this%mesh%ncell  
    
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
    
    ! Use a_flux_vol to update a_vof now
    do j = 1, this%mesh%ncell_onP
      a_vof(:,j) = a_vof(:,j) * abs(this%mesh%volume(j))    
      number_of_faces = this%mesh%xcface(j+1) - this%mesh%xcface(j)
      cell_volume_sum = this%mesh%volume(j)
      do f = 1, number_of_faces
        a_vof(:,j) = a_vof(:,j) - a_flux_vol(:, this%mesh%xcface(j)+f-1)
        cell_volume_sum = cell_volume_sum - sum(a_flux_vol(:, this%mesh%xcface(j)+f-1))
      end do     
      a_vof(:,j) = a_vof(:,j) / cell_volume_sum
    end do
    
    call gather_boundary(this%mesh%cell_ip, a_vof)
    
  end subroutine update_vof    
  
  subroutine clean_vof(this, a_vof)
      
    class(unsplit_geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(inout)  :: a_vof(:,:)
    
    integer :: j, k
    logical :: vof_modified
    
    do j = 1, this%mesh%ncell
      vof_modified = .false.
      do k = 1, this%nmat
        if(a_vof(k,j) < this%cutoff) then
          vof_modified = .true.
          a_vof(k,j) = 0.0_r8 
          ! TODO Also update IRL planar_separator plane to reflect this
        else if(a_vof(k,j) > 1.0_r8 - this%cutoff) then
          vof_modified = .true.
          a_vof(k,j) = 1.0_r8
          ! TODO Also update IRL planar_separator plane to reflect this
        end if
      end do
    
      if(vof_modified) then
        ! Rescale VOF to be valid (sum to 1)
        a_vof(:,j) = a_vof(:,j) / sum(a_vof(:,j))
      end if
    end do
    
    
  end subroutine clean_vof

end module unsplit_geometric_volume_tracker_type
