!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module geometric_volume_tracker_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use volume_tracker_class
  use truchas_logging_services
  use truchas_timers
  use unstr_mesh_type
  use index_partitioning
  implicit none
  private

  type, extends(volume_tracker), public :: geometric_volume_tracker
    private
    type(unstr_mesh), pointer :: mesh ! unowned reference
    integer :: location_iter_max ! maximum number of iterations to use in fitting interface
    integer :: subcycles
    logical :: nested_dissection
    real(r8), allocatable :: flux_vol_sub(:,:), normal(:,:,:)
    ! node/face/cell workspace
    real(r8), allocatable :: w_node(:,:), w_face(:,:), w_cell(:,:,:)
    integer, allocatable :: priority(:), bc_index(:), local_face(:), inflow_mat(:)
    integer :: nrealfluid, nfluid, nmat ! # of non-void fluids, # of fluids incl. void, # of materials
  contains
    procedure :: init
    procedure :: flux_volumes
    procedure :: set_inflow_material
    procedure :: write_interface
    procedure, private :: normals
    procedure, private :: donor_fluxes
    procedure, private :: donor_fluxes_nd_cell
    procedure, private :: donor_fluxes_os_cell
    procedure, private :: flux_renorm
    procedure, private :: flux_acceptor
    procedure, private :: flux_bc
    procedure, private :: accumulate_volume
    procedure, private :: enforce_bounded_vof
  end type geometric_volume_tracker

contains

  subroutine init(this, mesh, nrealfluid, nfluid, nmat, liq_matid, params)

    use parameter_list_type
    use property_module, only: get_truchas_material_id
    use f08_intrinsics, only: findloc

    class(geometric_volume_tracker), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    integer, intent(in) :: nrealfluid, nfluid, nmat, liq_matid(:)
    type(parameter_list), intent(inout) :: params
    integer :: i, j, k

    this%mesh => mesh
    this%nrealfluid = nrealfluid
    this%nfluid = nfluid
    this%nmat = nmat

    call params%get('location_iter_max', this%location_iter_max, default=40)
    call params%get('cutoff', this%cutoff, default=1.0e-8_r8)
    call params%get('subcycles', this%subcycles, default=2)
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

    allocate(this%flux_vol_sub(this%nmat,size(mesh%cface)))
    allocate(this%normal(3,this%nmat,mesh%ncell))
    allocate(this%w_node(2,mesh%nnode))
    allocate(this%w_face(this%nmat,mesh%nface))
    ! need this array so we can violate conservation in parallel
    allocate(this%w_cell(8,this%nfluid,mesh%ncell))

    ! list of boundary face ids
    j = count(this%mesh%fcell(2,:this%mesh%nface_onP) == 0)
    allocate(this%bc_index(j), this%local_face(j), this%inflow_mat(j))
    j = 1
    do i = 1, this%mesh%nface_onP
      if (this%mesh%fcell(2,i) == 0) then
        this%bc_index(j) = i
        k = this%mesh%fcell(1,i)
        this%local_face(j) = this%mesh%xcface(k) - 1 + &
            findloc(this%mesh%cface(this%mesh%xcface(k):this%mesh%xcface(k+1)-1), i)
        j = j + 1
      end if
    end do
    this%inflow_mat = 0

  end subroutine init

  ! flux volumes routine assuming vel/flux_vol is a cface-like array
  subroutine flux_volumes(this, vel, vel_cc, vof_n, vof, flux_vol, fluids, void, dt, a_interface_band)
    class(geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: vel(:), vel_cc(:,:), vof_n(:,:), dt
    real(r8), intent(out) :: flux_vol(:,:), vof(:,:)
    integer, intent(in) :: fluids, void
    integer, intent(in) :: a_interface_band(:)

    integer :: i
    real(r8) :: sub_dt

    flux_vol = 0.0_r8
    sub_dt = dt/real(this%subcycles, r8)
    vof = vof_n

    do i = 1, this%subcycles
      call this%normals(vof)

      call this%donor_fluxes(vel, vof, sub_dt)

      call this%flux_renorm(vel, vof_n, flux_vol, sub_dt)

      call this%flux_acceptor()

      call this%flux_bc(vel, vof_n, sub_dt)

      call this%accumulate_volume(vof, flux_vol)

      call this%enforce_bounded_vof(vof, flux_vol, fluids, void)
    end do

  end subroutine flux_volumes

  subroutine write_interface(this)
    class(geometric_volume_tracker), intent(in) :: this    
    return
  end subroutine write_interface
  
  !! Set the inflow material for the given boundary faces. A material index
  !! of 0 will result in materials being fluxed in proportion to the material
  !! fractions present in the cell. This is the preset default. The BC_INDEX
  !! array component is ordered, so we use a binary search to locate the faces
  !! in the the array, so we can set the corresponding inflow material index.
  !! TODO: If FACES is also ordered (likely) the search can be improved further.

  subroutine set_inflow_material(this, mat, faces)
    class(geometric_volume_tracker), intent(inout) :: this
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

  subroutine flux_bc(this, vel, vof_n, dt)

    class(geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: vel(:), vof_n(:,:), dt

    integer :: i, f, j, fl, m

    do i = 1, size(this%bc_index)
      f = this%bc_index(i)
      fl = this%local_face(i)
      j = this%mesh%fcell(1,f)
      if (j > this%mesh%ncell_onP) cycle
      if (vel(fl) < 0) then
        if (this%inflow_mat(i) > 0) then
          this%flux_vol_sub(this%inflow_mat(i),fl) = vel(fl) * dt * this%mesh%area(f)
        else ! flux material in proportion to starting fractions in cell
          this%flux_vol_sub(:,fl) = vof_n(:,j) * vel(fl) * dt * this%mesh%area(f)
        end if
      end if
    end do

  end subroutine flux_bc


  subroutine normals(this, vof)

    use flow_operators, only: gradient_cc
    use f08_intrinsics, only: findloc
    intrinsic :: norm2

    class(geometric_volume_tracker), intent(inout) :: this
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

      ! perform onion skin if requested
      ! normal vectors (gradients) include previous materials in the priority ordering
      if (.not.this%nested_dissection .and. c > 2) then
        do j = 2, this%nmat
          k = this%priority(j)
          this%normal(:,k,i) = this%normal(:,k,i) + this%normal(:,this%priority(j-1),i)
        end do
      end if

      ! normalize and remove smallish components due to robustness issues in nested disection
      do j = 1 , this%nmat
        !if (vof(j,i) <= this%cutoff) cycle ! should this be 0 or cutoff?
        ! remove small values
        do k = 1, 3
          if (abs(this%normal(k,j,i)) < epsilon(1.0_r8)) this%normal(k,j,i) = 0.0_r8
        end do
        ! normalize if possible
        mag = norm2(this%normal(:,j,i))
        if (mag > epsilon(1.0_r8)) this%normal(:,j,i) = this%normal(:,j,i)/mag
        ! remove slightly larger values
        do k = 1, 3
          if (abs(this%normal(k,j,i)) < 1.0e-6_r8) this%normal(k,j,i) = 0.0_r8
        end do
        ! normalize if possible
        mag = norm2(this%normal(:,j,i))
        if (mag > epsilon(1.0_r8)) then
          this%normal(:,j,i) = this%normal(:,j,i)/mag
        else
          this%normal(:,j,i) = 1.0_r8
        end if
      end do
    end do
    ! will need normals for vof reconstruction in ghost cells
    call gather_boundary(this%mesh%cell_ip, this%normal)

    call stop_timer('normals')

  end subroutine normals

  subroutine donor_fluxes_os_cell(this, i, vel, vof, dt)

    use cell_geom_type

    class(geometric_volume_tracker), intent(inout) :: this
    integer, intent(in) :: i
    real(r8), intent(in)  :: dt, vof(:,:), vel(:)

    real(r8) :: face_normal(3,6)
    integer :: j,k,ierr
    type(cell_geom) :: cell

    associate (cn => this%mesh%cnode(this%mesh%xcnode(i):this%mesh%xcnode(i+1)-1), &
        fi => this%mesh%cface(this%mesh%xcface(i):this%mesh%xcface(i+1)-1))

      do j = 1, size(fi)
        k = fi(j)
        if (btest(this%mesh%cfpar(i),pos=j)) then
          face_normal(:,j) = -this%mesh%normal(:,k)/this%mesh%area(k)
        else
          face_normal(:,j) = this%mesh%normal(:,k)/this%mesh%area(k)
        end if
      end do

      call cell%init(this%mesh%x(:,cn), this%mesh%volume(i), this%mesh%area(fi), &
          face_normal(:,1:size(fi)))

      call cell_volume_flux(dt, cell, vof(:,i), this%normal(:,:,i), &
          vel(this%mesh%xcface(i):this%mesh%xcface(i+1)-1), this%cutoff, &
          this%priority, this%nmat, this%location_iter_max, &
          this%flux_vol_sub(:,this%mesh%xcface(i):this%mesh%xcface(i+1)-1))
    end associate

  end subroutine donor_fluxes_os_cell

  ! get the volume flux for every material in the given cell
  subroutine cell_volume_flux(dt, cell, vof, int_norm, vel, cutoff, priority, nmat, maxiter, &
      flux_volume)

    use f08_intrinsics, only: findloc
    use locate_plane_os_function
    use plane_type
    use cell_geom_type

    real(r8), intent(in) :: dt, int_norm(:,:), vof(:), vel(:), cutoff
    integer, intent(in) :: priority(:), nmat, maxiter
    type(cell_geom), intent(in) :: cell
    real(r8), intent(out) :: flux_volume(:,:)

    real(r8) :: Vofint, dvol
    real(r8) :: flux_vol_sum(cell%nfc), flux_vol
    integer :: ni,f,nlast, nint, ierr,nmat_in_cell
    logical :: is_mixed_donor_cell
    type(plane) :: P

    flux_volume = 0.0_r8
    flux_vol_sum = 0.0_r8
    nmat_in_cell = count(vof > 0.0_r8)
    !nint = count(vof > 0.0_r8)
    ! Here, I am not certain the conversion from pri_ptr to direct material indices worked properly.
    ! This will be clear when trying 3 or more materials. -zjibben

    ! Loop over the interfaces in priority order
    do ni = 1,nmat-1
      ! check if this is a mixed material cell
      ! First accumulate the volume fraction of this material and materials with lower priorities.
      ! Force 0.0 <= Vofint <= 1.0
      Vofint = min(max(sum(vof(priority(:ni))), 0.0_r8), 1.0_r8)
      is_mixed_donor_cell = cutoff < Vofint .and. Vofint < 1 - cutoff
      ! locate each interface plane by computing the plane constant
      if (is_mixed_donor_cell) &
          P = locate_plane_os(int_norm(:,priority(ni)), vofint, cell%volume, cell%node, &
          cutoff, maxiter)

      ! calculate delta advection volumes for this material at each donor face and accumulate the sum
      ASSERT(priority(ni) <= size(flux_volume, dim=1))
      call compute_material_volume_flux(flux_volume(priority(ni),:), flux_vol_sum, P, cell, &
          is_mixed_donor_cell, vel, dt, vof(priority(ni)), cutoff)
    end do

    ! Compute the advection volume for the last material.
    nlast = priority(findloc(vof(priority) >= cutoff, .true., back=.true.))
    ASSERT(nlast <= size(flux_volume, dim=1))
    do f = 1,cell%nfc
      ! Recalculate the total flux volume for this face.
      flux_vol = dt*vel(f)*cell%face_area(f)
      if (abs(flux_vol) > 0.5_r8 * cell%volume) then
        write(*,*) dt,flux_vol,cell%volume,flux_vol/cell%volume
        call TLS_fatal('advection timestep too large')
      end if
      if (flux_vol <= cutoff*cell%volume) cycle

      ! For donor cells containing only one material, assign the total flux.
      if (nmat_in_cell==1) then
        flux_volume(nlast,f) = flux_vol
      else
        ! The volume flux of the last material shouldn't be less than
        ! zero nor greater than the volume of this material in the donor cell.
        dvol = min(max(abs(flux_vol - Flux_Vol_Sum(f)), 0.0_r8), Vof(nlast)*cell%volume)

        ! Store the last material's volume flux.
        if (dvol > cutoff*cell%volume) flux_volume(nlast,f) = dvol
      end if
    end do

  end subroutine cell_volume_flux

  ! calculate the flux of one material in a cell
  subroutine compute_material_volume_flux(material_volume_flux, flux_vol_sum, P, cell, &
      is_mixed_donor_cell, vel, dt, vof, cutoff)

    use cell_geom_type
    use truncation_volume_type
    use plane_type

    real(r8), intent(out) :: material_volume_flux(:)
    real(r8), intent(inout) :: flux_vol_sum(:)
    type(plane), intent(in) :: P
    type(cell_geom), intent(in) :: cell
    logical, intent(in) :: is_mixed_donor_cell
    real(r8), intent(in) :: vel(:), dt, vof, cutoff

    integer :: f
    real(r8) :: vp, flux_vol, flux_vol_node(3,8)
    type(truncation_volume) :: trunc_vol

    material_volume_flux = 0
    do f = 1,cell%nfc
      ! Flux volumes
      flux_vol = dt*vel(f)*cell%face_area(f)
      if (flux_vol <= cutoff*cell%volume) cycle

      if (is_mixed_donor_cell) then
        ! calculate the vertices describing the volume being truncated through the face
        call flux_vol_nodes(f, cell, vel(f)*dt, flux_vol, cutoff, flux_vol_node)

        ! compute the volume truncated by interface planes in each flux volumes
        call trunc_vol%init(flux_vol_node, P%normal)
        Vp = trunc_vol%volume(P%rho)
      else
        ! For clean donor cells, the entire Flux volume goes to the single donor material.
        Vp = merge(abs(flux_vol), 0.0_r8, vof >= 1-cutoff)
      end if

      ! If Vp is close to 0 set it to 0.  If it is close
      ! to 1 set it to 1. This will avoid numerical round-off.
      if (Vp > (1-cutoff)*abs(flux_vol)) Vp = abs(flux_vol)

      ! Make sure that the current material-integrated advection
      ! volume hasn't decreased from its previous value.  (This
      ! can happen if the interface significantly changed its
      ! orientation and now crosses previous interfaces.)  Also
      ! limit Volume_Flux_Sub to take no more than the material
      ! occupied volume in the donor cell.
      material_volume_flux(f) = min(max(Vp - flux_vol_sum(f), 0.0_r8), Vof*cell%volume)
      flux_vol_sum(f) = flux_vol_sum(f) + material_volume_flux(f)
    end do

  end subroutine compute_material_volume_flux

  ! Given the value of Flux_Vol (the volume of material that moves
  ! through the current advection cell face), find the vertices which
  ! describe this volume.  Four of these vertices will be the ones that
  ! describe the advection cell face. The other four vertices will lie
  ! approximately "DIST" away from the advection cell face. These are
  ! only approximately "DIST" away because the cell cross-sectional area
  ! may increase or decrease as one moves away from the advection cell
  ! face.  The value used is varied from "DIST" such that the vertices
  ! describe a hexagonal volume that matches the value of Flux_Vol.

  subroutine flux_vol_nodes(face, cell, dist, Flux_Vol, cutvof, flux_vol_node)

    use cell_geometry, only: hex_volume
    use cell_geom_type

    integer, intent(in) :: face
    type(cell_geom), intent(in) :: cell
    real(r8), intent(in) :: dist, cutvof, flux_vol
    real(r8), intent(out) :: flux_vol_node(:,:)

    integer, parameter :: flux_vol_iter_max = 10

    ! this array maps the nodes on a given cell's face
    ! to the nodes of a flux volume. The flux volume
    ! will either be a hex or a wedge, depending on
    ! if we are fluxing through a triangular or
    ! quadrilateral face. Example:
    !
    ! edge_ends_type(i,j,k,l) refers to a flux volume's
    ! node ID corresponding to the jth node of the kth
    ! face of a given cell of type l. i=1 refers to the
    ! node coinciding with the given node, and i=2
    ! refers to a node on the opposite end of the flux
    ! volume.
    integer, parameter :: &
        back_nodes(4,6,4) = reshape([&
        ! tet
        3, 3, 3, 3, &
        1, 1, 1, 1, &
        2, 2, 2, 2, &
        4, 4, 4, 4, &
        0, 0, 0, 0, &
        0, 0, 0, 0, &

        ! pyramid
        4, 3, 3, 4, &
        1, 4, 4, 1, &
        2, 1, 1, 2, &
        2, 2, 3, 3, &
        5, 5, 5, 5, &
        0, 0, 0, 0, &

        ! wedge
        3, 3, 6, 6, &
        1, 1, 4, 4, &
        2, 5, 5, 2, &
        4, 6, 5, 5, &
        1, 2, 3, 3, &
        0, 0, 0, 0, &

        ! hex
        4, 3, 7, 8, &
        1, 4, 8, 5, &
        2, 1, 5, 6, &
        2, 6, 7, 3, &
        5, 8, 7, 6, &
        1, 2, 3, 4],&
        [4,6,4])

    ! this is identical to the <shape>_faces arrays in cell_topology,
    ! except that triangular faces are treated like a degenerate quadrilaterals
    integer, parameter :: &
        front_nodes(4,6,4) = reshape([&
        ! tet
        1,2,4,4, &
        2,3,4,4, &
        1,4,3,3, &
        1,3,2,2, &
        0,0,0,0, &
        0,0,0,0, &

        ! pyramid
        1,2,5,5, &
        2,3,5,5, &
        3,4,5,5, &
        1,5,5,4, &
        1,4,3,2, &
        0,0,0,0, &

        ! wedge
        1,2,5,4, &
        2,3,6,5, &
        1,4,6,3, &
        1,3,2,2, &
        4,5,6,6, &
        0,0,0,0, &

        ! hex
        1,2,6,5, &
        2,3,7,6, &
        3,4,8,7, &
        1,5,8,4, &
        1,4,3,2, &
        5,6,7,8],&
        [4,6,4])

    integer, parameter :: nedge = 4, edge_ends(2,4) = reshape([1,5,  4,8,  3,7, 2,6], [2,4])

    integer :: ia, ib, e, iter
    real(r8) :: percnt(4), Uedge(3,4), volume, mult, ndotuedge

    ASSERT(cell%cell_type > 0 .and. cell%cell_type <= 4)

    ! initialize flux volume as entire given cell
    flux_vol_node = 0
    flux_vol_node(:,edge_ends(1,:)) = cell%node(:,front_nodes(:,face,cell%cell_type))
    flux_vol_node(:,edge_ends(2,:)) = cell%node(:,back_nodes(:,face,cell%cell_type))

    ! compute the edge unit vectors
    percnt = 0
    do e = 1, nedge
      ia  = edge_ends(1,e) ! front node
      ib  = edge_ends(2,e) ! back node
      Uedge(:,e) = flux_vol_node(:,ib) - flux_vol_node(:,ia)

      ndotuedge = dot_product(cell%face_normal(:,face), Uedge(:,e))
      Percnt(e) = -Dist/(ndotuedge+epsilon(1.0_r8))
    end do

    if (any(percnt < 0) .or. any(percnt > 1)) &
        call TLS_fatal('FLUX_VOL_NODE: invalid flux volume or inverted element')

    ! iterate to find the four vertices at the back end of the flux volume
    mult = 1
    do iter = 1, flux_vol_iter_max
      ! adjust the back nodes
      do e = 1, nedge
        ia = edge_ends(1,e)
        ib = edge_ends(2,e)
        flux_vol_node(:,ib) = flux_vol_node(:,ia) + Mult*Percnt(e)*Uedge(:,e)
      end do

      ! compute the new flux volume and increment multiplier for next iteration
      volume = hex_volume(flux_vol_node)
      Mult = Mult * flux_vol/volume

      if (abs(flux_vol - volume) < cutvof*cell%volume) exit
    end do

    if (iter > flux_vol_iter_max) &
        call TLS_fatal('Flux volume vertex iteration did not converge')

  end subroutine flux_vol_nodes

  subroutine donor_fluxes_nd_cell(this, i, vel, vof, dt)

    use cell_topology
    use multimat_cell_type

    class(geometric_volume_tracker), intent(inout) :: this
    integer, intent(in) :: i
    real(r8), intent(in)  :: dt, vof(:,:), vel(:)

    real(r8) :: face_normal(3,6)
    integer :: j,k,ierr, face_vid(4,6)
    type(multimat_cell) :: cell

    associate (cn => this%mesh%cnode(this%mesh%xcnode(i):this%mesh%xcnode(i+1)-1), &
        fi => this%mesh%cface(this%mesh%xcface(i):this%mesh%xcface(i+1)-1))

      do j = 1, size(fi)
        k = fi(j)
        if (btest(this%mesh%cfpar(i),pos=j)) then
          face_normal(:,j) = -this%mesh%normal(:,k)/this%mesh%area(k)
        else
          face_normal(:,j) = this%mesh%normal(:,k)/this%mesh%area(k)
        end if
      end do

      select case (size(cn))
      case (4)
        call cell%init(ierr, this%mesh%x(:,cn), face_normal(:,1:size(fi)), this%mesh%volume(i))

      case (5)
        ! zero treated as sentinel value in multimat_cell procedures
        face_vid = 0
        do j = 1, size(fi)
          face_vid(1:PYR5_FSIZE(j),j) = PYR5_FACES(PYR5_XFACE(j):PYR5_XFACE(j+1)-1)
        end do
        call cell%init(ierr, this%mesh%x(:,cn), face_vid(:,1:size(fi)), &
            PYR5_EDGES, face_normal(:,1:size(fi)), this%mesh%volume(i))

      case (6)
        ! zero treated as sentinel value in multimat_cell procedures
        face_vid = 0
        do j = 1, size(fi)
          face_vid(1:WED6_FSIZE(j),j) = WED6_FACES(WED6_XFACE(j):WED6_XFACE(j+1)-1)
        end do
        call cell%init(ierr, this%mesh%x(:,cn), face_vid(:,1:size(fi)), &
            WED6_EDGES, face_normal(:,1:size(fi)), this%mesh%volume(i))

      case (8)
        call cell%init(ierr, this%mesh%x(:,cn), reshape(source=HEX8_FACES,shape=[4,6]), &
            HEX8_EDGES, face_normal(:,1:size(fi)), this%mesh%volume(i))

      case default
        call TLS_fatal('unaccounted topology in donor_fluxes_nd')
      end select

      if (ierr /= 0) call TLS_fatal('cell_outward_volflux failed: could not initialize cell')

      call cell%partition(vof(:,i), this%normal(:,:,i), this%cutoff, this%priority, &
          this%location_iter_max)

      this%flux_vol_sub(:this%nfluid,this%mesh%xcface(i):this%mesh%xcface(i+1)-1) = &
          cell%outward_volflux(dt, vel(this%mesh%xcface(i):this%mesh%xcface(i+1)-1),&
          this%mesh%area(fi), this%cutoff, this%nfluid, ierr)
      if (ierr /= 0) call TLS_fatal('cell_outward_volflux failed')
    end associate

  end subroutine donor_fluxes_nd_cell

  subroutine donor_fluxes(this, vel, vof, dt)

    class(geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: dt, vof(:,:), vel(:)

    integer :: i, nmat

    ! calculate the flux volumes for each face
    call start_timer('reconstruct/advect')

    do i = 1, this%mesh%ncell
      nmat = count(vof(:,i) > this%cutoff)

      if (nmat > 2 .and. this%nested_dissection) then
        call this%donor_fluxes_nd_cell(i, vel, vof, dt)
      else
        call this%donor_fluxes_os_cell(i, vel, vof, dt)
      end if
    end do

    call stop_timer('reconstruct/advect')

  end subroutine donor_fluxes


  ! Scan all faces with an outward flux and determine if any material is over-exhausted from this
  ! cell.  If so lower the fluxes until the material is just exhausted. Then loop over the faces
  ! and balance the individual material fluxes with the total face flux. The sum of the material
  ! volume fluxes (Volume_Flux_Sub) for each face should sum to the total volume flux for that
  ! face.  This balancing is an iterative adjust procedure such that by the end, two criteria are
  ! satisfied: 1) the sum of material fluxes from each face of the cell equals the total face
  ! volume flux and 2) The cumulative sum of individual material fluxes (from this and previous
  ! volume_track_subcycles) does not exceed the volume of a particular material originally within a
  ! cell.  To make this happen, we decrease some fluxes to equal the volume of material still
  ! available to be fluxed, and increase other fluxes appropriately.
  !
  ! This routine assumes that the total fluxed volume is dt*v*face_area.  This is obviously not
  ! true for partially solid/void cells.  In these cases, the fluxes will be balanced upwards
  ! towards a ficticiously large flux volume.  This routine also permits the creation of volume
  ! fluxes for materials which may not be present in the cell.

  subroutine flux_renorm(this, vel, vof_n, flux_vol, dt)

    class(geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: vel(:), vof_n(:,:), flux_vol(:,:), dt

    integer :: i,j,o,m,f0,f1,navail,nmat, ierr
    logical  :: adjust_fluxes, avail(size(vof_n,dim=1))
    real(r8) :: mat_flux_cur, mat_flux_acc, mat_avail
    real(r8) :: face_flux_fixed, face_flux_adjustable, face_flux

    nmat = this%nfluid
    ierr = 0

    do i = 1, this%mesh%ncell_onP
      avail = .true.
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1

      do o = 1, nmat+1
        adjust_fluxes = .false.

        do m = 1, nmat
          mat_flux_cur = sum(this%flux_vol_sub(m,f0:f1))
          if (mat_flux_cur == 0.0_r8) cycle

          mat_flux_acc = 0.0_r8
          do j = f0, f1
            if (flux_vol(m,j) > 0.0_r8) mat_flux_acc = mat_flux_acc + flux_vol(m,j)
          end do
          mat_avail = vof_n(m,i)*this%mesh%volume(i)
          ! check/correct for overflow
          if ((mat_flux_acc + mat_flux_cur) > mat_avail) then
            adjust_fluxes = .true.
            avail(m) = .false.
            ! rescale so that mat_flux_acc + mat_flux_cur == mat_avail
            this%flux_vol_sub(m,f0:f1) = this%flux_vol_sub(m,f0:f1) * &
                (mat_avail-mat_flux_acc)/mat_flux_cur
          end if
        end do

        if (.not.adjust_fluxes) exit

        do j = f0, f1
          face_flux = dt*vel(j)*this%mesh%area(this%mesh%cface(j))
          if (face_flux > this%cutoff*this%mesh%volume(i)) then
            face_flux_fixed = 0.0_r8
            face_flux_adjustable = 0.0_r8
            do m = 1, nmat
              if (avail(m)) then
                face_flux_adjustable = face_flux_adjustable + this%flux_vol_sub(m,j)
              else
                face_flux_fixed = face_flux_fixed + this%flux_vol_sub(m,j)
              end if
            end do

            if (face_flux_adjustable > 0.0_r8) then
              ! rescale so face_flux_adj+face_flux_fixed == face_flux
              do m = 1, nmat
                if (avail(m)) this%flux_vol_sub(m,j) = this%flux_vol_sub(m,j) * &
                    (face_flux - face_flux_fixed)/face_flux_adjustable
              end do
            else
              navail = count(avail)
              if (navail == 0) then
                ierr = 1
                exit
              end if

              ! arbitrarily add volume flux to potentially non-existent material in cell to balance
              ! equations.  This seems really, really, really, really, really, really bad.
              do m = 1, nmat
                if (avail(m)) this%flux_vol_sub(m,j) = (face_flux-face_flux_fixed)/real(navail,r8)
              end do
            end if
          end if
        end do
      end do
    end do


    call TLS_fatal_if_any (ierr /= 0, 'FLUX_RENORM: cannot reassign face flux to any other material')

    ! Is there really not a better way?
    do i = 1, this%mesh%ncell_onP
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1
      do m = 1, nmat
        this%w_cell(1:f1-f0+1,m,i) = this%flux_vol_sub(m,f0:f1)
      end do
    end do
    call gather_boundary(this%mesh%cell_ip, this%w_cell)
    do i = this%mesh%ncell_onP+1, this%mesh%ncell
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1
      do m = 1, nmat
        this%flux_vol_sub(m,f0:f1) = this%w_cell(1:f1-f0+1,m,i)
      end do
    end do


  end subroutine flux_renorm

  ! On entrance, this%flux_vol_sub only contains outward (i.e. positive) flux volumes.  On exit,
  ! this%flux_vol_sub is made consistent with appropriate negative entries

  subroutine flux_acceptor(this)

    class(geometric_volume_tracker), intent(inout) :: this

    integer :: m,j,i,f0,f1,nmat

    nmat = size(this%flux_vol_sub, dim=1)
    this%w_face = 0

    do i = 1, this%mesh%ncell
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1
      do j = f0, f1
        do m = 1, nmat
          if (this%flux_vol_sub(m,j) > 0.0_r8) &
              this%w_face(m,this%mesh%cface(j)) = this%flux_vol_sub(m,j)
        end do
      end do
    end do

    do i = 1, this%mesh%ncell_onP
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1
      do j = f0, f1
        do m = 1, nmat
          if (this%flux_vol_sub(m,j) /= this%w_face(m,this%mesh%cface(j))) &
              this%flux_vol_sub(m,j) = -this%w_face(m,this%mesh%cface(j))
        end do
      end do
    end do

  end subroutine flux_acceptor


  subroutine accumulate_volume(this, vof, flux_vol)

    class(geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(inout) :: flux_vol(:,:), vof(:,:)

    integer :: i,j,f0,f1,m

    do i = 1, this%mesh%ncell_onP
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1

      do j = f0, f1
        do m = 1, this%nfluid
          flux_vol(m,j) = flux_vol(m,j) + this%flux_vol_sub(m,j)
          vof(m,i) = vof(m,i) - this%flux_vol_sub(m,j)/this%mesh%volume(i)
        end do
      end do
    end do

  end subroutine accumulate_volume

  ! Enforce boundedness by allowing _inconsistent_ material flux volumes at faces
  subroutine enforce_bounded_vof(this, vof, flux_vol, fluids, void)

    class(geometric_volume_tracker), intent(inout) :: this
    real(r8), intent(inout) :: vof(:,:), flux_vol(:,:)
    integer, intent(in) :: fluids, void

    integer :: i,m,f0,f1
    real(r8) :: a1, v, q, excess

    a1 = 1.0_r8 - this%cutoff

    do i = 1, this%mesh%ncell_onP
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1

      do m = 1, fluids
        if (vof(m,i) > a1 .and. vof(m,i) /= 1.0_r8) then
          call adjust_flux_matl(flux_vol(m,f0:f1), this%mesh%volume(i)*(1.0_r8-vof(m,i)))
          vof(m,i) = 1.0_r8
        else if (vof(m,i) < this%cutoff .and. vof(m,i) /= 0) then
          call adjust_flux_matl(flux_vol(m,f0:f1), -this%mesh%volume(i)*vof(m,i))
          vof(m,i) = 0.0_r8
        end if
      end do

      do m = fluids+1, fluids+void
        if (vof(m,i) > a1 .and. vof(m,i) /= 1.0_r8) then
          vof(m,i) = 1.0_r8
        else if (vof(m,i) < this%cutoff .and. vof(m,i) /= 0) then
          vof(m,i) = 0.0_r8
        end if
      end do

      v = sum(vof(:,i))
      excess = v - 1.0_r8
      if (excess == 0.0_r8) cycle

      q = sum(vof(fluids+1:fluids+void,i))

      if (q > 0.0_r8) then
         ! we can add or remove enough void from cell
        if (excess < 0.0_r8 .or. (excess > 0.0_r8 .and. q >= excess)) then
          do m = fluids+1, fluids+void
            vof(m,i) = vof(m,i) * (1.0_r8 - excess/q)
          end do
        else ! we cannot remove enough void from cell
          do m = fluids+1, fluids+void
            vof(m,i) = 0.0_r8
          end do
          call adjust_flux_all(flux_vol(:,f0:f1), vof(:,i), q-excess, this%mesh%volume(i), fluids)
        end if
        cycle
      end if

      ! There is no void to adjust in this cell.
      call adjust_flux_all(flux_vol(:,f0:f1), vof(:,i), -excess, this%mesh%volume(i), fluids)
    end do

    call gather_boundary(this%mesh%cell_ip, vof)

    ! Is there really not a better way?
    do i = 1, this%mesh%ncell_onP
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1
      do m = 1, fluids
        this%w_cell(1:f1-f0+1,m,i) = flux_vol(m,f0:f1)
      end do
    end do
    call gather_boundary(this%mesh%cell_ip, this%w_cell)
    do i = this%mesh%ncell_onP+1, this%mesh%ncell
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1
      do m = 1, fluids
        flux_vol(m,f0:f1) = this%w_cell(1:f1-f0+1,m,i)
      end do
    end do

  end subroutine enforce_bounded_vof

  ! Adjust the material volume fluxes on the faces of a single cell to match the evaluated material
  ! volume to a target value. vof_delta is desired_vof-actual_vof
  ! IGNORES BC's FOR NOW
  ! The result of this process is an _inconsitent_ view of material flux volumes

  subroutine adjust_flux_matl(flux_vol, vol_delta)

    real(r8), intent(inout) :: flux_vol(:)
    real(r8), intent(in) :: vol_delta

    integer        :: i
    real(r8)       :: flux_vol_adj, r

    ! both inflow and outflow will be rescaled so take abs
    flux_vol_adj = sum(abs(flux_vol))

    if (flux_vol_adj > 0.0_r8) then
      r = vol_delta/flux_vol_adj

      do i = 1, size(flux_vol)
        if (flux_vol(i) > 0.0_r8) then
          flux_vol(i) = (1.0_r8-r)*flux_vol(i)
        else
          flux_vol(i) = (1.0_r8+r)*flux_vol(i)
        end if
      end do

    else if (vol_delta < 0.0_r8) then
      ! jms Note:   If the material is to be removed from the cell
      ! look for a face that doesn't have incoming material, and
      ! flux it out through that face
      do i = 1, size(flux_vol)
        if (flux_vol(i) >= 0.0_r8) then
          flux_vol(i) = flux_vol(i) - vol_delta
          exit
        end if
      end do
      !
      ! WHY IS THERE NO CHECK THAT THE FLUXES HAVE ACTUALLY BEEN ADJUSTED?
      !
    else
      !
      ! WHY NOT ARBITRARILY FLUX VOLUME IN AS WELL IF NEED BE?
      !
    end if

  end subroutine adjust_flux_matl

  ! Adjust all the material volume fluxes on the faces of a single cell to match the evaluated
  ! material volume to a target value. vof_delta is desired_vof-actual_vof
  ! IGNORES BC's FOR NOW The
  ! result of this process is an _inconsitent_ view of material flux volumes

  subroutine adjust_flux_all(flux_vol, vof, vof_delta, vol, fluids)

    real(r8), intent(inout) :: flux_vol(:,:), vof(:)
    real(r8), intent(in) :: vof_delta, vol
    integer, intent(in) :: fluids

    integer :: i, j
    real(r8) :: flux_vol_adj, excess, r

    ! both inflow and outflow will be rescaled so take abs. Ignore BCs for now
    flux_vol_adj = sum(abs(flux_vol(1:fluids,:)))

    if (flux_vol_adj > 0.0_r8) then
      r = vof_delta/flux_vol_adj

      do j = 1, size(flux_vol,dim=2)
        do i = 1, fluids
          if (flux_vol(i,j) > 0.0_r8) then
            vof(i) = vof(i) + r*flux_vol(i,j)
            flux_vol(i,j) = (1.0_r8-r*vol)*flux_vol(i,j)
          else
            vof(i) = vof(i) - r*flux_vol(i,j)
            flux_vol(i,j) = (1.0_r8+r*vol)*flux_vol(i,j)
          end if
        end do
      end do
    else
      ! truchas prints a warning... why is the matl routine different?
      return
    end if

  end subroutine adjust_flux_all

end module geometric_volume_tracker_type
