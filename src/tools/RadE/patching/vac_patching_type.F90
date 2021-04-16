!!
!! VAC_PATCHING
!!
!! Concrete implementations of the abstract base class PATCHING that encapsulate
!! the Vertex Anchor Clustering (VAC) algorithm and the related Vertex Anchor
!! Cluster Paving algorithm (PAVE) for clustering faces of a radiation enclosure
!! mesh.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! 30 October 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! VAC_PATCHING and PAVE_PATCHING are extensions of the abstract base class
!! PATCHING.  See the base class comments for a description of the common type
!! bound procedures.  Since the parameters of each algorithm are nearly
!! identical, we only list them once, and highlight any differences.  Details
!! specific to each algorithm are discussed in the IMPLEMENTATION NOTES section
!! below.
!!
!!   INIT(E, PARAMS, STAT, ERRMSG) initializes the object.  The ENCL object E is
!!     the radiation enclosure.  The PARAMETER_LIST object PARAMS contains all
!!     other parameters needed to initialize the object.  If an error is
!!     encountered, STAT returns a non-zero value and ERRMSG an explanatory
!!     message.  The following parameters are recognized:
!!
!!       max-angle -- maximum allowable angle in degrees between normals of
!!           adjacent faces.  If two topologically adjacent faces are at an
!!           angle greater than max-angle degrees, they will not be considered
!!           adjacent during patch construction.
!!       verbosity-level -- specifies the verbosity level of messages printed by
!!           the object.  A verbosity level <= 0 suppresses all output.
!!       merge-level -- controls the aggressiveness of patch merging.  A merge
!!           level of 0 indicates no merging.  A merge level >= 3 indicates the
!!           most aggressive merging method.
!!       split-patch-size -- maximum size of patches that are 'split' into
!!           1-face patches during patch merging.  Actually, we queue these
!!           1-face patches along with the original patch.  The 1-face patches
!!           have a large weight, and are only used to 'fill-in the gaps' after
!!           merging.
!!       random-seed -- optional parameter. Only recognized by PAVE. If present,
!!           it is used to initialize the random number generator used to pick
!!           the initial seed patches.  If not present, the seed is taken from
!!           the system clock.
!!
!!
!!  IMPLEMENTATION NOTES
!!
!!  1. These algorithms are based on the idea that the set of all faces sharing a
!!  particular vertex (called the faces of a vertex) tend to form patches with
!!  desirable properties.  Such patches are connected, relatively small, and tend
!!  to be roughly circular.  These patches have an associated 'vertex anchor',
!!  namely the vertex shared by all their faces.  We use the term 'full patches'
!!  to refer to patches consisting of all the faces of their vertex anchor.
!!  'Partial patches', as expected, refer to patches whose faces are a strict
!!  subset of the faces of their vertex anchor.  These algorithms attempt to
!!  maximize the number of full patches generated.
!!
!!  2. VAC_PATCHING: In the VAC algorithm, we iterate through each node and add
!!  its corresponding full patch to a global priority queue.  We then pop the
!!  queue entries one by one until the queue is empty.  If all the faces of a
!!  queue entry are unassigned, we create a new patch from the entry.  Otherwise,
!!  we make a new entry from the subset of faces of the original entry that are
!!  still unassigned, and add the entry to the queue.  When the queue is empty,
!!  we have a valid patching of the enclosure.  Finally, we merge patches where
!!  possible, in accordance with the MERGE_LEVEL parameter
!!
!!  3. PAVE_PATCHING: The PAVE algorithm is very similar to VAC, with one key
!!  difference.  Instead of queueing all the nodes at the beginning, we instead
!!  queue a seed patch in each connected component of the enclosure.  When we pop
!!  an entry from the queue, we add new entries for each of the 1st and 2nd
!!  degree neighbors of the entry's vertex anchor.  This process continues until
!!  the queue is empty.  In this way, we 'pave' each component with full patches,
!!  starting from the component seeds.  For better results, the component seeds
!!  are chosen to be on the corners or edges of the component, if they exist.
!!


#include "f90_assert.fpp"

module vac_patching_type

  use kinds, only: r8
  use re_encl_type
  use vac_min_heap
  use patching_class
  use edge_neighbor_table_type
  use patching_tools, only: PI
  use parameter_list_type
  use graph_type
  implicit none

  !! Parameter defaults
  integer, parameter  :: VAC_SPLIT_PATCH_SIZE_DEFAULT = 3
  integer, parameter  :: VAC_MERGE_LEVEL_DEFAULT = 3

  type, private :: vac_patch
    integer :: node
    integer, allocatable :: face(:)
    real(r8) :: weight
  end type

  type, extends(patching), public :: vac_patching
    private
    type(encl), pointer :: e
    type(min_heap) :: heap
    type(edge_neighbor_table) :: enhbr
    type(graph) :: vgraph  ! Graph of mesh vertices and edges
    type(vac_patch), allocatable :: patch(:)
    real(r8), allocatable :: area(:), normal(:,:)
    integer, allocatable  :: vface(:,:)
    integer, allocatable  :: xfnhbr(:), fnhbr(:)  ! Faces adjacent to faces
    integer, allocatable  :: xvnhbr(:), vnhbr(:)  ! Vertices adjacent to vertices
    integer, allocatable  :: f2p_map(:)
    logical, allocatable  :: boundary_node(:)
    integer :: npatch
    integer :: split_patch_size  ! Maximum patch size to split
    integer :: merge_level  ! Aggressiveness of patch merging
    integer :: verbosity

    contains
      procedure, public  :: init => init_vac
      procedure, public  :: run => run_vac
      procedure, public  :: output
      procedure, private :: name => name_vac
      procedure, private :: print_stats
      procedure, private :: set_patches
      procedure, private :: split_patches
      procedure, private :: merge_patches
      procedure, private :: merge_patches_vertex_anchor
      procedure, private :: merge_patches_vertex_neighbors
      procedure, private :: add_patch
      procedure, private :: delete_patch
      procedure, private :: delete_patches
      procedure, private :: queue_connected_faces
      procedure, private :: compute_weight
  end type

  type, extends(vac_patching), public :: pave_patching
    contains
      procedure, public  :: init => init_pave
      procedure, public :: run => run_pave
      procedure, private :: name => name_pave
      procedure, private :: pick_seeds
      procedure, private :: queue_vertex_neighbors
  end type

contains

  function name_vac(this)
    class(vac_patching), intent(in) :: this
    character(:), allocatable :: name_vac
    name_vac = "VAC"
  end function name_vac

  function name_pave(this)
    class(pave_patching), intent(in) :: this
    character(:), allocatable :: name_pave
    name_pave = "PAVE"
  end function name_pave

  !! Allocate and initialize VAC_PATCHING data
  subroutine init_vac(this, e, params, stat, errmsg)

    use cell_topology, only: get_edge_nodes
    use cell_geometry, only: face_normal, vector_length, normalized
    use patching_tools, only: get_fnhbr_aux, faces_neighboring_vertices

    class(vac_patching), intent(out) :: this
    type(encl), intent(in), target :: e
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(edge_neighbor), allocatable :: nhbrs(:)
    integer, allocatable :: face_nodes(:), edge(:)
    integer :: i, j
    real(r8) :: max_angle  ! Maximum allowable angle for adjacent faces (in degrees)
    real(r8) :: normal(3), angle, max_angle_rad
    character(:), allocatable :: context

    !! Process the parameters.
    context = 'processing ' // params%name() // ': '
    call params%get('max-angle', max_angle, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    call params%get('verbosity-level', this%verbosity, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    call params%get('merge-level', this%merge_level, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    call params%get('split-patch-size', this%split_patch_size, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if

    this%e => e
    this%npatch = 0

    allocate(this%area(e%nface), this%normal(3,e%nface))
    allocate(this%f2p_map(e%nface))
    allocate(this%patch(e%nface))
    allocate(this%boundary_node(e%nnode))

    !! Compute face areas and normals
    do i = 1, e%nface
      face_nodes = e%fnode(e%xface(i):e%xface(i+1)-1)
      normal = face_normal(e%x(:,face_nodes))
      this%area(i) = vector_length(normal)
      this%normal(:,i) = normalized(normal)
    end do

    !! Get the faces of each vertex
    this%vface = faces_neighboring_vertices(e)

    !! Initialize inner data structures
    call this%heap%init(e%nnode * size(this%vface,dim=1))
    call this%enhbr%init(e%xface, e%fnode)
    call this%vgraph%init(e%nnode)

    !! Add edges to vertex graph
    do i = 1, e%nface
      face_nodes = e%fnode(e%xface(i):e%xface(i+1)-1)
      do j = 1, size(face_nodes)
        call get_edge_nodes(face_nodes, j, edge)
        call this%vgraph%add_edge(edge(1), edge(2))
      end do
    end do

    !! Get vertex adjacency matrix
    call this%vgraph%get_adjacency(this%xvnhbr, this%vnhbr)

    !! Get face adjacency matrix
    call get_fnhbr_aux(this%enhbr, this%e%xface, this%e%fnode, this%xfnhbr, this%fnhbr, stat, errmsg, this%normal, max_angle)
    if (stat/=0) return

    !! Mark vertices on mesh boundary
    max_angle_rad = PI*max_angle/180.0_r8
    if (.not. allocated(edge)) allocate(edge(2))
    this%boundary_node = .false.
    do i = 1, e%nnode
      edge(1) = i
      do j = this%xvnhbr(i), this%xvnhbr(i+1)-1
        edge(2) = this%vnhbr(j) !! neighbor of node i
        call this%enhbr%get_neighbors(edge, nhbrs)
        !! Boundary edges only neighbor one face
        if (size(nhbrs) <= 1) then
          this%boundary_node(edge) = .true.
          cycle
        end if
        !! Boundary edges neighbor faces with large angle
        angle = dot_product(this%normal(:,nhbrs(1)%j), this%normal(:,nhbrs(2)%j))
        if (angle > 1) angle = 1  ! Fix floating point errors
        angle = acos(angle)
        if (angle >= max_angle_rad) this%boundary_node(edge) = .true.
      end do
    end do

    if (this%verbosity > 1) then
      print '("INITIALIZING ", a, ":")', this%name()
      print '("  MAX ANGLE:", f6.2)', max_angle
      print '("  MERGE LEVEL:", i3)', this%merge_level
      print '("  SPLIT PATCH SIZE:", i3)', this%split_patch_size
    end if

  end subroutine init_vac


  !! Allocate and initialize PAVE_PATCHING data
  subroutine init_pave(this, e, params, stat, errmsg)

    use patching_tools, only: init_random_seed

    class(pave_patching), intent(out) :: this
    type(encl), intent(in), target :: e
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: context
    integer :: seed

    call init_vac(this, e, params, stat, errmsg)
    if (stat/=0) return

    !! Process the parameters.
    context = 'processing ' // params%name() // ': '
    if (params%is_parameter('random-seed')) then
      call params%get('random-seed', seed, stat=stat, errmsg=errmsg)
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      end if
    else
      call system_clock(count=seed)
    end if
    call init_random_seed(seed)

    if (this%verbosity > 1) then
      print '("  RANDOM SEED:", i10)', seed
    end if

  end subroutine init_pave


  !! The main VAC patching algorithm
  subroutine run_vac(this, stat, errmsg)

    class(vac_patching), target, intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer, allocatable :: faces(:)
    integer :: i

    stat = 0

    !! Reset face assignments
    this%f2p_map = -1

    !! Fill the heap
    do i = 1, this%e%nnode
      faces = pack(this%vface(:,i), this%vface(:,i) /= -1)
      call this%queue_connected_faces(i, faces)
    end do

    !! Assign patches
    call this%set_patches()
    ASSERT(all(this%f2p_map /= -1))

    !! Statistics after VAC
    if (this%verbosity > 1) then
      print '("AFTER ", a)', this%name()
      call this%print_stats()
    end if

    call this%merge_patches()

    !! Final statistics
    if (this%verbosity > 0) then
      if (this%verbosity > 1) print '("FINAL STATS")'
      call this%print_stats()
    end if

  end subroutine run_vac


  !! The main PAVE patching algorithm
  subroutine run_pave(this, stat, errmsg)

    class(pave_patching), target, intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(heap_entry) :: cur  ! Current patch being tested
    integer, allocatable :: faces(:)

    stat = 0

    !! Reset face assignments
    this%f2p_map = -1

    !! Pick seed patches
    call this%pick_seeds()

    !! Pave patches
    do while (.not. this%heap%empty())
      cur = this%heap%pop()

      if (all(this%f2p_map(cur%face) == -1)) then
        !! No conflicts, assign all faces
        call this%add_patch(cur%node, cur%face, cur%weight)

        !! Queue neighbors
        call this%queue_vertex_neighbors(cur%node)
      else
        !! Add connected subsets of faces to queue
        faces = pack(cur%face, this%f2p_map(cur%face) == -1)
        if (size(faces) < 1) cycle
        call this%queue_connected_faces(cur%node, faces)
      end if
    end do

    !! Statistics after paving
    if (this%verbosity > 1) then
      print '("AFTER ", a)', this%name()
      call this%print_stats()
    end if

    call this%merge_patches()

    !! Final statistics
    if (this%verbosity > 0) then
      if (this%verbosity > 1) print '("FINAL STATS")'
      call this%print_stats()
    end if

  end subroutine run_pave


  !! Write the patch data
  subroutine output(this, f2p_map, global_ids, npatch)

    class(vac_patching), intent(inout) :: this
    integer, allocatable, intent(out) :: f2p_map(:), global_ids(:)
    integer, intent(out) :: npatch

    integer :: i

    npatch = this%npatch
    call move_alloc(this%f2p_map, f2p_map)

    allocate(global_ids(npatch))

    do i = 1, npatch
      global_ids(i) = i
      ASSERT(any(f2p_map == i))
    end do

    ASSERT(maxval(f2p_map) == npatch)

  end subroutine output


  !! Report statistics for this run
  subroutine print_stats(this)

    class(vac_patching), intent(inout) :: this

    integer, allocatable :: nfp(:)  ! Number of faces in each patch
    integer, allocatable :: nps(:)  ! Number of patches of each size
    integer :: max_nfp, max_nfp_idx ! Max faces per patch
    real(r8) :: avg_nfp             ! Average faces per patch
    real(r8) :: std_nfp             ! Std. dev. of faces per patch
    integer :: i

    !! Initialize stats
    allocate(nfp(this%npatch))
    nfp = 0
    max_nfp = -1
    std_nfp = 0

    !! Collect face data
    avg_nfp = this%e%nface / REAL(this%npatch, kind=r8)
    do i = 1, this%e%nface
      nfp(this%f2p_map(i)) = nfp(this%f2p_map(i)) + 1
    end do

    !! Collect stats
    max_nfp_idx = maxloc(nfp, dim=1)
    max_nfp = nfp(max_nfp_idx)

    allocate(nps(max_nfp))
    nps = 0
    do i = 1, this%npatch
      nps(nfp(i)) = nps(nfp(i)) + 1
      std_nfp = std_nfp + (nfp(i) - avg_nfp)**2
    end do

    std_nfp = sqrt(std_nfp/this%npatch)

    print '("------------------------------------------------------------")'
    print '(a, " STATS:")', this%name()
    print '("  NFACE:  ", i8)', this%e%nface
    print '("  NPATCH: ", i8)', this%npatch
    print '("  AVG. FACES PER PATCH: ", es11.4)', avg_nfp
    print '("  S.D. FACES PER PATCH: ", es11.4)', std_nfp
    print '("  MAX FACES PER PATCH: ", i3)', max_nfp

    print '("  PATCHES BY SIZE:")'
    do i = 1, size(nps)
      print '("  ", i3, " : ", i8, "  (",  f6.2, "%)")', i, nps(i), nps(i)/REAL(this%npatch)*100
    end do
    print '("------------------------------------------------------------")'

  end subroutine print_stats


  !! Queues the initial seed patches for the PAVE algorithm.  For each connected
  !! component of the enclosure, PICK_SEEDS finds the nodes neighboring the
  !! highest number of boundary nodes, and picks one of those nodes at random.
  !! The seed patch for that component is then the full patch formed with that
  !! node as a vertex anchor.  Thus, for components with edges or corners, the
  !! seed patches are on the edge or in a corner, respectively.
  subroutine pick_seeds(this)

    class(pave_patching), intent(inout) :: this

    integer :: i, j, k, v, tmp, ncomp
    integer :: nbndry, max_nbndry      ! Number of boundary nodes neighboring vertex
    integer, allocatable :: cand(:,:)  ! Stores the node IDs of seed candidates
    integer, allocatable :: ncand(:)   ! ncand(j) is the number of nodes neighboring j boundary nodes
    integer, allocatable :: xcomp(:), comp(:)
    integer, allocatable :: faces(:)
    real(r8) :: r

    !! Get vertex graph connected components, excluding boundary nodes
    call this%vgraph%get_components(.not. this%boundary_node, ncomp, xcomp, comp)

    allocate(ncand(0:size(this%vface, dim=1)))  ! A node can only have "max degree" neighbors
    tmp = maxval(xcomp(2:ncomp+1) - xcomp(1:ncomp))  ! Maximum nodes in any component
    allocate(cand(tmp, 0:size(this%vface, dim=1)))

    if (this%verbosity > 1) then
      print '("PICK SEEDS")'
      print '("  Found ", i3, " component(s) in vertex graph.")', ncomp
      print '("  Max nodes per component: ", i5)', tmp
    end if

    !! Choose seeds for each component
    do i = 1, ncomp

      ncand = 0
      cand = -1

      !! Split component nodes by the number of adjacent boundary nodes
      do j = xcomp(i), xcomp(i+1)-1
        v = comp(j)
        nbndry = count(this%boundary_node(this%vnhbr(this%xvnhbr(v):this%xvnhbr(v+1)-1)))
        ncand(nbndry) = ncand(nbndry) + 1
        cand(ncand(nbndry), nbndry) = v
      end do

      !! Choose the set of nodes neighboring the most boundary nodes
      j = size(ncand)-1
      do while (ncand(j) == 0 .and. j > 0)
        j = j - 1
      end do
      max_nbndry = j

      if (this%verbosity > 1) then
        print '(/,"  Component", i3)', i
        print '("    Total nodes = ", i5)', xcomp(i+1)-xcomp(i)
        print '("    Total candidates = ", i4)', ncand(max_nbndry)
        print '("      # of adjacent boundary nodes = ", i2)', max_nbndry
      end if

      !! Pick random candidate
      call random_number(r)
      k = floor(r*(ncand(max_nbndry))) + 1
      v = cand(k, max_nbndry)
      faces = pack(this%vface(:,v), this%vface(:,v) /= -1)
      call this%queue_connected_faces(v, faces)

      if (this%verbosity > 1) then
        print '("      Chose random seed anchor:", i5)', v
      end if
    end do

    deallocate(cand, ncand)

  end subroutine pick_seeds


  !! Adds patches by clearing the priority queue.  QUEUE_SUBSET determines
  !! whether to queue subsets of the faces of a popped heap entry when one
  !! of the faces is already assigned.
  subroutine set_patches(this, queue_subset)

    class(vac_patching), intent(inout) :: this
    logical, intent(in), optional :: queue_subset

    type(heap_entry) :: cur  !! Current patch being tested
    integer, allocatable :: faces(:)
    logical :: qs

    qs = .true.
    if (present(queue_subset)) qs = queue_subset

    do while (.not. this%heap%empty())
      cur = this%heap%pop()

      !! No conflicts, assign all faces
      if (all(this%f2p_map(cur%face) == -1)) then
        call this%add_patch(cur%node, cur%face, cur%weight)
        cycle
      end if

      !! Add connected subsets of faces to queue
      if (qs) then
        faces = pack(cur%face, this%f2p_map(cur%face) == -1)
        if (size(faces) < 1) cycle
        call this%queue_connected_faces(cur%node, faces)
      end if
    end do

  end subroutine set_patches


  !! Deletes small patches and adds them to the queue.  The constituent faces
  !! of each patch are also queued as 1-face patches.
  subroutine split_patches(this)

    use integer_set_type

    class(vac_patching), intent(inout) :: this

    type(integer_set) :: del_patch_set  ! Set of patches to be deleted
    integer, allocatable :: del_patch(:)
    integer :: i, j, n, face(1)
    real(r8) :: weight

    !! Find patches to delete
    do i = 1, this%npatch
      associate (patch => this%patch(i))
        n = size(patch%face)
        if (n > 1 .and. n <= this%split_patch_size) then
          call del_patch_set%add(i)
          !! Queue old patches
          call this%heap%put(patch%node, patch%face, patch%weight)
          !! Queue 1-face patches
          do j = 1, n
            face = patch%face(j)
            weight = this%compute_weight(patch%node, face)
            call this%heap%put(patch%node, face, weight)
          end do
        end if
      end associate
    end do

    !! Delete patches
    del_patch = del_patch_set
    call this%delete_patches(del_patch)

    if (this%verbosity > 1) then
      print '("SPLIT PATCHES")'
      print '("  Found ",i5," patches of size <= ",i2," to split")', size(del_patch), this%split_patch_size
    end if

    deallocate(del_patch)

  end subroutine split_patches


  !! Merges patches using the various merging algorithms.
  subroutine merge_patches(this)

    class(vac_patching), intent(inout) :: this

    if (this%merge_level <= 0) return

    if (this%merge_level > 0) then
      !! Merge patches by vertex anchor
      call this%merge_patches_vertex_anchor()
      ASSERT(all(this%f2p_map /= -1))

      if (this%verbosity > 1) then
        print '("AFTER MERGE")'
        call this%print_stats()
      end if
    end if

    if (this%merge_level > 1) then
      !! Merge patches by vertex neighbors
      call this%merge_patches_vertex_neighbors(.false.)
      ASSERT(all(this%f2p_map /= -1))

      if (this%verbosity > 1) then
        print '("AFTER MERGE_VERTEX_NEIGHBORS WITHOUT OLD PATCH BIAS")'
        call this%print_stats()
      end if
    end if

    if (this%merge_level > 2) then
      !! Merge patches by vertex neighbors
      call this%merge_patches_vertex_neighbors(.true.)
      ASSERT(all(this%f2p_map /= -1))

      if (this%verbosity > 1) then
        print '("AFTER MERGE_VERTEX_NEIGHBORS WITH OLD PATCH BIAS")'
        call this%print_stats()
      end if
    end if

  end subroutine merge_patches


  !! Merge patches fully enclosed in the faces of a node
  subroutine merge_patches_vertex_anchor(this)

    use integer_set_type

    class(vac_patching), intent(inout) :: this

    integer, allocatable :: adj_faces(:), adj_patch(:), del_patch(:)
    integer :: i, j, p, cnt, nface, nadj_patch, npatch_old
    type(integer_set) :: del_patch_set  ! Set of patches to be deleted

    allocate(adj_patch(size(this%vface, dim=1)))
    npatch_old = this%npatch
    cnt = 0

    !! Remove small patches to encourage more merge candidates
    call this%split_patches()

    !! Find merge candidates
    do i = 1, this%e%nnode
      adj_faces = pack(this%vface(:,i), this%vface(:,i) /= -1)
      nadj_patch = 0
      adj_patch = -1
      nface = 0

      !! Find adjacent patches
      do j = 1, size(adj_faces)
        p = this%f2p_map(adj_faces(j))
        !! Ignore unassigned faces
        if (p == -1) then
          nface = nface + 1
          cycle
        end if

        if (.not. any(adj_patch(1:nadj_patch) == p)) then
          nadj_patch = nadj_patch + 1
          adj_patch(nadj_patch) = p
        end if
      end do

      !! Ignore optimal patches
      if (nadj_patch == 1 .and. nface == 0) cycle

      do j = 1, nadj_patch
        nface = nface + size(this%patch(adj_patch(j))%face)
      end do

      !! Ignore cases when patches exceed node range
      if (nface /= size(adj_faces)) cycle

      !! Record patches to delete
      if (nadj_patch > 0) call del_patch_set%add(adj_patch(1:nadj_patch))

      !! Queue merge candidates
      call this%queue_connected_faces(i, adj_faces)

      cnt = cnt + 1
    end do

    if (cnt > 0) then
      !! Delete patches
      del_patch = del_patch_set
      call this%delete_patches(del_patch)
    end if

    !! Add new patches or restore old ones.
    call this%set_patches()

    if (this%verbosity > 1) then
      print '("MERGE PATCHES BY VERTEX ANCHOR")'
      print '("  Found ", i5, " merge candidates")', cnt
      print '("  Reduced patch count by ", i5)', npatch_old-this%npatch
    end if

  end subroutine merge_patches_vertex_anchor


  !! Merge patches fully enclosed in the faces of pairs of adjacent nodes
  subroutine merge_patches_vertex_neighbors(this, old_patch_bias)

    use integer_set_type

    class(vac_patching), intent(inout) :: this
    logical, intent(in) :: old_patch_bias

    type(integer_set) :: del_patch_set  ! Set of patches to be deleted
    integer, allocatable :: adj_face(:), adj_patch(:), del_patch(:), faces(:)
    integer :: i, j, k, v, n, f, p, node, nodes(2), cnt, nface, nadj_patch, nadj_face, npatch_old
    real(r8) :: bias

    allocate(adj_patch(2*size(this%vface, dim=1)))
    allocate(adj_face(2*size(this%vface, dim=1)))
    npatch_old = this%npatch
    bias = merge(100_r8, 0_r8, old_patch_bias)
    cnt = 0

    !! Remove small patches to encourage more merge candidates
    call this%split_patches()

    !! Find merge candidates
    do v = 1, this%e%nnode

      nodes(1) = v

      do i = this%xvnhbr(v), this%xvnhbr(v+1)-1

        n = this%vnhbr(i)  ! neighbor of v
        if (v > n) cycle   ! Avoid duplicates
        nodes(2) = n

        nadj_patch = 0
        adj_patch = -1
        nadj_face = 0
        adj_face = -1
        nface = 0

        !! Find adjacent faces
        faces = pack(this%vface(:,nodes), this%vface(:,nodes) /= -1)
        do j = 1, size(faces)
          f = faces(j)
          if (.not. any(adj_face(1:nadj_face) == f)) then
            nadj_face = nadj_face + 1
            adj_face(nadj_face) = f
          end if
        end do

        !! Find adjacent patches
        do j = 1, nadj_face
          p = this%f2p_map(adj_face(j))
          !! Ignore unassigned faces
          if (p == -1) then
            nface = nface + 1
            cycle
          end if

          if (.not. any(adj_patch(1:nadj_patch) == p)) then
            nadj_patch = nadj_patch + 1
            adj_patch(nadj_patch) = p
          end if
        end do

        !! Ignore optimal patches
        if (nadj_patch == 1 .and. nface == 0) cycle

        do k = 1, nadj_patch
          nface = nface + size(this%patch(adj_patch(k))%face)
        end do

        !! Ignore cases when patches exceed node range
        if (nface /= nadj_face) cycle

        !! Record patches to delete
        if (nadj_patch > 0) call del_patch_set%add(adj_patch(1:nadj_patch))

        !! Queue merge candidate
        node = merge(v, n, this%boundary_node(n))  ! prefer non-boundary nodes as anchors
        call this%queue_connected_faces(node, adj_face(1:nadj_face))

        !! Queue old patches
        do j = 1, nadj_patch
          k = adj_patch(j)
          call this%heap%put(this%patch(k)%node, this%patch(k)%face, this%patch(k)%weight + bias)
        end do

        cnt = cnt + 1
      end do
    end do

    if (cnt > 0) then
      !! Delete patches
      del_patch = del_patch_set
      call this%delete_patches(del_patch)
    end if

    !! Add new patches or restore old ones.  To avoid incomplete patches, do
    !! not re-queue subsets of merge candidates.  The old patches will fill
    !! the spaces instead.
    call this%set_patches(.false.)

    if (this%verbosity > 1) then
      print '("MERGE PATCHES BY VERTEX NEIGHBORS")'
      print '("  Found ", i5, " merge candidates")', cnt
      print '("  Reduced patch count by ", i5)', npatch_old-this%npatch
    end if

  end subroutine merge_patches_vertex_neighbors


  !! Inserts a new patch
  subroutine add_patch(this, node, faces, weight)

    class(vac_patching), intent(inout) :: this
    integer, intent(in) :: node, faces(:)
    real(r8), intent(in) :: weight

    this%npatch = this%npatch + 1
    this%f2p_map(faces) = this%npatch

    associate (patch => this%patch(this%npatch))
      patch%node = node
      patch%face = faces
      patch%weight = weight
    end associate

  end subroutine add_patch


  !! Deletes a patch, resetting face assignments
  subroutine delete_patch(this, p)

    class(vac_patching), intent(inout) :: this
    integer, intent(in) :: p

    ASSERT(p <= this%npatch)

    !! Reset face assignments
    this%f2p_map(this%patch(this%npatch)%face) = p
    this%f2p_map(this%patch(p)%face) = -1

    this%patch(p) = this%patch(this%npatch)
    this%npatch = this%npatch - 1

  end subroutine delete_patch


  !! Deletes a list of patches
  subroutine delete_patches(this, patches)

    use sort_utilities

    class(vac_patching), intent(inout) :: this
    integer, intent(in) :: patches(:)

    integer, allocatable :: tmp(:)
    integer :: i, p, n

    ASSERT(size(patches) <= this%npatch)

    n = size(patches)
    allocate(tmp(n))
    call heap_sort(patches, tmp)

    do i = n, 1, -1
      !! Deleting in descending order avoids the issue of deleting a patch ID greater than this%npatch
      p = patches(tmp(i))
      call this%delete_patch(p)
    end do

    deallocate(tmp)

  end subroutine delete_patches


  !! Adds connected sets of faces to the global priority queue
  subroutine queue_connected_faces(this, node, faces)

    use patching_tools, only: get_connected_faces_subset

    class(vac_patching), intent(inout) :: this
    integer, intent(in) :: node, faces(:)

    integer, allocatable :: qfaces(:) ! Faces to be queued
    integer, allocatable :: tag(:)
    integer :: i, ncomp
    real(r8) :: weight

    !! Find connected components
    call get_connected_faces_subset(this%xfnhbr, this%fnhbr, faces, tag, ncomp)

    !! Queue connected components
    do i = 1, ncomp
      qfaces = pack(faces, tag == i)
      weight = this%compute_weight(node, qfaces)
      call this%heap%put(node, qfaces, weight)
    end do

  end subroutine queue_connected_faces


  !! Adds 1st and 2nd degree neighbors of a given node to the global priority queue
  subroutine queue_vertex_neighbors(this, node)

    class(pave_patching), intent(inout) :: this
    integer, intent(in) :: node

    integer, allocatable :: faces(:) ! Faces to be queued
    integer :: j, k, n1, n2, v

    !! Queue neighbors
    do j = this%xvnhbr(node), this%xvnhbr(node+1)-1
      n1 = this%vnhbr(j)  ! neighbor of v
      faces = pack(this%vface(:,n1), this%vface(:,n1) /= -1)
      call this%queue_connected_faces(n1, faces)

      !! Queue neighbors of neighbors
      do k = this%xvnhbr(n1), this%xvnhbr(n1+1)-1
        n2 = this%vnhbr(k)  ! neighbor of n1
        if (n2 == node) cycle
        !! Count n2 as a boundary node if we crossed a boundary.
        v = merge(n1, n2, this%boundary_node(n1))
        faces = pack(this%vface(:,n2), this%vface(:,n2) /= -1)
        call this%queue_connected_faces(v, faces)
      end do
    end do

  end subroutine queue_vertex_neighbors


  !! Computes the weight of a list of faces
  function compute_weight (this, node, faces) result(ret)

    use cell_geometry, only: vector_length, normalized
    use cell_topology, only: get_edge_nodes

    class(vac_patching), intent(in) :: this
    integer, intent(in) :: node, faces(:)
    real(r8) :: ret

    type(edge_neighbor), allocatable :: nhbrs(:)
    real(r8) :: area, perimeter, length, normal(3)
    real(r8) :: e_norm, e_shape, e_pos, e_nface
    integer, allocatable :: edge(:), adj_faces(:)
    integer :: i, j, k, f

    ASSERT(size(faces) > 0)

    area   = 0
    normal = 0
    perimeter = 0

    !! Compute face geometry
    do i = 1, size(faces)
      f = faces(i)
      area = area + this%area(f)
      normal = normal + this%normal(:,f)*this%area(f)

      !! Add length of exterior edges to perimeter
      associate (face_nodes => this%e%fnode(this%e%xface(f):this%e%xface(f+1)-1))
        do j = 1, size(face_nodes)
          call get_edge_nodes(face_nodes, j, edge)
          call this%enhbr%get_neighbors(edge, nhbrs)

          if (size(nhbrs) == 1) then
            !! Edge is on the border of the enclosure
            length = vector_length(this%e%x(:,edge(1)) - this%e%x(:,edge(2)))
            perimeter = perimeter + length
          else
            do k = 1, size(nhbrs)
              !! Ignore interior edges
              if (any(faces == nhbrs(k)%j)) cycle
              length = vector_length(this%e%x(:,edge(1)) - this%e%x(:,edge(2)))
              perimeter = perimeter + length
            end do
          end if

        end do
      end associate
    end do

    normal = normalized(normal)

    !! Compute normal error
    e_norm = 0
    do i = 1, size(faces)
      f = faces(i)
      e_norm = e_norm + sum((this%normal(:,f) - normal)**2)
    end do
    e_norm = e_norm / REAL(size(faces), kind=r8)

    !! Compute shape error
    e_shape = perimeter**2 / (4*PI*area)

    !! Bias against nodes on enclosure boundary
    e_pos = merge(100_r8, 0_r8, this%boundary_node(node))

    !! Compute patch size bias
    e_nface = merge(4_r8, 0_r8, size(faces) <= 1)
    adj_faces = pack(this%vface(:,node), this%vface(:,node) /= -1)
    do i = 1, size(adj_faces)
      f = adj_faces(i)
      if (.not. any(faces == f)) then
        e_nface = e_nface + 1_r8
        exit
      end if
    end do

    !! Total weight
    ret = e_norm + e_shape + e_pos + e_nface

  end function compute_weight


end module vac_patching_type
