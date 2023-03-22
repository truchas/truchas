!!
!! VSA_PATCHING_TYPE
!!
!! A concrete implementation of the abstract base class PATCHING that
!! encapsulates the Variational Shape Approximation (VSA) algorithm for
!! clustering faces of a radiation enclosure mesh.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! 2 May 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! VSA_PATCHING is an extension of the abstract base class PATCHING.  See the
!! base class comments for a description of the common type bound procedures.
!! Additional details specific to the VSA_PATCHING type are included below.
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
!!       face-patch-ratio -- ratio of total faces to total patches.  It defines
!!           the number of patch seeds with which to initialize the algorithm,
!!           and therefore the total number of patches in the output:
!!             [# SEEDS] = [# FACES] / face-patch-ratio
!!       max-iter -- maximum number of iterations allowed.  The procedure stops
!!           when max-iter is reached, regardless of other stopping conditions.
!!       min-delta -- threshold for the minimum change in patch proxies between
!!           successive iterations.  The procedure stops if the minimum change
!!           in patch proxies is less than min-delta.
!!       max-patch-radius -- maximum desired radius for a patch.  A face outside
!!           this radius will have a penalty added to its weight that is
!!           proportional to its distance from the patch center.
!!       normalize-dist -- determines whether the Voronoi distance bias should
!!           be normalized by the face radius.
!!       random-seed -- optional parameter. If present, it is used to initialize
!!           the random number generator used to pick the initial seeds.
!!           If not present, the seed is taken from the system clock.
!!


#include "f90_assert.fpp"

module vsa_patching_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use re_encl_type
  use vsa_min_heap
  use vsa_patch_type
  use patching_class
  use parameter_list_type
  implicit none

  !! Parameter defaults
  integer, parameter :: VSA_MAX_ITER_DEFAULT = 1000
  real(r8), parameter :: VSA_MIN_DELTA_DEFAULT = 1E-6_r8
  real(r8), parameter :: VSA_FACE_PATCH_RATIO_DEFAULT = 4.0_r8
  real(r8), parameter :: VSA_MAX_PATCH_RADIUS_DEFAULT = sqrt(huge(0.0_r8))
  logical, parameter :: VSA_NORMALIZE_DIST_DEFAULT = .true.

  integer, parameter :: TELEPORT_PATCH_EVERY_ITER = 5

  type, extends(patching), public :: vsa_patching
    private
    type(encl), pointer :: e
    type(min_heap) :: heap
    real(r8), allocatable :: area(:), normal(:,:), center(:,:), radius(:)
    type(vsa_patch), allocatable :: patch(:)
    integer, allocatable  :: xfnhbr(:), fnhbr(:)
    integer, allocatable  :: f2p_map(:), seeds(:)
    integer :: npatch, npatch_min
    integer :: max_iter     ! Maximum iterations
    real(r8) :: iter_delta  ! Change in proxies between successive iterations
    real(r8) :: min_delta   ! Stopping threshold
    logical :: dir          ! Determines direction of patch traversal
    real(r8) :: max_radius  ! Maximum desired patch radius
    logical :: normalize    ! Whether to normalize the Voronoi distance bias
    real(r8) :: max_face_err, max_patch_err
    integer :: max_face_idx, max_patch_idx
    integer :: verbosity

    contains
      procedure, public  :: init
      procedure, public  :: run
      procedure, public  :: output
      procedure, private :: print_stats
      procedure, private :: pick_start_seeds_spread
      procedure, private :: pick_start_seeds_iter
      procedure, private :: pick_seeds
      procedure, private :: partition
      procedure, private :: proxy_fit
      procedure, private :: patch_insert
      procedure, private :: patch_delete
      procedure, private :: patch_init
      procedure, private :: compute_weight
      procedure, private :: recompute_weights
  end type

contains


  !! Allocate and initialize VSA_PATCHING data
  subroutine init(this, e, params, stat, errmsg)

    use cell_geometry, only: face_normal, vector_length, normalized, polygon_center
    use patching_tools, only: init_random_seed, get_face_neighbor_array

    class(vsa_patching), intent(out) :: this
    type(encl), target, intent(in) :: e
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i, j, n, max_edges, seed
    real(r8) :: fp_ratio   ! Total faces / total patches
    real(r8) :: max_angle  ! Maximum allowable angle for adjacent faces (in degrees)
    real(r8) :: normal(3), rface
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
    call params%get('face-patch-ratio', fp_ratio, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    call params%get('max-iter', this%max_iter, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    call params%get('min-delta', this%min_delta, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    call params%get('max-patch-radius', this%max_radius, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    call params%get('normalize-dist', this%normalize, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if

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

    this%e => e
    this%dir = .true.
    this%npatch = 0
    this%npatch_min = e%nface / fp_ratio

    max_edges = 0
    do i = 1, e%nface
      n = e%xface(i+1)-e%xface(i)
      max_edges = max(n, max_edges)
    end do

    !! Allocate space
    call this%heap%init(max_edges*e%nface)
    allocate(this%area(e%nface), this%normal(3,e%nface), this%center(3,e%nface))
    allocate(this%radius(e%nface), this%f2p_map(e%nface))
    allocate(this%patch(e%nface), this%seeds(e%nface))

    !! Compute face areas, normals, and centers
    do i = 1, e%nface
      associate(face_nodes => e%fnode(e%xface(i):e%xface(i+1)-1))
        normal = face_normal(e%x(:,face_nodes))
        this%area(i) = vector_length(normal)
        this%normal(:,i) = normalized(normal)
        this%center(:,i) = polygon_center(e%x(:,face_nodes))
        rface = 0.0_r8
        do j = 1, size(face_nodes)
          rface = max(rface, norm2(e%x(:,face_nodes(j)) - this%center(:,i)))
        end do
        this%radius(i) = rface
      end associate
    end do

    !! Get face neighbors
    call get_face_neighbor_array(e%xface, e%fnode, this%xfnhbr, this%fnhbr, stat, errmsg, this%normal, max_angle)
    if (stat /= 0) then
      block
        integer :: n
        character(40) :: coord
        n = findloc(this%fnhbr, -1, dim=1) ! first side with invalid topology
        write(coord,'("(",es12.5,2(",",es12.5),")")') e%side_location(n)
        errmsg = errmsg // ' at ' // coord
      end block
      return
    end if

    if (this%verbosity > 1) then
      print '("INITIALIZING VSA:")'
      print '("  FACE-PATCH RATIO: ", f0.2)', fp_ratio
      print '("  MAX PATCH RADIUS: ", es11.3)', this%max_radius
      print '("  NORMALIZE DISTANCE: ", l1)', this%normalize
      print '("  MIN NPATCH: ", i0)', this%npatch_min
      print '("  MAX ANGLE: ", f0.2)', max_angle
      print '("  RANDOM SEED: ", i0)', seed
    end if

  end subroutine init


  !! The main VSA patching algorithm
  subroutine run(this, stat, errmsg)

    class(vsa_patching), target, intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i = 1

    stat = 0

    !! Skip first convergence check
    this%iter_delta = this%min_delta + 1

    !! Run first iteration
    call this%pick_start_seeds_iter()
    call this%partition()
    call this%proxy_fit()
    if (this%verbosity > 1) then
      print '("------------------------------------------------------------")'
      print '("VSA ITER ", i0)', i
      print '("  ITER_DELTA: ", es11.4)', this%iter_delta
      print '("  MAX_PATCH_ERR: ", es11.4, " | PATCH: ", i0)', this%max_patch_err, this%max_patch_idx
      print '("  MAX_FACE_ERR:  ", es11.4, " | FACE: ", i0)', this%max_face_err, this%max_face_idx
      print '("------------------------------------------------------------")'
    end if

    do while ((this%iter_delta > this%min_delta) .and.  (i < this%max_iter))
      i = i + 1
      call this%pick_seeds()
      call this%partition()

      !! Teleport patches
      if ((mod(i, TELEPORT_PATCH_EVERY_ITER) == 0) .and. (i < this%max_iter - 1)) then
        call this%patch_delete()
        call this%patch_insert()
      end if

      call this%proxy_fit()

      if (this%verbosity > 1) then
        print '("------------------------------------------------------------")'
        print '("VSA ITER ", i0)', i
        print '("  ITER_DELTA: ", es11.4)', this%iter_delta
        print '("  NPATCH: ", i0)', this%npatch
        print '("  MAX_PATCH_ERR: ", es11.4, " | PATCH: ", i0)', this%max_patch_err, this%max_patch_idx
        print '("  MAX_FACE_ERR:  ", es11.4, " | FACE: ", i0)', this%max_face_err, this%max_face_idx
        print '("------------------------------------------------------------")'
      end if
    end do

    ASSERT(all(this%f2p_map /= -1))

    if (this%verbosity > 0) then
      print '("------------------------------------------------------------")'
      print '("VSA STOPED:")'
      print '("  ITERATIONS: ", i0, " out of maximum ", i0)', i, this%max_iter
      print '("  ITER_DELTA: ", es11.4, " with MIN_DELTA of ", es11.4)', this%iter_delta, this%min_delta
      print '("------------------------------------------------------------")'
    end if

    !! Print statistics for this run
    call this%print_stats()

  end subroutine run


  !! Write the patch data
  subroutine output(this, f2p_map, global_ids, npatch)

    class(vsa_patching), intent(inout) :: this
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

    class(vsa_patching), intent(inout) :: this

    integer, allocatable :: nfp(:)  ! Number of faces in each patch
    integer, allocatable :: nps(:)  ! Number of patches of each size
    integer :: max_nfp, max_nfp_idx ! Max faces per patch
    real(r8) :: avg_nfp             ! Average faces per patch
    real(r8) :: std_nfp             ! Std. dev. of faces per patch
    real(r8) :: avg_fe              ! Average face error
    real(r8) :: avg_pe              ! Average patch error
    integer :: i, j

    !! Initialize stats
    allocate(nfp(this%npatch))
    nfp = 0
    max_nfp = -1
    std_nfp = 0
    avg_fe  = 0
    avg_pe  = 0

    !! Collect face data
    avg_nfp = this%e%nface / REAL(this%npatch, kind=r8)
    nfp = this%patch%nface

    !! Collect patch data
    max_nfp_idx = maxloc(nfp, dim=1)
    max_nfp = nfp(max_nfp_idx)

    allocate(nps(max_nfp))
    nps = 0

    do i = 1, this%npatch
      associate (patch => this%patch(i))
        !! Collect average error data
        do j = 1, patch%nface
          avg_fe = avg_fe + patch%weight(j)
          avg_pe = avg_pe + patch%weight(j)
        end do

        !! Collect face count data
        nps(nfp(i)) = nps(nfp(i)) + 1
        std_nfp = std_nfp + (patch%nface - avg_nfp)**2
      end associate
    end do

    std_nfp = sqrt( std_nfp/this%npatch )
    avg_pe = avg_pe / this%npatch
    avg_fe = avg_fe / this%e%nface

    print '("------------------------------------------------------------")'
    print '("VSA STATS:")'
    print '("  NFACE:  ", i0)', this%e%nface
    print '("  NPATCH: ", i0)', this%npatch
    print '("  AVG FACES PER PATCH: ", es11.4)', avg_nfp
    print '("  S.D. FACES PER PATCH: ", es11.4)', std_nfp
    print '("  PATCHES BY SIZE:")'
    do i = 1, size(nps)
      print '("  ", i0, " : ", i0, "  (",  f6.2, "%)")', i, nps(i), nps(i)/REAL(this%npatch)*100
    end do
    if (this%verbosity > 1) then
      print '("  AVG PATCH WEIGHT: ", es11.4)', avg_pe
      print '("  MAX PATCH WEIGHT: ", es11.4)', this%max_patch_err
      print '("    IN PATCH: ", i0)', this%max_patch_idx
      print '("  AVG FACE WEIGHT: ", es11.4)', avg_fe
      print '("  MAX FACE WEIGHT: ", es11.4)', this%max_face_err
      print '("    IN FACE, PATCH: ", i0, ", ", i0)', this%max_face_idx, this%f2p_map(this%max_face_idx)
    end if
    print '("------------------------------------------------------------")'

  end subroutine print_stats


  !! Selects this%npatch_min faces as initial seeds by spreading
  !! them out evenly amongst faces.
  subroutine pick_start_seeds_spread(this)

    class(vsa_patching), intent(inout) :: this
    integer :: p, f
    real(r8) :: fpr  ! face-patch ratio

    fpr = REAL(this%e%nface,kind=r8)/REAL(this%npatch_min,kind=r8)
    fpr = merge(1.0_r8, fpr, fpr < 1.0_r8)

    !! Set the number of patches
    this%npatch = this%npatch_min

    !! Reset patch assignments
    this%f2p_map = -1

    do p = 1,this%npatch_min
      f = (p-1)*fpr + 1
      this%seeds(p) = f   ! face f is seed for patch p
      this%f2p_map(f) = p ! face f belongs to patch p
      call this%patch_init(f, p)
    end do

    if (this%verbosity > 2) then
      do p = 1, this%npatch
        print '("Face ", i0, " seeds patch ", i0, " of ", i0)', this%seeds(p), p, this%npatch
      end do
    end if

  end subroutine pick_start_seeds_spread


  !! Choose start seeds by iteratively inserting patches at largest weight face
  subroutine pick_start_seeds_iter(this)

    use patching_tools, only: get_connected_faces

    class(vsa_patching), intent(inout) :: this

    integer, allocatable :: xcomp(:), comp(:)
    integer :: ncomp
    integer :: i, k, p, f, fmax
    real(r8) :: r

    !! Reset patch assignments
    this%f2p_map = -1

    !! Get connected components in face neighbor graph
    call get_connected_faces(this%e%nface, this%xfnhbr, this%fnhbr, ncomp, xcomp, comp)
    ASSERT(ncomp == size(xcomp) - 1)
    if (this%verbosity > 1) print '("Found ", i0, " connected components")', ncomp

    !! Pick random seed in each component
    do i = 1, ncomp
      !! Pick random face in this component
      call random_number(r)
      k = floor(r*(xcomp(i+1)-xcomp(i)))
      ASSERT(k < xcomp(i+1)-xcomp(i))
      f = comp(xcomp(i)+k)

      !! Assign seed to new patch
      this%npatch = this%npatch + 1
      call this%patch_init(f, this%npatch)
      this%seeds(this%npatch) = f   ! face f seeds last patch
      this%f2p_map(f) = this%npatch ! face f belongs to last patch
    end do

    if (this%verbosity > 1) then
      do i = 1, ncomp
        print '("Face ", i0, " seeds patch ", i0, " of ", i0, " in component ", i0)', this%seeds(i), i, this%npatch_min, i
      end do
    end if

    !! Iteratively pick start seeds
    do while (this%npatch < this%npatch_min)
      !! Grow patches and find maximum weight face
      call this%partition(fmax)
      ASSERT(.not. any(this%seeds(1:this%npatch)==fmax))

      !! Reset patch assignments
      this%f2p_map = -1

      !! Reset patches to seed faces
      do p = 1, this%npatch
        associate (patch => this%patch(p))
          f = this%seeds(p)
          this%f2p_map(f) = p
          call patch%reset(f, this%area(f), 0.0_r8)
        end associate
      end do

      !! Max weight face is the new seed
      this%npatch = this%npatch + 1
      this%seeds(this%npatch) = fmax
      this%f2p_map(fmax) = this%npatch
      call this%patch_init(fmax, this%npatch)

      !! Picking starts seeds takes a while, so report progress on each iteration
      if (this%verbosity > 2) print '("Face ", i0, " seeds patch ", i0, " of ", i0)', fmax, this%npatch, this%npatch_min
    end do

  end subroutine pick_start_seeds_iter


  !! Selects the minimum weight face from each patch as the new seeds
  subroutine pick_seeds(this)

    class(vsa_patching), intent(inout) :: this
    integer :: i, p, f, fidx
    real(r8) :: weight, min_weight

    !! Reset patch assignments
    this%f2p_map = -1

    do p = 1, this%npatch
      associate (patch => this%patch(p))
        !! Proxies were changed in last iteration, recalculate weight
        min_weight = huge(0_r8)
        do i = 1, patch%nface
          f = patch%face(i)
          weight = this%compute_weight(f, p)
          if (min_weight > weight) then
            min_weight = weight
            fidx = i
          end if
        end do

        !! Set minimum weight face as new seed
        f = patch%face(fidx)
        this%seeds(p) = f   ! face f is seed for patch p
        this%f2p_map(f) = p ! face f belongs to patch p
        call patch%reset(f, this%area(f), min_weight)
      end associate
    end do

  end subroutine pick_seeds


  !! Partition enclosure into connected patches. Optionally returns the index
  !! of the maximum weight face.
  subroutine partition(this, max_face)

    class(vsa_patching), intent(inout) :: this
    integer, optional, intent(out) :: max_face

    type(heap_entry) :: cur   ! Current face being tested
    real(r8) :: weight, max_weight
    integer :: p, f, k, n, max_face_l
    integer :: first, last, step

    !! Switch patch traversal order to reduce min-heap bias.
    if (this%dir) then
      first = 1
      last = this%npatch
      step = 1
    else
      first = this%npatch
      last = 1
      step = -1
    end if
    this%dir = (.not. this%dir)

    !! Add seed neighbors to the heap
    do p = first, last, step
      associate (patch => this%patch(p))
        f = this%seeds(p)
        do k = this%xfnhbr(f), this%xfnhbr(f+1)-1
          n = this%fnhbr(k)  ! Seed neighbor
          if (n <= 0) cycle  ! Skip missing neighbors
          if (this%f2p_map(n) /= -1) cycle  ! Skip assigned neighbors
          weight = this%compute_weight(n, p)
          call this%heap%put(n, p, weight)
        end do
      end associate
    end do

    !! Grow patches
    max_weight = -huge(0.0_r8)
    max_face_l = -1
    do while (.not. this%heap%empty())
      cur = this%heap%pop()
      p = cur%patch
      f = cur%faceid

      !! Skip assigned faces
      if (this%f2p_map(f) /= -1) cycle

      associate (patch => this%patch(p))
        !! Assign face to corresponding patch
        call patch%add_face(f, this%area(f), cur%weight)
        this%f2p_map(f) = p
        if (max_weight < cur%weight) then
          max_weight = cur%weight
          max_face_l = f
        end if
        !! Add neighbors to the heap
        do k = this%xfnhbr(f), this%xfnhbr(f+1)-1
          n = this%fnhbr(k)  ! Face neighbor
          if (n <= 0) cycle  ! Skip missing neighbors
          if (this%f2p_map(n) /= -1) cycle  ! Skip assigned neighbors
          weight = this%compute_weight(n, p)
          call this%heap%put(n, p, weight)
        end do
      end associate
    end do

    if (present(max_face)) max_face = max_face_l

  end subroutine partition


  !! Calculate proxy normals for each patch
  subroutine proxy_fit(this)

    use cell_geometry, only: vector_length

    class(vsa_patching), intent(inout) :: this

    real(r8) :: center(3), normal(3), delta, normal_mag
    integer :: p, f, k

    !! Reset max error
    this%iter_delta = -1
    this%max_face_err = -huge(0_r8)
    this%max_patch_err = -huge(0_r8)

    do p = 1, this%npatch
      associate (patch => this%patch(p))
        center = 0
        normal = 0
        do k = 1, patch%nface
          f = patch%face(k)
          center = center + this%center(:,f)*this%area(f)
          normal = normal + this%normal(:,f)*this%area(f)
        end do
        center = center / patch%area
        normal_mag = norm2(normal)
        normal = normal / normal_mag

        delta = vector_length(patch%normal-normal) + vector_length(patch%center-center)
        this%iter_delta = max(this%iter_delta, delta)

        !! Collect statistics
        if (this%verbosity > 1) then
          if (this%max_patch_err < patch%total_weight) then
            this%max_patch_err = patch%total_weight
            this%max_patch_idx = p
          end if

          k = maxloc(patch%weight(1:patch%nface), dim=1)
          if (this%max_face_err < patch%weight(k)) then
            this%max_face_err = patch%weight(k)
            this%max_face_idx = patch%face(k)
          end if
        end if

        !! Assign proxy
        patch%normal = normal
        patch%normal_mag = normal_mag
        patch%center = center
      end associate
    end do

  end subroutine proxy_fit


  !! Adds a patch with its seed at the largest weight face
  subroutine patch_insert(this)

    class(vsa_patching), intent(inout) :: this

    real(r8) :: max_weight
    integer :: i, p, fmax, pmax, imax

    ASSERT(this%npatch < size(this%patch))

    max_weight = -huge(0_r8)

    !! Find maximum weight face
    do p = 1, this%npatch
      associate (patch => this%patch(p))
        !! Skip patches with one face
        if (patch%nface <= 1) cycle
        do i = 1, patch%nface
          if (max_weight < patch%weight(i)) then
            max_weight = patch%weight(i)
            imax = i
            pmax = p
          end if
        end do
      end associate
    end do

    fmax = this%patch(pmax)%face(imax)
    call this%patch(pmax)%remove_face(imax, this%area(fmax))

    !! Insert patch at maximum weight face
    this%npatch = this%npatch + 1
    call this%patch_init(fmax, this%npatch)
    this%seeds(this%npatch) = fmax
    this%f2p_map(fmax) = this%npatch

    if (this%verbosity > 1) then
      print '("PATCH_INSERT:")'
      print '("  FOUND MAX_WEIGHT FACE ", i0.3, " IN PATCH ", i0.3)', fmax, pmax
      print '("  INSERTING PATCH ", i0.3, " AT FACE ", i0.3)', this%npatch, fmax
    end if

  end subroutine patch_insert


  !! Deletes a patch by merging the pair of adjacent patches with least distortion
  subroutine patch_delete(this)

    use cell_geometry, only: normalized

    class(vsa_patching), intent(inout) :: this

    integer :: p1, p2  ! The two patches being tested
    integer :: merge_ids(2) ! The two patches to be merged
    real(r8) :: center1(3), normal1(3)
    real(r8) :: center2(3), normal2(3)
    real(r8) :: centerT(3), normalT(3), normal_magT  ! Proxies of merge candidate
    real(r8) :: center_min(3), normal_min(3), normal_mag_min
    real(r8) :: delta, min_delta
    integer :: i, j, k, f

    !! Initialize error
    min_delta = huge(0_r8)

    ! Todo: if we maintained a patch connectivity graph, this would be a lot more efficient
    !! Find pair of patches with least distortion
    do p1 = 1, this%npatch
      associate (patch1 => this%patch(p1))
        center1 = patch1%center * patch1%area
        normal1 = patch1%normal * patch1%normal_mag

        !! Traverse face neighbors to find adjacent patches
        do i = 1, patch1%nface
          k = patch1%face(i)
          do j = this%xfnhbr(k), this%xfnhbr(k+1)-1
            f = this%fnhbr(j)
            !! Skip missing neighbors (i.e. k is on the mesh boundary)
            if (f <= 0) cycle
            p2 = this%f2p_map(f)
            !! Skip faces in this patch and avoid duplicate pairs
            if (p1 >= p2) cycle
            !! Compare the two patches
            associate (patch2 => this%patch(p2))
              center2 = patch2%center * patch2%area
              normal2 = patch2%normal * patch2%normal_mag

              centerT = (center1 + center2) / (patch1%area + patch2%area)
              normalT = normal1+normal2
              normal_magT = norm2(normalT)
              normalT = normalized(normalT)

              !! Change in error for this merge
              delta = sum((centerT - patch1%center)**2) + sum((centerT - patch2%center)**2)
              delta = delta + sum((normalT - patch1%normal)**2) + sum((normalT - patch2%normal)**2)

              if (min_delta > delta) then
                min_delta = delta
                merge_ids(1) = p1
                merge_ids(2) = p2
                center_min = centerT
                normal_min = normalT
                normal_mag_min = normal_magT
              end if

            end associate
          end do
        end do
      end associate
    end do

    if (this%verbosity > 1) then
      print '("PATCH_DELETE:")'
      print '("  FOUND MIN_WEIGHT PATCHES ", i0.3, " AND ", i0.3)', merge_ids(1), merge_ids(2)
      print '("    Replacing patch ", i0.3, " with patch ", i0.3)', merge_ids(2), this%npatch
    end if

    !! Merge minimum weight patches
    associate (patch1 => this%patch(merge_ids(1)), patch2 => this%patch(merge_ids(2)))
      !! Add faces to patch1
      do i = 1, patch2%nface
        f = patch2%face(i)
        call patch1%add_face(f, this%area(f), 0.0_r8)
        this%f2p_map(f) = merge_ids(1)
      end do

      !! Recompute face weights of new patch
      patch1%center = center_min
      patch1%normal = normal_min
      patch1%normal_mag = normal_mag_min
      call this%recompute_weights(merge_ids(1))

      !! Replace patch2 with last patch
      if (merge_ids(2) /= this%npatch) then
        call patch2%replace(this%patch(this%npatch))
        !! Fix f2p_map
        do i = 1, patch2%nface
          f = patch2%face(i)
          this%f2p_map(f) = merge_ids(2)
        end do
      end if

      this%npatch = this%npatch - 1
    end associate

  end subroutine patch_delete


  !! Initializes patch p with face f
  subroutine patch_init(this, f, p)
    class(vsa_patching), intent(inout) :: this
    integer, intent(in) :: f, p
    call this%patch(p)%init(f, this%area(f), this%center(:,f), this%normal(:,f))
  end subroutine patch_init


  !! Computes the weight of adding face f to patch p
  function compute_weight(this, f, p) result(ret)
    class(vsa_patching), intent(in) :: this
    integer, intent(in) :: f, p
    real(r8) :: ret
    ret = this%patch(p)%get_weight(this%center(:,f), this%normal(:,f), this%radius(f), this%max_radius, this%normalize)
  end function compute_weight


  !! Recomputes the weights of patch p
  subroutine recompute_weights(this, p)
    class(vsa_patching), intent(inout) :: this
    integer, intent(in) :: p
    real(r8) :: weight
    integer :: i, f
    associate (patch => this%patch(p))
      patch%total_weight = 0.0_r8
      do i = 1, patch%nface
        f = patch%face(i)
        weight = this%compute_weight(f, p)
        patch%weight(i) = weight
        patch%total_weight = patch%total_weight + weight
      end do
    end associate
  end subroutine recompute_weights


end module vsa_patching_type
