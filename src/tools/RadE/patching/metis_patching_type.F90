!!
!! METIS_PATCHING
!!
!! A concrete implementation of the abstract base class PATCHING that
!! encapsulates the METIS graph partitioner as a method of clustering
!! the faces of a radiation enclosure mesh.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! January 2021
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! METIS_PATCHING is an extension of the abstract base class PATCHING.  See the
!! base class comments for a description of the common type bound procedures.
!! Additional details specific to the METIS_PATCHING type are included below.
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
!!           the number of parts that METIS will partition the dual graph into,
!!           and therefore total the number of patches in the output:
!!             [# PARTS] = [# FACES] / FP_RATIO
!!           This is only an upper bound, since METIS is free to produce less
!!           partitions.  Consider changing the METIS options to ensure the
!!           desired number of partitions.
!!       face-weight -- determines whether the face areas should be used as
!!           vertex weights in the dual graph passed to METIS.
!!       edge-weight -- determines whether the lengths of enclosure edges should
!!           be used as edge weights in the dual graph passed to METIS.
!!       metis-options -- parameter list containing other options passed to
!!           METIS.  The following options are supported:
!!              ptype, ctype, iptype, objtype, no2hop, contig, minconn, ufactor,
!!              niter, ncuts, seed, dbglvl
!!           Refer to the METIS manual for details on these options.
!!


#include "f90_assert.fpp"

module metis_patching_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use re_encl_type
  use patching_class
  use parameter_list_type
  implicit none

  !! Parameter defaults
  real(r8), parameter :: METIS_FACE_PATCH_RATIO_DEFAULT = 4.0_r8
  logical, parameter :: METIS_FACE_WEIGHT_DEFAULT = .true.
  logical, parameter :: METIS_EDGE_WEIGHT_DEFAULT = .true.

  type, extends(patching), public :: metis_patching
    private
    type(encl), pointer :: e
    integer, allocatable  :: xfnhbr(:), fnhbr(:)  ! face adjacency graph
    integer, allocatable  :: vwgt(:)  ! vertex weights
    integer, allocatable  :: ewgt(:)  ! edge weights
    real(r8), allocatable :: area(:), normal(:,:)
    integer, allocatable  :: f2p_map(:)
    integer :: npatch, npatch_init
    integer :: verbosity
    type(parameter_list), pointer :: metis_opt => null()

    contains
      procedure, public  :: init
      procedure, public  :: run
      procedure, public  :: output
      procedure, private  :: print_stats
      procedure, private  :: split_disconnected
  end type

contains

  !! Allocate and initialize METIS_PATCHING data
  subroutine init(this, e, params, stat, errmsg)

    use cell_geometry, only: face_normal, vector_length, normalized
    use patching_tools, only: get_face_neighbor_array
    use parameter_list_json, only: parameter_list_to_json
    use, intrinsic :: iso_fortran_env, only : output_unit

    class(metis_patching), intent(out) :: this
    type(encl), intent(in), target :: e
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable :: length(:)
    integer, allocatable  :: xfnhbr(:), fnhbr(:)  ! Face adjacency graph
    integer, allocatable :: edge(:)
    logical, allocatable :: marker(:)
    real(r8) :: max_angle  ! maximum allowable angle for adjacent faces (in degrees)
    real(r8) :: fp_ratio   ! total faces / total patches
    real(r8) :: min_area, min_length, norm(3)
    integer :: i, j, n
    logical :: face_weight, edge_weight
    character(:), allocatable :: context

    !! Process the parameters.
    context = 'processing ' // params%path() // ': '
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
    call params%get('face-weight', face_weight, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    call params%get('edge-weight', edge_weight, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    this%metis_opt => params%sublist('metis-options', stat, errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if

    this%e => e
    this%npatch_init = e%nface / fp_ratio

    allocate(this%area(e%nface), this%normal(3,e%nface))
    allocate(this%f2p_map(e%nface))

    !! Compute face areas and normals
    min_area = huge(0.0_r8)
    do i = 1, e%nface
      associate (face_nodes => e%fnode(e%xface(i):e%xface(i+1)-1))
        norm = face_normal(e%x(:,face_nodes))
        this%area(i) = vector_length(norm)
        this%normal(:,i) = normalized(norm)
        min_area = min(min_area, this%area(i))
      end associate
    end do

    !! Get face adjacency matrix
    call get_face_neighbor_array(e%xface, e%fnode, xfnhbr, fnhbr, stat, errmsg, this%normal, max_angle)
    if (stat /= 0) then
      block
        integer :: n
        character(40) :: coord
        n = findloc(fnhbr, -1, dim=1) ! first side with invalid topology
        write(coord,'("(",es12.5,2(",",es12.5),")")') e%side_location(n)
        errmsg = errmsg // ' at ' // coord
      end block
      return
    end if

    !! Compress face adjacency matrix. We must remove the 0s, which denote "no neighbor".
    allocate(this%xfnhbr(e%nface+1))
    this%fnhbr = pack(fnhbr, mask=(fnhbr > 0))
    this%xfnhbr(1) = 1
    do i = 1, e%nface
      this%xfnhbr(i+1) = this%xfnhbr(i) + count(fnhbr(xfnhbr(i):xfnhbr(i+1)-1) > 0)
    end do
    ASSERT(this%xfnhbr(e%nface+1)-1 == size(this%fnhbr))

    !! Dual graph vertex weights
    if (face_weight) then
      allocate(this%vwgt(this%e%nface))
      this%vwgt = nint(this%area/min_area)
    end if

    !! Compute edge lengths
    if (edge_weight) then
      allocate(length(size(this%fnhbr)))
      allocate(marker(this%e%nnode), source=.false.)
      min_length = huge(0.0_r8)
      do i = 1, e%nface
        associate (face_nodes => e%fnode(e%xface(i):e%xface(i+1)-1))
          marker(face_nodes) = .true.
          do j = this%xfnhbr(i), this%xfnhbr(i+1)-1
            n = this%fnhbr(j)  ! neighbor of face i
            associate (nhbr_nodes => e%fnode(e%xface(n):e%xface(n+1)-1))
              edge = pack(nhbr_nodes, marker(nhbr_nodes))
              ASSERT(size(edge) == 2)
              length(j) = vector_length(e%x(:,edge(1)) - e%x(:,edge(2)))
              min_length = min(min_length, length(j))
            end associate
          end do
          marker(face_nodes) = .false.
        end associate
      end do

      !! Dual graph edge weights
      allocate(this%ewgt(size(this%fnhbr)))
      this%ewgt = nint(length/min_length)
    end if

    if (this%verbosity > 1) then
      print '("INITIALIZING METIS:")'
      print '("  FACE-PATCH RATIO: ", f0.2)', fp_ratio
      print '("  NFACE: ", i0)', e%nface
      print '("  NPATCH: ", i0)', this%npatch_init
      print '("  MAX ANGLE: ", f0.2)', max_angle
      print '("  FACE WEIGHT: ", l1)', face_weight
      print '("  EDGE WEIGHT: ", l1)', edge_weight
    end if
    if (this%verbosity > 2) then
      print '("  METIS OPTIONS: ")'
      call parameter_list_to_json(this%metis_opt, output_unit)
    end if

  end subroutine init


  !! The main METIS patching algorithm
  subroutine run(this, stat, errmsg)

    use metis_c_binding
    use string_utilities, only: i_to_c
    use,intrinsic :: iso_c_binding, only: c_loc, c_ptr, c_null_ptr

    class(metis_patching), target, intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer, allocatable :: map(:)
    integer(idx_t) :: ierr, objval
    integer(idx_t), target :: options(0:METIS_NOPTIONS-1)
    type(c_ptr) :: vwgt, ewgt
    integer :: i, ptype

    stat = 0

    vwgt = merge(c_loc(this%vwgt), c_null_ptr, allocated(this%vwgt))
    ewgt = merge(c_loc(this%ewgt), c_null_ptr, allocated(this%ewgt))

    !! Set METIS options
    ierr = METIS_SetDefaultOptions(options)
    INSIST(ierr == METIS_OK)  ! really this should never fail

    options(METIS_OPTION_NUMBERING) = 1 ! Fortran 1-based array indexing

    call this%metis_opt%get('ctype', options(METIS_OPTION_CTYPE), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%metis_opt%get('iptype', options(METIS_OPTION_IPTYPE), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%metis_opt%get('objtype', options(METIS_OPTION_OBJTYPE), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%metis_opt%get('no2hop', options(METIS_OPTION_NO2HOP), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%metis_opt%get('contig', options(METIS_OPTION_CONTIG), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%metis_opt%get('minconn', options(METIS_OPTION_MINCONN), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%metis_opt%get('ufactor', options(METIS_OPTION_UFACTOR), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%metis_opt%get('niter', options(METIS_OPTION_NITER), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%metis_opt%get('ncuts', options(METIS_OPTION_NCUTS), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%metis_opt%get('seed', options(METIS_OPTION_SEED), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%metis_opt%get('dbglvl', options(METIS_OPTION_DBGLVL), default=0, stat=stat, errmsg=errmsg)
    if (stat /= 0) return

    call this%metis_opt%get('ptype', ptype, default=METIS_PTYPE_RB, stat=stat, errmsg=errmsg)
    if (stat /= 0) return

    !! Run METIS
    select case (ptype)
    case (METIS_PTYPE_RB)
    ierr = METIS_PartGraphRecursive(this%e%nface, 1, this%xfnhbr, this%fnhbr, &
        vwgt, c_null_ptr, ewgt, this%npatch_init, c_null_ptr, c_null_ptr, &
        c_loc(options), objval, this%f2p_map)
    case (METIS_PTYPE_KWAY)
    ierr = METIS_PartGraphKway(this%e%nface, 1, this%xfnhbr, this%fnhbr, &
        vwgt, c_null_ptr, ewgt, this%npatch_init, c_null_ptr, c_null_ptr, &
        c_loc(options), objval, this%f2p_map)
    case default
      stat = 1
      errmsg = 'unknown metis partitioning method ptype=' // i_to_c(ptype)
      return
    end select

    if (ierr /= METIS_OK) then
      stat = 1
      select case (ierr)
      case (METIS_ERROR_INPUT)
        errmsg = 'metis input error'
      case (METIS_ERROR_MEMORY)
        errmsg = 'metis memory allocation error'
      case default
        errmsg = 'metis returned an error'
      end select
      return
    end if

    !! Compress patches.  METIS might produce less parts than requested.
    allocate(map(this%npatch_init), source=0)
    map(this%f2p_map) = 1
    do i = 1, this%npatch_init-1
      map(i+1) = map(i) + map(i+1)
    end do
    this%f2p_map = map(this%f2p_map)
    this%npatch = map(this%npatch_init)

    !! Ensure patches are connected
    call this%split_disconnected

    if (this%verbosity > 1) call this%print_stats
    if (this%verbosity > 0) then
      print '(/,"PATCHING COMPLETE:")'
      print '("  NFACE: ", i0)', this%e%nface
      print '("  INITIAL NPATCH: ",i0)', this%npatch_init
      print '("  ACTUAL NPATCH: ",i0," (",f0.2,"%)")', this%npatch, this%npatch/real(this%npatch_init)*100
      print '("  NFACE / MAX NPATCH: ",f0.2)', this%e%nface/real(this%npatch_init)
      print '("  NFACE / ACTUAL NPATCH: ",f0.2)', this%e%nface/real(this%npatch)
    end if

  end subroutine run


  !! Creates new patches for each connected component of a patch
  subroutine split_disconnected(this)

    use sort_utilities
    use patching_tools, only: get_connected_faces_subset

    class(metis_patching), intent(inout) :: this

    integer, allocatable :: xpface(:), pface(:)  ! Faces of a patch
    integer, allocatable :: tag(:)
    integer :: i, j, prev, offset, ncomp, cnt, npatch_old

    !! Sorting f2p_map yields the faces of a patch
    allocate(pface(this%e%nface))
    call heap_sort(this%f2p_map, pface)

    !! Get offset of each patch
    allocate(xpface(this%npatch+1), source=-1)
    prev = -1
    offset = 1
    do i = 1, this%e%nface
      !! Found next patch
      if (prev /= this%f2p_map(pface(i))) then
        xpface(offset) = i
        offset = offset + 1
        prev = this%f2p_map(pface(i))
      end if
    end do
    ASSERT(offset==this%npatch+1)
    xpface(this%npatch+1) = this%e%nface + 1

    !! Create new patches
    cnt = 0
    npatch_old = this%npatch
    do i = 1, npatch_old
      associate (faces => pface(xpface(i):xpface(i+1)-1))
        call get_connected_faces_subset(this%xfnhbr, this%fnhbr, faces, tag, ncomp)
        if (ncomp == 1) cycle  ! skip connected patches
        cnt = cnt + 1
        !! Assign new patch IDs. Faces of first component keep current ID.
        do j = 2, ncomp
          this%npatch = this%npatch + 1
          this%f2p_map(pack(faces, tag==j)) = this%npatch
        end do
      end associate
    end do

    if (this%verbosity > 1) then
      print '(/,"DISCONNECTED PATCHES")'
      print '("  Found ",i0," disconnected patches.")', cnt
      print '("  Split them into ",i0," new patches.")', this%npatch-npatch_old
    end if

  end subroutine split_disconnected


  !! Write the patch data
  subroutine output(this, f2p_map, global_ids, npatch)

    class(metis_patching), intent(inout) :: this
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

    use patching_tools, only: PI
    use cell_geometry, only: normalized
    class(metis_patching), intent(inout) :: this

    real(r8) :: nfp(this%npatch)  ! Number of faces in each patch
    real(r8) :: tpa(this%npatch)  ! Total patch area
    real(r8) :: fpna(this%e%nface)  ! Face to patch normal angle
    real(r8) :: normal(3,this%npatch)  ! Area-weighted normal
    integer, allocatable :: perm(:)
    real(r8) :: angle
    integer :: i

    nfp = 0.0_r8
    tpa = 0.0_r8
    normal = 0.0_r8

    !! Face-centered data
    do i = 1, this%e%nface
      nfp(this%f2p_map(i)) = nfp(this%f2p_map(i)) + 1
      tpa(this%f2p_map(i)) = tpa(this%f2p_map(i)) + this%area(i)
      normal(:,this%f2p_map(i)) = normal(:,this%f2p_map(i)) + this%normal(:,i)*this%area(i)
    end do

    !! Face-normal deviation
    do i = 1, this%npatch
      normal(:,i) = normalized(normal(:,i))
    end do
    do i = 1, this%e%nface
      angle = dot_product(normal(:,this%f2p_map(i)), this%normal(:,i))
      if (angle > 1) angle = 1  ! Fix floating point errors
      fpna(i) = acos(angle)*180.0_r8/PI
    end do

    print '(/,"PATCH STATS:")'
    print '("  NFACE:  ", i8)', this%e%nface
    print '("  NPATCH: ", i8)', this%npatch


    print '("  PATCH INDEXED:")'
    allocate(perm(this%npatch))
    call summary(nfp, perm, "FACES PER PATCH")
    call summary(tpa, perm, "PATCH AREA")
    deallocate(perm)

    print '("  FACE INDEXED:")'
    allocate(perm(this%e%nface))
    call summary(fpna, perm, "FACE VS. PATCH NORMAL ANGLE")

  contains
    subroutine summary(arr, perm, title)

      use sort_utilities

      real(r8), intent(in) :: arr(:)
      integer, intent(inout) :: perm(:)
      character(*), intent(in) :: title
      real(r8) :: avg, std
      integer :: idx(5)

      call heap_sort(arr, perm)
      avg = sum(arr) / size(arr)
      std = sqrt(sum((arr - avg)**2) / size(arr))

      idx(1) = 1
      idx(2:) = size(arr) * [0.25, 0.5, 0.75, 1.0]

      print '("    ", a, ": ")', title
      print '("      MEAN ± STD DEV:", es10.3, " ±", es10.3)', avg, std
      print '("      MIN, 25%, 50%, 75%, MAX: ",4(es10.3,", "),es10.3)', arr(perm(idx))
      print '("                        INDEX: ",4(i10,", "),i10)', perm(idx)

    end subroutine summary
  end subroutine print_stats

end module metis_patching_type
