!!
!! TOOLPATH_TYPE
!!
!! This module defines the derived type TOOLPATH that is intended to describe
!! the 3D path taken by a machine tool, such as a laser, cutting head, etc.
!! It is designed to emulate the paths that are typical of CNC processes.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! June 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The TOOLPATH type represents a 3D path as a sequence of continuous path
!!  segments on time intervals [-\infty,t1], [t1,t2], ..., [tn,\infty].  The
!!  individual segments are typically simple paths (e.g., no motion, linear
!!  motion), and associated with each segment is a set of flags (0 through 31)
!!  that are either set or clear.  A flag can be used to indicate whether the
!!  laser is on or off for the duration of the segment, for example.  While
!!  the overall path is typically continuous (though this isn't strictly
!!  required), the flags may change their settings from one segment to the
!!  next, and so the path is regarded as having potential discontinuities
!!  at the segment end points.  For this reason the toolpath has the notion
!!  of a current path segment.  Queries about position or flag settings are
!!  made using the info for this segment, including at the endpoints.  Client
!!  code makes the determination when to advance the toolpath to the next
!!  segment.
!!
!!  The type has the following type bound procedures:
!!
!!  GET_POSITION(T, R) returns the path position R at time T for the current
!!     path segment.  T must belong to the time interval for the segment.
!!
!!  IS_FLAG_SET(N) returns true if flag N is set for the current path segment.
!!
!!  SET_SEGMENT(T) sets the current path segment to the one whose time interval
!!    contains T, but is not bounded above by T.
!!
!!  NEXT_SEGMENT() advances the current path segment to the next segment.
!!
!!  GET_SEGMENT_STARTS(TIMES, DISCONT) returns the list of segment end-points
!!    t1, t2, ..., tn in the array TIMES.  The corresponding logical array
!!    DISCONT returns true if any of the flags changed their setting at the
!!    times.  Both arrays are allocatable and allocated by the subroutine.
!!
!!  IS_VALID() returns true if the toolpath data appears valid.  What this
!!    means is that there is a sequence of paths and that they properly tile
!!    the entire real line (time).  It does not examine the individual paths
!!    beyond their time end points.
!!
!!  The following partition subroutines are optional and have no impact on
!!  the behavior of the preceding procedures.  They assume (and check for)
!!  an overall continuous path.
!!
!!  SET_PARTITION(DS) generates an internal partition of the path that is
!!    subordinate to the path segments.  This is an ordered sequence of times
!!    and corresponding path coordinates that includes the end points of the
!!    path segments, and that also includes equally spaced points within each
!!    segment separated by an arclength approximately equal to, but no greater
!!    than, DS.  The generated partition data can be retrived with GET_PARTITION.
!!
!!  GET_PARTITION(TIME, COORD, HASH) retrieves a previously generated partition.
!!    all arguments are optional allocatable arrays, and must be specified with
!!    keywords.  TIME(:) is the ordered list of times, COORD(3,:) the list of
!!    path coordinates, and HASH(:) is a character array (deferred length)
!!    of a hash of the coordinates (first 7 characters of the sha1 hash), which
!!    are useful for naming files associated with the coordinates, for example.
!!    There is guaranteed to be no hash collision.
!!
!!  This module does not provide a method for creating a toolpath object (see
!!  the TOOLPATH_FACTORY module), but it does provide some basic tools for doing
!!  so.  The following functions return a pointer to a path segment:
!!
!!  NEW_PATH_SEGMENT(MOVE, FLAGS) creates a segment with path MOVE and flag
!!    settings given by the integer FLAGS (the flags are the individual bits
!!    of the integer). MOVE is an allocatable class XYZ_MOTION variable. Its
!!    allocation is taken by the segment result and returned unallocated.
!!
!!  NEW_START_SEGMENT(T, R, FLAGS) creates a segment with constant path at
!!    position R for t <= T, and flag settings given by FLAGS.  The initial
!!    toolpath segment must be one of these.
!!
!!  NEW_FINAL_SEGMENT(T, R, FLAGS) creates a segment with constant path at
!!    position R for t >= T, and flag settings given by FLAGS.  The final
!!    toolpath segment must be one of these.
!!
!!  Use of these path segment function expressions is limited to the actual
!!  argument of the TOOLPATH type bound subroutine APPEND_PATH_SEGMENT, as in
!!
!!    TYPE(TOOLPATH) :: TP
!!    CALL TP%APPEND_PATH_SEGMENT(NEW_START_SEGMENT(...))
!!
!!  See the TOOLPATH_FACTORY module for example usage.
!!
!! IMPLEMENTATION NOTES
!!
!! 1. The partition capability is driven by the moving radiation enclosure
!! use case.  Here a toolpath describes the motion of part of the enclosure.
!! Computing the changing view factor matrix at each time step is impractical.
!! Instead we compute (and save, typically offline) the view factor matrix
!! at a select set of configurations.  That is the partition.
!!
!! The current implementation feels a bit of a hack, being optional extra stuff
!! tacked onto the toolpath type. A more sensible design would simply generate
!! external partition data from a toolpath, but the problem is both enclosure
!! radiation and genre need this info and there is no good owner for the data.
!! Hence the expedient choice to let the toolpath hold it internally. A possible
!! alternative approach worth exploring is to define an extension of the
!! toolpath type that adds the extra data.
!!

#include "f90_assert.fpp"

module toolpath_type

  use xyz_motion_class
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, public :: toolpath
    private
#ifdef INTEL18_WORKAROUND
    type(path_segment), pointer :: first => null()
    type(path_segment), pointer :: last => null()
    type(path_segment), pointer :: curr => null()
#else
    type(path_segment), pointer :: first => null(), last => null(), curr => null()
#endif
    !! Optional partition of the path subordinate to the segments
    real(r8), allocatable :: time(:), coord(:,:)
    character(7), allocatable :: hash(:)
  contains
    procedure :: get_position
    procedure :: is_flag_set
    procedure :: start_time
    procedure :: final_time
    procedure :: prec_flag_flip_time
    procedure :: set_segment
    procedure :: next_segment
    procedure :: get_segment_starts
    procedure :: is_valid => valid_toolpath
    procedure :: write_plotfile
    procedure :: set_partition
    procedure :: get_partition
    procedure :: has_partition
    procedure :: append_path_segment
    final :: toolpath_delete
  end type

  type :: path_segment
    private
    type(path_segment), pointer :: next => null(), prev => null()
    class(xyz_motion), allocatable :: move
    integer :: flags = 0
  contains
    final :: path_segment_delete
  end type

  public :: new_path_segment, new_start_segment, new_final_segment

  !! Barrier used to force use of argument keywords
  type :: kwarg_barrier
  end type

contains

  !! Final subroutine for TOOLPATH objects.
  subroutine toolpath_delete (this)
    type(toolpath), intent(inout) :: this
    if (associated(this%first)) deallocate(this%first)
  end subroutine toolpath_delete

  !! Sets the current segment to the segment containing the given time.
  subroutine set_segment(this, t)
    class(toolpath), intent(inout) :: this
    real(r8), intent(in) :: t
    this%curr => this%first
    do while (associated(this%curr))
      if (t < this%curr%move%final_time()) exit
      this%curr => this%curr%next
    end do
    !INSIST(associated(this%curr))
    !INSIST(t >= this%curr%move%start_time())
  end subroutine set_segment

  !! Returns the toolpath position at the given time.
  subroutine get_position(this, t, r)
    class(toolpath), intent(in) :: this
    real(r8), intent(in)  :: t
    real(r8), intent(out) :: r(:)
    ASSERT(t >= this%curr%move%start_time())
    ASSERT(t <= this%curr%move%final_time())
    r = this%curr%move%coord(t)
  end subroutine get_position

  !! Returns true if flag N is set.
  logical function is_flag_set(this, n)
    class(toolpath), intent(in) :: this
    integer, intent(in) :: n
    ASSERT(n >= 0 .and. n < bit_size(this%curr%flags))
    is_flag_set = btest(this%curr%flags,n)
  end function is_flag_set

  !! Advance to the next path segment.
  subroutine next_segment(this)
    class(toolpath), intent(inout) :: this
    if (associated(this%curr%next)) this%curr => this%curr%next
  end subroutine next_segment

  !! Start time for the current path segment.
  function start_time(this) result(t)
    class(toolpath), intent(in) :: this
    real(r8) :: t
    t = this%curr%move%start_time()
  end function start_time

  !! Final time for the current path segment.
  function final_time(this) result(t)
    class(toolpath), intent(in) :: this
    real(r8) :: t
    t = this%curr%move%final_time()
  end function final_time

  !! Returns the most recent time preceding the current path segment when
  !! flag N flipped state, or -HUGE(1.0_r8) if flag N did not flip state.
  function prec_flag_flip_time(this, n) result(t)
    class(toolpath), intent(in) :: this
    integer, intent(in) :: n
    real(r8) :: t
    type(path_segment), pointer :: seg
    ASSERT(associated(this%curr))
    seg => this%curr
    do while (associated(seg%prev))
      if (btest(seg%prev%flags,n) .neqv. btest(this%curr%flags,n)) exit
      seg => seg%prev
    end do
    t = seg%move%start_time()
  end function prec_flag_flip_time

  !! Appends the given PATH_SEGMENT (as made by NEW_PATH_SEGMENT) to the list.
  subroutine append_path_segment(this, segment)
    class(toolpath), intent(inout) :: this
    type(path_segment), pointer, intent(in) :: segment
    if (associated(this%last)) then
      this%last%next => segment
      segment%prev => this%last
    else
      this%first => segment
    end if
    this%last => segment
  end subroutine append_path_segment

  subroutine get_segment_starts(this, times, discont)
    class(toolpath), intent(in) :: this
    real(r8), allocatable, intent(out) :: times(:)
    logical, allocatable, intent(out), optional :: discont(:)
    type(path_segment), pointer :: seg
    integer :: n
    integer :: flags
    n = 0; seg => this%first
    do while (associated(seg))
      n = n + 1
      seg => seg%next
    end do
    n = n - 1
    !INSIST(n > 0)
    allocate(times(n))
    if (present(discont)) allocate(discont(n))
    flags = this%first%flags
    n = 0
    seg => this%first%next
    do n = 1, size(times)
      times(n) = seg%move%start_time()
      if (present(discont)) discont(n) = (seg%flags /= flags)
      flags = seg%flags
      seg => seg%next
    end do
  end subroutine get_segment_starts

  logical function valid_toolpath(this)
    class(toolpath), intent(in) :: this
    real(r8) :: t
    type(path_segment), pointer :: seg
    valid_toolpath = .false.
    if (.not.associated(this%first)) return
    if (this%first%move%start_time() /= -huge(1.0_r8)) return
    if (.not.associated(this%last)) return
    if (this%last%move%final_time() /= huge(1.0_r8)) return
    t = this%first%move%final_time()
    seg => this%first%next
    do while (associated(seg))
      if (t /= seg%move%start_time()) return
      if (t >= seg%move%final_time()) return
      t = seg%move%final_time()
      seg => seg%next
    end do
    valid_toolpath = .true.
  end function valid_toolpath

  subroutine write_plotfile(this, unit, dt)

    class(toolpath), intent(in) :: this
    integer, intent(in) :: unit
    real(r8), intent(in) :: dt

    integer :: i, j, n, nseg
    real(r8) :: t0, t1, t, r(3)
    type(path_segment), pointer :: seg
    integer(kind(seg%flags)) :: bitmask
    integer, allocatable :: flags(:)

    n = 0; bitmask = 0; seg => this%first
    do while (associated(seg))
      n = n + 1
      bitmask = ior(bitmask,seg%flags)
      seg => seg%next
    end do
    nseg = n - 2 ! not counting first and last
    ASSERT(nseg >= 0)

    flags = [(n,n=0,bit_size(bitmask)-1)]
    flags = pack(flags, mask=btest(bitmask,pos=flags))

    write(unit,'(a,:,*("  F",i0,:))') '#  N  T  X  Y  Z', flags

    seg => this%first
    do i = 0, nseg+1
      if (i == 0) then  ! initial unbounded segment
        t1 = seg%move%final_time()
        n = 2
        t0 = t1 - n*dt
      else if (i == nseg+1) then ! final unbounded segment
        t0 = seg%move%start_time()
        n = 2
        t1 = t0 + n*dt
      else
        t0 = seg%move%start_time()
        t1 = seg%move%final_time()
        n = ceiling((t1-t0)/dt)
      end if
      do j = 0, n
        t = t0*(real(n-j,r8)/real(n,r8)) + t1*(real(j,r8)/real(n,r8))
        r = seg%move%coord(t)
        write(unit,'(i0,4es15.7,*(i2))') i, t, r, merge(1,0,btest(seg%flags,pos=flags))
      end do
      seg => seg%next
    end do

  end subroutine write_plotfile

  logical function is_continuous(this)
    class(toolpath), intent(in) :: this
    type(path_segment), pointer :: seg
    real(r8) :: r(3)
    ASSERT(associated(this%first))
    is_continuous = .false.
    r = this%first%move%final_coord()
    seg => this%first%next
    do while (associated(seg))
      if (any(seg%move%start_coord() /= r)) return
      r = seg%move%final_coord()
      seg => seg%next
    end do
    is_continuous = .true.
  end function is_continuous

  subroutine set_partition(this, ds)

    use sha1_hash_type

    class(toolpath), intent(inout) :: this
    real(r8), intent(in) :: ds

    integer :: n, j
    real(r8), allocatable :: times(:)
    type(path_segment), pointer :: seg
    type(sha1_hash) :: hash

    ASSERT(ds > 0)
    INSIST(this%is_valid())
    INSIST(is_continuous(this))

    n = 1; seg => this%first%next
    do while (associated(seg%next))
      times = seg%move%partition(ds)
      n = n + size(times) - 1
      seg => seg%next
    end do

    if (allocated(this%time)) then
      if (size(this%time) /= n) deallocate(this%time, this%coord, this%hash)
    end if
    if (.not.allocated(this%time)) allocate(this%time(n), this%coord(3,n), this%hash(n))

    n = 1; seg => this%first%next
    this%time(1) = seg%move%start_time()
    this%coord(:,1) = seg%move%start_coord()
    call hash%update(this%coord(:,1))
    this%hash(1) = hash%hexdigest()
    do while (associated(seg%next))
      times = seg%move%partition(ds)
      do j = 2, size(times)
        n = n + 1
        this%time(n) = times(j)
        this%coord(:,n) = seg%move%coord(times(j))
        call hash%update(this%coord(:,n))
        this%hash(n) = hash%hexdigest()
      end do
      seg => seg%next
    end do

    !! Ensure there are no hash collisions.
    do n = 2, size(this%hash)
      do j = 1, n-1
        if (this%hash(j) == this%hash(n)) then
          INSIST(all(this%coord(:,j) == this%coord(:,n)))
        end if
      end do
    end do

  end subroutine set_partition

  subroutine get_partition(this, kwarg, time, coord, hash)
    class(toolpath), intent(in) :: this
    type(kwarg_barrier), optional :: kwarg
    real(r8), allocatable, intent(out), optional :: time(:), coord(:,:)
    character(:), allocatable, intent(out), optional :: hash(:)
    INSIST(allocated(this%time))
    if (present(time)) time = this%time
    if (present(coord)) coord = this%coord
    if (present(hash)) hash = this%hash
  end subroutine get_partition

  logical function has_partition(this)
    class(toolpath), intent(in) :: this
    has_partition = allocated(this%time)
  end function has_partition

  !! Final subroutine for PATH_SEGMENT objects.  This recursively follows the
  !! NEXT pointer.  When deallocating a linked list, only the root needs to be
  !! explicitly deallocated. When the desire is to deallocate a single object,
  !! first nullify the NEXT pointer to prevent the recursive finalization from
  !! deallocating more than it should.
  recursive subroutine path_segment_delete (this)
    type(path_segment), intent(inout) :: this
    if (associated(this%next)) deallocate(this%next)
  end subroutine path_segment_delete

  !! Returns a pointer to a new initialized (but unlinked) PATH_SEGMENT.
  function new_path_segment(move, flags) result(segment)
    class(xyz_motion), allocatable, intent(inout) :: move
    integer, intent(in) :: flags
    type(path_segment), pointer :: segment
    allocate(segment)
    call move_alloc(move, segment%move)
    segment%flags = flags
  end function new_path_segment

  !! Returns a pointer to a new (but unlinked) start PATH_SEGMENT.
  function new_start_segment(t, r, flags) result(segment)
    use dwell_xyz_motion_type
    real(r8), intent(in) :: t, r(:)
    integer,  intent(in) :: flags
    type(path_segment), pointer :: segment
    allocate(segment)
#ifdef NO_2008_LHS_POLY_REALLOC
    allocate(segment%move, source=dwell_xyz_motion(r=r, t1=t))
#else
    segment%move = dwell_xyz_motion(r=r, t1=t)
#endif
    segment%flags = flags
  end function new_start_segment

  !! Returns a pointer to a new (but unlinked) final PATH_SEGMENT.
  function new_final_segment(t, r, flags) result(segment)
    use dwell_xyz_motion_type
    real(r8), intent(in) :: t, r(:)
    integer,  intent(in) :: flags
    type(path_segment), pointer :: segment
    allocate(segment)
#ifdef NO_2008_LHS_POLY_REALLOC
    allocate(segment%move, source=dwell_xyz_motion(r=r, t0=t))
#else
    segment%move = dwell_xyz_motion(r=r, t0=t)
#endif
    segment%flags = flags
  end function new_final_segment

end module toolpath_type
