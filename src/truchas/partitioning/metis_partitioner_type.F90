!!
!! METIS_PARTITIONER_TYPE
!!
!! A concrete implementation of the abstract GRAPH_PARTITIONER class that uses
!! the Metis library to compute a graph partition.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module metis_partitioner_type

  use graph_partitioner_class
  use parameter_list_type
  implicit none
  private

  type, extends(graph_partitioner), public :: metis_partitioner
    private
    type(parameter_list), pointer :: params => null()
    integer :: npart_major
  contains
    procedure :: compute
  end type

  !! Defined constructor
  interface metis_partitioner
    procedure :: metis_partitioner_params
  end interface

contains

  function metis_partitioner_params(params) result(this)
    type(parameter_list), intent(inout), target :: params
    type(metis_partitioner) :: this
    this%params => params
  end function

  subroutine compute(this, nvrtx, xadj, adjncy, ewgt, npart, part, stat, errmsg)

    use metis_c_binding
    use string_utilities, only: i_to_c
    use,intrinsic :: iso_c_binding, only: c_loc, c_null_ptr

    class(metis_partitioner), intent(inout) :: this
    integer(idx_t), intent(in)  :: nvrtx, xadj(:), adjncy(:) ! the graph
    real(real_t),   intent(in)  :: ewgt(:) ! edge weights
    integer(idx_t), intent(in)  :: npart   ! number of parts
    integer(idx_t), intent(out) :: part(:) ! graph vertex partition vector
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(idx_t) :: p, offset, npart_major, npart_minor, nvrtx_minor
    real(real_t), allocatable :: tpwgts(:), ewgt_minor(:)
    integer(idx_t), allocatable :: part_major(:), part_minor(:)
    integer(idx_t), allocatable :: sub_xadj(:), sub_adjncy(:)
    logical, allocatable :: mask(:), emask(:)
    integer :: j, n

    call this%params%get('major-partitions', npart_major, stat, errmsg, default=1_idx_t)
    if (stat /= 0) return
    if (npart_major < 1) then
      stat = 1
      errmsg = 'major-partitions must be >= 1'
      return
    end if

    if (npart_major == 1) then

      call compute_core(this, nvrtx, xadj, adjncy, ewgt, npart, null(), part, stat, errmsg)

    else

      if (npart_major > npart) then
        stat = 1
        errmsg = 'major-partitions exceeds the number of partitions'
        return
      end if

      !! First-level, major partitions

      ! Target major-partition weights
      tpwgts = npart / npart_major + merge(1, 0, [(p, p=1, npart_major)] <= modulo(npart, npart_major))
      tpwgts = tpwgts / npart ! must sum to 1

      allocate(part_major, mold=part)
      call compute_core(this, nvrtx, xadj, adjncy, ewgt, npart_major, tpwgts, part_major, stat, errmsg)
      if (stat /= 0) return
      INSIST(minval(part_major) == 1 .and. maxval(part_major) == npart_major)

      !! Second-level, minor partitions
      part = 0
      offset = 0
      do p = 1, npart_major
        mask = (part_major == p)
        npart_minor = npart / npart_major
        if (p <= modulo(npart, npart_major)) npart_minor = npart_minor + 1
        if (npart_minor > 1) then ! partition the major partition
          call get_subgraph(xadj, adjncy, mask, sub_xadj, sub_adjncy, emask)
          nvrtx_minor = size(sub_xadj) - 1
          ewgt_minor = pack(ewgt, emask)
          allocate(part_minor(nvrtx_minor))
          call compute_core(this, nvrtx_minor, sub_xadj, sub_adjncy, ewgt_minor, npart_minor, null(), part_minor, stat, errmsg)
          if (stat /= 0) return
          INSIST(minval(part_minor) == 1 .and. maxval(part_minor) == npart_minor)
          part = unpack(offset+part_minor, mask, part)
          deallocate(part_minor)
        else
          where (mask) part = offset + 1
        end if
        offset = offset + npart_minor
      end do
      INSIST(minval(part) == 1 .and. maxval(part) == npart)
    end if

  end subroutine compute

  subroutine compute_core(this, nvrtx, xadj, adjncy, ewgt, npart, tpwgts, part, stat, errmsg)

    use metis_c_binding
    use string_utilities, only: i_to_c
    use,intrinsic :: iso_c_binding, only: c_loc, c_null_ptr

    class(metis_partitioner), intent(inout) :: this
    integer(idx_t), intent(in)  :: nvrtx, xadj(:), adjncy(:) ! the graph
    real(real_t),   intent(in)  :: ewgt(:) ! edge weights
    real(real_t),   intent(in), optional  :: tpwgts(:) ! target partition weights
    integer(idx_t), intent(in)  :: npart   ! number of parts
    integer(idx_t), intent(out) :: part(:) ! graph vertex partition vector
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(idx_t) :: ierr, objval
    integer(idx_t), target :: options(0:METIS_NOPTIONS-1)
    integer :: ptype

    ASSERT(nvrtx >= 0)
    ASSERT(size(xadj) == nvrtx+1)
    ASSERT(size(adjncy) == xadj(nvrtx+1)-1)
    ASSERT(size(ewgt) == size(adjncy))
    ASSERT(npart >= 1)
    ASSERT(size(part) == nvrtx)
    !ASSERT(size(tpwgts) == npart)

    ierr = METIS_SetDefaultOptions(options)
    INSIST(ierr == METIS_OK)  ! really this should never fail

    options(METIS_OPTION_NUMBERING) = 1 ! Fortran 1-based array indexing

    !! One might think that the above call to METIS_SetDefaultOptions would set
    !! the options to their defaults, but no, it simply sets them to -1, which
    !! apparently indicates to METIS_PartGraphKway to use the builtin defaults,
    !! whatever they might be. Here, to simplify the code, we go ahead reassign
    !! that default value. To see what the defaults actually are you need to run
    !! with DBGLVL=1.

    !! Despite what the documentation suggests IPTYPE only applies to recursive
    !! bisection, RTYPE value is entirely ignored

    call this%params%get('iptype', options(METIS_OPTION_IPTYPE), stat, errmsg, default=-1)
    if (stat /= 0) return
    call this%params%get('ctype', options(METIS_OPTION_CTYPE), stat, errmsg, default=-1)
    if (stat /= 0) return
    call this%params%get('ncuts', options(METIS_OPTION_NCUTS), stat, errmsg, default=-1)
    if (stat /= 0) return
    call this%params%get('niter', options(METIS_OPTION_NITER), stat, errmsg, default=-1)
    if (stat /= 0) return
    call this%params%get('ufactor', options(METIS_OPTION_UFACTOR), stat, errmsg, default=-1)
    if (stat /= 0) return
    call this%params%get('minconn', options(METIS_OPTION_MINCONN), stat, errmsg, default=-1)
    if (stat /= 0) return
    call this%params%get('contig', options(METIS_OPTION_CONTIG), stat, errmsg, default=-1)
    if (stat /= 0) return
    call this%params%get('seed', options(METIS_OPTION_SEED), stat, errmsg, default=-1)
    if (stat /= 0) return
    call this%params%get('dbglvl', options(METIS_OPTION_DBGLVL), stat, errmsg, default=0)
    if (stat /= 0) return

    call this%params%get('ptype', ptype, stat, errmsg, default=METIS_PTYPE_RB)
    if (stat /= 0) return
    select case (ptype)
    case (METIS_PTYPE_RB)
      ierr = METIS_PartGraphRecursive(nvrtx, 1, xadj, adjncy, &
          c_null_ptr, c_null_ptr, c_null_ptr, npart, tpwgts, c_null_ptr, &
          c_loc(options), objval, part)
    case (METIS_PTYPE_KWAY)
      ierr = METIS_PartGraphKway(nvrtx, 1, xadj, adjncy, &
          c_null_ptr, c_null_ptr, c_null_ptr, npart, tpwgts, c_null_ptr, &
          c_loc(options), objval, part)
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

  end subroutine compute_core

  !! Return the subgraph corresponding to the nodes with true mask value.
  !! Also return the associated edge mask
  subroutine get_subgraph(xadj, adjncy, mask, sub_xadj, sub_adjncy, emask)

    use metis_c_binding

    integer(idx_t), intent(in) :: xadj(:), adjncy(:)
    logical, intent(in) :: mask(:)
    integer(idx_t), allocatable, intent(out) :: sub_xadj(:), sub_adjncy(:)
    logical, allocatable, intent(out) :: emask(:)

    integer(idx_t) :: i, j, n
    integer(idx_t), allocatable :: map(:)

    allocate(map(size(mask)), sub_xadj(count(mask)+1))

    n = 0
    map = 0
    sub_xadj(1) = 1
    do j = 1, size(mask)
      if (.not.mask(j)) cycle
      n = n + 1
      map(j) = n
      sub_xadj(n+1) = sub_xadj(n)
      do i = xadj(j), xadj(j+1)-1
        if (mask(adjncy(i))) sub_xadj(n+1) = sub_xadj(n+1) + 1
      end do
    end do
    allocate(sub_adjncy(sub_xadj(n+1)-1))

    n = 0
    allocate(emask(size(adjncy)), source=.false.)
    do j = 1, size(mask)
      if(.not.mask(j)) cycle
      do i = xadj(j), xadj(j+1)-1
        if (mask(adjncy(i))) then
          n = n + 1
          sub_adjncy(n) = map(adjncy(i))
          emask(i) = .true.
        end if
      end do
    end do
    INSIST(size(sub_adjncy) == count(emask))

  end subroutine get_subgraph

end module metis_partitioner_type
