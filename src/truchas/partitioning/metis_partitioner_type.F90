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

    integer(idx_t) :: ierr, objval
    integer(idx_t), target :: options(0:METIS_NOPTIONS-1)

    ASSERT(nvrtx >= 0)
    ASSERT(size(xadj) == nvrtx+1)
    ASSERT(size(adjncy) == xadj(nvrtx+1)-1)
    ASSERT(size(ewgt) == size(adjncy))
    ASSERT(npart >= 1)
    ASSERT(size(part) == nvrtx)

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

    call this%params%get('ctype', options(METIS_OPTION_CTYPE), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%params%get('ncuts', options(METIS_OPTION_NCUTS), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%params%get('niter', options(METIS_OPTION_NITER), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%params%get('ufactor', options(METIS_OPTION_UFACTOR), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%params%get('minconn', options(METIS_OPTION_MINCONN), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%params%get('contig', options(METIS_OPTION_CONTIG), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%params%get('seed', options(METIS_OPTION_SEED), default=-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call this%params%get('dbglvl', options(METIS_OPTION_DBGLVL), default=0, stat=stat, errmsg=errmsg)
    if (stat /= 0) return

    ierr = METIS_PartGraphKway(nvrtx, 1, xadj, adjncy, &
        c_null_ptr, c_null_ptr, c_null_ptr, npart, c_null_ptr, c_null_ptr, &
        c_loc(options), objval, part)
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

  end subroutine compute

end module metis_partitioner_type
