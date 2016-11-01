!!
!! CHACO_PARTITIONER_TYPE
!!
!! A concrete implementation of the abstract GRAPH_PARTITIONER class that uses
!! the Chaco library to compute a graph partition.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! May 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!!  This is simple repackaging of existing code that maintains the practice
!!  of hardwiring the values of all the partitioner parameters.  Some of these
!!  ought to be exposed as user-specified input.  For large meshes it is not
!!  uncommon to have Chaco warn it was unable to generate a fully acceptable
!!  partition with the current set of parameters.
!!

#include "f90_assert.fpp"

module chaco_partitioner_type

  use graph_partitioner_class
  implicit none
  private

  type, extends(graph_partitioner), public :: chaco_partitioner
  contains
    procedure :: compute
  end type chaco_partitioner

contains

  subroutine compute (this, nvrtx, xadj, adjncy, ewgt, npart, part, stat, errmsg)

    use chaco_c_binding, only: interface
    use,intrinsic :: iso_c_binding, only: c_short, c_long, c_null_ptr, c_double
    use string_utilities, only: i_to_c

    class(chaco_partitioner), intent(inout) :: this
    integer, intent(in)  :: nvrtx, xadj(:), adjncy(:) ! the graph
    real,    intent(in)  :: ewgt(:) ! edge weights
    integer, intent(in)  :: npart   ! number of parts
    integer, intent(out) :: part(:) ! graph vertex partition vector
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(c_short) :: part_short(nvrtx)

    ASSERT(nvrtx >= 0)
    ASSERT(size(xadj) == nvrtx+1)
    ASSERT(size(adjncy) == xadj(nvrtx+1)-1)
    ASSERT(size(ewgt) == size(adjncy))
    ASSERT(npart >= 1)
    ASSERT(size(part) == nvrtx)

    !! Chaco C arrays use 0-based indexing, so the values in XADJ, which point
    !! into ADJNCY, must be offset by 1.  However, graph vertices are numbered
    !! from one (as here), so we do not need to offset the ADJNCY values.
    stat = interface(nvrtx, xadj-1, adjncy, c_null_ptr, ewgt, &
        c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr, &
        part_short, 1, 0, [npart, 0, 0], c_null_ptr, &
        1, 1, 0, 1000, 1, 1.0e-4_c_double, 1234567_c_long)
    if (stat /= 0) errmsg = 'chaco returned an error ('//i_to_c(stat)//')'

    !! Chaco uses a 0-based part numbering; we want 1-based numbering.
    part = part_short + 1

  end subroutine compute

end module chaco_partitioner_type
