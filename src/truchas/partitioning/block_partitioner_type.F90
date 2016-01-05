!!
!! BLOCK_PARTITIONER_TYPE
!!
!! A concrete implementation of the abstract GRAPH_PARTITIONER class that
!! computes a simple block partition of the graph vertices into equally-sized
!! parts, ignoring the graph connectivity.
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

#include "f90_assert.fpp"

module block_partitioner_type

  use graph_partitioner_class
  implicit none
  private

  type, extends(graph_partitioner), public :: block_partitioner
  contains
    procedure :: compute
  end type block_partitioner

contains

  subroutine compute (this, nvrtx, xadj, adjncy, ewgt, npart, part, stat)

    class(block_partitioner), intent(inout) :: this
    integer, intent(in)  :: nvrtx, xadj(:), adjncy(:) ! the graph (not used here)
    real,    intent(in)  :: ewgt(:) ! edge weights (not used here)
    integer, intent(in)  :: npart   ! number of parts
    integer, intent(out) :: part(:) ! graph vertex partition vector
    integer, intent(out) :: stat

    integer :: j, n, xpart(npart+1)

    ASSERT(nvrtx >= 0)
    ASSERT(npart >= 1)
    ASSERT(size(part) == nvrtx)

    xpart(1) = 1
    do j = 1, npart
      n = nvrtx/npart
      if (j <= modulo(nvrtx, npart)) n = n + 1
      xpart(j+1) = xpart(j) + n
    end do
    ASSERT(xpart(npart+1) == nvrtx + 1)

    do j = 1, npart
      part(xpart(j):xpart(j+1)-1) = j
    end do

    stat = 0

  end subroutine compute

end module block_partitioner_type
