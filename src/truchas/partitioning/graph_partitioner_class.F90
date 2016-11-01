!!
!! GRAPH_PARTITIONER_CLASS
!!
!! A common interface to graph partitioners.
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

module graph_partitioner_class

  type, abstract, public :: graph_partitioner
  contains
    procedure(gpc), deferred :: compute
  end type graph_partitioner

  abstract interface
    subroutine gpc (this, nvrtx, xadj, adjncy, ewgt, npart, part, stat, errmsg)
      import graph_partitioner
      class(graph_partitioner), intent(inout) :: this
      integer, intent(in)  :: nvrtx, xadj(:), adjncy(:) ! the graph
      real,    intent(in)  :: ewgt(:) ! edge weights
      integer, intent(in)  :: npart   ! number of parts
      integer, intent(out) :: part(:) ! graph vertex partition vector
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
  end interface

end module graph_partitioner_class
