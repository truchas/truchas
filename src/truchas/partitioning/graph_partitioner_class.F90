!!
!! GRAPH_PARTITIONER_CLASS
!!
!! A common interface to graph partitioners.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! May 2015
!!

module graph_partitioner_class

  type, abstract, public :: graph_partitioner
  contains
    procedure(gpc), deferred :: compute
  end type graph_partitioner

  abstract interface
    subroutine gpc (this, nvrtx, xadj, adjncy, ewgt, npart, part, stat)
      import graph_partitioner
      class(graph_partitioner), intent(inout) :: this
      integer, intent(in)  :: nvrtx, xadj(:), adjncy(:) ! the graph
      real,    intent(in)  :: ewgt(:) ! edge weights
      integer, intent(in)  :: npart   ! number of parts
      integer, intent(out) :: part(:) ! graph vertex partition vector
      integer, intent(out) :: stat
    end subroutine
  end interface

end module graph_partitioner_class
