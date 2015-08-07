!!
!! GRAPH_PARTITIONER_FACTORY
!!
!! A procedure for instantiating a new GRAPH_PARTITIONER class object.
!! The dynamic type of the object is determined by parameter list input.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! May 2015
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL ALLOC_GRAPH_PARTITIONER (THIS, PARAMS) allocates the allocatable
!!    CLASS(GRAPH_PARTITIONER) argument THIS.  The parameter list PARAMS
!!    expects a single parameter 'partitioner' whose value is either 'block'
!!    or 'chaco'.  For 'block', the dynamic type of THIS is BLOCK_PARTITIONER
!!    which computes a simple block partition of the graph vertices (ignoring
!!    the graph connectivity).  For 'chaco', the dynamic type of THIS is
!!    CHACO_PARTITIONER, which uses the Chaco library to compute the partition.
!!
!! IMPLEMENTATION NOTES
!!
!!  Use of a parameter list argument was intended to prepare the way for adding
!!  partitioiner-specifc parameters.  Currently these parameters are hardwired
!!  in the chaco implementation.
!!

module graph_partitioner_factory

  use graph_partitioner_class
  use parameter_list_type
  implicit none
  private

  public :: graph_partitioner ! re-export
  public :: alloc_graph_partitioner

contains

  subroutine alloc_graph_partitioner (this, params)

    use block_partitioner_type
#ifdef USE_CHACO
    use chaco_partitioner_type
#endif
    use truchas_logging_services
    use string_utilities, only: raise_case

    class(graph_partitioner), allocatable, intent(out) :: this
    type(parameter_list) :: params

    integer :: stat
    character(:), allocatable :: partitioner, errmsg

    call params%get ('partitioner', partitioner, stat=stat, errmsg=errmsg)
    if (stat /= 0) call TLS_fatal ('ALLOC_GRAPH_PARTITIONER: ' // errmsg)
    select case (raise_case(partitioner))
    case ('BLOCK')
      allocate(this, source=block_partitioner())
    case ('CHACO')
#ifdef USE_CHACO
      allocate(this, source=chaco_partitioner())
#else
      call TLS_fatal ('ALLOC_GRAPH_PARTITIONER: chaco partitioning not available')
#endif
    case default
      call TLS_fatal ('ALLOC_GRAPH_PARTITIONER: unknown "partitioner": ' // partitioner)
    end select

  end subroutine alloc_graph_partitioner

end module graph_partitioner_factory