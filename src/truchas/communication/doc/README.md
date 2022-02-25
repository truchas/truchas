# Parallel Communication in Truchas

Truchas uses an application-specific layer over MPI to encapsulate (and
hide) its use of MPI. Two modules comprise this layer. The first is
`parallel_communication`. It handles the most basic parallel operations:
broadcasting data from one process to all processes, distributing data
from one process across all processes, collating distributed data onto
one process, and various reductions of distributed data. The second is
`index_map_type` which handles the complex parallel halo exchange
operations that are associated with the domain decomposition scheme
Truchas uses to parallelize computation.

* [`parallel_communication`](./parallel_communication.md)
* [`index_map_type`](./index_map_type.md)
