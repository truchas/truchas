module simplex_topology

  implicit none
  private

  integer, parameter, public :: TRI_EDGE_VERT(2,3) = &
      reshape(source=[2,3, 1,3, 1,2], shape=[2,3])

  integer, parameter, public :: TET_EDGE_VERT(2,6) = &
      reshape(source=[1,2, 1,3, 1,4, 2,3, 2,4, 3,4], shape=[2,6])
    
  integer, parameter, public :: TET_FACE_VERT(3,4) = &
      reshape(source=[2,3,4, 1,3,4, 1,2,4, 1,2,3], shape=[3,4])
    
  integer, parameter, public :: TET_FACE_EDGE(3,4) = &
      reshape(source=[6,5,4, 6,3,2, 5,3,1, 4,2,1], shape=[3,4])

end module simplex_topology
