!! NNC, October 2015.  Separate out parameters/variables from PARAMETER_MODULE
!! that belong specifically to the original mesh data structure.  This will be
!! jettisoned along with the rest of the original mesh.

MODULE MESH_PARAMETER_MODULE

  implicit none
  public
  
  ! Number of physical dimensions
  integer, parameter :: ndim = 3

  ! Number of faces per cell
  integer, parameter :: nfc = 2*ndim

  ! Number of vertices per cell
  integer, parameter :: nvc = 2**ndim

  ! Number of vertices per face
  integer, parameter :: nvf = 2**(ndim - 1)

  ! Number of faces per vertex
  integer, parameter :: nfv = ndim

  ! Number of edges per cell
  integer, parameter :: nec = ndim*2**(ndim - 1)

  ! Array sizes
  integer :: nnodes, nnodes_tot, ncells, ncells_tot

  ! Number of side sets
  integer :: nssets = 0

END MODULE MESH_PARAMETER_MODULE
