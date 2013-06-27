MODULE PARAMETER_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define all integer and scalar parameters.
  !
  ! Contains: None
  !
  ! Author(s): Douglas B. Kothe, LANL (dbk@lanl.gov)
  !
  !=======================================================================
  implicit none
  public
  
  save

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

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

  ! Number of rotation axes
  integer, parameter :: nrot = 2*ndim - 3

  ! Maximum number of outputs
  integer, parameter :: mops = 21

  ! Maximum output character string array size
  integer, parameter :: string_dim = 30, string_len = 256

  ! Maximum number of BC namelists and surfaces
  integer, parameter :: mbc_surfaces = 100
  integer, parameter :: mbcsrf       = 5
  integer, parameter :: mbc_nodes    = 50
  !  The radiation BC requires max_bc_dof be at least 4 greater than the number of time,reference-temperature pairs
  integer, parameter :: max_bc_dof   = 24 ! Maximum number of degrees of freedom any kind of BC
  

  ! Number of different BC variable character string forms. 
  integer, parameter :: bc_forms = 30

  ! Number of material relation character strings allowed
  integer, parameter :: max_relation_forms = 15

  ! Number of interface topology character strings allowed
  integer, parameter :: max_topology_models = 5

  ! Current number of BC types and variables
  integer, parameter :: nbcs = 32, nvar = 10

  ! Interface (body) initialization parameters
  integer, parameter :: mtype = 10
  integer, parameter :: msurf = 16
  integer, parameter :: mbody = 50
  integer, parameter :: mtab  = 50
  integer, parameter :: mcoef =  3
  integer, parameter :: mphi  =  5

  integer, parameter :: mregion =  50

  ! Maximum number of material constants          
  integer, parameter :: max_slots = 10
  integer, parameter :: maxmat = 64
  integer, parameter :: maxcon = 10
  integer, parameter :: maxform = 8

  ! Maximum number of mesh segments
  integer, parameter :: mseg = 21

  ! Maximum number of mesh domains
  integer, parameter :: max_domains = 4196

  ! Maximum number of long edit cells
  integer, parameter :: max_long_edit_cells = 500

  ! Maximum number of long edit boxes
  integer, parameter :: max_long_edit_boxes = 20

  ! Maximum number of alloy components allowed
  integer, parameter :: max_alloy_comps = 10

  ! Maximum number of phases allowed
  integer, parameter :: max_phases = 10

  ! Number of stress/strain components
  integer, Parameter :: ncomps = (ndim*(ndim + 1))/2

  ! Maximum number of chemical reactions allowed
  integer, parameter :: max_cr = 10

  ! Upper-to-lower character conversion constant
  integer, parameter :: upper_to_lower = IACHAR('a') -IACHAR('A')

  ! Array sizes
  integer, dimension(ndim) :: Nx, Mx, Nx_tot, Mx_tot
  integer :: nmat, mmat, nnodes, nnodes_tot, ncells,     &
                              ncells_tot, nicells, nicells_tot,           &
                              nfaces, boundary_faces, boundary_faces_tot, &
                              mat_slot = 0, mat_slot_new = 0,             &
                              mat_slot_tmp = 0, mat_slot_tmp_new = 0
  ! Number of side sets
  integer :: nssets = 0

  ! maximum number of probes allowed in input file
  integer, parameter :: MAX_PROBES = 20
  integer            :: nprobes

  ! Maximum time intervals for electromagnetics

  integer, parameter, public :: MAXSV = 32

END MODULE PARAMETER_MODULE
