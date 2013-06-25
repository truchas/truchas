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
  use kind_module, only: int_kind

  implicit none
  save

  ! Public Module
  public

  private :: int_kind

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Number of physical dimensions
  integer(int_kind), parameter :: ndim = 3

  ! Number of faces per cell
  integer(int_kind), parameter :: nfc = 2*ndim

  ! Number of vertices per cell
  integer(int_kind), parameter :: nvc = 2**ndim

  ! Number of vertices per face
  integer(int_kind), parameter :: nvf = 2**(ndim - 1)

  ! Number of faces per vertex
  integer(int_kind), parameter :: nfv = ndim

  ! Number of edges per cell
  integer(int_kind), parameter :: nec = ndim*2**(ndim - 1)

  ! Number of rotation axes
  integer(int_kind), parameter :: nrot = 2*ndim - 3

  ! Maximum number of outputs
  integer(int_kind), parameter :: mops = 21

  ! Maximum output character string array size
  integer(int_kind), parameter :: string_dim = 30, string_len = 256

  ! Maximum number of BC namelists and surfaces
  integer(int_kind), parameter :: mbc_surfaces = 100
  integer(int_kind), parameter :: mbcsrf       = 5
  integer(int_kind), parameter :: mbc_nodes    = 50
  !  The radiation BC requires max_bc_dof be at least 4 greater than the number of time,reference-temperature pairs
  integer(int_kind), parameter :: max_bc_dof   = 24 ! Maximum number of degrees of freedom any kind of BC
  

  ! Number of different BC variable character string forms. 
  integer(int_kind), parameter :: bc_forms = 30

  ! Number of material relation character strings allowed
  integer(int_kind), parameter :: max_relation_forms = 15

  ! Number of interface topology character strings allowed
  integer(int_kind), parameter :: max_topology_models = 5

  ! Current number of BC types and variables
  integer(int_kind), parameter :: nbcs = 32, nvar = 10

  ! Interface (body) initialization parameters
  integer(int_kind), parameter :: mtype = 10
  integer(int_kind), parameter :: msurf = 16
  integer(int_kind), parameter :: mbody = 50
  integer(int_kind), parameter :: mtab  = 50
  integer(int_kind), parameter :: mcoef =  3
  integer(int_kind), parameter :: mphi  =  5

  integer(int_kind), parameter :: mregion =  50

  ! Maximum number of material constants          
  integer(int_kind), parameter :: max_slots = 10
  integer(int_kind), parameter :: maxmat = 64
  integer(int_kind), parameter :: maxcon = 10
  integer(int_kind), parameter :: maxform = 8

  ! Maximum number of mesh segments
  integer(int_kind), parameter :: mseg = 21

  ! Maximum number of mesh domains
  integer(int_kind), parameter :: max_domains = 4196

  ! Maximum number of long edit cells
  integer(int_kind), parameter :: max_long_edit_cells = 500

  ! Maximum number of long edit boxes
  integer(int_kind), parameter :: max_long_edit_boxes = 20

  ! Maximum number of alloy components allowed
  integer(int_kind), parameter :: max_alloy_comps = 10

  ! Maximum number of phases allowed
  integer(int_kind), parameter :: max_phases = 10

  ! Number of stress/strain components
  Integer(KIND =  int_kind), Parameter :: ncomps = (ndim*(ndim + 1))/2

  ! Maximum number of chemical reactions allowed
  integer(KIND = int_kind), parameter :: max_cr = 10

  ! Upper-to-lower character conversion constant
  integer(int_kind), parameter :: upper_to_lower = IACHAR('a') -IACHAR('A')

  ! Array sizes
  integer(KIND = int_kind), dimension(ndim) :: Nx, Mx, Nx_tot, Mx_tot
  integer(KIND = int_kind) :: nmat, mmat, nnodes, nnodes_tot, ncells,     &
                              ncells_tot, nicells, nicells_tot,           &
                              nfaces, boundary_faces, boundary_faces_tot, &
                              mat_slot = 0, mat_slot_new = 0,             &
                              mat_slot_tmp = 0, mat_slot_tmp_new = 0
  ! Number of side sets
  integer(KIND = int_kind) :: nssets = 0

   ! maximum number of probes allowed in input file
   integer (int_kind), parameter :: MAX_PROBES = 20
   integer (int_kind)            :: nprobes

   ! Maximum time intervals for electromagnetics

   integer, parameter, public :: MAXSV = 32

END MODULE PARAMETER_MODULE
