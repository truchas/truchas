#include "f90_assert.fpp"

MODULE BC_DATA_MODULE
  !=======================================================================
  ! Purpose:
  !
  !   Defines data variables for boundary conditions. Note: This module is
  !   declared public and the variables are saved upon use of the module.
  !
  ! Author(s): Jerry S. Brock (jsbrock@lanl.gov)
  !            Douglas B. Kothe (dbk@lanl.gov)
  !            Bryan Lally (lally@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use bc_type_module,   only: BOUNDARY_CONDITION
  use parameter_module, only: bc_forms, ndim, nbcs, nvar, mbc_surfaces, &
                              string_len, mbcsrf, mbc_nodes, max_bc_dof
  use velocity_boundary_data_type

  implicit none

  ! Save all variables
  save

  ! Public module
  public

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Namelist input variables
  character(80), dimension(0:mbc_surfaces)            :: BC_Type
  character(80), dimension(0:mbc_surfaces)            :: BC_Name
  real(r8),      dimension(max_bc_dof,0:mbc_surfaces) :: BC_Value
  real(r8),      dimension(4,16,0:mbc_surfaces)       :: BC_Table
  character(80), dimension(0:mbc_surfaces)            :: BC_Variable

  integer,  dimension(0:mbc_surfaces)       :: Inflow_Material
  real(r8), dimension(0:mbc_surfaces)       :: Inflow_Temperature

  real(r8), dimension(mbcsrf,0:mbc_surfaces)      :: Conic_XX
  real(r8), dimension(mbcsrf,0:mbc_surfaces)      :: Conic_YY
  real(r8), dimension(mbcsrf,0:mbc_surfaces)      :: Conic_ZZ
  real(r8), dimension(mbcsrf,0:mbc_surfaces)      :: Conic_XY
  real(r8), dimension(mbcsrf,0:mbc_surfaces)      :: Conic_XZ
  real(r8), dimension(mbcsrf,0:mbc_surfaces)      :: Conic_YZ
  real(r8), dimension(mbcsrf,0:mbc_surfaces)      :: Conic_X
  real(r8), dimension(mbcsrf,0:mbc_surfaces)      :: Conic_Y
  real(r8), dimension(mbcsrf,0:mbc_surfaces)      :: Conic_Z
  real(r8), dimension(mbcsrf,0:mbc_surfaces)      :: Conic_Constant
  real(r8), dimension(mbcsrf,0:mbc_surfaces)      :: Conic_Tolerance
  character(string_len), dimension(mbcsrf,0:mbc_surfaces) :: Conic_Relation
  character(string_len), dimension(mbcsrf,0:mbc_surfaces) :: Surface_Name
  integer, dimension(2,mbcsrf,0:mbc_surfaces)  :: Surface_Materials
  real(r8), dimension(2,ndim,0:mbc_surfaces)      :: Bounding_Box


  integer,  dimension(mbcsrf,0:mbc_surfaces)  :: Mesh_Surface  

  ! Variables needed in processing the namelist input
  integer, dimension(mbcsrf,mbc_surfaces)      :: Srfmatl_Index
  integer,  dimension(mbc_surfaces)         :: Inflow_Index
  integer,  dimension(mbc_surfaces)         :: surfaces_in_this_bc
  integer                                   :: nbc_surfaces
  character(80),         dimension(nbcs,bc_forms)        :: Type_Forms
  character(80),         dimension(nvar,bc_forms)        :: Variable_Forms 
  character(80),         dimension(mbc_surfaces)         :: BC_Surface_Forms 
  ! Coordinates for individual node displacement BCs (read from input file)   
  real(r8), dimension(ndim,mbc_nodes,0:mbc_surfaces)         :: Node_Disp_Coords
  ! Mask arrays for face and node sets read from mesh file 
  integer, pointer, dimension(:,:,:)           :: Mesh_Face_Set     => NULL()
  integer, pointer, dimension(:,:,:)           :: Mesh_Face_Set_Tot => NULL()
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Define Arrays for Dirichlet BC .........................
  ! Applied concentration values (Each face of each cell)
  real(r8), pointer, dimension(:,:) :: BC_Conc => null()

  ! Applied pressure values (Each face of each cell)
  !     BC_Zero is used for pressure change solution
  !     BC_Pressure contains the actual pressures at the faces
  !     BC_Prs points to either BC_Zero or BC_Pressure as appropriate
  real(r8), pointer, dimension(:,:) :: BC_Pressure => null()
  real(r8), pointer, dimension(:,:) :: BC_Zero => null()
  real(r8), pointer, dimension(:,:) :: BC_Prs => NULL()

  ! Inflow temperature values.
  real(r8), pointer, dimension(:,:) :: BC_Temp => null()

  ! Applied velocity values (3 components on each face of each cell)
  integer, pointer, dimension(:,:)   :: BC_Mat => null()
  !! NNC, Jan 2014.  Time-dependent dirichlet velocity.
  !ORIG: real(r8),   pointer, dimension(:,:,:) :: BC_Vel => null()
  type(velocity_boundary_data) :: bndry_vel

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Concentration (C) BC flags, values, and offsets .......
  !
  ! The flag word is divided into 6 sections of 4 bits each, one for
  ! each face. The meaning of the bits in each section is:
  !
  ! bits 3 2 1 0 (value)   type
  !      ----------------------
  !        0 0 0    0      no boundary
  !        0 0 1    1      Dirichlet
  !        0 1 0    2      Homogeneous Neumann
  !        0 1 1    3      Neumann
  !        1 0 0    4      reserved
  !        1 0 1    5      reserved
  !        1 1 0    6      reserved
  !        1 1 1    7      reserved
  !      0                 external boundary
  !      1          8      internal boundary
  !
  ! The values necessary to specify the BCs are stored in
  ! BC_C_Value#(nfc,ncells)

  ! Concentration (C) Constants
  integer, parameter :: BC_C_NO_BC         = 0
  integer, parameter :: BC_C_DIRICHLET     = 1
  integer, parameter :: BC_C_HNEUMANN      = 2
  integer, parameter :: BC_C_NEUMANN       = 3
  integer, parameter :: BC_C_RESERVED1     = 4
  integer, parameter :: BC_C_RESERVED2     = 5
  integer, parameter :: BC_C_RESERVED3     = 6
  integer, parameter :: BC_C_RESERVED4     = 7

  integer, parameter :: BC_C_SHIFT = 4 ! Bit section size per face
  integer, parameter :: BC_C_EXIST = 7 ! Mask - zero means no BC
  integer, parameter :: BC_C_IFLAG = 8 ! External boundary bit

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Heat Transfer (T) BC flags, values, and offsets ........
  !
  ! The flag word is divided into 6 sections of 5 bits each, one for
  ! each face.  The meaning of the bits in each section is:
  !
  ! bits 4 3 2 1 0 (value)   type
  !      ----------------------
  !        0 0 0 0    0      no boundary
  !        0 0 0 1    1      Dirichlet
  !        0 0 1 0    2      HNeumann
  !        0 0 1 1    3      Neumann
  !        0 1 0 0    4      heat transfer coefficient
  !        0 1 0 1    5      radiation
  !        0 1 1 0    6      heat transfer coefficient + radiation
  !        0 1 1 1    7      reflective
  !        1 0 0 0    8      viewfactor based radiation
  !      0 1 1 1 1    15     internal boundary
  !      1 0 0 0 0    16     external boundary
  !
  ! The values necessary to specify the BCs are stored in
  ! BC_T_Value#(nfc,ncells)

  ! Heat Transfer (T) Constants
  integer, parameter :: BC_T_NO_BC         = 0
  integer, parameter :: BC_T_DIRICHLET     = 1
  integer, parameter :: BC_T_HNEUMANN      = 2
  integer, parameter :: BC_T_NEUMANN       = 3
  integer, parameter :: BC_T_HTC           = 4
  integer, parameter :: BC_T_RADIATION     = 5
  integer, parameter :: BC_T_HTC_RADIATION = 6
  integer, parameter :: BC_T_REFLECTIVE    = 7
  integer, parameter :: BC_T_VFRADIATION   = 8
  integer, parameter :: BC_T_HTC_GAP       = 9

  integer, parameter :: BC_T_SHIFT =  5 ! Bit section size per face
  integer, parameter :: BC_T_EXIST = 15 ! Mask - zero means no BC
  integer, parameter :: BC_T_IFLAG = 16 ! External boundary bit

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Pressure (P) BC flags, values, and offsets ............
  !
  ! The flag word is divided into 6 sections of 4 bits each, one for
  ! each face. The meaning of the bits in each section is:
  !
  ! bits 3 2 1 0 (value)   type
  !      ----------------------
  !        0 0 0    0      no boundary
  !        0 0 1    1      Dirichlet
  !        0 1 0    2      Homogeneous Neumann
  !        0 1 1    3      Neumann
  !        1 0 0    4      reserved
  !        1 0 1    5      reserved
  !        1 1 0    6      reserved
  !        1 1 1    7      reserved
  !      0                 external boundary
  !      1          8      internal boundary
  !
  ! The values necessary to specify the BCs are stored in
  ! BC_P_Value#(nfc,ncells)

  ! Pressure (P) Constants
  integer, parameter :: BC_P_NO_BC         = 0
  integer, parameter :: BC_P_DIRICHLET     = 1
  integer, parameter :: BC_P_HNEUMANN      = 2
  integer, parameter :: BC_P_NEUMANN       = 3
  integer, parameter :: BC_P_REFLECTIVE    = 4
  integer, parameter :: BC_P_RESERVED2     = 5
  integer, parameter :: BC_P_RESERVED3     = 6
  integer, parameter :: BC_P_RESERVED4     = 7

  integer, parameter :: BC_P_SHIFT = 4 ! Bit section size per face
  integer, parameter :: BC_P_EXIST = 7 ! Mask - zero means no BC
  integer, parameter :: BC_P_IFLAG = 8 ! External boundary bit

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Velocity (V) BC flags, values, and offsets ............
  !
  ! The flag word is divided into 6 sections of 4 bits each, one for
  ! each face. The meaning of the bits in each section is:
  !
  ! bits 3 2 1 0 (value)   type
  !      ----------------------
  !        0 0 0    0      no boundary
  !        0 0 1    1      Dirichlet
  !        0 1 0    2      Homogeneous Neumann
  !        0 1 1    3      Neumann
  !        1 0 0    4      reserved
  !        1 0 1    5      reserved
  !        1 1 0    6      reserved
  !        1 1 1    7      reserved
  !      0                 external boundary
  !      1          8      internal boundary
  !
  ! The values necessary to specify the BCs are stored in
  ! BC_V_Value#(nfc,ncells)

  ! Velocity (V) Constants
  integer, parameter :: BC_V_NO_BC         = 0
  integer, parameter :: BC_V_DIRICHLET     = 1
  integer, parameter :: BC_V_HNEUMANN      = 2
  integer, parameter :: BC_V_NEUMANN       = 3
  integer, parameter :: BC_V_RESERVED1     = 4
  integer, parameter :: BC_V_RESERVED2     = 5
  integer, parameter :: BC_V_RESERVED3     = 6
  integer, parameter :: BC_V_RESERVED4     = 7

  integer, parameter :: BC_V_SHIFT = 4 ! Bit section size per face
  integer, parameter :: BC_V_EXIST = 7 ! Mask - zero means no BC
  integer, parameter :: BC_V_IFLAG = 8 ! External boundary bit

  type(BOUNDARY_CONDITION), dimension(:), pointer :: BC => null()

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SKIP_SIDE_SET_DATA
 !!
 !! Neil N. Carlson <nnc@lanl.gov>
 !! 26 Mar 2010
 !!
 !! Version 2 restart files contain side set data (dropped in version 3).
 !! That data is now read from the mesh file and is no longer needed.  To be
 !! able to use old version 2 restart files, this routine will, if needed,
 !! skip over the side set data in a restart file opened (and pre-positioned)
 !! on UNIT.  VERSION is the restart file version.  For reference, the
 !! original code to read the side set data follows as a comment.
 !!

  subroutine skip_side_set_data (unit, version)

    use restart_utilities, only: read_var, skip_records

    integer, intent(in) :: unit, version
    
    integer :: nssets ! local value only

    if (version <= 2) then
      call read_var (unit, nssets, 'SKIP_SIDE_SET_DATA: error reading NSSETS record')
      nssets = max(0, nssets)
      if (nssets == 0) return ! no further side set records to skip
      call skip_records (unit, 6*nssets, 'SKIP_SIDE_SET_DATA: error skipping SIDESET data')
    end if

  end subroutine skip_side_set_data

!  subroutine read_side_set_data (unit, version)
!
!    use parameter_module, only: ncells, nfc, nssets
!    use mesh_module, only: pcell => unpermute_mesh_vector
!    use restart_utilities, only: read_var, read_dist_array
!
!    integer, intent(in) :: unit, version
!
!    ASSERT( .not.associated(mesh_face_set) )
!
!    !! Read the number of side sets.
!    call read_var (unit, nssets, 'READ_SIDE_SET_DATA: error reading NSSETS record')
!    nssets = max(0, nssets)
!    if (nssets == 0) return ! no further records to read
!
!    !! Allocate the MESH_FACE_SET array and read the array data.
!    allocate(mesh_face_set(nssets,nfc,ncells))
!    call read_dist_array (unit, mesh_face_set, pcell, 'READ_SIDE_SET_DATA: error reading SIDESET records')
!
!  end subroutine read_side_set_data

END MODULE BC_DATA_MODULE
