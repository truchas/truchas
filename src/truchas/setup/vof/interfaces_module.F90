MODULE INTERFACES_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define variables necessary for the generation of initial interface
  !   locations.
  !
  ! Public Interface(s): None
  !
  ! Contains: None
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use parameter_module, only: mbody, msurf, mtab, ndim, string_len, &
                              nrot, mcoef, mphi
  use kind_module,      only: int_kind, real_kind

  implicit none

  ! Private Module
  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Body surface character strings
  character(LEN = string_len), dimension(msurf,mbody), public, save :: Surface
  character(LEN = string_len), dimension(msurf), public, save :: Fill, Surface_Name, &
                                                                 Axis, Tab_Type, File_Name

  ! Body surface integers
  integer(KIND = int_kind), dimension(mbody), public, save :: Matnum, Nsurf, Mesh_Matnum
  integer(KIND = int_kind), dimension(msurf,mbody), public, save :: Isrftype
  integer(KIND = int_kind), public, save :: body_surfaces,        &
                                            material_number,      &
                                            mesh_material_number, &
                                            nbody,                &
                                            num_args,             &
                                            background_body
  
  ! Body surface reals
  real(KIND = real_kind), dimension(msurf), public, save :: Height
  real(KIND = real_kind), public, save :: old_total_mass, total_mass
  real(KIND = real_kind), dimension(ndim,mbody), public, save :: Body_Vel
  real(KIND = real_kind), dimension(0:mtab,mbody), public, save :: Rtab, Ztab
  real(KIND = real_kind), dimension(2,mtab), public, save :: RZ_Tabular_Pt, Rtheta_Tabular_Pt
  real(KIND = real_kind), dimension(mbody), public, save :: Body_Mass, Body_Enthalpy,         &
                                                            Body_Temp, Body_Conc, Old_Body_Mass
  real(KIND = real_kind), dimension(mbody,mphi), public, save :: Body_Phi
  real(KIND = real_kind), dimension(ndim,msurf), public, save :: Rotation_Pt, Translation_Pt, &
                                                                 Radius, Length
  real(KIND = real_kind), dimension(nrot,msurf), public, save :: Rotation_Angle
  real(KIND = real_kind), dimension(nrot,msurf,mbody), public, save :: Rotangl, Cosa, Sina
  real(KIND = real_kind), dimension(ndim,msurf,mbody), public, save :: Offset, Rotpt, Ar
  real(KIND = real_kind), dimension(mcoef,msurf,mbody), public, save :: Sgeom

  ! VOF initialization controls
  character(string_len), public, save :: vof_method         ! 'points' or 'divide'
  integer (int_kind),    public, save :: int_particles      ! 'points' control
  integer (int_kind),    public, save :: vof_particles      ! 'points' control
  real (real_kind),      public, save :: vof_tolerance      ! 'divide' control
  integer (int_kind),    public, save :: vof_max_recursion  ! 'divide' control

END MODULE INTERFACES_MODULE
