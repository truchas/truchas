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
  use kinds, only: r8
  use parameter_module, only: mbody, msurf, mtab, ndim, string_len, &
                              nrot, mcoef, mphi
  implicit none
  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Body surface character strings
  character(string_len), dimension(msurf,mbody), public, save :: Surface
  character(string_len), dimension(msurf), public, save :: Fill, Surface_Name, &
                                                           Axis, Tab_Type, File_Name

  ! Body surface integers
  integer, dimension(mbody), public, save :: Matnum, Nsurf, Mesh_Matnum
  integer, dimension(msurf,mbody), public, save :: Isrftype
  integer, public, save :: body_surfaces, material_number, mesh_material_number, &
                           nbody, num_args, background_body
  
  ! Body surface reals
  real(r8), dimension(msurf), public, save :: Height
  real(r8), public, save :: old_total_mass, total_mass
  real(r8), dimension(ndim,mbody), public, save :: Body_Vel
  real(r8), dimension(0:mtab,mbody), public, save :: Rtab, Ztab
  real(r8), dimension(2,mtab), public, save :: RZ_Tabular_Pt, Rtheta_Tabular_Pt
  real(r8), dimension(mbody), public, save :: Body_Mass, Body_Enthalpy,         &
                                                            Body_Temp, Body_Conc, Old_Body_Mass
  real(r8), dimension(mbody,mphi), public, save :: Body_Phi
  real(r8), dimension(ndim,msurf), public, save :: Rotation_Pt, Translation_Pt, &
                                                                 Radius, Length
  real(r8), dimension(nrot,msurf), public, save :: Rotation_Angle
  real(r8), dimension(nrot,msurf,mbody), public, save :: Rotangl, Cosa, Sina
  real(r8), dimension(ndim,msurf,mbody), public, save :: Offset, Rotpt, Ar
  real(r8), dimension(mcoef,msurf,mbody), public, save :: Sgeom

  ! VOF initialization controls
  character(string_len), public, save :: vof_method         ! 'points' or 'divide'
  integer,  public, save :: int_particles      ! 'points' control
  integer,  public, save :: vof_particles      ! 'points' control
  real(r8), public, save :: vof_tolerance      ! 'divide' control
  integer,  public, save :: vof_max_recursion  ! 'divide' control

END MODULE INTERFACES_MODULE
