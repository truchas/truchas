MODULE BODY_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Body force modeling variables and procedures
  !
  ! Author(s):  Jim Sicilian (sicilian@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module, only: ndim

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Body Force Vector
  real(r8), dimension(ndim), public, save :: Body_Force

  ! Time-weighting for treatment of body forces
  real(r8), public, save :: body_force_implicitness

  ! Turn on face-based treatment of body-forces
  logical, public, save :: body_force_face_method

  ! potential + kinetic energy bounds
  real(r8), public, save :: mechanical_energy_bound

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

END MODULE BODY_DATA_MODULE
