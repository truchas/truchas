MODULE BODY_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Body force modeling variables and procedures
  !
  ! Author(s):  Jim Sicilian (sicilian@lanl.gov)
  !
  !=======================================================================
  use kind_module,      only: real_kind, log_kind
  use parameter_module, only: ndim

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Body Force Vector
  real(real_kind), dimension(ndim), public, save :: Body_Force

  ! Time-weighting for treatment of body forces
  real(real_kind),   public, save :: body_force_implicitness

  ! Turn on face-based treatment of body-forces
  logical(log_kind), public, save :: body_force_face_method

  ! potential + kinetic energy bounds
  real(real_kind),   public, save :: mechanical_energy_bound

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

END MODULE BODY_DATA_MODULE
