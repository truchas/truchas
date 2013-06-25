module porous_drag_data

  use kind_module, only: int_kind, log_kind, real_kind

  ! PHYSICS namelist variables -
  ! Flag for enabling/disabling the porous drag model.
  logical(KIND = log_kind), public, save :: porous_flow

  ! NUMERICS namelist input
  ! Flag for weighting time-integrator
  real(KIND=real_kind),     public, save :: porous_implicitness


end module porous_drag_data
