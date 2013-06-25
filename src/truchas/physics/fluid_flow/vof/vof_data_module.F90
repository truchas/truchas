MODULE VOF_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define the variables associated with volume-tracking functions
  !
  !   Public Interface:
  !
  !     * None
  ! 
  ! Contains: No procedures
  !
  ! Author(s): The Telluridians (telluride-info@lanl.gov)
  !
  !=======================================================================
  use kind_module,      only: int_kind, log_kind, real_kind
  use parameter_module, only: string_len

  implicit none

  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  
  ! NUMERICS Namelist variables.

  real(real_kind),            public, save                 :: volume_track_iter_tol
  integer(int_kind),          public, save                 :: volume_track_iter_max,      &
                                                              volume_track_subcycles
  logical(log_kind),          public, save                 :: volume_track_interfaces,    &
                                                              volume_track_brents_method, &
                                                              interface_area
  character(LEN=string_len),  public, save                 :: interface_geometry

  ! Derived variables.

  real(real_kind),            public, save                 :: volume_track_iter_avg, adv_dt
  real(real_kind),            public, save, dimension(4)   :: Eps
  integer(int_kind),          public, save, dimension(5)   :: Cases
  integer(int_kind),          public, save                 :: id_dir, id_face, sweep, iface
  logical(log_kind),          public, save                 :: count_cases
  logical(log_kind), pointer, public, save, dimension(:,:) :: VT_Interface_Mask => NULL()

END MODULE VOF_DATA_MODULE
