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
  use kinds, only: r8
  use parameter_module, only: string_len
  implicit none
  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  
  ! NUMERICS Namelist variables.

  real(r8), public, save :: volume_track_iter_tol
  integer,  public, save :: volume_track_iter_max, volume_track_subcycles
  logical,  public, save :: volume_track_interfaces, volume_track_brents_method, interface_area
  character(string_len),  public, save :: interface_geometry

  ! Derived variables.

  real(r8), public, save :: volume_track_iter_avg, adv_dt
  real(r8), public, save, dimension(4) :: Eps
  integer,  public, save, dimension(5) :: Cases
  integer,  public, save :: id_dir, id_face, sweep, iface
  logical,  public, save :: count_cases
  logical, pointer, public, save, dimension(:,:) :: VT_Interface_Mask => NULL()

END MODULE VOF_DATA_MODULE
