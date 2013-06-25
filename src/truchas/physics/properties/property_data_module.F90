MODULE PROPERTY_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define material property variables and control parameters.
  !
  ! Contains: None
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !            Anand V. Reddy, LANL T-3 (jsbrock@lanl.gov)
  !
  !=======================================================================
  use kind_module,      only: int_kind, log_kind, real_kind
  use parameter_module, only: maxcon, maxmat, ndim

  implicit none

  ! Public Module
  public

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Material Namelist Variables
  integer(KIND = int_kind), save :: background_material

  character(LEN = 80), dimension(0:maxmat), save ::          &
                       Material_Name, Material_Feature, Viscoplastic_Model
  
  real(KIND = real_kind), dimension(0:maxmat), save :: &
                       density, Sound_Speed, Void_Temperature

  integer(KIND = int_kind), dimension(0:maxmat), save :: Priority

  ! Solid Mechanics input - not temperature dependent                     
  real(KIND = real_kind), dimension(0:maxmat), save :: &
                       MTS_k, MTS_mu_0, MTS_sig_a,     &
                       MTS_d, MTS_temp_0, MTS_b,       &
                       MTS_edot_0i, MTS_g_0i, MTS_q_i, &
                       MTS_p_i, MTS_sig_i, Pwr_Law_A,  &
                       Pwr_Law_N, Pwr_Law_Q, Pwr_Law_R

  real(KIND = real_kind), dimension(ndim,0:maxmat), save :: Permeability_Constant

  ! Interface-tracking advection priority control parameters
  integer(KIND = int_kind), dimension(0:maxmat), save :: Matpri

END MODULE PROPERTY_DATA_MODULE
