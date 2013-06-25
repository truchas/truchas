Module BC_MODULE
  !=======================================================================
  ! Purpose:
  !   Master module for boundary condition parameters, structures,
  !   functions, and subroutines. This is the only boundary condition
  !   related module to be accessed or 'used' outside of the bc directory.
  !
  ! Author: Jerry S. Brock, (jsbrock@lanl.gov)
  !         Douglas B. Kothe, (dbk@lanl.gov)
  !         Bryan Lally, (lally@lanl.gov)
  !=======================================================================

  use bc_data_module, only: bc_type, bc_value, bc_variable,                    &
                            Inflow_Material, Inflow_Temperature,               &
                            BC_Temp,                                           &
                            Inflow_Index, Type_Forms, variable_forms,          &
                            BC_Conc, BC_Prs, BC_Mat, BC_Vel,                   &
                            BC_T_DIRICHLET, BC_T_HNEUMANN, BC_T_NEUMANN,       &
                            BC_T_HTC, BC_T_RADIATION, BC_T_HTC_RADIATION,      &
                            BC_T_REFLECTIVE, BC_T_NO_BC, BC, BC_P_NO_BC,       &
                            BC_P_DIRICHLET, BC_P_HNEUMANN, BC_P_NEUMANN,       &
                            BC_P_REFLECTIVE, BC_T_VFRADIATION, BC_T_HTC_GAP

  use bc_flag_module, only: ASSIGN_BC_BITS, SET_DIRICHLET, SET_FREE_SLIP,      &
                            SET_DIRICHLET_VEL, SET_INTERNAL_BC, SET_NEUMANN,   &
                            SET_NEUMANN_VEL, SET_VELOCITY_BC, SET_NO_VEL_BC

  use bc_info_module, only: BCMatch, Boundary, ExternalBoundary,               &
                            InternalBoundary

  use bc_kind_module, only: DIRICHLET, FREE_SLIP, IN_FLOW, INTERNAL_BC,        &
                            NEUMANN, OUT_FLOW, DIRICHLET_VEL, NEUMANN_VEL

  use bc_set_module,  only: InitializeBoundaryConditions, SetBoundaryCondition

  use bc_type_module, only: BC_T, BC_P, Conc, Prs, Vel, BC_STRUCTURE,          &
                            BOUNDARY_CONDITION



  implicit none

  private

  ! bc_data_module
  public :: bc_type, bc_value, bc_variable,                  &
            Inflow_Material, Inflow_Temperature,             &
            BC_Temp,                                         &
            Inflow_Index, Type_Forms, variable_forms,        &
            BC_Conc, BC_Prs, BC_Mat, BC_Vel,                 &
            BC_T_DIRICHLET, BC_T_HNEUMANN, &
            BC_T_NEUMANN, BC_T_HTC, BC_T_RADIATION,          &
            BC_T_VFRADIATION, BC_T_HTC_GAP,                  &
            BC_T_HTC_RADIATION, BC_T_REFLECTIVE, BC_T_NO_BC, &
            BC, BC_P_NO_BC, BC_P_DIRICHLET, BC_P_HNEUMANN,   &
            BC_P_NEUMANN, BC_P_REFLECTIVE


  ! bc_flag_module
  public :: ASSIGN_BC_BITS, SET_DIRICHLET, SET_FREE_SLIP,    &
            SET_DIRICHLET_VEL, SET_INTERNAL_BC, SET_NEUMANN, &
            SET_NEUMANN_VEL, SET_VELOCITY_BC, SET_NO_VEL_BC

  ! bc_info_module
  public :: BCMatch, Boundary, ExternalBoundary,             &
            InternalBoundary

  ! bc_kind_module
  Public :: DIRICHLET, FREE_SLIP, IN_FLOW, INTERNAL_BC,      &
            NEUMANN, OUT_FLOW, DIRICHLET_VEL, NEUMANN_VEL

  ! bc_set_module
  Public :: InitializeBoundaryConditions,                    &
            SetBoundaryCondition

  ! bc_type_module
  Public :: BC_T, BC_P, Conc, Prs, Vel, BC_STRUCTURE,        &
            BOUNDARY_CONDITION

End Module BC_MODULE
