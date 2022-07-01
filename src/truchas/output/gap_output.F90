!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! The subroutine SET_GAP_ELEMENT_OUTPUT was moved from TBROOK_UTILITY to this
!! new module as a temporary measure to avoid the cyclic module dependency that
!! would arise if it were moved into one of the existing modules that would be
!! a natural home for it.
!!
!! This routine sets the values of various fields on gap cells to some graphics
!! friendly values; e.g. the field value at a neighboring real cell. 
!!

module gap_output

  use kinds, only: r8
  implicit none
  private
  
  public :: set_gap_element_output
  
contains
 
  SUBROUTINE SET_GAP_ELEMENT_OUTPUT () 
    !--------------------------------------------------------------------------- 
    ! Purpose: 
    ! 
    !    Set gap element output data to that of an adjacent cell 
    !---------------------------------------------------------------------------
    use parameter_module,            only: ncomps
    use legacy_mesh_api,             only: ncells, nfc, Mesh, EE_GATHER
    use legacy_mesh_api,             only: GAP_ELEMENT_1, GAP_ELEMENT_3, GAP_ELEMENT_5 
    use zone_module,                 only: Zone 
    use parallel_communication
    use physics_module, only: heat_transport

    integer :: icomp 
    real(r8) :: Tstemp(nfc, ncells), Estemp(nfc, ncells), Pstemp(nfc, ncells), &
                Prtemp(nfc, ncells), Rtemp(nfc, ncells) 
    real(r8), dimension(:), allocatable :: rotation_magnitude,                &
                                          plastic_strain_rate
    real(r8), dimension(:,:), allocatable :: total_strain, plastic_strain,     &
                                             elastic_stress
    !--------------------------------------------------------------------------- 
    if (.not. global_any(Mesh%Cell_Shape >= GAP_ELEMENT_1)) return 
    HEAT_COND: if (heat_transport) then 
       call EE_GATHER(Tstemp, Zone%Temp) 
       call EE_GATHER(Estemp, Zone%Enthalpy) 
       call EE_GATHER(Pstemp, Zone%Rho) 
       where (Mesh%Cell_Shape == GAP_ELEMENT_1) 
          Zone%Temp = Tstemp(1,:) 
          Zone%Enthalpy = Estemp(1,:) 
          Zone%Rho = Pstemp(1,:) 
       end where 
       where (Mesh%Cell_Shape == GAP_ELEMENT_3) 
          Zone%Temp = Tstemp(3,:) 
          Zone%Enthalpy = Estemp(3,:) 
          Zone%Rho = Pstemp(3,:) 
       end where 
       where (Mesh%Cell_Shape == GAP_ELEMENT_5) 
          Zone%Temp = Tstemp(5,:) 
          Zone%Enthalpy = Estemp(5,:) 
          Zone%Rho = Pstemp(5,:) 
       end where 
   end if HEAT_COND 
  END SUBROUTINE SET_GAP_ELEMENT_OUTPUT 
 
end module gap_output
