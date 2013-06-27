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
    use parameter_module,            only: ncells, ncomps, nfc 
    use solid_mechanics_data,        only: SMech_Cell, Rotation_Magnitude, solid_mechanics 
    use zone_module,                 only: Zone 
    use mesh_module,                 only: Mesh, GAP_ELEMENT_1, GAP_ELEMENT_3, GAP_ELEMENT_5 
    use pgslib_module,               only: PGSLib_GLOBAL_ANY 
    use gs_module,                   only: EE_GATHER 
    use physics_module, only: heat_transport, heat_species_transport

    integer :: icomp 
    real(r8) :: Tstemp(nfc, ncells), Estemp(nfc, ncells), Pstemp(nfc, ncells), &
                Prtemp(nfc, ncells), Rtemp(nfc, ncells) 
    !--------------------------------------------------------------------------- 
    if (.not. PGSLib_GLOBAL_ANY(Mesh%Cell_Shape >= GAP_ELEMENT_1)) return 
    SOLID_MECH: if (solid_mechanics) then 
       do icomp = 1,ncomps 
          call EE_GATHER(Tstemp, SMech_Cell%Total_Strain(icomp,:)) 
          call EE_GATHER(Estemp, SMech_Cell%Elastic_Stress(icomp,:)) 
          call EE_GATHER(Pstemp, SMech_Cell%Plastic_Strain(icomp,:)) 
          where (Mesh%Cell_Shape == GAP_ELEMENT_1) 
             SMech_Cell%Total_Strain(icomp,:) = Tstemp(1,:) 
             SMech_Cell%Elastic_Stress(icomp,:) = Estemp(1,:) 
             SMech_Cell%Plastic_Strain(icomp,:) = Pstemp(1,:) 
          end where 
          where (Mesh%Cell_Shape == GAP_ELEMENT_3) 
             SMech_Cell%Total_Strain(icomp,:) = Tstemp(3,:) 
             SMech_Cell%Elastic_Stress(icomp,:) = Estemp(3,:) 
             SMech_Cell%Plastic_Strain(icomp,:) = Pstemp(3,:) 
          end where 
          where (Mesh%Cell_Shape == GAP_ELEMENT_5) 
             SMech_Cell%Total_Strain(icomp,:) = Tstemp(5,:) 
             SMech_Cell%Elastic_Stress(icomp,:) = Estemp(5,:) 
             SMech_Cell%Plastic_Strain(icomp,:) = Pstemp(5,:) 
          end where 
       end do 
       call EE_GATHER(Prtemp, SMech_Cell%Plastic_Strain_Rate) 
       call EE_GATHER(Rtemp, Rotation_Magnitude) 
       where (Mesh%Cell_Shape == GAP_ELEMENT_1) 
          SMech_Cell%Plastic_Strain_Rate(:) = Prtemp(1,:) 
          Rotation_Magnitude(:) = Rtemp(1,:) 
       end where 
       where (Mesh%Cell_Shape == GAP_ELEMENT_3) 
          SMech_Cell%Plastic_Strain_Rate(:) = Prtemp(3,:) 
          Rotation_Magnitude(:) = Rtemp(3,:) 
       end where 
       where (Mesh%Cell_Shape == GAP_ELEMENT_5) 
          SMech_Cell%Plastic_Strain_Rate(:) = Prtemp(5,:) 
          Rotation_Magnitude(:) = Rtemp(5,:) 
       end where 
    end if SOLID_MECH 
    HEAT_COND: if (heat_transport .or. heat_species_transport) then 
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
