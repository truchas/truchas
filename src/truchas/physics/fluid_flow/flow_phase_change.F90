!!
!! FLOW_PHASE_CHANGE
!!
!! Neil Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This module implements functionality originally implemented in the
!! heat_transfer module.  It has been separated out here for two reasons:
!! it actually has little to do with heat transfer, and it is independent
!! of a particular heat transfer implementation, working for both old and
!! new heat transfer models.
!!
!! The Navier-Stokes update step (fluid_driver) needs to account for any
!! change in the amount of fluid material that has taken place between the
!! end of the mass advection step and the beginning of the Navier-Stokes
!! step.  Currently, any such change will be due to phase change that may
!! occur (implicitly) from the heat transfer step that is interposed between
!! the those two 'flow' steps.  This module provides the following methods
!! for communicating this info to the Navier-Stokes step:
!!
!!  CALL SET_REFERENCE_FLUID_DENSITY () should be called after the the mass
!!    advection step (advect_mass) and before the heat transfer step to set
!!    the reference value of the fluid density.
!!
!!  CALL SET_SOLIDIFIED_DENSITY () should be called after the heat transfer
!!    step and before the Navier-Stokes step (fluid_flow_driver).  This
!!    computes the density of fluid that has solidified since the call to
!!    SET_REFERENCE_FLUID_DENSITY.  This data is then available in the public
!!    module array SOLIDIFIED_RHO(:) for use by flow.  Note that the values
!!    are always non-negative; if the amount of fluid has increased (melting)
!!    the value is set to 0.
!!
!! If fluid flow is not enabled, both of these methods reduce to no-ops, and
!! so they may be called without worrying about what physics are enabled.
!!
!! The function HAVE_SOLIDIFYING_FLOW() returns true if SOLIDIFIED_RHO(:) has
!! been defined by a call to SET_SOLIDIFIED_DENSITY; any use of the array,
!! which may not be defined, should be protected inside if-blocks guarded by
!! the return value of this function.
!!

#include "f90_assert.fpp"

module flow_phase_change

  use kinds, only: r8
  implicit none
  private
  
  public :: have_solidifying_flow
  public :: set_reference_fluid_density, set_solidified_density
  public :: flow_pc_deallocate

  real(r8), allocatable, save :: fluid_rho(:)
  real(r8), allocatable, save, public :: solidified_rho(:)

contains

  logical function have_solidifying_flow ()
    have_solidifying_flow = allocated(solidified_rho)
  end function have_solidifying_flow

  subroutine set_reference_fluid_density ()
    use fluid_data_module, only: fluid_flow
    use legacy_mesh_api, only: ncells
    if (fluid_flow) then
      if (.not.allocated(fluid_rho)) allocate (fluid_rho(ncells))
      call get_fluid_density (fluid_rho)
    end if
  end subroutine set_reference_fluid_density
  
  subroutine set_solidified_density ()
    use fluid_data_module, only: fluid_flow
    if (fluid_flow) then
      INSIST(allocated(fluid_rho))
      if (.not.allocated(solidified_rho)) allocate(solidified_rho(size(fluid_rho)))
      call get_fluid_density (solidified_rho)
      solidified_rho = max(fluid_rho - solidified_rho, 0.0_r8)
    end if
  end subroutine set_solidified_density
  
  subroutine flow_pc_deallocate ()
    if (allocated(fluid_rho)) deallocate(fluid_rho)
    if (allocated(solidified_rho)) deallocate(solidified_rho)
  end subroutine flow_pc_deallocate
  
  subroutine get_fluid_density (rho)

    use legacy_mesh_api, only: ncells
    use fluid_data_module, only: isImmobile
    use property_module, only: density_material
    use legacy_matl_api, only: Matl, mat_slot

    real(r8), intent(out) :: rho(:)

    integer :: j, s, m
    real(r8) :: sum
    
    ASSERT(size(rho) == ncells)

    do j = 1, ncells 
      sum = 0.0_r8
      do s = 1, mat_slot
        m = matl(s)%cell(j)%id
        if(m <= 0) cycle  ! unused slot
        if(isImmobile(m)) cycle ! not a fluid material
        sum = sum + Matl(s)%cell(j)%vof * density_material(m)
      end do
      rho(j) = sum
    end do

  end subroutine get_fluid_density

end module flow_phase_change
