!!
!! NS_HT_2D_ENTHALPY_ADVECTOR_TYPE
!!
!! This module defines NS_HT_2D_ENTHALPY_ADVECTOR, the explicit donor-cell
!! enthalpy-advection helper used by coupled two-dimensional
!! Navier--Stokes/thermal transport. It converts material-resolved signed
!! cell-face flux volumes to a conservative on-process cell energy increment.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module ns_ht_2d_enthalpy_advector_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use matl_prop_class
  use bndry_func1_class
  implicit none
  private

  type :: matl_prop_box
    class(matl_prop), allocatable :: prop
  end type

  type, public :: ns_ht_2d_enthalpy_advector
    private
    type(unstr_2d_mesh), pointer :: mesh => null() ! unowned reference
    type(matl_prop_box), allocatable :: enthalpy(:)
    integer, allocatable :: flow_material_ids(:)
    real(r8), allocatable :: temp(:) ! on- and off-process cell workspace
    class(bndry_func1), pointer :: inflow_temperature => null()
  contains
    procedure :: init
    procedure :: get_advected_enthalpy
  end type

contains

  subroutine init(this, mesh, matl_model, flow_material_ids, stat, errmsg, inflow_temperature)

    use material_model_type

    class(ns_ht_2d_enthalpy_advector), intent(out) :: this
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(material_model), intent(in) :: matl_model
    integer, intent(in) :: flow_material_ids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    class(bndry_func1), target, intent(in), optional :: inflow_temperature

    integer :: m

    if (size(flow_material_ids) == 0 .or. any(flow_material_ids < 1) .or. &
        any(flow_material_ids > matl_model%nmatl_real)) then
      stat = 1
      errmsg = 'invalid flow material IDs for enthalpy advection'
      return
    end if
    this%mesh => mesh
    this%flow_material_ids = flow_material_ids
    allocate(this%enthalpy(size(flow_material_ids)), this%temp(mesh%ncell))
    do m = 1, size(flow_material_ids)
      call matl_model%get_matl_prop(flow_material_ids(m), 'enthalpy', this%enthalpy(m)%prop, errmsg)
      if (.not.allocated(this%enthalpy(m)%prop)) then
        stat = 1
        return
      end if
    end do
    if (present(inflow_temperature)) this%inflow_temperature => inflow_temperature

  end subroutine


  !! Return the enthalpy increment due to the given pending-step material flux
  !! volumes. Positive fluxes leave their associated cell. At a boundary
  !! inflow face, an optional thermal inflow boundary condition supplies the
  !! donor temperature; absent that condition, the cell temperature is
  !! retained.
  subroutine get_advected_enthalpy(this, time, cell_temp, flux_volumes, dQ)

    class(ns_ht_2d_enthalpy_advector), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: cell_temp(:)
    real(r8), intent(in) :: flux_volumes(:,:)
    real(r8), intent(out) :: dQ(:)

    integer :: c, i, f, m, neighbor, n
    real(r8) :: donor_temp, flux_enthalpy

    ASSERT(associated(this%mesh))
    ASSERT(size(cell_temp) == this%mesh%ncell_onP)
    ASSERT(size(dQ) == this%mesh%ncell_onP)
    ASSERT(size(flux_volumes,1) >= size(this%flow_material_ids))
    ASSERT(size(flux_volumes,2) == size(this%mesh%cface))

    this%temp(:this%mesh%ncell_onP) = cell_temp
    call this%mesh%cell_imap%gather_offp(this%temp)
    if (associated(this%inflow_temperature)) call this%inflow_temperature%compute(time)

    do c = 1, this%mesh%ncell_onP
      dQ(c) = 0.0_r8
      do i = this%mesh%cstart(c), this%mesh%cstart(c+1)-1
        if (any(flux_volumes(:size(this%flow_material_ids),i) > 0.0_r8)) then
          donor_temp = this%temp(c)
        else
          neighbor = this%mesh%cnhbr(i)
          if (neighbor > 0) then
            donor_temp = this%temp(neighbor)
          else
            donor_temp = this%temp(c)
            if (associated(this%inflow_temperature)) then
              f = this%mesh%cface(i)
              n = findloc(this%inflow_temperature%index, f, dim=1)
              if (n > 0) donor_temp = this%inflow_temperature%value(n)
            end if
          end if
        end if
        do m = 1, size(this%flow_material_ids)
          call this%enthalpy(m)%prop%compute_value([donor_temp], flux_enthalpy)
          dQ(c) = dQ(c) - flux_volumes(m,i)*flux_enthalpy
        end do
      end do
    end do

  end subroutine


end module ns_ht_2d_enthalpy_advector_type
