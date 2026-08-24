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
  use parameter_list_type
  implicit none
  private

  type :: matl_prop_box
    class(matl_prop), allocatable :: prop
  end type

  type, public :: ns_ht_2d_enthalpy_advector
    private
    type(unstr_2d_mesh), pointer :: mesh => null() ! unowned reference
    type(matl_prop_box), allocatable :: enthalpy(:)
    integer, allocatable :: material_index(:)
    real(r8), allocatable :: temp(:) ! on- and off-process cell workspace
    class(bndry_func1), allocatable :: inflow_temperature
  contains
    procedure :: init
    procedure :: get_advected_enthalpy
  end type

contains

  subroutine init(this, mesh, matl_model, material_index, flow_bc_params, stat, errmsg)

    use material_model_type

    class(ns_ht_2d_enthalpy_advector), intent(out) :: this
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(material_model), intent(in) :: matl_model
    integer, intent(in) :: material_index(:)
    type(parameter_list), target, intent(inout) :: flow_bc_params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: m

    if (size(material_index) == 0 .or. any(material_index < 1) .or. &
        any(material_index > matl_model%nmatl_real)) then
      stat = 1
      errmsg = 'invalid material indices for enthalpy advection'
      return
    end if
    this%mesh => mesh
    this%material_index = material_index
    allocate(this%enthalpy(size(material_index)), this%temp(mesh%ncell))
    do m = 1, size(material_index)
      call matl_model%get_matl_prop(material_index(m), 'enthalpy', this%enthalpy(m)%prop, errmsg)
      if (.not.allocated(this%enthalpy(m)%prop)) then
        stat = 1
        return
      end if
    end do
    call init_inflow_temperature(this, flow_bc_params, stat, errmsg)

  end subroutine


  !! Return the enthalpy increment due to the given pending-step material flux
  !! volumes. Positive fluxes leave their associated cell. At a boundary
  !! inflow face, an optional INFLOW-TEMPERATURE defined in the flow BC list
  !! supplies the donor temperature; absent that value, the cell temperature
  !! is retained, matching mainline's historical fallback.
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
    ASSERT(size(flux_volumes,1) >= size(this%material_index))
    ASSERT(size(flux_volumes,2) == size(this%mesh%cface))

    this%temp(:this%mesh%ncell_onP) = cell_temp
    call this%mesh%cell_imap%gather_offp(this%temp)
    if (allocated(this%inflow_temperature)) call this%inflow_temperature%compute(time)

    do c = 1, this%mesh%ncell_onP
      dQ(c) = 0.0_r8
      do i = this%mesh%cstart(c), this%mesh%cstart(c+1)-1
        if (any(flux_volumes(:size(this%material_index),i) > 0.0_r8)) then
          donor_temp = this%temp(c)
        else
          neighbor = this%mesh%cnhbr(i)
          if (neighbor > 0) then
            donor_temp = this%temp(neighbor)
          else
            donor_temp = this%temp(c)
            if (allocated(this%inflow_temperature)) then
              f = this%mesh%cface(i)
              n = findloc(this%inflow_temperature%index, f, dim=1)
              if (n > 0) donor_temp = this%inflow_temperature%value(n)
            end if
          end if
        end if
        do m = 1, size(this%material_index)
          call this%enthalpy(m)%prop%compute_value([donor_temp], flux_enthalpy)
          dQ(c) = dQ(c) - flux_volumes(m,i)*flux_enthalpy
        end do
      end do
    end do

  end subroutine


  subroutine init_inflow_temperature(this, params, stat, errmsg)

    use bndry_face_func_type
    use scalar_func_class
    use scalar_func_factories, only: alloc_scalar_func

    class(ns_ht_2d_enthalpy_advector), intent(inout) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(bndry_face_func), allocatable :: bff
    class(scalar_func), allocatable :: func
    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    integer, allocatable :: setids(:)

    stat = 0
    piter = parameter_list_iterator(params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      if (plist%is_parameter('inflow-temperature')) then
        call plist%get('face-set-ids', setids, stat, errmsg)
        if (stat /= 0) exit
        call alloc_scalar_func(plist, 'inflow-temperature', func, stat, errmsg)
        if (stat /= 0) exit
        if (.not.allocated(bff)) then
          allocate(bff)
          call bff%init(this%mesh)
        end if
        call bff%add(func, setids, stat, errmsg)
        if (stat /= 0) exit
      end if
      call piter%next()
    end do
    if (stat /= 0) then
      errmsg = 'processing flow BC [' // piter%name() // ']: ' // errmsg
      return
    end if
    if (allocated(bff)) then
      call bff%add_complete()
      call move_alloc(bff, this%inflow_temperature)
    end if

  end subroutine

end module ns_ht_2d_enthalpy_advector_type
