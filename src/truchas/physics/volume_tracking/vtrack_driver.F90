!!
!! VOLUME_TRACKING_DRIVER
!!
!! This code is meant to be deleted and subsumed by high-level cycle driver.
!! At the moment, this driver serves as a mediator between the volume tracking
!! code with the desired api and the rest of truchas
!!
!! Peter Brady <ptb@lanl.gov>
!! 2017
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROGRAMMING INTERFACE
!!
!!  This module provides procedures for driving the microstructure modeling
!!  kernel.  It serves as an adapter between top-level Truchas drivers (input,
!!  initialization, cycle, and output) and the microstructure modeling object.
!!  The existing top-level drivers manage no data; they merely provide for
!!  orchestrating the sequence of procedural steps.  Thus a major role of
!!  these subroutines is to assemble the data from various Truchas components
!!  that must be passed the the microstructure modeling object.
!!
!!  CALL READ_VOLUMETRACKING_NAMELIST (LUN) reads the first instance of the
!!    VOLUMETRACKING namelist from the file opened on logical unit LUN.  This
!!    is a collective procedure; the file is read on the I/O process and the
!!    data replicated to all other processes.  If an error occurs (I/O error
!!    or invalid data) a message is written and execution is halted.  The
!!    presence of this namelist serves to enable the the volumetracker; it
!!    is permissible for there to be no instance of this namelist.
!!    NB: this should be called ONLY when flow is active.
!!
!!  CALL VTRACK_DRIVER_INIT (T) initializes the driver.  T is the initial time.
!!    This should be called only if heat transfer physics is enabled, and after
!!    its initialization.  If microstructure analysis has not been enabled this
!!    subroutine does nothing, and so may always be called.
!!
!!  CALL VTRACK_UPDATE (T) updates or advances the microstructure model to the
!!    next time level.  The subroutine handles collecting all the necessary
!!    mesh-based state arrays required by the model; only the current time T
!!    needs to be passed.
!!
!! NOTES
!!
!! The allocatation status of the private module variable THIS serves to
!! indicate whether the microstructure modeling kernel is enabled.  It is
!! allocated by READ_MICROSTRUCTURE_NAMELIST if the microstructure namelist
!! is present in the input file.
!!
#include "f90_assert.fpp"

module vtrack_driver

  use kinds, only: r8
  use unstr_mesh_type
  use volume_tracker_type
  use parameter_list_type
  use truchas_logging_services
  use truchas_timers
  use index_partitioning
  implicit none
  private

  public :: vtrack_driver_init, vtrack_update, vtrack_driver_final
  public :: vtrack_enabled
  public :: vtrack_vof_view, vtrack_flux_vol_view, vtrack_liq_matid_view

  !! Bundle up all the driver state data as a singleton THIS of private
  !! derived type.  All procedures use/modify this object.
  type :: vtrack_driver_data
    type(unstr_mesh), pointer :: mesh => null()  ! reference only -- do not own
    integer, allocatable :: liq_matid(:), sol_matid(:)
    type(volume_tracker) :: vt
    real(r8), allocatable  :: vof(:,:) ! volume fractions for all materials
    real(r8), allocatable :: fvof_i(:,:) ! fluid/void volume fractions at start of update
    real(r8), allocatable :: fvof_o(:,:) ! fluid/void volume fractions at end of update
    real(r8), allocatable :: flux_vol(:,:) ! flux volumes
    real(r8), allocatable :: flux_vel(:) ! fluxing velocity - ragged form
    integer :: fluids ! number of fluids (not including void) to advance
    integer :: void ! 1 if simulation includes void else 0
    integer :: solid ! 1 if simulation includes solid else 0
    logical :: active
  end type vtrack_driver_data
  type(vtrack_driver_data), allocatable, target :: this

contains

  subroutine vtrack_driver_final
    if (allocated(this)) deallocate(this)
  end subroutine vtrack_driver_final

  logical function vtrack_enabled()
    vtrack_enabled = allocated(this)
  end function vtrack_enabled

  function vtrack_vof_view() result(p)
    real(r8), pointer :: p(:,:)
    ASSERT(vtrack_enabled())
    p => this%fvof_o(:this%fluids+this%void,:)
  end function vtrack_vof_view

  function vtrack_flux_vol_view() result(p)
    real(r8), pointer :: p(:,:)
    ASSERT(vtrack_enabled())
    p => this%flux_vol
  end function vtrack_flux_vol_view

  ! material ids corresponding to flux_vol columns
  function vtrack_liq_matid_view() result(p)
    integer, pointer :: p(:)
    p => this%liq_matid(:this%fluids)
  end function


  subroutine vtrack_driver_init(params)

    use mesh_manager, only: unstr_mesh_ptr
    use material_interop, only: void_material_index
    use property_data_module, only: isImmobile
    use parameter_module, only: nmat

    type(parameter_list), intent(inout) :: params

    integer :: i,j,k

    allocate(this)

    call params%get('track_interfaces', this%active)

    call TLS_info('')
    call TLS_info('Configuring volume tracking ...')

    this%mesh => unstr_mesh_ptr('MAIN')
    INSIST(associated(this%mesh))

    this%solid = merge(1, 0, any(isImmobile(:nmat)))
    this%void = merge(1, 0, void_material_index > 0)
    this%fluids = count(.not.isImmobile(:nmat)) - this%void

    if (this%fluids == 0) then
      print *, 'no non-void fluids'
      return ! hmmm
    end if

    allocate(this%liq_matid(this%fluids+this%void))
    allocate(this%sol_matid(nmat-(this%fluids+this%void)))
    allocate(this%fvof_i(this%fluids+this%void+this%solid, this%mesh%ncell))
    allocate(this%fvof_o(this%fluids+this%void+this%solid, this%mesh%ncell))
    allocate(this%flux_vol(this%fluids+this%void,size(this%mesh%cface)))
    allocate(this%flux_vel(size(this%mesh%cface)))
    allocate(this%vof(nmat, this%mesh%ncell_onP))

    j = 1
    k = 1
    do i = 1, nmat
      if (i == void_material_index) then
        this%liq_matid(this%fluids+1) = i
      else if (.not.isImmobile(i)) then
        this%liq_matid(j) = i
        j = j + 1
      else
        this%sol_matid(k) = i
        k = k + 1
      end if
    end do

    call start_timer('Vof Initialization')
    call get_vof_from_matl()
    this%fvof_o = this%fvof_i
    call this%vt%init(this%mesh, this%fluids, this%fluids+this%void, &
        this%fluids+this%void+this%solid, this%liq_matid, params)
    call stop_timer('Vof Initialization')

  end subroutine vtrack_driver_init

  subroutine get_vof_from_matl()

    use matl_utilities, only: matl_get_vof

    integer :: i, n

    n = this%mesh%ncell_onP
    call matl_get_vof(this%vof)
    do i = 1, size(this%liq_matid)
      this%fvof_i(i,:n) = this%vof(this%liq_matid(i),:)
    end do
    call gather_boundary(this%mesh%cell_ip, this%fvof_i)

    if (this%solid > 0) then
      do i = 1, n
        this%fvof_i(this%fluids+this%void+this%solid,i) = sum(this%vof(this%sol_matid,i))
      end do
    end if
    ! don't need to communicate svof

  end subroutine get_vof_from_matl

  subroutine put_vof_into_matl()
    use matl_utilities, only: matl_set_vof
    integer :: i
    do i = 1, size(this%liq_matid)
      this%vof(this%liq_matid(i),:) = this%fvof_o(i,:this%mesh%ncell_onP)
    end do
    call matl_set_vof(this%vof)
  end subroutine put_vof_into_matl

  ! vel_fn is the outward oriented face-normal velocity
  subroutine vtrack_update(t, dt, vel_fn, initial)
    use constants_module
    real(r8), intent(in) :: t, dt
    real(r8), intent(in) :: vel_fn(:)
    logical, intent(in), optional :: initial

    real(r8) :: args(0:3), vel(3)
    integer :: i, j, k, f0, f1, fi

    if (.not.allocated(this)) return
    if (.not.this%active) then
      ! compute flux volumes for transport
      associate (m => this%mesh)
        do i = 1, m%ncell
          f0 = m%xcface(i)
          f1 = m%xcface(i+1)-1
          do j = f0, f1
            k = m%cface(j)
            if (btest(m%cfpar(i),pos=1+j-f0)) then ! normal points inward
              this%flux_vel(j) = -vel_fn(k)
            else
              this%flux_vel(j) = vel_fn(k)
            end if
            this%flux_vol(1,j) = this%flux_vel(j)*m%area(k)*dt
          end do
        end do
      end associate
      return
    end if

    if (this%fluids == 0) return

    call start_timer('Volumetracking')

    ! copy face velocities into cell-oriented array
    do i = 1,this%mesh%ncell
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1
      do j = f0, f1
        k = this%mesh%cface(j)

        if (btest(this%mesh%cfpar(i),pos=1+j-f0)) then ! normal points inward
          this%flux_vel(j) = -vel_fn(k)
        else
          this%flux_vel(j) = vel_fn(k)
        end if

      end do
    end do

    call get_vof_from_matl()

    call this%vt%flux_volumes(this%flux_vel, this%fvof_i, this%fvof_o, this%flux_vol, &
        this%fluids, this%void, dt)

    ! update the matl structure if this isn't the initial pass
    if (present(initial)) then
      if (.not.initial) &
          call put_vof_into_matl()
    else
      call put_vof_into_matl()
    end if

    call stop_timer('Volumetracking')
  end subroutine vtrack_update

end module vtrack_driver
