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

  public :: read_volumetracking_namelist
  public :: vtrack_driver_init, vtrack_update, vtrack_driver_final
  public :: vtrack_enabled
  public :: vtrack_vof_view, vtrack_flux_vol_view

  !! Bundle up all the driver state data as a singleton THIS of private
  !! derived type.  All procedures use/modify this object.
  type :: vtrack_driver_data
    type(unstr_mesh), pointer :: mesh => null()  ! reference only -- do not own
    integer, allocatable :: liq_matid(:), sol_matid(:)
    integer, allocatable :: liq_pri(:)
    type(volume_tracker) :: vt
    real(r8), allocatable  :: vof(:,:) ! volume fractions for all materials
    real(r8), allocatable :: svof(:) ! sum of solid volume fractions
    real(r8), allocatable :: fvof_i(:,:) ! fluid/void volume fractions at start of update
    real(r8), allocatable :: fvof_o(:,:) ! fluid/void volume fractions at end of update
    real(r8), allocatable :: flux_vol(:,:) ! flux volumes
    real(r8), allocatable :: flux_vel(:) ! fluxing velocity - ragged form
    integer :: fluids ! number of fluids (not including void) to advance
    integer :: void ! 1 if simulation includes void else 0
    logical :: active
  end type vtrack_driver_data
  type(vtrack_driver_data), allocatable, target :: this

  !! Input data cached in a private parameter list.
  type(parameter_list), save :: params

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
    p => this%fvof_o
  end function vtrack_vof_view

  function vtrack_flux_vol_view() result(p)
    real(r8), pointer :: p(:,:)
    ASSERT(vtrack_enabled())
    p => this%flux_vol
  end function vtrack_flux_vol_view

  !! Current Truchas design requires that parameter input and object
  !! initialization be separated and occur at distinct execution stages.
  !! The following procedure reads the VOLUMETRACKING namelist and stores
  !! the data in private module data.  The data pertaining to VOF initialization
  !! is still read in using the interfaces/body namelists.

  subroutine read_volumetracking_namelist(lun)

    use string_utilities, only: i_to_c
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist

    integer, intent(in) :: lun

    integer :: ios, location_iter_max, subcycles
    real(r8) :: location_tol, cutoff
    logical :: found, use_brents_method, nested_dissection, active
    character(128) :: iom

    namelist /volumetracking/ use_brents_method, location_iter_max, subcycles, location_tol, &
        cutoff, nested_dissection, active

    !! Locate the MICROSTRUCTURE namelist (first occurrence)
    if (is_IOP) then
      rewind lun
      call seek_to_namelist(lun, 'VOLUMETRACKING', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('Error reading input file: iostat=' // i_to_c(ios))

    call broadcast(found)
    if (.not.found) return  ! the namelist is optional

    call TLS_info('')
    call TLS_info('Reading VOLUMETRACKING namelist ...')

    !! Read the namelist.
    if (is_IOP) then
      active = .true.
      nested_dissection = .true.
      use_brents_method = .true.
      location_iter_max = 20
      location_tol = 1.0e-8_r8
      subcycles = 2
      cutoff = 1.0e-8_r8
      read(lun,nml=volumetracking,iostat=ios,iomsg=iom)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading VOLUMETRACKING namelist: ' // trim(iom))

    !! Broadcast the namelist variables
    call broadcast(active)
    call broadcast(use_brents_method)
    call broadcast(location_tol)
    call broadcast(location_iter_max)
    call broadcast(subcycles)
    call broadcast(cutoff)
    call broadcast(nested_dissection)

    call params%set('use_brents_method', use_brents_method)
    call params%set('location_tol', location_tol)
    call params%set('location_iter_max', location_iter_max)
    call params%set('subcycles', subcycles)
    call params%set('cutoff', cutoff)
    call params%set('nested_dissection', nested_dissection)

    allocate(this)
    this%active = active

    if (this%active) call this%vt%read_params(params)

  end subroutine read_volumetracking_namelist


  subroutine vtrack_driver_init(mesh)

    use material_interop, only: phase_to_material, material_to_phase, void_material_index
    use phase_property_table
    use parameter_module, only: nmat

    type(unstr_mesh), pointer, intent(in) :: mesh
    integer :: i,j,v,k
    real(r8) :: cent(3)

    if (.not.allocated(this)) return

    call TLS_info('')
    call TLS_info('Configuring volume tracking ...')

    this%mesh => mesh
    INSIST(associated(this%mesh))

    this%fluids = 0
    this%void = 0
    if (void_material_index > 0) this%void = 1

    if (ppt_has_property("viscosity")) then
      v = ppt_property_id("viscosity")
    else
      print *, "no viscosity found in table"
      return ! hmmm
    end if

    do i = 1, nmat
      if (i == void_material_index) cycle
      if (ppt_has_phase_property(material_to_phase(i), v)) this%fluids = this%fluids + 1
    end do

    print *, "NMAT: ", nmat
    allocate(this%liq_matid(this%fluids+this%void))
    allocate(this%sol_matid(nmat-(this%fluids+this%void)))
    allocate(this%liq_pri(this%fluids+this%void))
    allocate(this%fvof_i(this%fluids+this%void, mesh%ncell))
    allocate(this%fvof_o(this%fluids+this%void, mesh%ncell))
    allocate(this%flux_vol(this%fluids,size(mesh%cface)))
    allocate(this%flux_vel(size(mesh%cface)))
    allocate(this%vof(nmat, mesh%ncell))
    allocate(this%svof(mesh%ncell))

    j = 1
    k = 1
    do i = 1, nmat
      if (i == void_material_index) then
        this%liq_matid(this%fluids+1) = i
        this%liq_pri(this%fluids+1) = 0
      else if (ppt_has_phase_property(material_to_phase(i), v)) then
        this%liq_matid(j) = i
        j = j + 1
      else
        this%sol_matid(k) = i
        k = k + 1
      end if
    end do

    call start_timer('Vof Initialization')
    call get_vof_from_matl()
    this%fvof_o(:,:) = this%fvof_i(:,:)
    call this%vt%init(mesh, this%fluids+this%void)
    call stop_timer('Vof Initialization')

  end subroutine vtrack_driver_init

  subroutine get_vof_from_matl()

    use matl_module, only: gather_vof
    integer :: i

    associate (n => this%mesh%ncell_onP)

      do i = 1, size(this%liq_matid) + size(this%sol_matid)
        call gather_vof(i, this%vof(i,:n))
      end do

      do i = 1, size(this%liq_matid)
        this%fvof_i(i,:n) = this%vof(this%liq_matid(i),:n)
      end do

      this%svof = 0.0_r8

      do i = 1, size(this%sol_matid)
        this%svof(:n) = this%svof(:n) + this%vof(this%sol_matid(i),:n)
      end do
      call gather_boundary(this%mesh%cell_ip, this%fvof_i)
      ! don't need to communicate svof
    end associate
  end subroutine get_vof_from_matl

  ! vel_fn is the outward oriented face-normal velocity
  subroutine vtrack_update(t, dt, vel_fn)
    use constants_module
    use matl_utilities, only: MATL_SET_VOF
    real(r8), intent(in) :: t, dt
    real(r8), intent(in) :: vel_fn(:)

    real(r8) :: args(0:3), vel(3)
    integer :: i, n, j, ncell, nfc, k, f0, f1, fi

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
        this%fluids, this%void, dt, this%svof)


    call MATL_SET_VOF(this%fvof_o(:,:n))
    call stop_timer('Volumetracking')
  end subroutine vtrack_update

end module vtrack_driver
