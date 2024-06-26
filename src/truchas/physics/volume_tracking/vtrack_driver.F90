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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
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

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use volume_tracker_class
  use parameter_list_type
  use truchas_logging_services
  use truchas_timers
  use material_model_driver, only: matl_model
  implicit none
  private

  public :: vtrack_driver_init, vtrack_update, vtrack_driver_final
  public :: vtrack_enabled
  public :: vtrack_vof_view, vtrack_flux_vol_view, vtrack_liq_matid_view
  public :: get_vof_from_matl
  public :: vtrack_set_inflow_bc, vtrack_set_inflow_material
  public :: put_orig_vof_into_matl

  !! Bundle up all the driver state data as a singleton THIS of private
  !! derived type.  All procedures use/modify this object.
  type :: vtrack_driver_data
    type(unstr_mesh), pointer :: mesh => null()  ! reference only -- do not own
    integer, allocatable :: liq_matid(:), sol_matid(:)
    class(volume_tracker), allocatable :: vt
    real(r8), allocatable  :: vof(:,:) ! volume fractions for all materials
    real(r8), allocatable :: fvof_i(:,:) ! fluid/void volume fractions at start of update
    real(r8), allocatable :: fvof_o(:,:) ! fluid/void volume fractions at end of update
    real(r8), allocatable :: flux_vol(:,:) ! flux volumes
    real(r8), allocatable :: flux_vel(:) ! fluxing velocity - ragged form
    integer :: fluids ! number of fluids (not including void) to advance
    integer :: void ! 1 if simulation includes void else 0
    integer :: solid ! 1 if simulation includes solid else 0
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
    p => this%fvof_o
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
    use geometric_volume_tracker_type
    use simple_volume_tracker_type

    type(parameter_list), intent(inout) :: params

    integer :: i,j,k
    logical :: track_interfaces

    allocate(this)

    call TLS_info('')
    call TLS_info('Initializing volume tracker ...')

    this%mesh => unstr_mesh_ptr('MAIN')
    INSIST(associated(this%mesh))

    this%solid = 0; this%fluids = 0
    do i = 1, matl_model%nphase_real
      if (matl_model%is_fluid(i)) then
        this%fluids = this%fluids + 1
      else
        this%solid = 1
      end if
    end do
    this%void = merge(1, 0, matl_model%have_void)

    if (this%fluids == 0) then
      print *, 'no non-void fluids'
      return ! hmmm
    end if

    allocate(this%liq_matid(this%fluids+this%void))
    allocate(this%sol_matid(matl_model%nphase-(this%fluids+this%void)))
    allocate(this%fvof_i(this%fluids+this%void+this%solid, this%mesh%ncell))
    allocate(this%fvof_o(this%fluids+this%void+this%solid, this%mesh%ncell))
    allocate(this%flux_vol(this%fluids+this%void,size(this%mesh%cface)))
    allocate(this%flux_vel(size(this%mesh%cface)))
    allocate(this%vof(matl_model%nphase, this%mesh%ncell_onP))

    j = 1
    k = 1
    do i = 1, matl_model%nphase_real
      if (matl_model%is_fluid(i)) then
        this%liq_matid(j) = i
        j = j + 1
      else
        this%sol_matid(k) = i
        k = k + 1
      end if
    end do
    if (matl_model%have_void) this%liq_matid(this%fluids+1) = matl_model%nphase

    call params%get('track_interfaces', track_interfaces)
    if (track_interfaces) then
      allocate(geometric_volume_tracker :: this%vt)
    else
      allocate(simple_volume_tracker :: this%vt)
    end if

    call start_timer('Vof Initialization')
    call get_vof_from_matl(this%fvof_i)
    this%fvof_o = this%fvof_i
    call this%vt%init(this%mesh, this%fluids, this%fluids+this%void, &
        this%fluids+this%void+this%solid, this%liq_matid, params)
    call stop_timer('Vof Initialization')

  end subroutine vtrack_driver_init

  subroutine get_vof_from_matl(vof)

    use legacy_matl_api, only: matl_get_vof

    real(r8), intent(out) :: vof(:,:)

    integer :: i, n

    n = this%mesh%ncell_onP
    call matl_get_vof(this%vof)
    do i = 1, size(this%liq_matid)
      vof(i,:n) = this%vof(this%liq_matid(i),:)
    end do

    if (this%solid > 0) then
      do i = 1, n
        vof(this%fluids+this%void+this%solid,i) = sum(this%vof(this%sol_matid,i))
      end do
    end if

    call this%mesh%cell_imap%gather_offp(vof)

  end subroutine get_vof_from_matl

  subroutine put_vof_into_matl()
    use legacy_matl_api, only: define_matl
    integer :: i
    do i = 1, size(this%liq_matid)
      this%vof(this%liq_matid(i),:) = this%fvof_o(i,:this%mesh%ncell_onP)
    end do
    call define_matl(this%vof)
  end subroutine put_vof_into_matl

  subroutine put_orig_vof_into_matl
    use legacy_matl_api, only: define_matl
    integer :: i
    do i = 1, size(this%liq_matid)
      this%vof(this%liq_matid(i),:) = this%fvof_i(i,:this%mesh%ncell_onP)
    end do
    call define_matl(this%vof)
  end subroutine

  ! vel_fn is the outward oriented face-normal velocity
  subroutine vtrack_update(t, dt, vel_fn, initial)

    real(r8), intent(in) :: t, dt
    real(r8), intent(in) :: vel_fn(:)
    logical, intent(in), optional :: initial

    integer :: i, j, k, f0, f1

    if (.not.allocated(this)) return
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

    call get_vof_from_matl(this%fvof_i)

    call this%vt%flux_volumes(this%flux_vel, this%fvof_i, this%fvof_o, this%flux_vol, &
        this%fluids, this%void, dt, t)

    ! update the matl structure if this isn't the initial pass
    if (present(initial)) then
      if (.not.initial) &
          call put_vof_into_matl()
    else
      call put_vof_into_matl()
    end if

    call stop_timer('Volumetracking')

  end subroutine vtrack_update

  !! Sets the inflow material for the given list of boundary faces. MATID is
  !! the Truchas material number which must be a fluid. MATID can also be the
  !! special value 0 which indicates that materials should be fluxed into the
  !! adjacent cell in proportion to the material fractions currently present
  !! in the cell. The initial default for all boundary faces is the latter.

  subroutine vtrack_set_inflow_material(matid, faces)
    use vector_func_class
    class(vector_func), intent(in) :: matid
    integer, intent(in) :: faces(:)
    call this%vt%set_inflow_material(matid, faces)
  end subroutine vtrack_set_inflow_material

  !! Sets the inflow material boundary conditions specified by the given
  !! parameter list. The structure of the list is a list of sublists. Only
  !! those that define an "inflow-material" parameter are significant. Its
  !! value and the corresponding value of the required "face-set-ids"
  !! parameter are used to set the inflow material on a portion of the
  !! boundary. Other sublists are ignored.

  subroutine vtrack_set_inflow_bc(params, stat, errmsg)

    use bndry_face_group_builder_type
    use vector_func_containers
    use vector_func_factories

    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    type(bndry_face_group_builder) :: builder
    integer, allocatable :: setids(:), xgroup(:), index(:)
    type(vector_func_box), allocatable :: mlist(:)
    character(:), allocatable :: name
    integer :: j, n, matid, vofid
    real(r8) :: vof(this%fluids+this%void)

    call builder%init(this%mesh, omit_offp=.true.)

    piter = parameter_list_iterator(params, sublists_only=.true.)
    n = piter%count()
    allocate(mlist(n))
    n = 0
    do while (.not.piter%at_end())
      plist => piter%sublist()
      if (plist%is_parameter('inflow-material') &
          .or. plist%is_parameter('inflow-material-func')) then

        call plist%get('face-set-ids', setids, stat, errmsg)
        if (stat /= 0) return
        call builder%add_face_group(setids, stat, errmsg)
        if (stat /= 0) return

        if (plist%is_parameter('inflow-material')) then
          call plist%get('inflow-material', name, stat, errmsg)
          if (stat /= 0) return
          n = n + 1
          matid = matl_model%phase_index(name)
          do vofid = 1, size(this%liq_matid)
            if (this%liq_matid(vofid) == matid) exit
          end do
          ASSERT(vofid <= size(this%liq_matid))
          vof = 0
          vof(vofid) = 1
          call alloc_const_vector_func(mlist(n)%f, vof)
        else
          n = n + 1
          call alloc_vector_func(plist, 'inflow-material-func', mlist(n)%f, stat, errmsg)
          if (stat /= 0) return
        end if
      end if
      call piter%next
    end do

    call builder%get_face_groups(n, xgroup, index)

    do j = 1, n
      call vtrack_set_inflow_material(mlist(j)%f, index(xgroup(j):xgroup(j+1)-1))
    end do

  end subroutine vtrack_set_inflow_bc

end module vtrack_driver
