!!
!! USTRUC_DRIVER
!!
!! Driver for the microstructure modeling kernel.
!!
!! Neil N. Carlson <nnc@lanl.gov>
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
!!  CALL READ_MICROSTRUCTURE_NAMELIST (LUN) reads the first instance of the
!!    MICROSTRUCTURE namelist from the file opened on logical unit LUN.  This
!!    is a collective procedure; the file is read on the I/O process and the
!!    data replicated to all other processes.  If an error occurs (I/O error
!!    or invalid data) a message is written and execution is halted.  The
!!    presence of this namelist serves to enable the the microstructure
!!    modeling kernel; it is permissible for there to be no instance of this
!!    namelist.  NB: this should be called ONLY when heat transfer is active.
!!
!!  CALL USTRUC_DRIVER_INIT (T) initializes the driver.  T is the initial time.
!!    This should be called only if heat transfer physics is enabled, and after
!!    its initialization.  If microstructure analysis has not been enabled this
!!    subroutine does nothing, and so may always be called.
!!
!!  CALL USTRUC_UPDATE (T) updates or advances the microstructure model to the
!!    next time level.  The subroutine handles collecting all the necessary
!!    mesh-based state arrays required by the model; only the current time T
!!    needs to be passed.
!!
!!  CALL USTRUC_OUTPUT (SEQ_ID) outputs the computed microstructure analysis
!!    quantities.  SEQ_ID is an HDF5 handle that is available in the
!!    TDO_WRITE_TIMESTEP subroutine that will call this.
!!
!! NOTES
!!
!! The allocatation status of the private module variable THIS serves to
!! indicate whether the microstructure modeling kernel is enabled.  It is
!! allocated by READ_MICROSTRUCTURE_NAMELIST if the microstructure namelist
!! is present in the input file.
!!

#include "f90_assert.fpp"

module ustruc_driver

  use kinds, only: r8
  use unstr_mesh_type
  use ustruc_model_type
  use parameter_list_type
  use truchas_logging_services
  use timing_tree
  implicit none
  private

  public :: read_microstructure_namelist
  public :: ustruc_driver_init, ustruc_update, ustruc_output

  !! Bundle up all the driver state data as a singleton THIS of private
  !! derived type.  All procedures use/modify this object.
  type :: ustruc_driver_data
    type(unstr_mesh), pointer :: mesh => null()  ! reference only -- do not own
    integer, allocatable :: sol_matid(:), liq_matid(:)
    type(ustruc_model) :: model
  end type ustruc_driver_data
  type(ustruc_driver_data), allocatable, save :: this

  !! Input data cached in a private parameter list.
  type(parameter_list), save :: params

contains

  !! Current Truchas design requires that parameter input and object
  !! initialization be separated and occur at distinct execution stages.
  !! The following procedure reads the MICROSTRUCTURE namelist and stores
  !! the data in private module data.

  subroutine read_microstructure_namelist (lun)

    use string_utilities, only: i_to_c
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist, NULL_I, NULL_R, NULL_C
    use material_table, only: MT_MAX_NAME_LEN, mt_has_material

    integer, intent(in) :: lun

    integer :: ios, cell_set_ids(32), symmetry_face_sets(32)
    real(r8) :: vel_max, vel_lo_solid_frac, vel_hi_solid_frac
    real(r8) :: theta1, theta1p, theta2, theta2p
    logical :: found
    integer, allocatable :: setids(:)
    character(MT_MAX_NAME_LEN) :: material
    character(128) :: iom

    namelist /microstructure/ material, cell_set_ids, symmetry_face_sets, &
        vel_max, vel_lo_solid_frac, vel_hi_solid_frac, &
        theta1, theta1p, theta2, theta2p

    !! Locate the MICROSTRUCTURE namelist (first occurrence)
    if (is_IOP) then
      rewind lun
      call seek_to_namelist (lun, 'MICROSTRUCTURE', found, iostat=ios)
    end if
    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('Error reading input file: iostat=' // i_to_c(ios))

    call broadcast (found)
    if (.not.found) return  ! the namelist is optional

    call TLS_info ('')
    call TLS_info ('Reading MICROSTRUCTURE namelist ...')

    !! Read the namelist.
    if (is_IOP) then
      material = NULL_C
      cell_set_ids = NULL_I
      symmetry_face_sets = NULL_I
      vel_max = NULL_R
      vel_lo_solid_frac = NULL_R
      vel_hi_solid_frac = NULL_R
      theta1 = NULL_R
      theta1p = NULL_R
      theta2 = NULL_R
      theta2p = NULL_R
      read(lun,nml=microstructure,iostat=ios,iomsg=iom)
    end if
    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading MICROSTRUCTURE namelist: ' // trim(iom))

    !! Broadcast the namelist variables
    call broadcast (material)
    call broadcast (cell_set_ids)
    call broadcast (symmetry_face_sets)
    call broadcast (vel_max)
    call broadcast (vel_lo_solid_frac)
    call broadcast (vel_hi_solid_frac)
    call broadcast (theta1)
    call broadcast (theta1p)
    call broadcast (theta2)
    call broadcast (theta2p)

    !! Check the values of the namelist variables as best we can before
    !! stuffing them into a parameter list.  The parameters are checked
    !! again when finally used, but only with assertions; we need proper
    !! error handling for this type of input.

    !! Check the MATERIAL value.
    if (material == NULL_C) then
      call TLS_fatal ('no value assigned to MATERIAL')
    else if (mt_has_material(material)) then
      call params%set ('material', trim(material))
    else
      call TLS_fatal ('unknown MATERIAL: "' // trim(material) // '"')
    end if

    !! Check CELL-SET-IDS.
    setids = pack(cell_set_ids, mask=(cell_set_ids /= NULL_I))
    if (size(setids) == 0) then
      call TLS_fatal ('no values assigned to CELL-SET-IDS')
    else
      call params%set ('cell-set-ids', setids)
    end if

    !! Check SYMMETRY-FACE-SETS.
    setids = pack(symmetry_face_sets, mask=(symmetry_face_sets /= NULL_I))
    call params%set ('symmetry-face-sets', setids)

    !! Check VEL-MAX
    if (vel_max == NULL_R) then
      call TLS_fatal ('no value assigned to VEL-MAX')
    else if (vel_max < 0.0_r8) then
      call TLS_fatal ('VEL-MAX must be >= 0')
    else
      call params%set ('vel-max', vel_max)
    end if

    !! Check VEL-LO-SOLID-FRAC and VEL-HI-SOLID-FRAC.
    if (vel_lo_solid_frac == NULL_R) then
      call TLS_fatal ('no value assigned to VEL-LO-SOLID-FRAC')
    else if (vel_lo_solid_frac <= 0.0_r8 .or. vel_lo_solid_frac >= 1.0_r8) then
      call TLS_fatal ('VEL-LO-SOLID-FRAC must be > 0.0 and < 1.0')
    end if
    if (vel_hi_solid_frac == NULL_R) then
      call TLS_fatal ('no value assigned to VEL-HI-SOLID-FRAC')
    else if (vel_hi_solid_frac <= 0.0_r8 .or. vel_hi_solid_frac >= 1.0_r8) then
      call TLS_fatal ('VEL-HI-SOLID-FRAC must be > 0.0 and < 1.0')
    end if
    if (vel_lo_solid_frac <= vel_hi_solid_frac) then
      call params%set ('vel-lo-solid-frac', vel_lo_solid_frac)
      call params%set ('vel-hi-solid-frac', vel_hi_solid_frac)
    else
      call TLS_fatal ('require VEL-LO-SOLID-FRAC <= VEL-HI-SOLID-FRAC')
    end if

    !! Check THETA1 and THETA2.
    if (theta1 == NULL_R) then
      call TLS_fatal ('no value assigned to THETA1')
    else if (theta1 < 0.0_r8 .or. theta1 > 1.0_r8) then
      call TLS_fatal ('THETA1 must be >= 0.0 and <= 1.0')
    end if
    if (theta2 == NULL_R) then
      call TLS_fatal ('no value assigned to THETA2')
    else if (theta2 < 0.0_r8 .or. theta2 > 1.0_r8) then
      call TLS_fatal ('THETA2 must be >= 0.0 and <= 1.0')
    end if
    if (theta1 <= theta2) then
      call params%set ('theta1', theta1)
      call params%set ('theta2', theta2)
    else
      call TLS_fatal ('require THETA1 <= THETA2')
    end if

    !! Check THETA1P.
    if (theta1p == NULL_R) then
      ! do not set a parameter value -- accept its default value
    else if (theta1p < 0.0_r8 .or. theta1p > theta1) then
      call TLS_fatal ('THETA1P must be >= 0.0 and <= THETA1')
    else
      call params%set ('theta1p', theta1p)
    end if

    !! Check THETA2P.
    if (theta2p == NULL_R) then
      ! do not set a parameter value -- accept its default value
    else if (theta2p < theta1 .or. theta2p > theta2) then
      call TLS_fatal ('THETA2P must be >= THETA1 and <= THETA2')
    else
      call params%set ('theta2p', theta2p)
    end if

    !! Enable microstructure modeling by allocating the data object THIS.
    allocate(this)

  end subroutine read_microstructure_namelist

  !! This initializes the driver.  It should only be called if heat transfer is
  !! enabled, and after its initialization (mesh_interop data).

  subroutine ustruc_driver_init (t)

    use mesh_manager, only: unstr_mesh_ptr
    use diffusion_solver_data, only: mesh_name
    use material_table, only: mt_material_id, mt_get_material
    use material_system, only: mat_system, ms_get_phase_id
    use material_interop, only: phase_to_material
    use fluid_data_module, only: isImmobile

    real(r8), intent(in) :: t

    integer :: n
    logical :: valid_mat
    character(:), allocatable :: material
    type(mat_system), pointer :: ms
    integer, pointer :: phase_id(:)
    integer, allocatable :: matid(:)

    if (.not.allocated(this)) return ! microstructure modeling not enabled

    call TLS_info ('Configuring microstructure modeling ...')
    call start_timer ('Microstructure')

    !! We should be using the same mesh as the heat transfer physics;
    !! it should eventually be provided by the MPC.
    this%mesh => unstr_mesh_ptr(mesh_name)
    INSIST(associated(this%mesh))

    !! Form the list of old material ids that are the phases of this material;
    !! they are linearly ordered from low to high temperature phases.
    call params%get ('material', material)
    ms => mt_get_material(mt_material_id(material))
    call ms_get_phase_id (ms, phase_id) ! new phase ids
    matid = phase_to_material(phase_id) ! translate to old material index
    deallocate(phase_id)

    !! Split the material id list between solid and liquid phases.
    do n = 1, size(matid) ! find the lowest temperature liquid phase
      if (.not.isImmobile(matid(n))) exit
    end do
    valid_mat = (n > 1) .and. (n <= size(matid))
    if (valid_mat) valid_mat = all(isImmobile(matid(:n-1))) .and. all(.not.isImmobile(matid(n:)))
    if (valid_mat) then
      this%sol_matid = matid(:n-1)
      this%liq_matid = matid(n:)
    else
      call TLS_fatal ('no compatible liquid-solid transformation found for MATERIAL: "' // trim(material) // '"')
    end if

    call this%model%init (this%mesh, params)
    call ustruc_set_initial_state (t)

    call stop_timer ('Microstructure')

  end subroutine ustruc_driver_init

  subroutine ustruc_set_initial_state (t)
    real(r8), intent(in) :: t
    real(r8), allocatable :: tcell(:), tface(:), liq_vf(:), sol_vf(:)
    call ustruct_update_aux (tcell, tface, liq_vf, sol_vf)
    call this%model%set_state (t, tcell, tface, liq_vf, sol_vf)
  end subroutine ustruc_set_initial_state

  subroutine ustruc_update (t)
    real(r8), intent(in) :: t
    real(r8), allocatable :: tcell(:), tface(:), liq_vf(:), sol_vf(:)
    if (.not.allocated(this)) return ! microstructure modeling not enabled
    call start_timer ('Microstructure')
    call start_timer ('collect input')
    call ustruct_update_aux (tcell, tface, liq_vf, sol_vf)
    call stop_timer ('collect input')
    call this%model%update_state (t, tcell, tface, liq_vf, sol_vf)
    call stop_timer ('Microstructure')
  end subroutine ustruc_update

  !! Output the microstructure analysis data to the HDF output file.
  !! No control over what gets written is provided; we write most everything
  !! that is available.  NB: The choice of things to write must be consistent
  !! with the analysis modules instantiated by USTRUC_COMP_FACTORY.

  subroutine ustruc_output (seq_id)

    use parameter_module, only: ncells
    use mesh_interop, only: pcell_t_to_ds
    use parallel_permutations, only: rearrange
    use truchas_danu_output_tools
    use,intrinsic :: iso_c_binding, only: c_ptr

    type(c_ptr), intent(in) :: seq_id

    real(r8), allocatable :: ds_scalar(:), t_scalar(:), ds_vector(:,:), t_vector(:,:)

    if (.not.allocated(this)) return
    call start_timer ('Microstructure')

    allocate(ds_scalar(this%mesh%ncell_onP),   t_scalar(ncells))
    allocate(ds_vector(3,this%mesh%ncell_onP), t_vector(3,ncells))

    !! Core module: temperature and gradient -- modeled cells only
    call write_scalar_field (data_name='temp',      hdf_name='uStruc-T',     viz_name='T')
    call write_vector_field (data_name='temp-grad', hdf_name='uStruc-gradT', viz_name=['dT/dx','dT/dy','dT/dz'])

    !! Core module: solid fraction, gradient, and rate of change
    call write_scalar_field (data_name='frac',      hdf_name='uStruc-F',     viz_name='F')
    call write_vector_field (data_name='frac-grad', hdf_name='uStruc-gradF', viz_name=['dF/dx','dF/dy','dF/dz'])
    call write_scalar_field (data_name='frac-rate', hdf_name='uStruc-Fdot',  viz_name='dF/dt')

    !! VEL1 analysis module: solidification front velocity and speed
    call write_vector_field (data_name='velocity',  hdf_name='uStruc-veloc', viz_name=['Vx','Vy','Vz'])
    call write_scalar_field (data_name='speed',     hdf_name='uStruc-speed', viz_name='solid-speed')

    !! TIME or GV1 analysis modules: time to solidify
    call write_scalar_field (data_name='solid-time', hdf_name='uStruc-solid-time', viz_name='solid-time')

    !! GV1 analysis module: temp gradient and solidification front speed at onset
    call write_scalar_field (data_name='g', hdf_name='uStruc-G', viz_name='G')
    call write_scalar_field (data_name='v', hdf_name='uStruc-V', viz_name='V')

    !! GV1 analysis module: count of steps in mushy zone
    !call write_scalar_field (data_name='count', hdf_name='uStruc-count', viz_name='count')

    call stop_timer ('Microstructure')

  contains

    subroutine write_scalar_field (data_name, hdf_name, viz_name)
      character(*), intent(in) :: data_name, hdf_name, viz_name
      call this%model%get (data_name, ds_scalar)
      call rearrange (pcell_t_to_ds, t_scalar, ds_scalar)
      call write_seq_cell_field (seq_id, t_scalar, hdf_name, for_viz=.true., viz_name=viz_name)
    end subroutine

    subroutine write_vector_field (data_name, hdf_name, viz_name)
      character(*), intent(in) :: data_name, hdf_name, viz_name(:)
      call this%model%get (data_name, ds_vector)
      call rearrange (pcell_t_to_ds, t_vector, ds_vector)
      call write_seq_cell_field (seq_id, t_vector, hdf_name, for_viz=.true., viz_name=viz_name)
    end subroutine

  end subroutine ustruc_output

!!!!!! AUXILIARY PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! This routine gathers the temperature and volume fraction data that will
  !! be passed to certain USTRUCT_MODEL procedures.  Cell and face temperature
  !! data is obtained directly from the diffusion solver, and the liquid and
  !! solid volume fraction data is assembled from MATL and mapped to the mesh
  !! used by the microstructure modeling component.

  subroutine ustruct_update_aux (tcell, tface, liq_vf, sol_vf)

    use diffusion_solver, only: ds_get_cell_temp, ds_get_face_temp

    real(r8), allocatable, intent(out) :: tcell(:), tface(:), liq_vf(:), sol_vf(:)

    allocate(tcell(this%mesh%ncell_onP), tface(this%mesh%nface_onP))
    call ds_get_cell_temp (tcell)
    call ds_get_face_temp (tface)

    allocate(liq_vf(this%mesh%ncell_onP), sol_vf(this%mesh%ncell_onP))
    call get_vol_frac (this%liq_matid, liq_vf)
    call get_vol_frac (this%sol_matid, sol_vf)

  end subroutine ustruct_update_aux

  !! This routine returns the total volume fraction VF occupied by the phases
  !! corresponding to the old material indices in the list MATID. The returned
  !! cell-centered volume fractions are based on the new mesh.

  subroutine get_vol_frac (matid, vf)

    use matl_module, only: gather_vof
    use parameter_module, only: ncells
    use mesh_interop, only: pcell_ds_to_t
    use parallel_permutations, only: rearrange

    integer, intent(in) :: matid(:)
    real(r8), intent(out) :: vf(:)

    integer :: n
    real(r8), allocatable :: vf1(:), vf2(:)

    ASSERT(size(vf) == this%mesh%ncell_onP)

    allocate(vf1(ncells))
    call gather_vof (matid(1), vf1)
    if (size(matid) > 1) then
      allocate(vf2(ncells))
      do n = 2, size(matid)
        call gather_vof (matid(n), vf2)
        vf1 = vf1 + vf2
      end do
    end if
    call rearrange (pcell_ds_to_t, vf, vf1)

  end subroutine get_vol_frac

end module ustruc_driver
