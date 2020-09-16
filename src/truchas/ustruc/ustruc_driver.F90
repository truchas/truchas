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
  use truchas_timers
  use material_model_driver, only: matl_model
  implicit none
  private

  public :: read_microstructure_namelist
  public :: ustruc_driver_init, ustruc_update, ustruc_output, ustruc_driver_final
  public :: ustruc_read_checkpoint, ustruc_skip_checkpoint
  public :: ustruc_enabled

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

  subroutine ustruc_driver_final
    if (allocated(this)) deallocate(this)
  end subroutine ustruc_driver_final

  logical function ustruc_enabled ()
    ustruc_enabled = allocated(this)
  end function ustruc_enabled

  !! Current Truchas design requires that parameter input and object
  !! initialization be separated and occur at distinct execution stages.
  !! The following procedure reads the MICROSTRUCTURE namelist and stores
  !! the data in private module data.

  subroutine read_microstructure_namelist (lun)

    use string_utilities, only: i_to_c
    use truchas_env, only: input_dir
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist, NULL_I, NULL_R, NULL_C

    integer, intent(in) :: lun

    integer :: ios, cell_set_ids(32)
    real(r8) :: vel_max, vel_lo_solid_frac, vel_hi_solid_frac
    real(r8) :: theta1, theta1p, theta2, theta2p, theta_gv
    logical :: found
    integer, allocatable :: setids(:)
    character(64) :: material
    character(128) :: iom, gv_model_file

    namelist /microstructure/ material, cell_set_ids, &
        vel_max, vel_lo_solid_frac, vel_hi_solid_frac, &
        theta1, theta1p, theta2, theta2p, theta_gv, gv_model_file

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
      vel_max = NULL_R
      vel_lo_solid_frac = NULL_R
      vel_hi_solid_frac = NULL_R
      theta1 = NULL_R
      theta1p = NULL_R
      theta2 = NULL_R
      theta2p = NULL_R
      theta_gv = NULL_R
      gv_model_file = NULL_C
      read(lun,nml=microstructure,iostat=ios,iomsg=iom)
    end if
    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading MICROSTRUCTURE namelist: ' // trim(iom))

    !! Broadcast the namelist variables
    call broadcast (material)
    call broadcast (cell_set_ids)
    call broadcast (vel_max)
    call broadcast (vel_lo_solid_frac)
    call broadcast (vel_hi_solid_frac)
    call broadcast (theta1)
    call broadcast (theta1p)
    call broadcast (theta2)
    call broadcast (theta2p)
    call broadcast (theta_gv)
    call broadcast (gv_model_file)

    !! Check the values of the namelist variables as best we can before
    !! stuffing them into a parameter list.  The parameters are checked
    !! again when finally used, but only with assertions; we need proper
    !! error handling for this type of input.

    !! Check the MATERIAL value.
    if (material == NULL_C) then
      call TLS_fatal ('no value assigned to MATERIAL')
    else if (material == 'VOID') then
      call TLS_fatal ('VOID is an invalid value for MATERIAL')
    else if (matl_model%has_matl(material)) then
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

    !! Check VEL-MAX
    if (vel_max == NULL_R) then
      call TLS_fatal ('no value assigned to VEL_MAX')
    else if (vel_max < 0.0_r8) then
      call TLS_fatal ('VEL_MAX must be >= 0')
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

    if (gv_model_file /= NULL_C) then ! file mode

      !! Check GV_MODEL_FILE.
      if (gv_model_file(1:1) /= '/') gv_model_file = trim(input_dir) // trim(gv_model_file)
      inquire(file=trim(gv_model_file), exist=found)  ! NB: all processes will read
      if (.not.found) call TLS_fatal (' GV_MODEL_FILE not found: "' // trim(gv_model_file) // '"')
      call params%set ('gv-model-file', trim(gv_model_file))

    else  ! basic GV fall back mode

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

      !! Check THETA_GV.
      if (theta_gv == NULL_R) then
        ! do not set a parameter value -- accept its default value
      else if (theta_gv < theta1 .or. theta_gv >= theta2) then
        call TLS_fatal ('THETA_GV must be >= THETA1 and < THETA2')
      else
        call params%set ('theta-gv', theta_gv)
      end if

    end if

    !! Enable microstructure modeling by allocating the data object THIS.
    allocate(this)

  end subroutine read_microstructure_namelist

  !! This initializes the driver.  It should only be called if heat transfer is
  !! enabled, and after its initialization (mesh_interop data).

  subroutine ustruc_driver_init (t)

    use mesh_manager, only: unstr_mesh_ptr
    use diffusion_solver_data, only: mesh_name

    real(r8), intent(in) :: t

    integer :: m, n, p, p1, p2
    logical :: valid_mat
    character(:), allocatable :: name

    if (.not.allocated(this)) return ! microstructure modeling not enabled

    call TLS_info ('')
    call TLS_info ('Configuring microstructure modeling ...')
    call start_timer ('Microstructure')

    !! We should be using the same mesh as the heat transfer physics;
    !! it should eventually be provided by the MPC.
    this%mesh => unstr_mesh_ptr(mesh_name)
    INSIST(associated(this%mesh))

    !! Split the phase id list between solid and liquid phases.
    !! They should be ordered from low to high temperature phases.
    call params%get ('material', name)
    m = matl_model%matl_index(name)
    INSIST(m > 0)
    call matl_model%get_matl_phase_index_range(m, p1, p2)
    do n = p1, p2 ! find the lowest temperature liquid phase
      if (matl_model%is_fluid(n)) exit
    end do
    valid_mat = (n > p1) .and. (n <= p2)
    if (valid_mat) valid_mat = all(matl_model%is_fluid(n:p2))
    if (valid_mat) then
      this%sol_matid = [(p, p=p1,n-1)]
      this%liq_matid = [(p, p=n,p2)]
    else
      call TLS_fatal ('no compatible liquid-solid transformation found for MATERIAL: "' // name // '"')
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
    call start_timer ('analysis')
    call this%model%update_state (t, tcell, tface, liq_vf, sol_vf)
    call stop_timer ('analysis')
    call stop_timer ('Microstructure')
  end subroutine ustruc_update

  !! Output the microstructure analysis data to the HDF output file.
  !! No control over what gets written is provided; we write most everything
  !! that is available.  NB: The choice of things to write must be consistent
  !! with the analysis modules instantiated by USTRUC_COMP_FACTORY.

  subroutine ustruc_output (seq)

    use legacy_mesh_api, only: ncells
    use truchas_danu_output_tools
    use truchas_h5_outfile, only: th5_seq_group

    class(th5_seq_group), intent(in) :: seq

    real(r8), allocatable :: scalar_out(:), vector_out(:,:)

    if (.not.allocated(this)) return
    call start_timer ('Microstructure')

    allocate(scalar_out(ncells), vector_out(3,ncells))

    !! Core module: temperature and gradient -- modeled cells only
    call write_scalar_field (data_name='temp',      hdf_name='uStruc-T',     viz_name='T')
    call write_vector_field (data_name='temp-grad', hdf_name='uStruc-gradT', viz_name=['dT/dx','dT/dy','dT/dz'])
    call write_scalar_field (data_name='temp-rate', hdf_name='uStruc-Tdot',  viz_name='dT/dt')

    !! Core module: solid fraction, gradient, and rate of change
    call write_scalar_field (data_name='frac',      hdf_name='uStruc-Fs',     viz_name='Fs')

    !! VEL1 analysis module: solidification front velocity and speed
    call write_vector_field (data_name='velocity',  hdf_name='uStruc-veloc', viz_name=['Vx','Vy','Vz'])
    call write_scalar_field (data_name='speed',     hdf_name='uStruc-speed', viz_name='solid-speed')

    !! GV0 or GV1 analysis modules: time to solidify
    call write_scalar_field (data_name='solid-time', hdf_name='uStruc-solid-time', viz_name='solid-time')

    !! GV0 analysis module: temp gradient and solidification front speed at onset
    call write_scalar_field (data_name='g', hdf_name='uStruc-G', viz_name='G')
    call write_scalar_field (data_name='v', hdf_name='uStruc-V', viz_name='V')

    !! GV0 analysis module: count of steps in mushy zone
    call write_scalar_field (data_name='count', hdf_name='uStruc-count', viz_name='count')

    !! GV1 analysis module
    call write_scalar_field (data_name='ustruc',  hdf_name='uStruc-gv1-ustruc',  viz_name='gv1-ustruc')
    call write_scalar_field (data_name='lambda1', hdf_name='uStruc-gv1-lambda1', viz_name='gv1-lambda1')
    call write_scalar_field (data_name='lambda2', hdf_name='uStruc-gv1-lambda2', viz_name='gv1-lambda2')

    call stop_timer ('Microstructure')

    !! Additional checkpoint data gets written at the same time.
    call ustruc_write_checkpoint (seq)

  contains

    subroutine write_scalar_field (data_name, hdf_name, viz_name)
      character(*), intent(in) :: data_name, hdf_name, viz_name
      if (this%model%has(data_name)) then
        call this%model%get (data_name, scalar_out)
        scalar_out(this%mesh%ncell_onP+1:) = 0.0_r8 ! gap elements, if any
        call write_seq_cell_field (seq, scalar_out, hdf_name, for_viz=.true., viz_name=viz_name)
      end if
    end subroutine

    subroutine write_vector_field (data_name, hdf_name, viz_name)
      character(*), intent(in) :: data_name, hdf_name, viz_name(:)
      if (this%model%has(data_name)) then
        call this%model%get (data_name, vector_out)
        vector_out(:,this%mesh%ncell_onP+1:) = 0.0_r8 ! gap elements, if any
        call write_seq_cell_field (seq, vector_out, hdf_name, for_viz=.true., viz_name=viz_name)
      end if
    end subroutine

  end subroutine ustruc_output

  !! Output the integrated internal state of the analysis components to the HDF
  !! file needed for restarts.  This is additional internal state that is not set
  !! by USTRUC_SET_INITIAL_STATE, and it corresponds to data recorded during the
  !! course of integration.

  subroutine ustruc_write_checkpoint (seq)

    use,intrinsic :: iso_fortran_env, only: int8
    use truchas_h5_outfile, only: th5_seq_group
    use parallel_communication, only: global_sum
    use string_utilities, only: i_to_c

    class(th5_seq_group), intent(in) :: seq

    integer :: j, n
    integer, allocatable :: cid(:), lmap(:)
    integer(int8), allocatable :: lar(:,:)
    character(:), allocatable :: name

    if (.not.allocated(this)) return
    call start_timer ('Microstructure')

    !! Write the external cell indices that correspond to the state data.
    call this%model%get_map (lmap)
    lmap = this%mesh%xcell(lmap)
    name = 'CP-USTRUC-MAP'
    n = global_sum(size(lmap))
    call seq%write_dist_array(name, n, lmap)

    !! Get the list of analysis components ...
    call this%model%get_comp_list (cid)
    do j = 1, size(cid)
      !! and for each one, fetch its state and write it.
      call this%model%serialize (cid(j), lar)
      if (.not.allocated(lar))  cycle ! no checkpoint data for this analysis component
      if (size(lar,dim=1) == 0) cycle ! no checkpoint data for this analysis component
      n = global_sum(size(lar,dim=2))
      name = 'CP-USTRUC-COMP-' // i_to_c(j)
      call seq%write_dist_array(name, n, lar)
      call seq%write_dataset_attr(name, 'COMP-ID', cid(j))
    end do

    call stop_timer ('Microstructure')

  end subroutine ustruc_write_checkpoint

  !! Reads the microstructure checkpoint data from the restart file, and pushes
  !! it to the microstructure analysis components (deserialize).  Note that this
  !! is internal integrated state that is *additional* to the state that is set
  !! via USTRUC_SET_INITIAL_STATE, and if this procedure is called (optional) it
  !! must be after the call to USTRUC_DRIVER_INIT.

  subroutine ustruc_read_checkpoint (lun)

    use,intrinsic :: iso_fortran_env, only: int8
    use parallel_communication, only: is_IOP, global_sum, global_any, broadcast, collate, distribute
    use sort_utilities, only: heap_sort
    use permutations, only: reorder, invert_perm
    use restart_utilities, only: read_var

    integer, intent(in) :: lun

    integer :: j, n, ios, ncell, ncomp, ncell_tot, cid, nbyte
    integer, allocatable :: lmap(:), gmap(:), map(:), perm(:), perm1(:)
    integer(int8), allocatable :: garray(:,:), larray(:,:)

    call TLS_info ('  Reading microstructure state data from the restart file.')

    !! Get the cell list and translate to external cell numbers.  This needs to
    !! exactly match the cell list read from the restart file, modulo the order.
    call this%model%get_map (lmap)
    ncell = size(lmap)
    do j = 1, size(lmap)
      lmap(j) = this%mesh%xcell(lmap(j))
    end do
    n = global_sum(size(lmap))
    allocate(gmap(merge(n,0,is_IOP)))
    call collate (gmap, lmap)

    !! Read the number of microstructure cells and verify the value.
    call read_var (lun, ncell_tot, 'USTRUC_READ_CHECKPOINT: error reading USTRUC-NCELL record')
    if (ncell_tot /= n) call TLS_fatal ('USTRUC_READ_CHECKPOINT: incorrect number of cells')

    !! Read the cell list; these are external cell numbers.
    allocate(map(merge(ncell_tot,0,is_IOP)))
    if (is_IOP) read(lun,iostat=ios) map
    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('USTRUC_READ_CHECKPOINT: error reading USTRUC-CELL record')

    !! We need to verify that GMAP and MAP specify the same cells and also determine the
    !! pull-back permutation that takes MAP to GMAP.  This will be used to permute the
    !! checkpoint data to match the current cell order.  This is accomplished by generating
    !! the permutations that sort both GMAP and MAP.
    if (is_IOP) then
      allocate(perm(ncell_tot), perm1(ncell_tot))
      call heap_sort (map, perm1)
      call reorder (map, perm1) ! map is now sorted
      call heap_sort (gmap, perm)
      call reorder (gmap, perm) ! gmap is now sorted
    end if
    if (global_any(map /= gmap)) call TLS_fatal ('USTRUC_READ_CHECKPOINT: incorrect cell list')
    deallocate(lmap, gmap, map)
    if (is_IOP) then
      call invert_perm (perm)
      do j = 1, size(perm)
        perm(j) = perm1(perm(j))
      end do
      deallocate(perm1)
    end if

    !! Read the number of component sections to follow.
    call read_var (lun, ncomp, 'USTRUC_READ_CHECKPOINT: error reading USTRUC-NCOMP record')

    !! For each component ...
    do n = 1, ncomp
      !! Read the component ID and the number of data bytes (per cell)
      call read_var (lun, cid, 'USTRUC_READ_CHECKPOINT: error reading USTRUC-COMP-ID record')
      call read_var (lun, nbyte, 'USTRUC_READ_CHECKPOINT: error reading USTRUC-COMP-NBYTE record')
      !! Read the component checkpoint data and reorder it according to PERM
      allocate(garray(nbyte,merge(ncell_tot,0,is_IOP)))
      if (is_IOP) read(lun,iostat=ios) garray
      call broadcast (ios)
      if (ios /= 0) call TLS_fatal ('USTRUC_READ_CHECKPOINT: error reading USTRUC-COMP-DATA record')
      if (is_IOP) call reorder (garray, perm)
      !! Distribute the checkpoint data and push it into the analysis component to consume.
      allocate(larray(nbyte,ncell))
      call distribute (larray, garray)
      call this%model%deserialize (cid, larray)
      deallocate(garray,larray)
    end do

  end subroutine ustruc_read_checkpoint

  !! In the event the restart file contains microstructure checkpoint data but
  !! microstructure modeling is not enabled, this subroutine should be called
  !! to skip over the data so that the file is left correctly positioned.

  subroutine ustruc_skip_checkpoint (lun)
    use restart_utilities, only: skip_records, read_var
    integer, intent(in) :: lun
    integer :: n
    call skip_records (lun, 2, 'USTRUC_SKIP_CHECKPOINT: error skipping the USTRUC data')
    call read_var (lun, n, 'USTRUC_SKIP_CHECKPOINT: error skipping the USTRUC data')
    call skip_records (lun, 3*n,  'USTRUC_SKIP_CHECKPOINT: error skipping the USTRUC data')
  end subroutine ustruc_skip_checkpoint

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
    use legacy_mesh_api, only: ncells

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
    vf = vf1(:this%mesh%ncell_onP)

  end subroutine get_vol_frac

end module ustruc_driver
