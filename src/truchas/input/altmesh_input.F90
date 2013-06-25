MODULE ALTMESH_INPUT
  !=======================================================================
  ! Purpose(s):
  !
  !   Input parameters for alternative mesh (used for electromagnetics)
  !
  !   Public Interface:
  !
  !     * call READ_ALTMESH_INPUT ()
  !
  !       Reads the AltMesh namelist, stores input data
  !       in the altmesh_data_store.
  !
  ! Author(s): Andrew Kuprat
  !
  !=======================================================================
  use parameter_module,       only: string_len
  use kind_module,            only: real_kind, log_kind
  use truchas_logging_services
  implicit none

  PRIVATE
  PUBLIC :: Read_AltMesh_Input

  !! Magic values used to detect variables not initialized by input
  character,       parameter :: NULL_C = char(0)
  integer,         parameter :: NULL_I = huge(1)
  real(real_kind), parameter :: NULL_R = huge(1.0_real_kind)
 
  ! Allowed Input parameters (listed in namelist, below)
  ! These are all local variables.

  logical(log_kind),         PUBLIC, SAVE :: altmesh_exists = .false.
  character(len=string_len), PUBLIC, SAVE :: altmesh_file = NULL_C
  real(real_kind),           PUBLIC, SAVE :: altmesh_coordinate_scale_factor = 1.d0
  character(len=string_len), PUBLIC, SAVE :: grid_transfer_file = NULL_C
  
  ! AltMesh namelist
  namelist /AltMesh/ altmesh_file, altmesh_coordinate_scale_factor, grid_transfer_file

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS
  SUBROUTINE Read_AltMesh_Input (lun)
    !=======================================================================
    ! Purpose(s):
    !
    !   Read AltMesh namelist.
    !   After namelist read, all variables are broadcast.
    !
    !=======================================================================
    use input_utilities, only: seek_to_namelist
    use parallel_info_module,   only: p_info
    
    integer, intent(in) :: lun
    
    ! Local variables
    logical :: fatal
    integer :: ioerror

    fatal = .FALSE.
    ! Input is done only on the IO PE.  Results of read are broadcast to all PE's
    ! For initialization, need to rewind input file.  
    altmesh_exists = .false.
    READ_IO_PE_ONLY: if (p_info%IOP) then
       ! Look for an ALTMESH namelist, if found read it
       ! Flag error if found namelist but had trouble reading it
       rewind lun
       call seek_to_namelist (lun, 'ALTMESH', found=altmesh_exists)
       if (altmesh_exists) then
          read (lun, NML = ALTMESH, IOSTAT = ioerror)
          fatal = (ioerror /= 0)
       end if
    end if READ_IO_PE_ONLY

    ! Check for fatal errors before broadcasting data to other PE's.
    call TLS_fatal_if_any (fatal, 'error reading ALTMESH namelist')

    ! Broadcast all variables in the ALTMESH namelist.
    call ALTMESH_INPUT_BROADCAST ()

    if (altmesh_exists) then
      ! Check input for obvious errors.
      call check_altmesh_input (fatal)
      if (fatal) call TLS_fatal ('terminating execution due to previous input errors')
    end if

  END SUBROUTINE Read_AltMesh_Input
  
  subroutine check_altmesh_input (fatal)

    use truchas_env, only: input_dir

    logical, intent(out) :: fatal
    
    fatal = .false.
    
    if (altmesh_file == NULL_C) then
      call input_error ('ALTMESH_FILE must be assigned a value')
    end if

    !! If not an absolute path, make it relative to the input directory.
    altmesh_file = adjustl(altmesh_file)
    if (altmesh_file(1:1) /= '/') altmesh_file = trim(input_dir) // trim(altmesh_file)

    !! If specified and not absolute, make the grid-transfer file relative to the input dir.
    if (grid_transfer_file /= NULL_C) then
      grid_transfer_file = adjustl(grid_transfer_file)
      if (grid_transfer_file /= '/') grid_transfer_file = trim(input_dir) // trim(grid_transfer_file)
    end if

  contains
                                                                                
    subroutine input_error (message)
      character(*), intent(in) :: message
      call TLS_error (message)
      fatal = .true.
    end subroutine input_error
                                                                                
  end subroutine check_altmesh_input

  Subroutine AltMesh_Input_Broadcast()
    use pgslib_module, only: pgslib_bcast
    ! Broadcast all the input variables
    
    call pgslib_bcast(altmesh_exists)
    call pgslib_bcast(altmesh_file)
    call pgslib_bcast(altmesh_coordinate_scale_factor)
    call pgslib_bcast(grid_transfer_file)

  end subroutine AltMesh_Input_Broadcast
    
END MODULE ALTMESH_INPUT
    
