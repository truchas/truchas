!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!
!!  The EM_input Module
!!
!!    Robert Ferrell <ferrell@diablotech.com>, original version
!!    Neil N. Carlson <nnc@newmexico.com>
!!
!!  Input parameters for electromagnetic simulation.
!!
!!  The module provides the following procedure:
!!
!!    * call read_EM_input ()
!!
!!      Reads the electromagnetics namelist, checks the validity of the
!!      input, and broadcasts the namelist variables.
!!
!!  NB: Only EM_data_proxy is allowed direct access to the data stored in
!!  this module.  All other direct access is strictly forbidden!  Access
!!  must go through EM_data_proxy.
!!

#include "f90_assert.fpp"

module EM_input

  use kinds, only: r8
  use parameter_module, only: string_len, MAXSV
  use string_utilities, only: i_to_c
  use truchas_logging_services
  implicit none
  private

  public :: read_EM_input

  !! Magic values used to detect variables not initialized by input
  character, parameter :: NULL_C = char(0)
  integer,   parameter :: NULL_I = huge(1)
  real(r8),  parameter :: NULL_R = huge(1.0_r8)

  !! Domain type: 'FULL_CYLINDER', 'HALF_CYLINDER', or 'QUARTER_CYLINDER'
  character(string_len), public, save :: EM_Domain_Type = NULL_C
  character, public, save :: Symmetry_Axis  = NULL_C
  
  !! Container for the INDUCTION_COIL namelist data.
  type, public :: coil_data
    real(r8) :: center(3)
    real(r8) :: radius
    real(r8) :: length
    integer       :: nturns
    real(r8), pointer :: current(:)   => null()
    !real(r8), pointer :: frequency(:) => null()
    !real(r8), pointer :: times(:)     => null()
  end type coil_data

  !! Magnetic source field parameters

  real(r8), pointer, public, save :: src_time(:)   => null()
  real(r8), pointer, public, save :: src_freq(:)   => null()
  real(r8), pointer, public, save :: unif_src(:)   => null()
  type(coil_data), pointer, public, save :: coil_array(:) => null()
  
  !! EM solver control parameters
  integer,  public, save :: Steps_Per_Cycle       = NULL_I
  integer,  public, save :: Maximum_Source_Cycles = NULL_I
  real(r8), public, save :: SS_Stopping_Tolerance = NULL_R
  integer,  public, save :: Maximum_CG_Iterations = NULL_I
  real(r8), public, save :: CG_Stopping_Tolerance = NULL_R
  real(r8), public, save :: Num_Etasq = NULL_R
  real(r8), public, save :: Material_Change_Threshold = NULL_R
  
  !! EM output control parameters
  integer,  public, save :: Output_Level = NULL_I
  logical,  public, save :: Graphics_Output = .false.
  real(r8), public, save :: Probe_Points(3,10) = NULL_R
  integer,  public, save :: num_probes = 0
  
contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ THE ELECTROMAGNETIC NAMELIST, CHECK INPUT, AND BROADCAST
 !!
    
  subroutine read_EM_input (lun)
  
    integer, intent(in) :: lun
    integer :: stat

    call TLS_info ('Reading ELECTROMAGNETICS and INDUCTION_COIL Namelists ...')
    
    !! Read the ELECTROMAGNETICS namelist from the input file.
    call read_electromagnetics (lun, stat)
    if (stat /= 0) call TLS_fatal ('Error reading ELECTROMAGNETICS namelist')
    
    !! Read the INDUCTION_COIL namelists from the input file.
    call read_induction_coil (lun, stat)
    if (stat /= 0) call TLS_fatal ('Error reading INDUCTION_COIL namelist ' // i_to_c(stat))
    
    !! Broadcast the input data to all processors
    call broadcast_EM_input ()

    !! Assign default values and check input for obvious errors (takes place on all PE :-/)
    call check_EM_input (stat)
    if (stat /= 0) call TLS_fatal ('terminating execution due to previous input errors')

  end subroutine read_EM_input
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_ELECTROMAGNETICS
 !!
 !! This procedure reads the ELECTROMAGNETICS namelist.  Only the first
 !! occurrence of the namelist is read; any others are ignored.  It is an
 !! error if the namelist is not found.  A global nonzero status value is
 !! returned in STAT if any error occurs (I/O error or a missing namelist).
 !!
 !! NOTES
 !!
 !! o This is a single-use procedure; it (or rather CHECK_EM_INPUT) relies
 !!   on the default initialization of the module variables being read for
 !!   correct behavior.  To make this a multiple-use procedure, the module
 !!   variables must be explicitly assigned their default values prior to
 !!   reading the namelist.
 !!
  
  subroutine read_electromagnetics (lun, stat)
  
    use input_utilities, only: seek_to_namelist
    use parallel_info_module, only: p_info
    use pgslib_module, only: pgslib_bcast
    
    integer, intent(in)  :: lun
    integer, intent(out) :: stat

    logical :: found
    real(r8) :: Source_Times(MAXSV-1) = NULL_R
    real(r8) :: Source_Frequency(MAXSV) = NULL_R
    real(r8) :: Uniform_Source(MAXSV) = NULL_R
    
    namelist /electromagnetics/ EM_Domain_Type, Symmetry_Axis, &
      Source_Times, Source_Frequency, Uniform_Source, &
      Steps_Per_Cycle, Maximum_Source_Cycles, SS_Stopping_Tolerance, &
      Maximum_CG_Iterations, CG_Stopping_Tolerance, Material_Change_Threshold, &
      Num_Etasq, Output_Level, Graphics_Output, Probe_Points
    
    !! Read the namelist on the IO processor; it is an error if it is missing.
    if (p_info%IOP) then
      stat = 1
      rewind lun
      call seek_to_namelist (lun, 'ELECTROMAGNETICS', found)
      if (found) read (lun, nml=electromagnetics, iostat=stat)
      if (stat == 0) then
        call copy_to_packed_array (Source_Times, src_time)
        call copy_to_packed_array (Source_Frequency, src_freq)
        call copy_to_packed_array (Uniform_Source, unif_src)
      end if
    end if
    
    !! Return a global status value.
    call pgslib_bcast (stat)
    
  end subroutine read_electromagnetics
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_INDUCTION_COIL
 !!
 !! This procedure reads the INDUCTION_COIL namelists.  Each occurrence of the
 !! namelist will generate a new instance of a coil data structure.  The
 !! namelist need not appear at all.  A global zero status value is returned
 !! if no error occurred; otherwise a global nonzero status value is returned
 !! that corresponds to the instance (first, second, etc) of the namelist in
 !! which the error occurred.
 !!
  
  subroutine read_induction_coil (lun, stat)
  
    use input_utilities, only: seek_to_namelist
    use parallel_info_module, only: p_info
    use pgslib_module, only: pgslib_bcast
    
    integer, intent(in)  :: lun
    integer, intent(out) :: stat

    logical :: found
    integer :: nturns
    real(r8) :: center(3), radius, length, current(MAXSV)
    
    namelist /induction_coil/ center, radius, length, nturns, current
    
    type :: list_node
      type(coil_data) :: coil
      type(list_node), pointer :: next => null()
    end type list_node
    type(list_node), pointer :: list, old_list
    
    integer :: ncoil, n
    
    !! Read every instance of the namelist (if any) on the IO processor.
    if (p_info%IOP) then
      stat = 0
      ncoil = 0
      list => null()
      rewind lun
      do
        call seek_to_namelist (lun, 'INDUCTION_COIL', found)
        if (.not.found) exit
        ncoil = ncoil + 1
        !! Set default values for the namelist variables.
        center = NULL_R
        radius = NULL_R
        length = NULL_R
        nturns = NULL_I
        current = NULL_R
        !! Read the namelist.
        read(lun, nml=induction_coil, iostat=stat)
        if (stat /= 0) exit
        !! Prepend a new instance of the coil data structure to the list.
        old_list => list
        allocate(list)
        list%next => old_list
        !! Copy the namelist variables into the coil data structure.
        list%coil%center = center
        list%coil%radius = radius
        list%coil%length = length
        list%coil%nturns = nturns
        !! Copy the significant part of the current vector.
        call copy_to_packed_array (current, list%coil%current)
      end do
      if (stat == 0) then
        !! Copy the coil data into the output array and deallocate the list.
        allocate(coil_array(ncoil))
        do n = ncoil, 1, -1
          coil_array(n) = list%coil
          old_list => list%next
          deallocate(list)
          list => old_list
        end do
      else
        !! Deallocate the accumulated list.
        do while (associated(list))
          old_list => list%next
          deallocate(list)
          list => old_list
        end do
        stat = ncoil  ! The offending instance of the namelist.
      end if
    end if
    
    !! Return a global status value
    call pgslib_bcast (stat)

  end subroutine read_induction_coil
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! COPY_TO_PACKED_ARRAY
 !!
 !! This auxillary procedure allocates the array pointer ARRAY with size equal
 !! to the number of input values in the array SOURCE, and copies those values
 !! (in order) to ARRAY.
 !!
    
    subroutine copy_to_packed_array (source, array)
      real(r8), intent(in) :: source(:)
      real(r8), pointer :: array(:)
      allocate(array(count(source /= NULL_R)))
      array = pack(source, mask=(source /= NULL_R))
    end subroutine copy_to_packed_array
    
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ASSIGN DEFAULT VALUES AND CHECK INPUT FOR OBVIOUS ERRORS
 !!
  
  subroutine check_EM_input (stat)
  
    use string_utilities, only: raise_case, i_to_c
    use mesh_manager, only: enable_mesh
    
    integer, intent(out) :: stat
  
    integer :: n
    logical :: exists
    
    stat = 0
    
    !! Ensure string input is in standard form.
    EM_Domain_Type = raise_case(trim(adjustl(EM_Domain_Type)))
    Symmetry_Axis  = raise_case(trim(adjustl(Symmetry_Axis)))
    
    select case (EM_Domain_Type)
      case ('FULL_CYLINDER')
      case ('HALF_CYLINDER')
      case ('QUARTER_CYLINDER')
      case ('CYLINDER')
      case ('FRUSTUM')
      case ('VERIFICATION1')
      
      case (NULL_C)
        call input_error ('EM_Domain_Type must be assigned a value')
        
      case default
        call input_error ('"'//trim(EM_Domain_Type)//'" is not a valid value for EM_Domain_Type')
    end select
    
    select case (Symmetry_Axis)
      case ('X', 'Y', 'Z')
      
      case (NULL_C)
        call input_info ('Using default value "Z" for Symmetry_Axis')
        Symmetry_Axis = 'Z'
        
      case default
        call input_error ('"'//trim(Symmetry_Axis)//'" is not a valid value for Symmetry_Axis')
    end select
    
    if (size(src_time) > 1) then
      n = size(src_time)
      if (any(src_time(2:n) <= src_time(:n-1))) then
        call input_error ('Source_Times values must be strictly increasing')
      end if
    end if
    
    if (size(src_freq) == 0) then
      call input_error ('Source_Frequency must be assigned a value')
    else if (any(src_freq <= 0.0_r8)) then
      call input_error ('Source_Frequency values must be > 0.0')
    else if (size(src_freq) /= size(src_time) + 1) then
      call input_error ('Wrong number of values provided for Source_Frequency')
    end if
    
    if (size(unif_src) > 0) then
      if (size(unif_src) /= size(src_freq)) then
        call input_error ('Wrong number of values provided for Uniform_Source')
      end if
    end if
          
    CHECK_INDUCTION_COILS: do n = 1, size(coil_array)
    
      if (any(coil_array(n)%center == NULL_R)) then
        if (all(coil_array(n)%center == NULL_R)) then
          call input_info ('Using default value (0,0,0) for Center of INDUCTION_COIL ' // i_to_c(n))
          coil_array(n)%center = 0.0_r8
        else
          call input_error ('Center of INDUCTION_COIL ' // i_to_c(n) // ' requires 3 values')
        end if
      end if
      
      if (coil_array(n)%nturns == NULL_I) then
        call input_error ('NTurns must be assigned a value for INDUCTION_COIL ' // i_to_c(n))
      else if (coil_array(n)%nturns <= 0) then
        call input_error ('NTurns must be > 0 for INDUCTION_COIL ' // i_to_c(n))
      end if
      
      if (coil_array(n)%nturns == 1) then
        if (coil_array(n)%length /= NULL_R) then
          call input_info ('Ignoring Length value for INDUCTION_COIL ' // i_to_c(n))
        end if
        coil_array(n)%length = 0.0_r8
      end if
      
      if (coil_array(n)%length == NULL_R) then
        call input_error ('Length must be assigned a value for INDUCTION_COIL ' // i_to_c(n))
      else if (coil_array(n)%length < 0.0_r8) then
        call input_error ('Length must be >= 0.0 for INDUCTION_COIL ' // i_to_c(n))
      end if

      if (coil_array(n)%radius == NULL_R) then
        call input_error ('Radius must be assigned a value for INDUCTION_COIL ' // i_to_c(n))
      else if (coil_array(n)%radius <= 0.0_r8) then
        call input_error ('Radius must be > 0.0 for INDUCTION_COIL ' // i_to_c(n))
      end if

      if (size(coil_array(n)%current) /= size(src_freq)) then
        call input_error ('Inconsistent number of Current values for INDUCTION_COIL ' // i_to_c(n))
      end if
        
    end do CHECK_INDUCTION_COILS
    
    if (Steps_Per_Cycle == NULL_I) then
      Steps_Per_Cycle = 20
      call input_info ('Using default value 20 for Steps_Per_Cycle')
    else if (Steps_Per_Cycle < 1) then
      call input_error ('Steps_Per_Cycle must be > 0')
    else if (Steps_Per_Cycle < 10) then
      call input_warn ('For decent accuracy, Steps_Per_Cycle should be at least 10')
    end if
    
    if (Maximum_Source_Cycles == NULL_I) then
      Maximum_Source_Cycles = 10
      call input_info ('Using default value 10 for Maximum_Source_Cycles')
    else if (Maximum_Source_Cycles < 1) then
      call input_error ('Maximum_Source_Cycles must be > 0')
    end if
    
    if (SS_Stopping_Tolerance == NULL_R) then
      SS_Stopping_Tolerance = 1.0e-2_r8
      call input_info ('Using default value 0.01 for SS_Stopping_Tolerance')
    else if (SS_Stopping_Tolerance <= 0.0_r8) then
      call input_error ('SS_Stopping_Tolerance must be > 0.0')
    else if (SS_Stopping_Tolerance > 0.1_r8) then
      call input_warn ('SS_Stopping_Tolerance is very loose; consider decreasing the value')
    end if
    
    if (Maximum_CG_Iterations == NULL_I) then
      call input_info ('Using default value 500 for Maximum_CG_Iterations')
      Maximum_CG_Iterations = 500
    else if (Maximum_CG_Iterations < 1) then
      call input_error ('Maximum_CG_Iterations must be > 0')
    end if
    
    if (CG_Stopping_Tolerance == NULL_R) then
      CG_Stopping_Tolerance = 1.0e-5_r8
      call input_info ('Using default value 1.0e-5 for CG_Stopping_Tolerance')
    else if (CG_Stopping_Tolerance <= 0.0_r8 .or. CG_Stopping_Tolerance >= 0.1_r8) then
      call input_error ('CG_Stopping_Tolerance must be > 0.0 and < 0.1')
    else if (CG_Stopping_Tolerance > 1.0e-4_r8) then
      call input_warn ('CG_Stopping_Tolerance is very loose and may lead to the build up of errors.')
    else if (CG_Stopping_Tolerance < 1.e-3_r8 * epsilon(1.0_r8)) then
      call input_warn ('CG_Stopping_Tolerance is too tight; CG iterations are unlikely to converge.')
    end if
    
    if (Material_Change_Threshold == NULL_R) then
      Material_Change_Threshold = 0.3_r8
      call input_info ('Using default value 0.3 for Material_Change_Threshold')
    else if (Material_Change_Threshold <= 0.0_r8) then
      call input_error ('Material_Change_Threshold must be > 0.0')
    end if
    
    if (Num_Etasq == NULL_R) then
      Num_Etasq = 0.0_r8
    else if (Num_Etasq < 0.0_r8) then
      call input_error ('Num_Etasq must be >= 0.0')
    end if
    
    if (Output_Level == NULL_I) then
      Output_Level = 1
    else if (Output_Level < 1 .or. Output_Level > 4) then
      call input_error ('Output_Level must be >= 1 and <= 4')
    end if
    
    num_probes = 0
    do n = 1, size(Probe_Points,dim=2)
      if (any(Probe_Points(:,n) == NULL_R)) exit
      num_probes = n
    end do
    
    !! Logical variable Graphics_Output is default initialized.
    
    call enable_mesh ('alt', exists)
    if (.not.exists) call input_error ('ALTMESH namelist was not specified.')
    
  contains
  
    subroutine input_error (message)
      use truchas_logging_services
      character(*), intent(in) :: message
      call TLS_error (message)
      stat = 1
    end subroutine input_error
    
    subroutine input_warn (message)
      use truchas_logging_services
      character(*), intent(in) :: message
      call TLS_warn (message)
    end subroutine input_warn
    
    subroutine input_info (message)
      use truchas_logging_services
      character(*), intent(in) :: message
      call TLS_info (message)
    end subroutine input_info
    
  end subroutine check_EM_input
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! BROADCAST THE ELECTROMAGNETICS NAMELIST VARIABLES
 !!

  subroutine broadcast_EM_input ()
  
    use pgslib_module, only: pgslib_bcast
    use parallel_info_module, only: p_info
    
    integer :: n
    
    call pgslib_bcast (EM_Domain_Type)
    call pgslib_bcast (Symmetry_Axis)
    
    call pgslib_bcast_pointer (src_time)
    call pgslib_bcast_pointer (src_freq)
    call pgslib_bcast_pointer (unif_src)
    
    ASSERT( p_info%IOP .eqv. associated(coil_array) )
    if (associated(coil_array)) n = size(coil_array)
    call pgslib_bcast (n)
    if (.not.associated(coil_array)) allocate(coil_array(n))
    do n = 1, size(coil_array)
      call pgslib_bcast (coil_array(n)%center)
      call pgslib_bcast (coil_array(n)%radius)
      call pgslib_bcast (coil_array(n)%length)
      call pgslib_bcast (coil_array(n)%nturns)
      call pgslib_bcast_pointer (coil_array(n)%current)
    end do
    
    call pgslib_bcast (Steps_Per_Cycle)
    call pgslib_bcast (Maximum_Source_Cycles)
    call pgslib_bcast (SS_Stopping_Tolerance)
    call pgslib_bcast (Maximum_CG_Iterations)
    call pgslib_bcast (CG_Stopping_Tolerance)
    call pgslib_bcast (Material_Change_Threshold)
    
    call pgslib_bcast (Output_Level)
    call pgslib_bcast (Graphics_Output)
    call pgslib_bcast (Probe_Points)
    call pgslib_bcast (num_probes)
    
    call pgslib_bcast (Num_Etasq)
    
  contains
  
    subroutine pgslib_bcast_pointer (ptr)
      real(r8), pointer :: ptr(:)
      integer :: n
      ASSERT( p_info%IOP .eqv. associated(ptr) )
      if (associated(ptr)) n = size(ptr)
      call pgslib_bcast (n)
      if (.not.associated(ptr)) allocate(ptr(n))
      call pgslib_bcast (ptr)
    end subroutine pgslib_bcast_pointer
    
  end subroutine broadcast_EM_input

end module EM_input
