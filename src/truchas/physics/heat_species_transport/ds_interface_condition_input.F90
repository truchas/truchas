!!
!! DS_INTERFACE_CONDITION_INPUT
!!
!! This module provides procedures for reading and processing the diffusion
!! solver's interface condition namelists.
!!
!! Neil Carlson <nnc@lanl.gov>
!! Markus Berndt <berndt@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! The following routines read multiple instances of the DS_INTERFACE_CONDITION
!! namelist from a file and process the data from selected instances to define
!! IF_DATA-type objects that encapsulate the data associated with various types
!! of interface conditions.  The namelist has the form
!!
!!    &DS_INTERFACE_CONDITION
!!      VARIABLE  = string
!!      CONDITION = string
!!      FACE_SET_IDS = list-of-integers
!!      DATA_CONSTANT = vector-of-reals
!!      DATA_FUNCTION = vector-of-strings
!!    /
!!
!! Namelist selection is based on the value of the VARIABLE and CONDITION
!! strings (case insensitive).  Each type of condition will have 0 or more
!! parameters associated with it.  A heat transfer condition, for example,
!! in which the heat flux is proportional to the temperature difference
!! across the interface, would have the proportionality constant as the only
!! parameter.  Other conditions may have more or less parameters.  The
!! namelist arrays DATA_CONSTANT and DATA_FUNCTION are used to specify this
!! parameter data: for each array index (up to the number of parameters) the
!! data is specified either as a constant via DATA_CONSTANT or as the name of
!! a known function (to the FUNCTION_TABLE module) via DATA_FUNCTION, but not
!! both.  The face set IDs will be matched to link set IDs defined in the
!! distributed mesh, and specify the portion of the interfaces where the
!! condition is to be applied.
!!
!! Note that this module does not codify a particular set of conditions or
!! variables.  The using application code determines the condition names,
!! the number of associated parameters and their interpretation, as well as
!! the variable names through its use of GET_INTERFACE_DATA, and the input
!! file data is expected to be consistent with it.
!!
!!  CALL READ_DS_INTERFACE_CONDITION (LUN) reads all the DS_INTERFACE_CONDITION
!!    namelists from the file opened on logical unit LUN.  The actual reading
!!    of the file occurs only on the I/O process and the data is replicated on
!!    all other processors.  If an error occurs (I/O error or invalid data) a
!!    message is written and execution is halted.  It is permissible for there
!!    to be no instance of the namelist.
!!
!!  CALL GET_INTERFACE_DATA (MESH, VARIABLE, CONDITION, NPAR, IDATA) initializes
!!    the IF_DATA type variable IDATA using data from the namelist instances
!!    whose VARIABLE and CONDITION strings match the value of the corresponding
!!    arguments (case insensitive).  NPAR specifies the number of parameters
!!    associated with the condition, and MESH is a pointer to the distributed
!!    mesh object over which the interface condition is being defined.  If an
!!    error is encountered a message is written and execution is halted.
!!
!! IMPLEMENTATION NOTES
!!
!! The face set IDs specified in the namelist are matched to link set IDs in
!! the distributed mesh object.  The reason for this seeming incongruity is
!! that interfaces (represented as face-to-face links) inherit their link set
!! IDs from the face set IDs of the internal surfaces from which they were
!! created, and it is these face set IDs that are meaningful to the user.
!!
!! This is a two stage input implementation in keeping with Truchas' input
!! design.  All data is read and copied into memory in its raw (or lightly
!! processed) form, and then in a later stage accessed to perform the final
!! initialization of the code.  An alternative, and simpler, implementation
!! would directly read the input file for data as needed for initialization,
!! avoiding the creation of a temporary in-memory copy of the input file.
!!
!! To make this module as flexible as possible, it was decided not to codify
!! an explicit set of known conditions or variables.  The downside is that
!! it is not possible to verify the input as thoroughly as could be done
!! otherwise.  It is possible to specify variable and condition names that
!! the application code will not use without receiving any warning.  To help
!! mitigate this issue, READ_DS_INTERFACE_CONDITION echos the specified
!! condition and variable for each namelist and calls to GET_INTERFACE_DATA
!! will echo the variable and condition being accessed.  The user can use
!! this info to confirm what they supplied and what the application code is
!! using.
!!
!! The private module data is only intended to hold the raw input data
!! temporarily until it can be processed into IF_DATA objects and acquired
!! by the diffusion solver.  However this private data is never deallocated.
!!

#include "f90_assert.fpp"

module ds_interface_condition_input

  use kinds
  use parallel_communication
  use unstr_mesh_type
  use interface_data
  use scalar_func_containers
  use truchas_logging_services
  use string_utilities, only: raise_case, i_to_c
  implicit none
  private

  public :: read_ds_interface_condition, get_interface_data

  integer, parameter :: MAX_NAME_LEN = 31, MAX_FACE_SET_IDS = 32, MAX_COND_PAR = 4

  integer,   parameter :: NULL_I = HUGE(1)
  real(r8),  parameter :: NULL_R = HUGE(1.0_r8)
  character, parameter :: NULL_C = char(0)

  type :: list_node
    integer :: seq = 0
    character(len=MAX_NAME_LEN) :: name
    character(len=MAX_NAME_LEN) :: variable
    character(len=MAX_NAME_LEN) :: condition
    integer, pointer :: face_set_ids(:) => null()
    type(scalar_func_box), allocatable :: farray(:)
    type(list_node), pointer :: next => null()
  end type list_node
  type(list_node), pointer, save :: list => null(), last => null()

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_DS_INTERFACE_CONDITION
 !!

  subroutine read_ds_interface_condition (lun)

    use input_utilities, only: seek_to_namelist
    use scalar_func_factories, only: alloc_const_scalar_func
    use function_namelist, only: lookup_func

    integer, intent(in) :: lun

    logical :: found
    integer :: i, n, npar, stat
    type(scalar_func_box), allocatable :: farray(:)
    character(len=8+MAX_NAME_LEN) :: label
    character(len=127) :: errmsg

    !! Namelist variables
    integer :: face_set_ids(MAX_FACE_SET_IDS)
    real(r8) :: data_constant(MAX_COND_PAR)
    character(len=MAX_NAME_LEN) :: name, variable, condition, data_function(MAX_COND_PAR)
    namelist /ds_interface_condition/ name, variable, condition, face_set_ids, &
                                      data_constant, data_function

    call TLS_info ('')
    call TLS_info ('Reading DS_INTERFACE_CONDITION namelists ...')

    if (is_IOP) rewind(lun)
    n = 0 ! namelist counter
    stat = 0

    do ! until all DS_INTERFACE_CONDITION namelists have been read or an error occurs.

      n = n + 1
      write(label,'(a,i0,a)') '  [', n, ']'

      if (is_IOP) call seek_to_namelist (lun, 'DS_INTERFACE_CONDITION', found, iostat=stat)

      call broadcast (stat)
      if (stat /= 0) then
        errmsg = trim(label) // ' seek error, iostat=' // i_to_c(stat)
        exit
      end if

      call broadcast (found)
      if (.not.found) exit

      !! Read the namelist variables, assigning default values first.
      if (is_IOP) then
        name          = NULL_C
        variable      = NULL_C
        condition     = NULL_C
        face_set_ids  = NULL_I
        data_constant = NULL_R
        data_function = NULL_C
        read(lun,nml=ds_interface_condition,iostat=stat)
      end if

      call broadcast (stat)
      if (stat /= 0) then
        errmsg = trim(label) // ' read error, iostat=' // i_to_c(stat)
        exit
      end if

      !! Replicate the namelist variables on all processes.
      call broadcast (name)
      call broadcast (variable)
      call broadcast (condition)
      call broadcast (face_set_ids)
      call broadcast (data_constant)
      call broadcast (data_function)

      !! Verify that NAME was assigned a unique value.
      if (name == NULL_C) then
        stat = -1
        errmsg = trim(label) // ' error: no value assigned to NAME'
        exit
      else if (name_exists(name)) then
        stat = -1
        errmsg = trim(label) // ' error: another DS_INTERFACE_CONDITION namelist has this NAME: "' &
                             // trim(name) // '"'
        exit
      end if

      !! Use NAME for the diagnostics label now.
      label = '  [' // trim(name) // ']'

      !! Verify VARIABLE was assigned a value.
      if (variable == NULL_C) then
        stat = -1
        errmsg = trim(label) // ' error: no value assigned to VARIABLE'
        exit
      end if

      !! Verify CONDITION was assigned a value.
      if (condition == NULL_C) then
        stat = -1
        errmsg = trim(label) // ' error: no value assigned to CONDITION'
        exit
      end if

      !! Check for a non-empty FACE_SET_IDS.
      if (count(face_set_ids /= NULL_I) == 0) then
        stat = -1
        errmsg = trim(label) // ' error: no values assigned to FACE_SET_IDS'
        exit
      endif

      !! Infer the number of condition parameters.
      do npar = size(data_constant), 1, -1
        if (data_constant(npar) /= NULL_R .or. data_function(npar) /= NULL_C) exit
      end do

      allocate(farray(npar))
      do i = 1, npar
        !! Verify that only one of DATA_CONSTANT and DATA_FUNCTION were specified.
        if (data_constant(i) == NULL_R .eqv. data_function(i) == NULL_C) then
          stat = -1
          errmsg = trim(label) // ' error: exactly one of DATA_CONSTANT(' // i_to_c(i) // &
                        ') and DATA_FUNCTION(' // i_to_c(i) // ') must be assigned a value'
          exit
        end if
        if (data_constant(i) /= NULL_R) then
          call alloc_const_scalar_func (farray(i)%f, data_constant(i))
        else
          call lookup_func (data_function(i), farray(i)%f)
          if (.not.allocated(farray(i)%f)) then
            stat = -1
            errmsg = trim(label) // ' error: unknown function name: ' // trim(data_function(i))
            exit
          end if
        end if
      end do
      if (stat /= 0) exit

      !! Append the data for this namelist to the list of namelist data.
      if (associated(last)) then
        allocate(last%next)
        last => last%next
      else  ! list is empty
        allocate(list)
        last => list
      end if
      last%seq = n
      last%name = name
      last%variable = variable
      last%condition = condition
      call move_alloc (farray, last%farray)
      allocate(last%face_set_ids(count(face_set_ids /= NULL_I)))
      last%face_set_ids = pack(face_set_ids, mask=(face_set_ids /= NULL_I))

      call TLS_info (trim(label) // ' read "' // trim(condition) // '" condition for "' // trim(variable) // '" variable')

    end do

    if (stat /= 0) then
      call TLS_info (trim(errmsg))
      call TLS_fatal ('error reading DS_INTERFACE_CONDITION namelists')
    end if

  end subroutine read_ds_interface_condition

  logical function name_exists (name)
    character(len=*), intent(in) :: name
    type(list_node), pointer :: l
    name_exists = .true.
    l => list
    do while (associated(l))
      if (l%name == name) return
      l => l%next
    end do
    name_exists = .false.
  end function name_exists

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! GET_INTERFACE_DATA
 !!

  subroutine get_interface_data (mesh, variable, condition, npar, idata)

    type(unstr_mesh), intent(in), target :: mesh
    character(len=*), intent(in)  :: variable
    character(len=*), intent(in)  :: condition
    integer,          intent(in)  :: npar
    type(if_data),    intent(out) :: idata

    logical :: error
    integer :: j, stat, stat_array(nPE)
    character(len=127) :: errmsg, errmsg_array(nPE)
    type(list_node), pointer :: l

    call TLS_info ('  Generating "' // trim(condition) // '" interface condition for "' &
                                   // trim(variable)  // '" variable')
    call if_data_prep (idata, mesh, npar)

    l => list
    stat = 0
    errmsg = ''
    do while (associated(l))
      if (raise_case(l%variable) == raise_case(variable)) then
        if (raise_case(l%condition) == raise_case(condition)) then
          write(errmsg,'(4x,a,i0,2a)') 'using DS_INTERFACE_CONDITION[', l%seq, ']: ', trim(l%name)
          call TLS_info (trim(errmsg)) ! not an error message!
          if (size(l%farray) /= npar) then
            stat = -1
            write(errmsg,'(6x,2(a,i0))') 'condition requires ', npar, ' data values but read ', size(l%farray)
            exit
          end if
          call if_data_add (idata, l%farray, l%face_set_ids, stat, errmsg)
          if (stat /= 0) exit
        end if
      end if
      l => l%next
    end do

    if (global_any(stat /= 0)) then
      call collate (stat_array, stat)
      call broadcast (stat_array)
      call collate (errmsg_array, errmsg)
      call broadcast (errmsg_array)
      do j = 1, nPE
        if (stat_array(j) /= 0) then
          call TLS_info ('Error[' // i_to_c(j) // ']: ' // trim(errmsg_array(j)))
        end if
      end do
      call TLS_fatal ('error generating interface condition')
    end if

    call if_data_done (idata)

  end subroutine get_interface_data

end module ds_interface_condition_input
