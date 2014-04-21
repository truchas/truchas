!!
!! DS_SOURCE_INPUT
!!
!! This module provides procedures for reading the DS_SOURCE namelist and
!! for processing its data into a source mesh function object used by the
!! diffusion solver.
!!
!! Neil Carlson <nnc@lanl.gov>
!! 4 Apr 2009, revised 5 Jun 2010
!!
!! PROGRAMMING INTERFACE
!!
!! The following routines read multiple instances of the DS_SOURCE namelist
!! from a file and process the data from selected instances to define SOURCE_MF
!! objects that define an external source function over the mesh.  The namelist
!! has the form
!!
!!    &DS_SOURCE
!!      EQUATION  = string
!!      CELL_SET_IDS = list-of-integers
!!      SOURCE_CONSTANT = real
!!      SOURCE_FUNCTION = string
!!    /
!!
!! Namelist selection is based on the value of the EQUATION string (case
!! insensitive).  The namelist arrays SOURCE_CONSTANT and SOURCE_FUNCTION
!! are used to specify the source: either a constant via SOURCE_CONSTANT
!! or the name of a known function (to the FUNCTION_TABLE module) via
!! SOURCE_FUNCTION, but not both.  The specified cell set IDs will be matched
!! to cell set IDs defined in the distributed mesh, and specify that part
!! of the domain where the source applies.  Different instances of the
!! namelist for a given EQUATION value must apply to disjoint parts of
!! the domain.
!!
!! Note that this module does not codify a particular set of equations.  The
!! using application code determines the equation names through its particular
!! use of DEFINE_EXTERNAL_SOURCE, and the input file data is expected to be
!! consistent with it.
!!
!! The following routines should be called once each, in the order presented.
!!
!!  CALL READ_DS_SOURCE (LUN) reads all the DS_SOURCE namelists from the file
!!    opened on logical unit LUN.  The actual reading of the file occurs only
!!    on the I/O process and the data is replicated on all other processes.
!!    If an error occurs (I/O error or invalid data) a message is written and
!!    execution is halted.  It is permissible for there to be no instance of
!!    the namelist.
!!
!!  CALL DEFINE_EXTERNAL_SOURCE (MESH, EQUATION, Q) initializes the SOURCE_MF
!!    type argument Q using data from the namelist instances whose EQUATION
!!    string match the value of the argument (case insensitive).  MESH is a
!!    pointer to the distributed mesh object over which the source is being
!!    defined.  If an error is encountered a message is written and execution
!!    is halted.  A 0-valued source is assigned to parts of the domain not
!!    referenced by one of the namelist instances.
!!
!! IMPLEMENTATION NOTES
!!
!! This is a two stage input implementation in keeping with Truchas' input
!! design.  All data is read and copied into memory in its raw (or lightly
!! processed) form, and then in a later stage accessed to perform the final
!! initialization of the code.  An alternative, and simpler, implementation
!! would directly read the input file for data as needed for initialization,
!! avoiding the creation of a temporary in-memory copy of the input file.
!!
!! To make this module as flexible as possible, it was decided not to codify
!! an explicit set of known equations. The downside is that it is not possible
!! to verify the input as thoroughly as could be done otherwise.  It is
!! possible to specify equation names that the application code will not use
!! without receiving any warning.  To help mitigate this issue, READ_DS_SOURCE
!! echos the specified equation for each namelist and calls to
!! DEFINE_EXTERNAL_SOURCE will echo the equation being accessed.  The user can
!! use this info to confirm what they supplied and what the application code
!! is using.
!!
!! The private module data is only intended to hold the raw input data
!! temporarily until it can be processed into SOURCE_MF objects and acquired
!! by the diffusion solver.  However this private data is never deallocated.
!!

#include "f90_assert.fpp"

module ds_source_input

  use kinds
  use scalar_func_class
  use source_mesh_function
  use distributed_mesh
  use truchas_logging_services
  use parallel_communication
  use string_utilities, only: raise_case, i_to_c
  implicit none
  private

  public :: read_ds_source, define_external_source

  integer, parameter :: MAX_NAME_LEN = 31, MAX_CELL_SET_IDS = 32

  integer,   parameter :: NULL_I = HUGE(1)
  real(r8),  parameter :: NULL_R = HUGE(1.0_r8)
  character, parameter :: NULL_C = char(0)

  !! Private module variables defined by READ_DS_SOURCE.
  type :: list_node
    integer :: seq
    character(len=MAX_NAME_LEN) :: equation
    integer, pointer :: cell_set_ids(:) => null()
    class(scalar_func), allocatable :: srcf
    type(list_node), pointer :: next => null()
  end type list_node
  type(list_node), pointer, save :: list => null(), last => null()

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_DS_SOURCE
 !!

  subroutine read_ds_source (lun)

    use input_utilities, only: seek_to_namelist
    use scalar_func_factories, only: alloc_const_scalar_func
    use function_namelist, only: lookup_func

    integer, intent(in) :: lun

    logical :: found
    integer :: n, stat
    class(scalar_func), allocatable :: srcf
    character(len=7) :: label
    character(len=127) :: errmsg

    !! Namelist variables
    integer  :: cell_set_ids(MAX_CELL_SET_IDS)
    real(r8) :: source_constant
    character(len=MAX_NAME_LEN) :: equation, source_function
    namelist /ds_source/ equation, cell_set_ids, source_constant, source_function

    call TLS_info ('')
    call TLS_info ('Reading DS_SOURCE namelists ...')

    if (is_IOP) rewind(lun)
    n = 0 ! namelist counter
    stat = 0

    do ! until all DS_SOURCE namelists have been read or an error occurs.

      n = n + 1
      write(label,'(a,i0,a)') '  [', n, ']'

      if (is_IOP) call seek_to_namelist (lun, 'DS_SOURCE', found, iostat=stat)

      call broadcast (stat)
      if (stat /= 0) then
        errmsg = trim(label) // ' seek error, iostat=' // i_to_c(stat)
        exit
      end if

      call broadcast (found)
      if (.not.found) exit

      !! Read the namelist variables, assigning default values first.
      if (is_IOP) then
        equation      = NULL_C
        cell_set_ids  = NULL_I
        source_constant = NULL_R
        source_function = NULL_C
        read(lun,nml=ds_source,iostat=stat)
      end if

      call broadcast (stat)
      if (stat /= 0) then
        errmsg = trim(label) // ' read error, iostat=' // i_to_c(stat)
        exit
      end if

      !! Replicate the namelist variables on all processes.
      call broadcast (equation)
      call broadcast (cell_set_ids)
      call broadcast (source_constant)
      call broadcast (source_function)

      !! Verify EQUATION was assigned a value.
      if (equation == NULL_C) then
        stat = -1
        errmsg = trim(label) // ' error: no value assigned to EQUATION'
        exit
      end if

      !! Check for a non-empty CELL_SET_IDS.
      if (count(cell_set_ids /= NULL_I) == 0) then
        stat = -1
        errmsg = trim(label) // 'error: no values assigned to CELL_SET_IDS'
        exit
      endif

      !! Verify that only one of SOURCE_CONSTANT and SOURCE_FUNCTION were specified.
      if (source_constant == NULL_R .eqv. source_function == NULL_C) then
        stat = -1
        errmsg = trim(label) // ' error: exactly one of SOURCE_CONSTANT and SOURCE_FUNCTION must be assigned a value'
        exit
      end if

      !! Create or get the source function.
      if (source_constant /= NULL_R) then
        call alloc_const_scalar_func (srcf, source_constant)
      else
        call lookup_func (source_function, srcf)
        if (.not.allocated(srcf)) then
          stat = -1
          errmsg = trim(label) // ' error: unknown function name: ' // trim(source_function)
          exit
        end if
      end if

      !! Append the data for this namelist to the list of namelist data.
      if (associated(last)) then
        allocate(last%next)
        last => last%next
      else  ! list is empty
        allocate(list)
        last => list
      end if
      last%seq = n
      last%equation = equation
      call move_alloc (srcf, last%srcf)
      allocate(last%cell_set_ids(count(cell_set_ids /= NULL_I)))
      last%cell_set_ids = pack(cell_set_ids, mask=(cell_set_ids /= NULL_I))

      call TLS_info (trim(label) // ' read source for "' // trim(equation) // '" equation')

    end do

    if (stat /= 0) then
      call TLS_info (trim(errmsg))
      call TLS_fatal ('error reading DS_SOURCE namelists')
    else if (n == 1) then ! if no errors, N is one more than the number of namelists found
      call TLS_info ('  No DS_SOURCE namelists found.')
    end if

  end subroutine read_ds_source

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DEFINE_EXTERNAL_SOURCE
 !!

  subroutine define_external_source (mesh, equation, q)

    type(dist_mesh),  intent(in), target :: mesh
    character(len=*), intent(in)  :: equation
    type(source_mf),  intent(out) :: q

    integer :: stat
    character(len=31) :: label
    character(len=127) :: errmsg
    type(list_node), pointer :: l

    call TLS_info ('  Generating external source for "' // trim(equation) // '" equation')
    call smf_prep (q, mesh)

    l => list
    do while (associated(l))
      if (raise_case(l%equation) == raise_case(equation)) then
        write(label,'(a,i0,a)') 'DS_SOURCE[', l%seq, ']'
        call smf_add_function (q, l%srcf, l%cell_set_ids, stat, errmsg)
        call TLS_fatal_if_any (stat /= 0, trim(label) // ' error: ' // trim(errmsg))
      end if
      l => l%next
    end do

    call smf_set_default (q, 0.0_r8)  ! 0-source for all other cells
    call smf_done (q, stat, errmsg)
    call TLS_fatal_if_any (stat /= 0, 'Error defining external source: ' // trim(errmsg))

  end subroutine define_external_source

end module ds_source_input
