!!
!!  The DATA_EXPLORER Module
!!
!!  Neil N. Carlson <nnc@newmexico.com>
!!  Initial version prior to 9/96 with many subsequent revisions
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  This module provides an API to easily create native format DX files
!!  used by the OpenDX visualization tool (previously IBM Visualization
!!  Data Explorer).
!!
!!  22 Aug 2003, NNC
!!  Unless the fpp variable NAG_4_2 is defined, binary append mode is
!!  not available (a message is printed, and execution stops).  This
!!  mode requires system access (FLUSH and SYSTEM) and I haven't yet
!!  worked out how to do this cleanly across platforms.  Note that the
!!  default value for APPEND remains .true..
!!
!!  30 Aug 2003, NNC
!!  Because of the broken Truchas build procedure, I've had to simply
!!  comment out the code protected by the NAG_4_2 ifdefs.  The problem:
!!  The automatic dependency generator looks at the original code
!!  instead of the output of the preprocessor as it should.  As a
!!  result the build process believes it needs a NAG system module.

module data_explorer

  use string_utilities, i2c => i_to_c
  implicit none
  private
  
  public :: dx_open, dx_close, dx_export_array, dx_export_connections, dx_export_field
  public :: dx_new_series, dx_append_to_series, dx_export_series, dx_delete_series
  public :: dx_export_regular_grid
  public :: dx_new_multigrid, dx_append_to_multigrid, dx_export_multigrid, dx_delete_multigrid
  
  !! Maximum character-type length for file names, object names, etc.
  integer, parameter, private :: CHAR_LEN = 128
  
  type, public :: dx_file
    private
    integer :: hunit = -1
    integer :: dunit = -1
    logical :: embed_data = .false.
    logical :: append_data = .true.
    integer :: object = 0 ! next object number
    integer :: offset = 0 ! current (byte) offset in the data file
    character(len=CHAR_LEN) :: path       ! directory path
    character(len=CHAR_LEN) :: file       ! DX header file name
    character(len=CHAR_LEN) :: dfile      ! DX data file name
  end type dx_file
  
  type, public :: dx_object
    private
    integer :: number
    character(len=CHAR_LEN) :: name
    character(len=CHAR_LEN), pointer :: file => null()   ! DX header file name
  end type dx_object
  
  type, public :: dx_series
    private
    integer :: n = 0          ! number of objects in series
    integer :: max = 0        ! size of allocated arrays
    integer :: n_init = 100   ! initial size to allocate
    integer :: n_incr = 50    ! incremental size when reallocating
    real,            dimension(:), pointer :: time   => null()
    type(dx_object), dimension(:), pointer :: object => null()
  end type dx_series
  
  type, public :: dx_multigrid
    private
    integer :: n = 0          ! number of members in group
    integer :: max = 0        ! size of allocated arrays
    integer :: n_init = 10    ! initial size to allocate
    integer :: n_incr = 5     ! incremental size when reallocating
    type(dx_object), dimension(:), pointer :: object => null()
  end type dx_multigrid
  
  interface dx_export_field
    module procedure export_field, export_real_field_data, export_real_field_data_1, export_int_field_data
  end interface
  private :: export_field, export_real_field_data, export_real_field_data_1, export_int_field_data
  
  interface dx_export_array
    module procedure export_integer_array_data_0, export_integer_array_data_1
    module procedure export_integer_array_0_0, export_integer_array_1_0
    module procedure export_integer_array_0_1, export_integer_array_1_1
    module procedure export_real_array_data_0, export_real_array_data_1
    module procedure export_real_array_0_0, export_real_array_1_0
    module procedure export_real_array_0_1, export_real_array_1_1
  end interface
  private :: export_integer_array_data_0, export_integer_array_data_1
  private :: export_integer_array_0_0, export_integer_array_1_0
  private :: export_integer_array_0_1, export_integer_array_1_1
  private :: export_real_array_data_0, export_real_array_data_1
  private :: export_real_array_0_0, export_real_array_1_0
  private :: export_real_array_0_1, export_real_array_1_1
  
  interface dx_is_defined
    module procedure is_defined_file, is_defined_object
  end interface
  private :: is_defined_file, is_defined_object
  
  interface write_attributes
    module procedure write_string_attribute_one, write_string_attribute_many
  end interface
  private :: write_attributes
  private :: write_string_attribute_one, write_string_attribute_many
  
  
  character(len=20), parameter, private :: DATA_MODE = 'data mode lsb binary'
  integer, parameter, private :: SIZE_OF_REC_HEAD = 4, SIZE_OF_REC_TAIL = 4
  integer, parameter, private :: SIZE_OF_REC_SEP = SIZE_OF_REC_HEAD + SIZE_OF_REC_TAIL
  integer, parameter, private :: SIZE_OF_INTEGER = 4, SIZE_OF_REAL = 4
  
contains

  elemental logical function is_defined_object (this)
    type(dx_object), intent(in) :: this
    is_defined_object = associated(this % file)
  end function is_defined_object

  elemental logical function is_defined_file (this)
    type(dx_file), intent(in) :: this
    is_defined_file = (this % object > 0)
  end function is_defined_file
  
  subroutine dx_open (dxf, file, append, embed)

    use,intrinsic :: iso_fortran_env, only: stdout => output_unit

    type(dx_file), intent(out) :: dxf
    character(len=*), intent(in), optional :: file
    logical, intent(in), optional :: append, embed
    
    character(len=CHAR_LEN) :: base
    integer :: n
    
    if (present(embed))  dxf%embed_data  = embed
    if (present(append)) dxf%append_data = append
    
#ifndef NAG_4_2
    if (dxf%append_data) then
      print *, 'Sorry, binary append mode is not available on this platform;'
      print *, 'use the the optional argument APPEND=.FALSE. to DX_OPEN.'
      stop
    end if
#endif
    
    if (present(file)) then ! We'll be writing to a named file
    
      if (len_trim(file) > CHAR_LEN) then
        print *, 'PANIC: file name is too long: "' // trim(file) // '"'
        stop
      end if
      
      n = scan(file, '/', back=.true.)
      if (n > 0) then
        dxf%path = file(:n)
        base = file(n+1:)
      else
        dxf%path = './'
        base = file
      end if

      n = scan(file,'.',back=.true.)
      if (n > 0) then
        if (file(n:) == '.dx') then
          base = file(1:n-1)
        end if
      end if

      if (len_trim(base) + 4 > CHAR_LEN) then
        print *, 'PANIC: file name is too long: "' // trim(file) // '"'
        stop
      end if

      dxf%file  = trim(base) // '.dx'
      dxf%dfile = trim(base) // '.bin'

      open(newunit=dxf%hunit, file=trim(dxf%path)//trim(dxf%file), status='replace', action='write', position='rewind')

      open(newunit=dxf%dunit, file=trim(dxf%path)//trim(dxf%dfile), &
        status='replace', action='write', position='rewind', form='unformatted')

    else  ! We'll be writing to stdout
      if (.not.dxf%embed_data) then
        print *, 'PANIC: must use embedded data when writing to stdout'
        stop
      end if
      dxf%hunit = STDOUT
    end if
        
    dxf % object = 1
    dxf % offset = SIZE_OF_REC_HEAD
    
    if (.not.dxf%embed_data) write(unit=dxf % hunit, fmt='(a)') DATA_MODE
    
  end subroutine dx_open
  
  subroutine dx_close (dxf)
  
#ifdef NAG_4_2
!    use f90_unix_io, only: flush
!    use f90_unix_proc, only: system
#endif
    
    type(dx_file), intent(inout) :: dxf
    
    if (dxf%dunit >= 0) then  ! we're not writing to stdout
      !! Close the header file and flush any buffered output to the data file.
      write(unit=dxf % hunit, fmt='(/,a)') 'end'
      close(unit=dxf % hunit)
#ifdef NAG_4_2
      if (dxf%append_data) then
        !! Append the binary data file to the end of the header file
        call flush (dxf % dunit)
        call system ('cat '//trim(dxf%path)//trim(dxf%dfile)//' >> '//trim(dxf%path)//trim(dxf%file))
        !close(unit=dxf % dunit)
        close(unit=dxf % dunit, status='delete')
      else
        close(unit=dxf%dunit)
      end if
#else
      close(unit=dxf%dunit)
#endif

    end if
    
    dxf % hunit = -1
    dxf % dunit = -1
    dxf % object = 0
    dxf % offset = 0
    
  end subroutine dx_close
  
  subroutine export_field (dxf, dxo, pos_obj, dat_obj, con_obj, name)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    type(dx_object), intent(in) :: pos_obj, dat_obj
    type(dx_object), intent(in), optional :: con_obj
    character(len=*), intent(in), optional :: name
    
    if (present(name)) then
      write(unit=dxf % hunit, fmt='(/,a)') 'object "' // trim(name) // '" class field'
      dxo % file => dxf % file
      dxo % number = 0
      dxo % name = name
    else
      write(unit=dxf % hunit, fmt='(/,a)') 'object ' // i2c(dxf % object) // ' class field'
      dxo % file => dxf % file
      dxo % number = dxf % object
      dxf % object = 1 + dxf % object
    end if
    if (present(con_obj)) call write_component (dxf, 'connections', con_obj)
    call write_component (dxf, 'positions', pos_obj)
    call write_component (dxf, 'data', dat_obj)

  end subroutine export_field
  
  subroutine export_real_field_data (dxf, dxo, con_obj, pos_obj, data, name, cc)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    type(dx_object), intent(in) :: con_obj, pos_obj
    real, dimension(:), intent(in) :: data
    character(len=*), intent(in), optional :: name
    logical, intent(in), optional :: cc
    
    logical :: cell_centered
    type(dx_object) :: dat_obj
    
    if (present(cc)) then
      cell_centered = cc
    else
      cell_centered = .false.
    end if
    
    if (cell_centered) then
      call dx_export_array (dxf, dat_obj, data, 'dep', 'connections')
    else
      call dx_export_array (dxf, dat_obj, data, 'dep', 'positions')
    end if
    call dx_export_field (dxf, dxo, pos_obj, dat_obj, con_obj, name)
    
  end subroutine export_real_field_data
        
  subroutine export_real_field_data_1 (dxf, dxo, con_obj, pos_obj, data, name, cc)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    type(dx_object), intent(in) :: con_obj, pos_obj
    real, dimension(:,:), intent(in) :: data
    character(len=*), intent(in), optional :: name
    logical, intent(in), optional :: cc
    
    logical :: cell_centered
    type(dx_object) :: dat_obj
    
    if (present(cc)) then
      cell_centered = cc
    else
      cell_centered = .false.
    end if
    
    if (cell_centered) then
      call dx_export_array (dxf, dat_obj, data, 'dep', 'connections')
    else
      call dx_export_array (dxf, dat_obj, data, 'dep', 'positions')
    end if
    call dx_export_field (dxf, dxo, pos_obj, dat_obj, con_obj, name)
    
  end subroutine export_real_field_data_1
        
  subroutine export_int_field_data (dxf, dxo, con_obj, pos_obj, data, name, cc)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    type(dx_object), intent(in) :: con_obj, pos_obj
    integer, dimension(:), intent(in) :: data
    character(len=*), intent(in), optional :: name
    logical, intent(in), optional :: cc
    
    logical :: cell_centered
    type(dx_object) :: dat_obj
    
    if (present(cc)) then
      cell_centered = cc
    else
      cell_centered = .false.
    end if
    
    if (cell_centered) then
      call dx_export_array (dxf, dat_obj, data, 'dep', 'connections')
    else
      call dx_export_array (dxf, dat_obj, data, 'dep', 'positions')
    end if
    call dx_export_field (dxf, dxo, pos_obj, dat_obj, con_obj, name)
    
  end subroutine export_int_field_data
        
  subroutine write_component (dxf, comp, obj)
  
    type(dx_file), intent(in), target :: dxf
    character(len=*), intent(in) :: comp
    type(dx_object), intent(in) :: obj
  
    write(unit=dxf%hunit, fmt='(t3,a)', advance='no') 'component "' // trim(comp) // '" value '
    if (.not.associated(obj % file, dxf % file)) then ! object belongs to a different file
      write(unit=dxf%hunit, fmt='(a)', advance='no') 'file ' // trim(obj % file) // ','
    end if
  
    if (obj % number > 0) then ! this is a numbered object
      write(unit=dxf%hunit, fmt='(a)') i2c(obj % number)
    else  ! this is a named object
      write(unit=dxf%hunit, fmt='(a)') '"' // trim(obj % name) // '"'
    end if
    
  end subroutine write_component
    
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! SPECIFIC ROUTINES FOR EXPORTING ARRAYS
!!

  subroutine export_integer_array_data_0 (dxf, dxo, array)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    integer, dimension(:), intent(in) :: array
    
    !! Write array descriptor to header file
    write(unit=dxf%hunit, fmt='(/,a)',advance='no') 'object ' // i2c(dxf%object) // &
      ' class array type int items ' // i2c(size(array)) // ' data '
    if (dxf%embed_data) then
      write(unit=dxf%hunit,fmt='(a)') 'follows'
      write(unit=dxf%hunit,fmt='(t3,i0)') array
    else
      if(.not.dxf%append_data) then
        write(unit=dxf%hunit,fmt='(a)',advance='no') 'file ' // trim(dxf%dfile) // ','
      end if
      write(unit=dxf%hunit,fmt='(a)') i2c(dxf%offset)
      write(unit=dxf % dunit) array
      dxf % offset = dxf % offset + SIZE_OF_REC_SEP + size(array) * SIZE_OF_INTEGER
    end if
    
    !! DX object handle
    dxo % file => dxf % file
    dxo % number = dxf % object
    
    !! Update the counters in the DX file structure
    dxf % object = dxf % object + 1
    dxf % offset = dxf % offset + SIZE_OF_REC_SEP + size(array) * SIZE_OF_INTEGER
    
  end subroutine export_integer_array_data_0
  
  subroutine export_integer_array_data_1 (dxf, dxo, array)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    integer, dimension(:,:), intent(in) :: array
    
    !! Write array descriptor to header file
    write(unit=dxf % hunit, fmt='(/,a)',advance='no') 'object ' // i2c(dxf % object) // &
      ' class array type int rank 1 shape ' // i2c(size(array,1)) // &
      ' items ' // i2c(size(array,2)) // ' data '
    if (dxf%embed_data) then
      write(unit=dxf%hunit,fmt='(a)') 'follows'
      write(unit=dxf%hunit,fmt='((t3,'//i2c(size(array,1))//'(i0,:,1x)))') array
    else
      if(.not.dxf%append_data) then
        write(unit=dxf%hunit,fmt='(a)',advance='no') 'file ' // trim(dxf%dfile) // ','
      end if
      write(unit=dxf%hunit,fmt='(a)') i2c(dxf%offset)
      write(unit=dxf % dunit) array
      dxf % offset = dxf % offset + SIZE_OF_REC_SEP + size(array) * SIZE_OF_INTEGER
    end if
    
    !! DX object handle
    dxo % file => dxf % file
    dxo % number = dxf % object
    
    !! Update the counters in the DX file structure
    dxf % object = dxf % object + 1
    
  end subroutine export_integer_array_data_1
  
  subroutine export_real_array_data_0 (dxf, dxo, array)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    real, dimension(:), intent(in) :: array

    !! Write array descriptor to header file
    write(unit=dxf % hunit, fmt='(/,a)',advance='no') 'object ' // i2c(dxf % object) // &
      ' class array type float items ' // i2c(size(array)) // ' data '
    if (dxf%embed_data) then
      write(unit=dxf%hunit,fmt='(a)') 'follows'
      write(unit=dxf%hunit,fmt='(t3,es16.8)') array
    else
      if(.not.dxf%append_data) then
        write(unit=dxf%hunit,fmt='(a)',advance='no') 'file ' // trim(dxf%dfile) // ','
      end if
      write(unit=dxf%hunit,fmt='(a)') i2c(dxf%offset)
      write(unit=dxf % dunit) array
      dxf % offset = dxf % offset + SIZE_OF_REC_SEP + size(array) * SIZE_OF_REAL
    end if
    
    !! DX object handle
    dxo % file => dxf % file
    dxo % number = dxf % object
    
    !! Update the counters in the DX file structure
    dxf % object = dxf % object + 1
    
  end subroutine export_real_array_data_0
  
  subroutine export_real_array_data_1 (dxf, dxo, array)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    real, dimension(:,:), intent(in) :: array
    
    !! Write array descriptor to header file
    write(unit=dxf % hunit, fmt='(/,a)',advance='no') 'object ' // i2c(dxf % object) // &
      ' class array type float rank 1 shape ' // i2c(size(array,1)) // &
      ' items ' // i2c(size(array,2)) // ' data '
    if (dxf%embed_data) then
      write(unit=dxf%hunit,fmt='(a)') 'follows'
      write(unit=dxf%hunit,fmt='(t3,'//i2c(size(array,1))//'es16.8)') array
    else
      if(.not.dxf%append_data) then
        write(unit=dxf%hunit,fmt='(a)',advance='no') 'file ' // trim(dxf%dfile) // ','
      end if
      write(unit=dxf%hunit,fmt='(a)') i2c(dxf%offset)
      write(unit=dxf % dunit) array
      dxf % offset = dxf % offset + SIZE_OF_REC_SEP + size(array) * SIZE_OF_REAL
    end if
    
    !! DX object handle
    dxo % file => dxf % file
    dxo % number = dxf % object
    
    !! Update the object counter in the DX file structure
    dxf % object = dxf % object + 1
    
  end subroutine export_real_array_data_1
  
  subroutine export_integer_array_0_0 (dxf, dxo, array, attr, value)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    integer, dimension(:), intent(in) :: array
    character(len=*), intent(in) :: attr
    character(len=*), intent(in) :: value
    
    call dx_export_array (dxf, dxo, array)
    call write_attributes (dxf, attr, value)

  end subroutine export_integer_array_0_0
  
  subroutine export_integer_array_1_0 (dxf, dxo, array, attr, value)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    integer, dimension(:,:), intent(in) :: array
    character(len=*), intent(in) :: attr
    character(len=*), intent(in) :: value
    
    call dx_export_array (dxf, dxo, array)
    call write_attributes (dxf, attr, value)

  end subroutine export_integer_array_1_0
  
  subroutine export_integer_array_0_1 (dxf, dxo, array, attr, value)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    integer, dimension(:), intent(in) :: array
    character(len=*), dimension(:), intent(in) :: attr
    character(len=*), dimension(:), intent(in) :: value
    
    call dx_export_array (dxf, dxo, array)
    call write_attributes (dxf, attr, value)

  end subroutine export_integer_array_0_1
  
  subroutine export_integer_array_1_1 (dxf, dxo, array, attr, value)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    integer, dimension(:,:), intent(in) :: array
    character(len=*), dimension(:), intent(in) :: attr
    character(len=*), dimension(:), intent(in) :: value
    
    call dx_export_array (dxf, dxo, array)
    call write_attributes (dxf, attr, value)

  end subroutine export_integer_array_1_1
  
  subroutine export_real_array_0_0 (dxf, dxo, array, attr, value)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    real, dimension(:), intent(in) :: array
    character(len=*), intent(in) :: attr
    character(len=*), intent(in) :: value
    
    call dx_export_array (dxf, dxo, array)
    call write_attributes (dxf, attr, value)

  end subroutine export_real_array_0_0
  
  subroutine export_real_array_1_0 (dxf, dxo, array, attr, value)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    real, dimension(:,:), intent(in) :: array
    character(len=*), intent(in) :: attr
    character(len=*), intent(in) :: value
    
    call dx_export_array (dxf, dxo, array)
    call write_attributes (dxf, attr, value)

  end subroutine export_real_array_1_0
  
  subroutine export_real_array_0_1 (dxf, dxo, array, attr, value)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    real, dimension(:), intent(in) :: array
    character(len=*), dimension(:), intent(in) :: attr
    character(len=*), dimension(:), intent(in) :: value
    
    call dx_export_array (dxf, dxo, array)
    call write_attributes (dxf, attr, value)

  end subroutine export_real_array_0_1
  
  subroutine export_real_array_1_1 (dxf, dxo, array, attr, value)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    real, dimension(:,:), intent(in) :: array
    character(len=*), dimension(:), intent(in) :: attr
    character(len=*), dimension(:), intent(in) :: value
    
    call dx_export_array (dxf, dxo, array)
    call write_attributes (dxf, attr, value)

  end subroutine export_real_array_1_1
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! AUXILLARY ROUTINES FOR WRITING NAMED ATTRIBUTE CLAUSES
!!
  
  subroutine write_string_attribute_one (dxf, attr, value)
    type(dx_file), intent(in) :: dxf
    character(len=*), intent(in) :: attr, value
    write(unit=dxf % hunit, fmt='(t3,a)') 'attribute "' // trim(attr) // '" string "' // trim(value) // '"'
  end subroutine write_string_attribute_one
  
  subroutine write_string_attribute_many (dxf, attr, value)
    type(dx_file), intent(in) :: dxf
    character(len=*), dimension(:), intent(in) :: attr, value
    integer :: j
    !! ASSERT( size(attr) == size(value) )
    do j = 1, size(attr)
      call write_string_attribute_one (dxf, attr(j), value(j))
    end do
  end subroutine write_string_attribute_many
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!

  subroutine dx_export_connections (dxf, dxo, array, type)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxo
    integer, dimension(:,:), intent(in) :: array
    character(len=*), intent(in) :: type
    
    character(len=16), dimension(2) :: attr, value
    
    attr(1) = 'element type'
    attr(2) = 'ref'
    
    value(1) = type
    value(2) = 'positions'
    
    select case (type)
    case ('triangles')
      !! ASSERT( size(array,dim=1) == 3 )
    case ('tetrahedra')
      !! ASSERT( size(array,dim=1) == 4 )
    case ('lines')
      !! ASSERT( size(array,dim=1) == 2 )
    case ('quads')
      !! ASSERT( size(array,dim=1) == 4 )
    case ('cubes')
      !! ASSERT( size(array,dim=1) == 6 ) ! Hey! this has to be wrong!
    case default
      print *, 'PANIC!  dx_export_connections: unknown element type "' // trim(type) // '"'
    end select
    
    !! NB: DX does 0-based indexing, so we shift the connections array.
    call dx_export_array (dxf, dxo, array - 1, attr, value)
    
  end subroutine dx_export_connections
  
  subroutine dx_export_regular_grid (dxf, dxcon, dxpos, counts, origin, deltas)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(out) :: dxpos, dxcon
    integer, dimension(:), intent(in) :: counts
    real, dimension(:), intent(in) :: origin
    real, dimension(:,:), intent(in) :: deltas
    
    integer :: j
    
    !! Need some error checking here
    
    !! Write gridposition data to header file
    write(unit=dxf % hunit, fmt='(/,4a)') 'object ' // i2c(dxf % object) // &
      ' class gridpositions counts', ( ' ' // i2c(counts(j)), j = size(counts), 1, -1)
    write(unit=dxf % hunit, fmt='(t3,a,3es13.5)') 'origin', origin
    do j = size(deltas,dim=2), 1, -1
      write(unit=dxf % hunit, fmt='(t3,a,3es13.5)') 'delta', deltas(:,j)
    end do
    
    !! DX object handle
    dxpos % file => dxf % file
    dxpos % number = dxf % object
    
    !! Update the counters in the DX file structure
    dxf % object = dxf % object + 1
    
    !! Write gridconnection data to header file
    write(unit=dxf % hunit, fmt='(/,4a)') 'object ' // i2c(dxf % object) // &
      ' class gridconnections counts', ( ' ' // i2c(counts(j)), j = size(counts), 1, -1)
    
    !! DX object handle
    dxcon % file => dxf % file
    dxcon % number = dxf % object
    
    !! Update the counters in the DX file structure
    dxf % object = dxf % object + 1
    
  end subroutine dx_export_regular_grid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!

  subroutine dx_new_series (dxs, max)
    type(dx_series), intent(inout) :: dxs
    integer, intent(in), optional :: max
    call dx_delete_series (dxs)  ! in case it has been previously used
    if (present(max)) then
      dxs % max = max
    else
      dxs % max = dxs % n_init
    end if
    allocate(dxs % time(dxs % max), dxs % object(dxs % max))
  end subroutine dx_new_series
  
  subroutine dx_delete_series (dxs)
    type(dx_series), intent(inout) :: dxs
    if (associated(dxs % time))   deallocate(dxs % time)
    if (associated(dxs % object)) deallocate(dxs % object)
    dxs % n = 0
    dxs % max = 0
  end subroutine dx_delete_series
  
  subroutine dx_append_to_series (dxs, time, object)
  
    type(dx_series), intent(inout) :: dxs
    real, intent(in) :: time
    type(dx_object), intent(in) :: object
    
    integer :: old_size
    real, dimension(:), pointer :: old_time
    type(dx_object), dimension(:), pointer :: old_object
    
    dxs % n = 1 + dxs % n
    
    !! Verify space exists; (re)allocate if necessary
    if (dxs % n > dxs % max) then  ! out of space; allocate more
      if (associated(dxs % object)) then  ! reallocate, and copy
        old_size = size(dxs % object)
        dxs % max = old_size + dxs % n_incr
        old_object => dxs % object
        old_time   => dxs % time
        allocate(dxs % time(dxs % max), dxs % object(dxs % max))
        dxs % object(1:old_size) = old_object
        dxs % time  (1:old_size) = old_time
        deallocate(old_time, old_object)
      else  ! initial allocation
        dxs % max = dxs % n_init
        allocate(dxs % time(dxs % max), dxs % object(dxs % max))
      end if
    end if
    
    dxs % time(dxs % n) = time
    dxs % object(dxs % n) = object
    
  end subroutine dx_append_to_series
  
  subroutine dx_export_series (dxf, dxo, dxs, name)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(inout) :: dxo
    type(dx_series), intent(in) :: dxs
    character(len=*), intent(in), optional :: name
    
    integer :: j
    
    if (present(name)) then
      write(unit=dxf%hunit, fmt='(/,a)') 'object "' // trim(name) // '" class series'
      dxo % file => dxf % file
      dxo % number = 0
      dxo % name = name
    else
      write(unit=dxf%hunit, fmt='(/,a)') 'object ' // i2c(dxf % object) // ' class series'
      dxo % file => dxf % file
      dxo % number = dxf % object
      dxf % object = 1 + dxf % object
    end if
    
    do j = 1, dxs % n
      call write_series_member (dxf, j-1, dxs % time(j), dxs % object(j))
    end do
    
  end subroutine dx_export_series
  
  subroutine write_series_member (dxf, n, time, object)
  
    type(dx_file), intent(in), target :: dxf
    integer, intent(in) :: n
    real, intent(in) :: time
    type(dx_object), intent(in) :: object
    
    write(unit=dxf%hunit, fmt='(t3,a,es13.5,a)', advance='no') 'member ' // i2c(n) // &
        ' position ', time, ' value '
    if (.not.associated(object % file, dxf % file)) then ! object belongs to a different file
      write(unit=dxf%hunit, fmt='(a)', advance='no') 'file ' // trim(object % file) // ','
    end if
  
    if (object % number > 0) then ! this is a numbered object
      write(unit=dxf%hunit, fmt='(a)') i2c(object % number)
    else  ! this is a named object
      write(unit=dxf%hunit, fmt='(a)') '"' // trim(object % name) // '"'
    end if
    
  end subroutine write_series_member

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!

  subroutine dx_new_multigrid (dxmg, max)
    type(dx_multigrid), intent(inout) :: dxmg
    integer, intent(in), optional :: max
    call dx_delete_multigrid (dxmg)  ! in case it has been previously used
    if (present(max)) then
      dxmg% max = max
    else
      dxmg%max = dxmg%n_init
    end if
    allocate(dxmg%object(dxmg%max))
  end subroutine dx_new_multigrid
  
  subroutine dx_delete_multigrid (dxmg)
    type(dx_multigrid), intent(inout) :: dxmg
    if (associated(dxmg%object)) deallocate(dxmg%object)
    dxmg%n = 0
    dxmg%max = 0
  end subroutine dx_delete_multigrid
  
  subroutine dx_append_to_multigrid (dxmg, object)
  
    type(dx_multigrid), intent(inout) :: dxmg
    type(dx_object), intent(in) :: object
    
    integer :: old_size
    type(dx_object), dimension(:), pointer :: old_object
    
    dxmg%n = 1 + dxmg%n
    
    !! Verify space exists; (re)allocate if necessary
    if (dxmg%n > dxmg%max) then  ! out of space; allocate more
      if (associated(dxmg%object)) then  ! reallocate, and copy
        old_size = size(dxmg%object)
        dxmg%max = old_size + dxmg%n_incr
        old_object => dxmg%object
        allocate(dxmg%object(dxmg%max))
        dxmg%object(1:old_size) = old_object
        deallocate(old_object)
      else  ! initial allocation
        dxmg%max = dxmg%n_init
        allocate(dxmg%object(dxmg%max))
      end if
    end if
    
    dxmg%object(dxmg%n) = object
    
  end subroutine dx_append_to_multigrid
  
  subroutine dx_export_multigrid (dxf, dxo, dxmg, name)
  
    type(dx_file), intent(inout), target :: dxf
    type(dx_object), intent(inout) :: dxo
    type(dx_multigrid), intent(in) :: dxmg
    character(len=*), intent(in), optional :: name
    
    integer :: j
    
    if (present(name)) then
      write(unit=dxf%hunit, fmt='(/,a)') 'object "' // trim(name) // '" class multigrid'
      dxo%file => dxf%file
      dxo%number = 0
      dxo%name = name
    else
      write(unit=dxf%hunit, fmt='(/,a)') 'object ' // i2c(dxf%object) // ' class multigrid'
      dxo%file => dxf%file
      dxo%number = dxf%object
      dxf%object = 1 + dxf%object
    end if
    
    do j = 1, dxmg%n
      call write_multigrid_member (dxf, j-1, dxmg%object(j))
    end do
    
  end subroutine dx_export_multigrid
  
  subroutine write_multigrid_member (dxf, n, object)
  
    type(dx_file), intent(in), target :: dxf
    integer, intent(in) :: n
    type(dx_object), intent(in) :: object
    
    write(unit=dxf%hunit, fmt='(t3,a)', advance='no') 'member ' // i2c(n) // ' value '
    if (.not.associated(object%file, dxf%file)) then ! object belongs to a different file
      write(unit=dxf%hunit, fmt='(a)', advance='no') 'file ' // trim(object%file) // ','
    end if
  
    if (object%number > 0) then ! this is a numbered object
      write(unit=dxf%hunit, fmt='(a)') i2c(object%number)
    else  ! this is a named object
      write(unit=dxf%hunit, fmt='(a)') '"' // trim(object%name) // '"'
    end if
    
  end subroutine write_multigrid_member

end module data_explorer
