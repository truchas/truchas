!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module write_probes_command_line

  implicit none
  private
  
  public :: parse_command_line
  
  logical, public :: list_only = .false.
  integer, public :: probe_num = -1
  character(256), public :: h5file
  character(32),  public :: prog
  logical :: found_n = .false.
  
contains

  subroutine parse_command_line

    integer :: n, num_arg, ios
    character(256) :: arg

    !! Get the name of the invoked program.
    call get_command_argument (0, arg)
    n = scan(arg, '/', back=.true.)
    prog = arg(n+1:)  ! remove the leading path component, if any

    n = 0
    num_arg = command_argument_count()

    do while (n < num_arg)

      n = n + 1
      call get_command_argument (n, arg)

      select case (arg)
      case ('-h', '--help')
      
        call usage_summary
        
      case ('-n')

        n = n + 1
        if (n > num_arg) call usage_halt ('option requires an argument: ' // trim(arg))
        call get_command_argument (n, arg)
        read(arg,*,iostat=ios) probe_num
        if (ios /= 0) call usage_halt ('invalid value for -n: ' // trim(arg))
        found_n = .true.

      case ('-l')
      
        list_only = .true.
        
      case default

        if (arg(1:1) == '-') call usage_halt ('invalid option: ' // trim(arg))
        n = n - 1 ! back-up
        exit

      end select

    end do
    
    !! Final argument is the Truchas h5 output file.
    if (num_arg - n /= 1) call usage_halt ('expecting Truchas h5 file argument')
    call get_command_argument (n+1, h5file)
    
    !! Some further option checking.
    if (.not.(found_n .or. list_only)) call usage_halt ('-n option not specified')
    
  end subroutine parse_command_line
  
  subroutine usage_summary
    use,intrinsic :: iso_fortran_env, only: output_unit
    write(output_unit,'(a)') 'Usage: ' // trim(prog) // ' [options] H5FILE'
    write(output_unit,'(a)') 'Writes probe data from a Truchas h5 output file to stdout.'
    write(output_unit,'(a)') 'Options:'
    write(output_unit,'(a)') '  -h, --help  Display this help and exit.'
    write(output_unit,'(a)') '  -l          Print a list of the available probes.'
    write(output_unit,'(a)') '  -n N        Data for probe index N is written to stdout.'
    stop
  end subroutine usage_summary
  
  subroutine usage_halt (message)
    use,intrinsic :: iso_fortran_env, only: error_unit
    character(*), intent(in) :: message
    write(error_unit,'(a)') trim(prog) // ': ' // message
    write(error_unit,'(a)') 'Try `' // trim(prog) // ' -h'' for more information.'
    stop
  end subroutine usage_halt

end module write_probes_command_line

program write_probes

  use h5out_type
  use write_probes_command_line
  use string_utilities, only: i_to_c
  use,intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none
  
  integer :: stat, n, seq_num
  character(256) :: errmsg
  type(h5out) :: ofile
  
  call parse_command_line
  
  !! Access the Truchas h5 output file.
  call h5out_init (ofile, h5file, stat, errmsg)
  if (stat /= 0) then
    write(error_unit,'(a)') trim(errmsg)
    stop
  end if
  
  if (list_only) then
  
    call h5out_write_probe_list (ofile, output_unit)
    
  else
    
    if (probe_num < 1 .or. probe_num > h5out_num_probes(ofile)) then
      write(error_unit,'(a,i0)') 'Invalid probe index: ', probe_num
      write(error_unit,'(a)') 'Use `' // trim(prog) // ' -l H5FILE'' to get a list of the available probes.'
      stop
    end if
    
    call h5out_write_probe_data (ofile, probe_num, output_unit)
    
  end if
  
  call h5out_delete (ofile)
  
end program write_probes
