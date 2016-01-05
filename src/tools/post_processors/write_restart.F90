!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module write_restart_command_line

  use kinds, only: r8
  implicit none
  private
  
  public :: parse_command_line
  
  logical, public :: list_only = .false., mapped = .false.
  integer, public :: cycle_num = -1
  real(r8), public :: coord_scale_factor = 1.0_r8
  character(256), public :: rfile = '', mfile = ''
  character(256), public :: h5file
  character(32),  public :: prog
  
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
        read(arg,*,iostat=ios) cycle_num
        if (ios /= 0) call usage_halt ('invalid value for -n: ' // trim(arg))

      case ('-l')
      
        list_only = .true.
        
      case ('-o')
      
        n = n + 1
        if (n > num_arg) call usage_halt ('option requires an argument: ' // trim(arg))
        call get_command_argument (n, rfile)
        if (len_trim(rfile) == 0) call usage_halt ('invalid value for -o: ' // trim(arg))
        
      case ('-m')
      
        n = n + 1
        if (n > num_arg) call usage_halt ('option requires an argument: ' // trim(arg))
        call get_command_argument (n, mfile)
        mapped = .true.
        
      case ('-s')
      
        n = n + 1
        if (n > num_arg) call usage_halt ('option requires an argument: ' // trim(arg))
        call get_command_argument (n, arg)
        read(arg,*,iostat=ios) coord_scale_factor
        if (ios /= 0) call usage_halt ('invalid value for -s: ' // trim(arg))
        if (coord_scale_factor <= 0.0) call usage_halt ('invalid value for -s: ' // trim(arg))

      case default

        if (arg(1:1) == '-') call usage_halt ('invalid option: ' // trim(arg))
        n = n - 1 ! back-up
        exit

      end select

    end do
    
    !! Final argument is the Truchas h5 output file.
    if (num_arg - n /= 1) call usage_halt ('expecting Truchas h5 file argument')
    call get_command_argument (n+1, h5file)
    
  end subroutine parse_command_line
  
  subroutine usage_summary
    use,intrinsic :: iso_fortran_env, only: output_unit
    write(output_unit,'(a)') 'Usage: ' // trim(prog) // ' [options] H5FILE'
    write(output_unit,'(a)') 'Creates a Truchas restart file using data from an h5 output file.'
    write(output_unit,'(a)') 'Options:'
    write(output_unit,'(a)') '  -h, --help  Display this help and exit.'
    write(output_unit,'(a)') '  -l          Print a list of the available cycles from which the restart'
    write(output_unit,'(a)') '                file can be created.  No restart file is written.'
    write(output_unit,'(a)') '  -n N        Data from cycle N is used to generate the restart file;'
    write(output_unit,'(a)') '                if not specified the last cycle is used.'
    write(output_unit,'(a)') '  -o FILE     Write restart data to FILE.  If not specified, FILE is taken'
    write(output_unit,'(a)') '                to be the H5FILE name with the .h5 suffix replaced by'
    write(output_unit,'(a)') '                .restart.N, where N is the cycle number.'
    write(output_unit,'(a)') '  -m FILE     Create a mapped restart file using the specified ExodusII mesh'
    write(output_unit,'(a)') '                FILE as the target mesh.'
    write(output_unit,'(a)') '  -s FLOAT    Scale the mapped restart mesh by the factor FLOAT.'
    stop
  end subroutine usage_summary
  
  subroutine usage_halt (message)
    use,intrinsic :: iso_fortran_env, only: error_unit
    character(*), intent(in) :: message
    write(error_unit,'(a)') trim(prog) // ': ' // message
    write(error_unit,'(a)') 'Try `' // trim(prog) // ' -h'' for more information.'
    stop
  end subroutine usage_halt

end module write_restart_command_line

program write_restart

  use h5out_type
  use write_restart_command_line
  use string_utilities, only: i_to_c
  use,intrinsic :: iso_fortran_env, only: output_unit, error_unit
  use exodus_mesh_type
  use exodus_mesh_io, only: read_exodus_mesh
  use mapped_restart
  implicit none
  
  integer :: stat, n, seq_num
  character(256) :: errmsg
  type(h5out) :: ofile
  type(exodus_mesh) :: exomesh
  character(:), allocatable :: errstr
  
  call parse_command_line
  
  !! Access the Truchas h5 output file.
  call h5out_init (ofile, h5file, stat, errmsg)
  if (stat /= 0) then
    write(error_unit,'(a)') trim(errmsg)
    stop
  end if
  
  if (list_only) then
  
    call h5out_write_cycle_list (ofile, output_unit)
    
  else
    
    !! Determine the cycle sequence number.
    if (cycle_num == -1) then
      seq_num = h5out_last_seq_num(ofile)
      cycle_num = h5out_cycle_num(ofile, seq_num)
      if (seq_num == 0) then
        write(error_unit,'(a)') 'No cycle output exists!'
        stop
      end if
    else
      seq_num = h5out_seq_num(ofile, cycle_num)
      if (seq_num == 0) then
        write(error_unit,'(a)') 'No output for cycle ' // i_to_c(cycle_num) // ' exists!'
        write(error_unit,'(a)') 'Use `' // trim(prog) // ' -l H5FILE'' to get a list of available cycles.'
        stop
      end if
    end if
    INSIST(seq_num > 0)
  
    !! Generate the restart file name if not specified.
    if (rfile == '') then
      n = scan(h5file, '/', back=.true.)
      rfile = h5file(n+1:)  ! remove any leading path component
      n = index(rfile, '.h5', back=.true.)
      if (n > 0) rfile(n:) = ''  ! drop the .h5 suffix
      rfile = trim(rfile) // '.restart.' // i_to_c(cycle_num)
    end if
    
    !! Write the restart file.
    open(10,file=trim(rfile),form='unformatted',status='replace')
    if (mapped) then
    
      !! Read the ExodusII mesh file.
      call read_exodus_mesh (mfile, exomesh, stat, errstr)
      if (stat /= 0) then
        write(error_unit,'(3a)') 'Error reading ', trim(mfile), ': ' // errstr
        stop
      end if
      
      if (coord_scale_factor /= 1.0) exomesh%coord = coord_scale_factor * exomesh%coord

      call write_mapped_restart (ofile, seq_num, exomesh, 10)
    
    else
    
      call h5out_write_restart_file (ofile, seq_num, 10)
      
    end if
    close(10)
    
  end if
  
  call h5out_delete (ofile)
  
end program write_restart
