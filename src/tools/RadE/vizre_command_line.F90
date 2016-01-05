!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vizre_command_line

  implicit none
  private

  public  :: parse_command_line
  private :: usage_halt, usage_summary

  !! Public variables initialized from the command line.
  character(len=32), public :: prog

contains

  subroutine parse_command_line (infile, outfile, row, col, sym)

    character(len=*), intent(out) :: infile, outfile
    integer, pointer :: row(:), col(:)
    logical, intent(out) :: sym

    integer :: i, n, num_arg, ios
    character(len=128) :: arg

    !! Default values for the optional arguments.
    row => null()
    col => null()
    sym = .false.

    call get_command_argument (0, arg)
    i = scan(arg, '/', back=.true.)
    prog = arg(i+1:)  ! Remove the leading path component, if any.

   !!!
   !!! PARSE COMMAND-LINE OPTIONS

    n = 0
    num_arg = command_argument_count()

    do while (n < num_arg)

      n = n + 1
      call get_command_argument (n, arg)

      select case (arg)

      case ('-h','--help')

        call usage_summary ()

      case ('-s')

        sym = .true.

      case ('-r')

        n = n + 1
        if (n > num_arg) call usage_halt ('option requires an argument: ' // trim(arg))
        call get_command_argument (n, arg)
        if (associated(row)) deallocate(row)
        call parse_range_list (arg, row, ios)
        if (ios /= 0) call usage_halt ('invalid argument for -r: ' // trim(arg))

      case ('-c')

        n = n + 1
        if (n > num_arg) call usage_halt ('option requires an argument: ' // trim(arg))
        call get_command_argument (n, arg)
        if (associated(col)) deallocate(col)
        call parse_range_list (arg, col, ios)
        if (ios /= 0) call usage_halt ('invalid argument for -c: ' // trim(arg))

      !!
      !! Example option with following numeric value
      !!
      !! case ('-n')
      !!
      !!   n = n + 1
      !!   if (n > num_arg) call usage_halt ('option requires an argument: ' // trim(arg))
      !!   call get_command_argument (n, arg)
      !!   read(unit=arg,fmt=*,iostat=ios) n_value
      !!   if (ios /= 0) call usage_halt ('invalid value for -n: ' // trim(arg))
      !!
      !! Continue to do sanity checking of the value.
      !!

      case default

        if (arg(1:1) == '-') call usage_halt ('invalid option: ' // trim(arg))
        n = n - 1 ! Back-up.
        exit

      end select

    end do

    if (num_arg - n /= 2) call usage_halt ('exactly two file arguments are required')

    !! Next two arguments are interpreted as file paths.
    call get_command_argument (n+1, infile)
    call get_command_argument (n+2, outfile)

  end subroutine parse_command_line

  subroutine usage_halt (message)

    character(len=*), intent(in) :: message

    write(unit=*, fmt='(a)') trim(prog) // ': ' // message
    write(unit=*, fmt='(a)') 'Try `' // trim(prog) // ' -h'' for more information.'
    stop

  end subroutine usage_halt

  subroutine usage_summary ()

    write(unit=*,fmt='(a)') 'Usage: ' // trim(prog) // ' [options] enclosure_file gmv_file'
    write(unit=*,fmt='(a)') ' '
    write(unit=*,fmt='(a)') 'Writes a GMV-format visualization file for the specified enclosure.'
    write(unit=*,fmt='(a)') 'If the enclosure includes view factor data, the ambient view factor'
    write(unit=*,fmt='(a)') 'and view factor matrix row sums are written as face variables.'
    write(unit=*,fmt='(a)') ' '
    write(unit=*,fmt='(a)') 'Options:'
    write(unit=*,fmt='(a)') '  -r list      Write the specified rows or columns of the view factor matrix'
    write(unit=*,fmt='(a)') '  -c list      as face variables.  List is a comma-separated list of ranges;'
    write(unit=*,fmt='(a)') '               a range is index, first:last, or first:last:stride.'
    write(unit=*,fmt='(a)') '  -s           Write the fully-developed enclosure surface defined by the'
    write(unit=*,fmt='(a)') '               enclosure''s symmetries.  The default is to write just the'
    write(unit=*,fmt='(a)') '               generating surface.'
    write(unit=*,fmt='(a)') '  -h, --help   Display this help and exit.'
    stop

  end subroutine usage_summary

  subroutine split_on_char (string, char, substrings)
    character(len=*), intent(in) :: string
    character(len=1), intent(in) :: char
    character(len=*), pointer :: substrings(:)
    integer :: l, n
    character(len=len(string)) :: tmp
    !! Count the number of substrings; one more than the number of char characters.
    l = 1 + len(string)
    n = 0
    do while (l > 0)
      n = n + 1
      l = scan(string(:l-1), set=char, back=.true.)
    end do
    allocate(substrings(n))
    !! Extract the substrings off the front and shift the string to the left.
    tmp = string
    do n = 1, size(substrings)
      l = scan(tmp, set=char)
      if (l == 0) then
        substrings(n) = tmp
        exit
      else
        substrings(n) = tmp(:l-1)
        tmp = tmp(l+1:)
      end if
    end do
  end subroutine split_on_char

  subroutine parse_range_list (string, list, stat)

    character(len=*), intent(in) :: string
    integer, pointer :: list(:)
    integer, intent(out) :: stat

    integer :: j, k, n, ios
    character(len=15), pointer :: substring1(:), substring2(:)
    integer, allocatable :: range(:,:)

    !! Get the list of range tokens.
    call split_on_char (string, ',', substring1)
    allocate(range(3,size(substring1)))

    !! Get the range values.
    n = 0 ! running count of integer list size
    PARSE: do j = 1, size(substring1)
      stat = -1
      !! Split the range token into components.
      call split_on_char (substring1(j), ':', substring2)
      select case (size(substring2))
      case (1:3)
        !! Read the integer value of each range component.
        do k = 1, size(substring2)
          read(unit=substring2(k),fmt=*,iostat=ios) range(k,j)
          if (ios /= 0) exit PARSE
        end do
        !! Fill in the assumed range components.
        select case (size(substring2))
        case (1)
          range(2,j) = range(1,j)
          range(3,j) = 1
        case (2)
          range(3,j) = 1
        case (3)
          if (range(3,j) == 0) exit PARSE
        end select
        n = n + max(0, 1 + (range(2,j)-range(1,j))/range(3,j))
      case default
        exit PARSE
      end select
      deallocate(substring2)
      stat = 0
    end do PARSE
    deallocate(substring1)

    if (stat == 0) then
      !! Generate the specified list of integers
      allocate(list(n))
      n = 0
      do j = 1, size(range,2)
        do k = range(1,j), range(2,j), range(3,j)
          n = n + 1
          list(n) = k
        end do
      end do
    else
      deallocate(substring2)
      list => null()
    end if

    deallocate(range)

  end subroutine

end module vizre_command_line
