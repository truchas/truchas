!!
!! SIMULATION_PROVENANCE
!!
!! This module writes compact simulation-provenance records to a simulation
!! log. These identify the code, build, run environment, and staged input.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes:
!!   The simulation_build_info module is configured by CMake with build-time
!!   compiler and host information.
!!

module simulation_provenance

  use mpi_f08, only: MPI_Bcast, MPI_CHARACTER, MPI_Comm, MPI_INTEGER, MPI_SUCCESS
  use simulation_build_info
  use simulation_environment_type
  use version_info, only: version
  implicit none
  private

  public :: write_simulation_prologue

contains

  subroutine write_simulation_prologue(env, program, simulation, input_file, stat, errmsg)

    type(simulation_environment), intent(in) :: env
    character(*), intent(in) :: program, simulation, input_file
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(8) :: date
    character(10) :: time, zone
    character(:), allocatable :: code_version, cpu, run_host, sha256, timestamp

    call version(code_version)
    cpu = 'unknown'
    run_host = 'unknown'
    if (env%rank == 0) then
      call run_cpu_model(env%output_dir, cpu)
      call get_run_host(env%output_dir, run_host)
    end if
    call date_and_time(date, time, zone)
    timestamp = date(:4) // '-' // date(5:6) // '-' // date(7:8) // 'T' // &
        time(:2) // ':' // time(3:4) // ':' // time(5:6) // trim(zone)

    call env%simlog%info('program=' // log_quoted(program) // ' simulation=' // &
        log_quoted(simulation) // ' version=' // log_quoted(code_version))
    call env%simlog%info('build compiler=' // log_quoted(compiler_id) // &
        ' compiler_version=' // log_quoted(compiler_version) // ' host=' // log_quoted(host) // &
        ' flags=' // log_quoted(build_flags))
    call env%simlog%info('run timestamp=' // log_quoted(timestamp) // ' host=' // log_quoted(run_host) // &
        ' cpu=' // log_quoted(cpu) // ' mpi_processes=' // integer_string(env%nproc))
    call stage_input_file(env, input_file, 'input.json', sha256, stat, errmsg)
    if (stat /= 0) return
    if (env%rank == 0) call env%simlog%info('input file=' // log_quoted(input_file) // &
        ' sha256=' // log_quoted(sha256))
    call env%simlog%info('')

  end subroutine write_simulation_prologue


  subroutine run_cpu_model(output_dir, cpu)

    character(*), intent(in) :: output_dir
    character(:), allocatable, intent(inout) :: cpu

    character(256) :: line
    character(:), allocatable :: command, temp_file
    integer :: exitstat, ios, unit

    line = ''
    temp_file = join_path(output_dir, '.run_cpu.tmp')
    command = 'case "$(uname -s)" in Linux) lscpu | awk -F: ' // &
        '''/^Model name:/ {sub(/^[[:space:]]*/, "", $2); print $2; exit}'' ;; ' // &
        'Darwin) sysctl -n machdep.cpu.brand_string ;; *) uname -m ;; esac > ' // shell_quoted(temp_file) // ' 2>/dev/null'
    call execute_command_line('sh -c ' // shell_quoted(command), wait=.true., exitstat=exitstat)
    if (exitstat /= 0) then
      call delete_file(temp_file)
      return
    end if

    open(newunit=unit, file=temp_file, action='read', status='old', iostat=ios)
    if (ios == 0) then
      read(unit, '(a)', iostat=ios) line
      close(unit, status='delete')
    end if
    if (ios == 0 .and. len_trim(line) > 0) cpu = trim(adjustl(line))

  end subroutine run_cpu_model


  subroutine get_run_host(output_dir, run_host)

    character(*), intent(in) :: output_dir
    character(:), allocatable, intent(inout) :: run_host

    character(256) :: line
    character(:), allocatable :: command, temp_file
    integer :: exitstat, ios, unit

    line = ''
    temp_file = join_path(output_dir, '.run_host.tmp')
    command = 'hostname > ' // shell_quoted(temp_file) // ' 2>/dev/null'
    call execute_command_line('sh -c ' // shell_quoted(command), wait=.true., exitstat=exitstat)
    if (exitstat /= 0) then
      call delete_file(temp_file)
      return
    end if

    open(newunit=unit, file=temp_file, action='read', status='old', iostat=ios)
    if (ios == 0) then
      read(unit, '(a)', iostat=ios) line
      close(unit, status='delete')
    end if
    if (ios == 0 .and. len_trim(line) > 0) run_host = trim(adjustl(line))

  end subroutine get_run_host


  subroutine broadcast_errmsg(comm, io_process, errmsg)

    type(MPI_Comm), intent(in) :: comm
    logical, intent(in) :: io_process
    character(:), allocatable, intent(inout) :: errmsg

    integer :: ierr, length

    length = 0
    if (io_process) length = len(errmsg)
    call MPI_Bcast(length, 1, MPI_INTEGER, 0, comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    if (.not.io_process) allocate(character(length) :: errmsg)
    call MPI_Bcast(errmsg, length, MPI_CHARACTER, 0, comm, ierr)

  end subroutine broadcast_errmsg


  function integer_string(value) result(string)

    integer, intent(in) :: value
    character(:), allocatable :: string
    character(32) :: buffer

    write(buffer, '(i0)') value
    string = trim(buffer)

  end function integer_string


  function log_quoted(value) result(string)

    character(*), intent(in) :: value
    character(:), allocatable :: string

    integer :: i

    string = achar(34)
    do i = 1, len_trim(value)
      if (value(i:i) == achar(34) .or. value(i:i) == achar(92)) string = string // achar(92)
      string = string // value(i:i)
    end do
    string = string // achar(34)

  end function log_quoted


  function join_path(directory, name) result(path)

    character(*), intent(in) :: directory, name
    character(:), allocatable :: path

    if (trim(directory) == '.') then
      path = trim(name)
    else if (directory(len_trim(directory):len_trim(directory)) == '/') then
      path = trim(directory) // trim(name)
    else
      path = trim(directory) // '/' // trim(name)
    end if

  end function join_path


  subroutine delete_file(filename)

    character(*), intent(in) :: filename

    logical :: exists
    integer :: ios, unit

    inquire(file=filename, exist=exists)
    if (.not.exists) return
    open(newunit=unit, file=filename, status='old', action='readwrite', iostat=ios)
    if (ios == 0) close(unit, status='delete')

  end subroutine delete_file


  function shell_quoted(value) result(string)

    character(*), intent(in) :: value
    character(:), allocatable :: string

    integer :: i

    string = "'"
    do i = 1, len_trim(value)
      if (value(i:i) == "'") then
        string = string // achar(39) // achar(34) // achar(39) // achar(34) // achar(39)
      else
        string = string // value(i:i)
      end if
    end do
    string = string // "'"

  end function shell_quoted


  logical function same_path(path1, path2)

    character(*), intent(in) :: path1, path2
    character(:), allocatable :: normalized1, normalized2

    normalized1 = trim(path1)
    normalized2 = trim(path2)
    do while (index(normalized1, './') == 1)
      normalized1 = normalized1(3:)
    end do
    do while (index(normalized2, './') == 1)
      normalized2 = normalized2(3:)
    end do
    same_path = normalized1 == normalized2

  end function same_path


  subroutine stage_input_file(env, source_file, copy_name, sha256, stat, errmsg)

    type(simulation_environment), intent(in) :: env
    character(*), intent(in) :: source_file, copy_name
    character(:), allocatable, intent(out) :: sha256
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: ierr

    stat = 0
    if (env%rank == 0) then
      if (scan(copy_name, '/\') /= 0) then
        stat = 1
        errmsg = 'input-copy filename must be a logical filename'
      else
        call copy_file(source_file, join_path(env%output_dir, copy_name), stat, errmsg)
        if (stat == 0) call file_sha256(join_path(env%output_dir, copy_name), sha256, stat, errmsg)
      end if
    end if
    call MPI_Bcast(stat, 1, MPI_INTEGER, 0, env%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      stat = 1
      errmsg = 'could not communicate input-staging status'
      return
    end if
    if (stat /= 0) call broadcast_errmsg(env%comm, env%rank == 0, errmsg)

  end subroutine stage_input_file


  subroutine copy_file(source_file, destination_file, stat, errmsg)

    character(*), intent(in) :: source_file, destination_file
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: command
    integer :: exitstat

    stat = 0
    if (len_trim(source_file) == 0) then
      stat = 1
      errmsg = 'input filename must not be empty'
      return
    end if
    if (len_trim(destination_file) == 0) then
      stat = 1
      errmsg = 'input-copy filename must not be empty'
      return
    end if
    if (same_path(source_file, destination_file)) return

    command = 'cp ' // shell_quoted(source_file) // ' ' // shell_quoted(destination_file) // ' >/dev/null 2>&1'
    call execute_command_line('sh -c ' // shell_quoted(command), wait=.true., exitstat=exitstat)
    if (exitstat /= 0) then
      stat = 1
      errmsg = 'could not stage input file for provenance'
    end if

  end subroutine copy_file


  subroutine file_sha256(filename, sha256, stat, errmsg)

    character(*), intent(in) :: filename
    character(:), allocatable, intent(out) :: sha256
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(256) :: line
    character(:), allocatable :: command, temp_file
    integer :: exitstat, ios, unit

    line = ''
    temp_file = trim(filename) // '.sha256.tmp'
    command = 'if command -v sha256sum >/dev/null 2>&1; then sha256sum < ' // shell_quoted(filename) // &
        ' > ' // shell_quoted(temp_file) // ' 2>/dev/null; else shasum -a 256 < ' // shell_quoted(filename) // &
        ' > ' // shell_quoted(temp_file) // ' 2>/dev/null; fi'
    call execute_command_line('sh -c ' // shell_quoted(command), wait=.true., exitstat=exitstat)
    if (exitstat /= 0) then
      call delete_file(temp_file)
      stat = 1
      errmsg = 'could not calculate input SHA-256'
      return
    end if

    open(newunit=unit, file=temp_file, action='read', status='old', iostat=ios)
    if (ios == 0) then
      read(unit, '(a)', iostat=ios) line
      close(unit, status='delete')
    end if
    if (ios /= 0 .or. .not.is_sha256(line)) then
      if (ios /= 0) call delete_file(temp_file)
      stat = 1
      errmsg = 'could not read input SHA-256'
      return
    end if

    sha256 = line(:64)
    stat = 0

  end subroutine file_sha256


  logical function is_sha256(string)

    character(*), intent(in) :: string

    integer :: i

    is_sha256 = len_trim(string) >= 64
    if (.not.is_sha256) return
    do i = 1, 64
      if (index('0123456789abcdefABCDEF', string(i:i)) == 0) then
        is_sha256 = .false.
        return
      end if
    end do

  end function is_sha256

end module simulation_provenance
