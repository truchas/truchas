!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program genre

  use re_utilities
  use re_dist_vf_type
  use re_encl_type
  use re_exodus_encl
  use re_chaparral_vf
  use re_patch_type
  use re_toolpath
  use genre_command_line
  use parameter_list_type
  use scl
  use timer_tree_type, only: start_timer, stop_timer, write_timer_tree
  use,intrinsic :: iso_fortran_env, only: output_unit
  implicit none

  integer :: lun, ios, n, num_encl, stat
  logical :: is_IOP, found, overwrite, exist
  type(encl_list) :: e
  type(re_patch) :: ep
  type(parameter_list) :: encl_params, chap_params, patch_params
  type(dist_vf) :: vf
  character(:), allocatable :: infile, outfile, ext, basename, errmsg
  character(255) :: iom

  call scl_init
  is_IOP = (scl_rank()==1)

  call start_timer('Total genre execution time')

  if (is_IOP) call parse_command_line(infile, outfile, overwrite)
  call scl_bcast_alloc(infile)
  call scl_bcast_alloc(outfile)
  call scl_bcast(overwrite)

  !! When generating a single radiation enclosure OUTFILE will be used as is.
  !! For multiple enclosures, unique output file names are generated using
  !! OUTFILE by inserting distinguishing labels before the file extension.
  ext = file_extension(outfile)
  n = index(outfile,ext,back=.true.) - 1 ! strip file extension
  basename = outfile(:n) // '.'

  !! If OUTFILE is a path, ensure the directory exists.
  if (is_IOP) call execute_command_line('mkdir -p `dirname ' // outfile // '`', exitstat=stat)
  call scl_bcast(stat)
  if (stat /= 0) call re_halt('unable to create directory for output file: ' // outfile)

  if (is_IOP) open(newunit=lun,file=infile,status='old',action='read',iostat=ios,iomsg=iom)
  call scl_bcast(ios)
  if (ios /= 0) call re_halt('unable to open input file: ' // infile // ': ' // trim(iom))

  !! Construct enclosure
  call start_timer('enclosure mesh creation')
  call read_toolpath_namelists(lun)
  call read_enclosure_namelist(lun, encl_params)
  call init_encl_list(e, encl_params, stat, errmsg)
  if (stat /= 0) call re_halt('error creating the enclosure: ' // errmsg)
  call stop_timer('enclosure mesh creation')

  !! Generate patches
  call read_patches_namelist(lun, patch_params, found)
  if (found) call start_timer('  mesh patch generation')
  call ep%generate_patches(e, patch_params)
  if (found) call stop_timer('  mesh patch generation')

  call read_chaparral_namelist(lun, chap_params, found)

  !! Compute view factors
  num_encl = e%num_encl()
  do n = 1, num_encl
    if (num_encl > 1) outfile = basename // e%this_label() // ext
    if (.not.overwrite) then
      if (is_IOP) inquire(file=outfile,exist=exist)
      call scl_bcast(exist)
      if (exist) then
        call re_info('refusing to overwrite ' // outfile // '; use "-f" to overwrite')
        call e%next_encl
        cycle
      end if
    end if
    call start_timer('  enclosure mesh output')
    call e%write(outfile)
    call ep%write(outfile)
    call stop_timer('  enclosure mesh output')
    if (found) then
      call start_timer('view factor computation')
      call calculate_vf(e, chap_params, ep, vf)
      call stop_timer('view factor computation')
      call start_timer('     view factor output')
      call vf%write(outfile)
      call stop_timer('     view factor output')
    end if
    call re_info('wrote ' // outfile)
    call e%next_encl
  end do

  call stop_timer('Total genre execution time')
  if (is_IOP) then
    call re_info('')
    call write_timer_tree(unit=output_unit, indent=3)
  end if

  call scl_finalize

contains

  !! Return the file extension, or empty string if none
  function file_extension(filename) result(ext)
    character(*), intent(in) :: filename
    character(:), allocatable :: ext
    integer :: n
    n = scan(filename, '/', back=.true.)
    ext = filename(n+1:)
    n = scan(ext, '.', back=.true.)
    if (n > 1) then
      ext = ext(n:)
    else
      ext = ''
    end if
  end function

end program genre
