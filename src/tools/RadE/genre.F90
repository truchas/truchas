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
  implicit none

  integer :: lun, ios, n, num_encl, stat
  logical :: is_IOP, found
  type(encl_list) :: e
  type(re_patch) :: ep
  type(parameter_list) :: encl_params, chap_params, patch_params
  type(dist_vf) :: vf
  character(:), allocatable :: infile, outfile, ext, basename, msg, errmsg
  character(255) :: iom

  call scl_init
  is_IOP = (scl_rank()==1)

  if (is_IOP) call parse_command_line(infile, outfile)
  call scl_bcast_alloc(infile)
  call scl_bcast_alloc(outfile)

  !! When generating a single radiation enclosure OUTFILE will be used as is.
  !! For multiple enclosures, unique output file names are generated using
  !! OUTFILE by inserting distinguishing labels before the file extension.
  ext = file_extension(outfile)
  n = index(outfile,ext,back=.true.) - 1 ! strip file extension
  basename = outfile(:n) // '.'

  if (is_IOP) open(newunit=lun,file=infile,status='old',action='read',iostat=ios,iomsg=iom)
  call scl_bcast(ios)
  if (ios /= 0) call re_halt('unable to open input file: ' // infile // ': ' // trim(iom))

  !! Construct enclosure
  call read_toolpath_namelists(lun)
  call read_enclosure_namelist(lun, encl_params)
  call init_encl_list(e, encl_params, stat, errmsg)
  if (stat /= 0) call re_halt('error creating the enclosure: ' // errmsg)

  !! Generate patches
  call read_patches_namelist(lun, patch_params, found)
  call ep%generate_patches(e, patch_params)

  call read_chaparral_namelist(lun, chap_params, found)

  !! Compute view factors
  num_encl = e%num_encl()
  do n = 1, num_encl
    if (num_encl > 1) outfile = basename // e%this_label() // ext
    call e%write(outfile)
    call ep%write_patch_data(outfile)
    if (found) then
      call calculate_vf(e, chap_params, ep, vf)
      call write_dist_vf(vf, outfile)
    end if
    call re_info('wrote ' // outfile)
    call e%next_encl
  end do

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
