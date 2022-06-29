!!
!! TRUCHAS_ENV
!!
!! This evolving module is intended to collect data and procedures that are
!! associated with the run-time environment of Truchas -- things like file
!! and directory paths.  Much of this is found in output_module, which is
!! slated for removal, and other places.  I expect this to eventually morph
!! into something reasonable.
!!
!! Neil Carlson <nnc@lanl.gov>
!! October 2012
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module truchas_env

  implicit none
  private
  
  public :: output_file_name
  
  ! Stuff formerly in output_data_dir
  character(256), public :: prefix, input_dir, output_dir, title, input_file

  logical, public :: overwrite_output

contains

  !! A simplified replacement for make_file_name.  Given a string EXT, it
  !! returns the path of the output file in the output directory having
  !! EXT as the extention.  The base name of the file is determined by the
  !! input file name.  This does no more than concatenate strings -- I don't
  !! know how necessary this really is.

  function output_file_name (ext) result (ofile)
    character(*), intent(in) :: ext
    character(len_trim(prefix)+len_trim(ext)+1) :: ofile
    ofile = trim(prefix) // '.' // trim(ext)
  end function output_file_name

end module truchas_env
