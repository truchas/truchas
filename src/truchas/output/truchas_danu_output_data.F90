!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module truchas_danu_output_data
  use truchas_h5_outfile, only: th5_file, th5_sim_group
  implicit none
  private
  integer, public :: io_group_size
  type(th5_file), target, public :: outfile
  type(th5_sim_group), public :: sim
end module truchas_danu_output_data
