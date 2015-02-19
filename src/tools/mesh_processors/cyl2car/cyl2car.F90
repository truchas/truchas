!!
!! CYL2CAR
!!
!! An Exodus mesh utility program that maps a mesh with cylindrical-coordinate
!! data to a mesh with Cartesian coordinate data, transforming element types
!! and side sets as necessary.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 10 November 2004
!!

program cyl2car

  use exodus_mesh_type
  use exodus_mesh_io
  use command_line
  use cyl2car_proc
  implicit none

  character(len=128) :: infile, outfile
  type(exodus_mesh)  :: inmesh, outmesh

  character(len=8), parameter :: CREATOR = 'cyl2car', VERSION = '0.1'

  call parse_command_line (infile, outfile)
  call read_exodus_mesh (infile, inmesh)
  call transform_mesh (inmesh, outmesh)
  call write_exodus_mesh (outfile, outmesh, CREATOR, VERSION)

end program cyl2car
