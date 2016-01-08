!!
!! SCRAMBLE
!!
!! Scramble is a Genesis/Exodus mesh utility that replaces the connectivity
!! data for each element with an equivalent, randomly-selected, node ordering.
!! It can be used to generate families of equivalent (but different) meshes for
!! testing purposes.  Since a 'scrambled' mesh is equivalent to the original,
!! numerical results using the mesh should be identical, save for differences
!! arising from order-dependence.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 2 December 2006
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program scramble

  use exodus_mesh_type
  use exodus_mesh_io
  use command_line
  use scramble_proc
  implicit none

  character(len=128) :: infile, outfile
  type(exodus_mesh)  :: mesh
  integer :: seed

  character(len=8), parameter :: CREATOR = 'scramble', VERSION = '1.0'

  call parse_command_line (infile, outfile, seed)
  call read_exodus_mesh (infile, mesh)
  call scramble_elem_nodes (seed, mesh)
  call write_exodus_mesh (outfile, mesh, CREATOR, VERSION)

end program scramble
