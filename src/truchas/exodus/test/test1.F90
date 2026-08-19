!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test1

  use exodus_mesh_type
  use exodus_mesh_io
  implicit none

  integer :: status
  logical :: exist
  character(len=64) :: file, file_copy
  type(exodus_mesh) :: mesh, mesh_copy

  call get_command_argument (1, file, status=status)
  if (status /= 0) stop 'FAIL'
  file_copy = trim(file)//'.copy'
print *, '@'
  
  !! Read the mesh file...
  call read_exodus_mesh (file, mesh, status)
  if (status /= 0) stop 'FAIL'
print *, '@@'
  
  !! Write it to another file...
  call write_exodus_mesh (file_copy, mesh, 'test', '1', status)
  if (status /= 0) stop 'FAIL'
print *, '@@@'
  
  !! Read the mesh file just written...
  call read_exodus_mesh (file_copy, mesh_copy, status)
  if (status /= 0) stop 'FAIL'
print *, '@@@@'
  
  !! Compare the two meshes
  if (mesh /= mesh_copy) stop 'FAIL'
  
  stop 'PASS'
  
end program test1
  
  
