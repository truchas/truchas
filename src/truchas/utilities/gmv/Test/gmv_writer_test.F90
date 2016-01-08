!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program gmv_writer_test

#ifdef NAG
  use f90_unix_proc, only: system, exit
#endif
  use mesh_manager
  use gmv_writer
  use parallel_communication
  implicit none

  integer, parameter :: r8 = selected_real_kind(10,50)  ! 8-byte IEEE float
  integer :: status, j
  
  type(dist_mesh), pointer :: mesh
  real(kind=r8), allocatable :: u(:)
  real(kind=r8) :: xc(3)
  
  call init_parallel_communication ()
  
  if (is_IOP) open(unit=10,file='gmv_writer_test.inp',position='rewind',action='read',status='old')
  call read_mesh_namelists (10)
  call init_mesh_manager ()
  
  mesh => named_mesh_ptr('mesh1')
  
  allocate(u(mesh%ncell))
  do j = 1, mesh%ncell
    xc = sum(mesh%x(:,mesh%cnode(:,j)),dim=2) / size(mesh%cnode,dim=1)  ! cell centroid
    u(j) = sum(xc)  ! linear function x + y + z
  end do
  
  call gmv_open ('mesh1.gmv')
  call gmv_write_dist_mesh (mesh)
  call gmv_begin_variables(time=3.14_r8, seq=666)
  call gmv_write_dist_cell_var (mesh, u, name='linear')
  call gmv_end_variables ()
  call gmv_close ()

  if (is_IOP) call system ('cmp mesh1.gmv mesh1.gmv.ref', status)
  call broadcast (status)
  if (status /= 0) status = 1
  
  deallocate(u)
  
  call halt_parallel_communication ()
  
  call exit (status)
  
end program gmv_writer_test
