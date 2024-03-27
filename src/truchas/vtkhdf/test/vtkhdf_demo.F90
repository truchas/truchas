program vtkhdf_demo

  use,intrinsic :: iso_fortran_env, only: r8 => real64, int8
  !use hdf5_c_binding
  !use hl_hdf5
  use vtkhdf_file_type
  implicit none

  integer :: stat
!  integer(hid_t) :: file_id, grp_id, steps_id, pgrp_id, pogrp_id, cgrp_id, cogrp_id
  real(r8), allocatable :: x(:,:), scalar_cell_data(:), vector_cell_data(:,:)
  real(r8), allocatable :: scalar_point_data(:), vector_point_data(:,:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  character(:), allocatable :: errmsg

  type(vtkhdf_file) :: vizfile
  
  call vizfile%create('demo.vtkhdf', stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if

  call get_mesh1_data(x, cnode, xcnode, types)
  call vizfile%write_mesh(x, cnode, xcnode, types, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if

  call get_scalar_cell_data(x, cnode, xcnode, scalar_cell_data)
  call get_vector_cell_data(x, cnode, xcnode, vector_cell_data)

  call get_scalar_point_data(x, scalar_point_data)
  call get_vector_point_data(x, vector_point_data)
  
  !! Register the datasets that evolve with time. At this stage the data arrays
  !! are only used to glean their types and shapes.

  call vizfile%register_temporal_cell_dataset('cell-radius', scalar_cell_data, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if
  
  call vizfile%register_temporal_cell_dataset('cell-velocity', vector_cell_data, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if
  
  call vizfile%register_temporal_point_dataset('point-radius', scalar_point_data, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if
  
  call vizfile%register_temporal_point_dataset('point-velocity', vector_point_data, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if

  !!!! Write the datasets for the first time step !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call vizfile%write_time_step(0.0_r8)

  call vizfile%write_temporal_cell_dataset('cell-radius', scalar_cell_data, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if

  call vizfile%write_temporal_cell_dataset('cell-velocity', vector_cell_data, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if

  call vizfile%write_temporal_point_dataset('point-radius', scalar_point_data, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if

  call vizfile%write_temporal_point_dataset('point-velocity', vector_point_data, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if

  !!!! Write the datasets for the second time step !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call vizfile%write_time_step(10.0_r8)

  call vizfile%write_temporal_cell_dataset('cell-radius', scalar_cell_data+1, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if

  call vizfile%write_temporal_cell_dataset('cell-velocity', vector_cell_data+1, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if

  call vizfile%write_temporal_point_dataset('point-radius', scalar_point_data+1, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if

  call vizfile%write_temporal_point_dataset('point-velocity', vector_point_data+1, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if

  !! At any point you can write a dataset that isn't time dependent, but its name must
  !! be unique from any other dataset temporal or not of the same type (cell or point).

  call vizfile%write_cell_dataset('static-cell-scalar', -scalar_cell_data, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if

  call vizfile%write_cell_dataset('static-cell-vector', -vector_cell_data, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if

  call vizfile%write_point_dataset('static-point-scalar', -scalar_point_data, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if

  call vizfile%write_point_dataset('static-point-vector', -vector_point_data, stat, errmsg)
  if (stat /= 0) then
    print *, errmsg
    stop 1
  end if

  call vizfile%close

contains

  ! A 5-tet subdivision of the squished unit cube.
  subroutine get_mesh1_data(x, cnode, xcnode, types)
    real(r8), allocatable, intent(out) :: x(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable :: types(:)
    x = 0.5_r8*reshape([0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1], shape=[3,8])
    ! distort to catch C/Fortran index ordering errors
    x(2,:) = 0.75_r8*x(2,:)
    x(3,:) = 1.25_r8*x(3,:)
    cnode = [1,2,4,5, 2,3,4,7, 2,5,6,7, 4,5,7,8, 2,4,5,7]
    xcnode = [1,5,9,13,17,21]
    types = [10, 10, 10, 10, 10] ! all tets
  end subroutine

  subroutine get_scalar_point_data(x, pdata)
    real(r8), intent(in) :: x(:,:)
    real(r8), allocatable, intent(out) :: pdata(:)
    integer :: j
    allocate(pdata(size(x,dim=2)))
    do j = 1, size(x,dim=2)
      pdata(j) = norm2(x(:,j))
    end do
  end subroutine

  subroutine get_vector_point_data(x, pdata)
    real(r8), intent(in) :: x(:,:)
    real(r8), allocatable, intent(out) :: pdata(:,:)
    integer :: j
    pdata = x
  end subroutine

  subroutine get_scalar_cell_data(x, cnode, xcnode, cdata)
    real(r8), intent(in) :: x(:,:)
    integer, intent(in) :: cnode(:), xcnode(:)
    real(r8), allocatable, intent(out) :: cdata(:)
    integer :: j
    allocate(cdata(size(xcnode)-1))
    do j = 1, size(cdata)
      associate(pid => cnode(xcnode(j):xcnode(j+1)-1))
        cdata(j) = norm2(sum(x(:,pid),dim=2)/size(pid))
      end associate
    end do
  end subroutine

  subroutine get_vector_cell_data(x, cnode, xcnode, cdata)
    real(r8), intent(in) :: x(:,:)
    integer, intent(in) :: cnode(:), xcnode(:)
    real(r8), allocatable, intent(out) :: cdata(:,:)
    integer :: j
    allocate(cdata(size(x,dim=1),size(xcnode)-1))
    do j = 1, size(cdata,dim=2)
      associate(pid => cnode(xcnode(j):xcnode(j+1)-1))
        cdata(:,j) = sum(x(:,pid),dim=2)/size(pid)
      end associate
    end do
  end subroutine

end program
