!!
!! VTKHDF_FILE_TYPE
!!
!! This module defines a derived type for exporting mesh-based solution data
!! to a VTKHDF format file that can be read by the ParaView visualizaton tool.
!! The format uses HDF5 for on-disk storage.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! March 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! The up-to-date specification for VTKHDF is at
!! https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#vtkhdf-file-format
!! and is in development and evolving. The VTKHDF development road map is at
!! https://discourse.vtk.org/t/vtkhdf-roadmap/13257
!!
!! This module was written for version 2.1 of the format.
!!
!! This module uses the "UnstructuredGrid" type of format, and supports static
!! and time dependent datasets (node and cell based), but both assume a single
!! static mesh.
!!
!! This module is currently serial only: data must be collated onto a single
!! process and written from that process using this object.
!!

#include "f90_assert.fpp"

module vtkhdf_file_type

  use,intrinsic :: iso_fortran_env
  use hdf5_c_binding
  use hl_hdf5
  implicit none
  private

  type, public :: vtkhdf_file
    private
    integer(hid_t) :: file_id = -1
    integer(hid_t) :: vtk_id=-1, cgrp_id=-1, pgrp_id=-1, steps_id=-1, cogrp_id=-1, pogrp_id=-1
    integer :: nnode, ncell
    logical :: temporal = .true.
    integer :: nsteps = -1
    type(temporal_dataset), pointer :: temporal_point_dsets => null()
    type(temporal_dataset), pointer :: temporal_cell_dsets => null()
  contains
    procedure :: create
    procedure :: write_mesh
    procedure :: write_time_step
    generic :: write_cell_dataset  => write_cell_dataset_real64
    generic :: write_point_dataset => write_point_dataset_real64
    generic :: register_temporal_cell_dataset  => register_temporal_cell_dataset_real64
    generic :: register_temporal_point_dataset => register_temporal_point_dataset_real64
    generic :: write_temporal_cell_dataset  => write_temporal_cell_dataset_real64
    generic :: write_temporal_point_dataset => write_temporal_point_dataset_real64
    procedure, private :: write_cell_dataset_real64
    procedure, private :: write_point_dataset_real64
    procedure, private :: register_temporal_cell_dataset_real64
    procedure, private :: register_temporal_point_dataset_real64
    procedure, private :: write_temporal_cell_dataset_real64
    procedure, private :: write_temporal_point_dataset_real64
    procedure :: close
  end type

  type :: temporal_dataset
    character(:), allocatable :: name
    integer :: next_offset = 0
    logical :: flag = .false.
    type(temporal_dataset), pointer :: next => null()
  contains
    final :: temporal_dataset_delete
  end type

  !! The VTK cell types that are relevant to Truchas.
  !! NB: The local node ordering for a VTK wedge cell may differ from Truchas
  integer(int8), parameter, public :: VTK_TETRA = 10
  integer(int8), parameter, public :: VTK_HEXAHEDRON = 12
  integer(int8), parameter, public :: VTK_WEDGE = 13
  integer(int8), parameter, public :: VTK_PYRAMID = 14

  integer, parameter, public :: vtkhdf_version(*) = [2,1]

contains

  !! Create a new VTKHDF format file, overwriting any existing file. A file
  !! supporting time-dependent data is created by default. Set the optional
  !! TEMPORAL argument to false for a static file.

  subroutine create(this, filename, stat, errmsg, temporal)

    class(vtkhdf_file), intent(out) :: this
    character(*), intent(in) :: filename
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    logical, intent(in), optional :: temporal

    call init_hdf5

    this%file_id = H5Fcreate(filename, H5F_ACC_TRUNC)
    if (this%file_id < 0) then
      stat = 1
      errmsg = 'h5fcreate error'  !TODO: refine msg
      return
    end if

    this%vtk_id = H5Gcreate(this%file_id, 'VTKHDF')
    INSIST(this%vtk_id > 0)

    call h5_write_attr(this%vtk_id, 'Version', vtkhdf_version, stat, errmsg)
    INSIST(stat == 0)
    call h5_write_attr(this%vtk_id, 'Type', 'UnstructuredGrid', stat, errmsg)
    INSIST(stat == 0)

    this%cgrp_id = H5Gcreate(this%file_id, 'VTKHDF/CellData')
    INSIST(this%cgrp_id > 0)

    this%pgrp_id = H5Gcreate(this%file_id, 'VTKHDF/PointData')
    INSIST(this%pgrp_id > 0)

    if (present(temporal)) this%temporal = temporal
    if (this%temporal) then

      this%steps_id = H5Gcreate(this%file_id, 'VTKHDF/Steps')
      INSIST(this%steps_id > 0)

      this%cogrp_id = H5Gcreate(this%file_id, 'VTKHDF/Steps/CellDataOffsets')
      INSIST(this%cogrp_id > 0)

      this%pogrp_id = H5Gcreate(this%file_id, 'VTKHDF/Steps/PointDataOffsets')
      INSIST(this%pogrp_id > 0)

      this%nsteps = 0
      call h5_write_attr(this%steps_id, 'NSteps', this%nsteps, stat, errmsg)
      INSIST(stat == 0)

      associate (imold => [1], rmold => [1.0_real64], chunk_size => 100)
        call h5_create_unlimited_dataset(this%steps_id, 'Values', rmold, chunk_size, stat, errmsg)
        INSIST(stat == 0)
        call h5_create_unlimited_dataset(this%steps_id, 'PointOffsets', imold, chunk_size, stat, errmsg)
        INSIST(stat == 0)
        call h5_create_unlimited_dataset(this%steps_id, 'CellOffsets', imold, chunk_size, stat, errmsg)
        INSIST(stat == 0)
        call h5_create_unlimited_dataset(this%steps_id, 'ConnectivityIdOffsets', imold, chunk_size, stat, errmsg)
        INSIST(stat == 0)
        call h5_create_unlimited_dataset(this%steps_id, 'NumberOfParts', imold, chunk_size, stat, errmsg)
        INSIST(stat == 0)
        call h5_create_unlimited_dataset(this%steps_id, 'PartOffsets', imold, chunk_size, stat, errmsg)
        INSIST(stat == 0)
      end associate

    end if

  end subroutine

  subroutine close(this)
    class(vtkhdf_file), intent(inout) :: this
    integer :: istat
    type(temporal_dataset), pointer :: dset
    if (this%cogrp_id > 0) istat = H5Gclose(this%cogrp_id)
    if (this%pogrp_id > 0) istat = H5Gclose(this%pogrp_id)
    if (this%steps_id > 0) istat = H5Gclose(this%steps_id)
    if (this%cgrp_id > 0) istat = H5Gclose(this%cgrp_id)
    if (this%pgrp_id > 0) istat = H5Gclose(this%pgrp_id)
    if (this%vtk_id > 0) istat = H5Gclose(this%vtk_id)
    if (this%file_id > 0) istat = H5Fclose(this%file_id)
    if (associated(this%temporal_cell_dsets)) deallocate(this%temporal_cell_dsets)
    if (associated(this%temporal_point_dsets)) deallocate(this%temporal_point_dsets)
    call default_initialize(this)
  contains
    subroutine default_initialize(this)
      class(vtkhdf_file), intent(out) :: this
    end subroutine
  end subroutine

  recursive subroutine temporal_dataset_delete(this)
    type(temporal_dataset), intent(inout) :: this
    if (associated(this%next)) deallocate(this%next)
  end subroutine

  !! Write the unstructured mesh to the file. The mesh is described in the
  !! conventional manner by the X, CNODE, and XCNODE arrays. The additional
  !! array TYPES provides an unambiguous specification of the cell types,
  !! that we would otherwise infer from the OFFSETS array. This procedure
  !! must be called before any of the following procedures.

  subroutine write_mesh(this, x, cnode, xcnode, types, stat, errmsg)

    class(vtkhdf_file), intent(inout) :: this
    real(real64), intent(in) :: x(:,:)
    integer, intent(in) :: cnode(:), xcnode(:)
    integer(int8), intent(in) :: types(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    INSIST(this%file_id > 0)

    this%nnode = size(x, dim=2)
    this%ncell = size(types)

    ASSERT(size(xcnode) == this%ncell+1)
    ASSERT(size(cnode) == xcnode(size(xcnode))-1)
    ASSERT(minval(cnode) >= 1)
    ASSERT(maxval(cnode) <= this%nnode)

    call h5_write_dataset(this%vtk_id, 'NumberOfPoints', [this%nnode], stat, errmsg)
    INSIST(stat == 0)
    call h5_write_dataset(this%vtk_id, 'NumberOfCells',  [this%ncell], stat, errmsg)
    INSIST(stat == 0)
    call h5_write_dataset(this%vtk_id, 'NumberOfConnectivityIds', [size(cnode)], stat, errmsg)
    INSIST(stat == 0)

    call h5_write_dataset(this%vtk_id, 'Connectivity', cnode-1, stat, errmsg)  ! 0-based indexing
    INSIST(stat == 0)
    call h5_write_dataset(this%vtk_id, 'Offsets', xcnode-1, stat, errmsg) ! offsets instead of starting indices
    INSIST(stat == 0)
    call h5_write_dataset(this%vtk_id, 'Types', types, stat, errmsg)
    INSIST(stat == 0)
    call h5_write_dataset(this%vtk_id, 'Points', x, stat, errmsg)
    INSIST(stat == 0)

  end subroutine

  !! Writes the cell-based data ARRAY to a new named cell dataset. Scalar,
  !! vector, and tensor cell-based data are supported. In the case of a
  !! temporal file supporting time-dependent datasets, this dataset is
  !! static and not associated with any time step.

  subroutine write_cell_dataset_real64(this, name, array, stat, errmsg)
    class(vtkhdf_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, allocatable :: dims(:)
    dims = shape(array)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%ncell)
    call h5_write_dataset(this%cgrp_id, name, array, stat, errmsg)
  end subroutine

  !! Writes the point-based data ARRAY to a new named point dataset. Scalar,
  !! vector, and tensor point-based data are supported. In the case of a
  !! temporal file supporting time-dependent datasets, this dataset is
  !! static and not associated with any time step.

  subroutine write_point_dataset_real64(this, name, array, stat, errmsg)
    class(vtkhdf_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, allocatable :: dims(:)
    dims = shape(array)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%nnode)
    call h5_write_dataset(this%pgrp_id, name, array, stat, errmsg)
  end subroutine

  !! Register the specified name as a time-dependent point dataset. This writes
  !! no data, but only configures some necessary internal metadata. The MOLD
  !  array argument shall have the same type, kind, and rank as the actual
  !! dataset, and the same extent in all but the last dimension, but the array
  !! values themselves are not accessed. Scalar, vector, and tensor-valued mesh
  !! data are supported (rank-1, 2, and 3 MOLD).

  subroutine register_temporal_point_dataset_real64(this, name, mold, stat, errmsg)

    class(vtkhdf_file), intent(inout) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: mold(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(temporal_dataset), pointer :: new

    INSIST(this%nsteps == 0)

    allocate(new)
    new%name = name
    associate (chunk_size => this%nnode)
      call h5_create_unlimited_dataset(this%pgrp_id, name, mold, chunk_size, stat, errmsg)
    end associate
    new%next => this%temporal_point_dsets
    this%temporal_point_dsets => new

    associate (mold => [1], chunk_size => 100)
      call h5_create_unlimited_dataset(this%pogrp_id, name, mold, chunk_size, stat, errmsg)
    end associate

  end subroutine

  !! Register the specified name as a time-dependent cell dataset. This writes
  !! no data, but only configures some necessary internal metadata. The MOLD
  !  array argument shall have the same type, kind, and rank as the actual
  !! dataset, and the same extent in all but the last dimension, but the array
  !! values themselves are not accessed. Scalar, vector, and tensor-valued mesh
  !! data are supported (rank-1, 2, and 3 MOLD).

  subroutine register_temporal_cell_dataset_real64(this, name, mold, stat, errmsg)

    class(vtkhdf_file), intent(inout) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: mold(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(temporal_dataset), pointer :: new

    INSIST(this%nsteps == 0)

    allocate(new)
    new%name = name
    associate (chunk_size => this%ncell)
      call h5_create_unlimited_dataset(this%cgrp_id, name, mold, chunk_size, stat, errmsg)
    end associate
    new%next => this%temporal_cell_dsets
    this%temporal_cell_dsets => new

    associate (mold => [1], chunk_size => 100)
      call h5_create_unlimited_dataset(this%cogrp_id, name, mold, chunk_size, stat, errmsg)
    end associate

  end subroutine

  !! Mark the start of a new time step with time value TIME. Subsequent output
  !! of time-dependent datasets will be associated with this time step.

  subroutine write_time_step(this, time)

    class(vtkhdf_file), intent(inout) :: this
    real(real64), intent(in) :: time

    type(temporal_dataset), pointer :: tmp
    integer :: stat
    character(:), allocatable :: errmsg

    INSIST(this%nsteps >= 0)

    this%nsteps = this%nsteps + 1
    call h5_write_attr(this%steps_id, 'NSteps', this%nsteps, stat, errmsg)
    INSIST(stat == 0)

    !! A single mesh is used for all time steps so there are no offsets
    call h5_append_to_dataset(this%steps_id, 'Values', time, stat, errmsg)
    INSIST(stat == 0)
    call h5_append_to_dataset(this%steps_id, 'PointOffsets', 0, stat, errmsg)
    INSIST(stat == 0)
    call h5_append_to_dataset(this%steps_id, 'CellOffsets', 0, stat, errmsg)
    INSIST(stat == 0)
    call h5_append_to_dataset(this%steps_id, 'ConnectivityIdOffsets', 0, stat, errmsg)
    INSIST(stat == 0)
    call h5_append_to_dataset(this%steps_id, 'NumberOfParts', 1, stat, errmsg)
    INSIST(stat == 0)
    call h5_append_to_dataset(this%steps_id, 'PartOffsets', 0, stat, errmsg)
    INSIST(stat == 0)

    !! Set default offsets into point and cell datasets for each of the temporal
    !! datasets. The default points to the dataset for the previous time step
    !! (except initially). This will be overwritten when the dataset is written
    !! to for this time step.

    tmp => this%temporal_point_dsets
    do while (associated(tmp))
      tmp%flag = .false.  ! dataset not yet written for this time step
      call h5_append_to_dataset(this%pogrp_id, tmp%name, tmp%next_offset, stat, errmsg)
      INSIST(stat == 0)
      tmp => tmp%next
    end do

    tmp => this%temporal_cell_dsets
    do while (associated(tmp))
      tmp%flag = .false.  ! dataset not yet written for this time step
      call h5_append_to_dataset(this%cogrp_id, tmp%name, tmp%next_offset, stat, errmsg)
      INSIST(stat == 0)
      tmp => tmp%next
    end do

  end subroutine write_time_step

  !! Write the cell-based data ARRAY for the named time-dependent cell dataset.
  !! The data is associated with the current time step.

  subroutine write_temporal_cell_dataset_real64(this, name, array, stat, errmsg)

    class(vtkhdf_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(..)
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg

    integer, allocatable :: dims(:)
    type(temporal_dataset), pointer :: dset

    INSIST(this%nsteps >= 0)

    dims = shape(array)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%ncell)

    dset => this%temporal_cell_dsets
    do while (associated(dset))
      if (dset%name == name) exit
      dset => dset%next
    end do

    if (associated(dset)) then
      if (dset%flag) then
        stat = 1
        errmsg = 'dataset "' // name // '" already written for this time step'
        return
      else
        call h5_append_to_dataset(this%cgrp_id, name, array, stat, errmsg)
        if (stat /= 0) return
        call h5_write_dataset_element(this%cogrp_id, name, this%nsteps, dset%next_offset, stat, errmsg)
        if (stat /= 0) return
        dset%next_offset = dset%next_offset + this%ncell
        dset%flag = .true. ! dataset has been written for this time step
      end if
    else
      stat = 1
      errmsg = '"' // name // '" is not a valid temporal cell dataset'
      return
    end if

  end subroutine

  !! Write the node-based data ARRAY for the named time-dependent cell dataset.
  !! The data is associated with the current time step.

  subroutine write_temporal_point_dataset_real64(this, name, array, stat, errmsg)

    class(vtkhdf_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(..)
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg

    integer, allocatable :: dims(:)
    type(temporal_dataset), pointer :: dset

    INSIST(this%nsteps >= 0)

    dims = shape(array)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%nnode)

    dset => this%temporal_point_dsets
    do while (associated(dset))
      if (dset%name == name) exit
      dset => dset%next
    end do

    if (associated(dset)) then
      if (dset%flag) then
        stat = 1
        errmsg = 'dataset "' // name // '" already written for this time step'
        return
      else
        call h5_append_to_dataset(this%pgrp_id, name, array, stat, errmsg)
        if (stat /= 0) return
        call h5_write_dataset_element(this%pogrp_id, name, this%nsteps, dset%next_offset, stat, errmsg)
        if (stat /= 0) return
        dset%next_offset = dset%next_offset + this%nnode
        dset%flag = .true. ! dataset has been written for this time step
      end if
    else
      stat = 1
      errmsg = '"' // name // '" is not a valid temporal point dataset'
      return
    end if

  end subroutine

end module vtkhdf_file_type
