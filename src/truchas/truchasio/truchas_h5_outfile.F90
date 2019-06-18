!!
!! TRUCHAS_H5_OUTFILE
!!
!! This module defines the high-level interface for generating the Truchas
!! HDF5 output file.
!!
!! Ondrej Certik <certik@lanl.gov>
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2017
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! The Truchas HDF5 output file has the following layout, as originally
!! established by the Danu library:
!!
!!  /Meshes/<mesh1_name>
!!         /<mesh2_name>
!!          ...
!!  /Simulations/<sim1_name>/Mesh (a link to one of the preceding mesh groups)
!!                          /Non-series Data
!!                          /Probes
!!                          /Series Data/Series 1
!!                                      /Series 2
!!                                       ...
!!              /<sim2_name>/Mesh
!!                          /Non-series Data
!!                          /Probes
!!                          /Series Data
!!               ...
!!
!! The data hierarchy is represented by objects of a number of derived types
!! defined by this module: TH5_FILE, TH5_MESH_GROUP, TH5_SIM_GROUP,
!! TH5_SEQ_GROUP, and TH5_PROBE.
!!
!! An instance of the TH5_FILE type describes the HDF5 file, with methods to
!! open and close the file, and methods to add a named mesh group and a named
!! simulation group (with subtree).
!!
!! An instance of the TH5_MESH_GROUP type describes a mesh group, with methods
!! to write mesh data into the group.
!!
!! An instance of the TH5_SIM_GROUP type describes a simulation group. It has
!! a method for creating the "Mesh" link and writing time-independent data to
!! the "Non-series Data" group. It also has methods for creating probe datasets
!! within the "Probes" group, and "Series <n>" groups within the "Series Data"
!! group.
!!
!! An instance of the TH5_SEQ_GROUP type describes a "Series <n>" group, with
!! methods for writing time-dependent data to the group.
!!
!! An instance of the TH5_PROBE type describes a probe dataset within the
!! "Probes" group, with methods for incrementally writing probe data to the
!! dataset
!!
!! IMPLEMENTATION NOTES
!!
!! The current organization of the code largely mirrors the Danu API in order
!! to minimize required changes to Truchas.  Improvements should be made when
!! a further refactoring of the Truchas output driver is undertaken.
!!
!! The backend Scorpio library really only knows how to write distributed arrays
!! in parallel.  Truchas needs to write some data that is replicated on each
!! process.  We handle that by presenting Scorpio a distributed array with zero
!! sizes on all but the Truchas IO process (rank 0), and the actual data on the
!! IO process.  A Scorpio extension that only writes from rank 0 should be
!! possible and preferable (FIXME).
!!
!! A similar situation holds for probe data output.  That data is currently
!! replicated on all processes (FIXME).  Our own Scorpio extension also does
!! a parallel write (zero-sized on ranks other than 0) but without doing any
!! (needless) mpi gather calls.  Really should be providing the probe data
!! on the single process where it exists and let Scorpio migrate it to the
!! process doing io for the group and then write in serial (FIXME).
!!

#include "f90_assert.fpp"

module truchas_h5_outfile

  use,intrinsic :: iso_fortran_env, only: int8, int32, int64, real64
  use pgslib_module, only: broadcast  => PGSLib_BCast
  use scorpio_file_type
  implicit none
  private

  type, public :: th5_file
    private
    type(scorpio_file) :: file
    logical :: is_IOP
  contains
    procedure :: open  => th5_file_open
    procedure :: close => th5_file_close
    procedure :: add_unstr_mesh_group
    procedure :: add_interface_mesh_group
    procedure :: add_sim_group
  end type th5_file


  type, public :: th5_mesh_group
    private
    type(scorpio_file), pointer :: file => null() ! reference only -- do not own
    character(:), allocatable :: path !
  contains
    procedure :: write_coordinates
    procedure :: write_connectivity
    procedure, private :: mesh_write_attr_real64
    procedure, private :: mesh_write_attr_int32    
    generic, public :: write_attr => mesh_write_attr_real64, mesh_write_attr_int32
  end type th5_mesh_group


  type, public :: th5_sim_group
    private
    type(scorpio_file), pointer :: file => null() ! reference only -- do not own
    character(:), allocatable :: path ! simulation group path
    integer(int64) :: groupid = -1 ! simulation group hdf5 id
    integer :: seqno = 0  ! sequence counter
    logical :: is_IOP
  contains
    private
    procedure, public :: add_mesh_link
    procedure, public :: next_seq_group
    procedure, public :: create_probe
    generic, public :: write_attr => sim_write_attr_int32, sim_write_attr_real64
    generic, public :: write_dist_array => sim_write_dist_array_int32_r1, &
        sim_write_dist_array_real64_r1, sim_write_dist_array_real64_r2
    generic, public :: write_repl_data => sim_write_repl_data_int32_r0, &
        sim_write_repl_data_int32_r1, sim_write_repl_data_reall64_r0, &
        sim_write_repl_data_reall64_r1, sim_write_repl_data_reall64_r2
    procedure, private :: sim_write_repl_data_int32_r0
    procedure, private :: sim_write_repl_data_int32_r1
    procedure, private :: sim_write_repl_data_reall64_r0
    procedure, private :: sim_write_repl_data_reall64_r1
    procedure, private :: sim_write_repl_data_reall64_r2
    procedure, private :: sim_write_dist_array_int32_r1
    procedure, private :: sim_write_dist_array_real64_r1
    procedure, private :: sim_write_dist_array_real64_r2
    procedure, private :: sim_write_attr_int32
    procedure, private :: sim_write_attr_real64
  end type th5_sim_group


  type, public :: th5_seq_group
    private
    type(scorpio_file), pointer :: file => null() ! reference only -- do not own
    character(:), allocatable :: path ! path of the series group in the HDF5 file
  contains
    private
    generic, public :: write_dist_array => &
        seq_write_dist_array_int8_r2, seq_write_dist_array_int32_r1, &
        seq_write_dist_array_real64_r1, seq_write_dist_array_real64_r2
    generic, public :: write_attr => &
        seq_write_attr_int32, seq_write_attr_real64, seq_write_attr_real64_r1
    generic, public :: write_dataset_attr => &
        seq_write_dataset_attr_int32, seq_write_dataset_attr_real64, seq_write_dataset_attr_string
    procedure :: seq_write_dist_array_int8_r2
    procedure :: seq_write_dist_array_int32_r1
    procedure :: seq_write_dist_array_real64_r1
    procedure :: seq_write_dist_array_real64_r2
    procedure :: seq_write_attr_int32
    procedure :: seq_write_attr_real64
    procedure :: seq_write_attr_real64_r1
    procedure :: seq_write_dataset_attr_int32
    procedure :: seq_write_dataset_attr_real64
    procedure :: seq_write_dataset_attr_string
  end type th5_seq_group


  type, public :: th5_probe
    private
    type(scorpio_file), pointer :: file => null() ! reference only -- do not own
    character(:), allocatable :: path ! dataset path in the HDF5 file
    integer(int64) :: datasetid  ! HDF5 dataset ID
  contains
    private
    generic, public :: write_data => probe_write_data_real64_r2
    generic, public :: write_attr => &
        probe_write_attr_int32, probe_write_attr_real64, probe_write_attr_string
    procedure :: probe_write_data_real64_r2
    procedure :: probe_write_attr_int32
    procedure :: probe_write_attr_real64
    procedure :: probe_write_attr_string
   end type th5_probe

contains

!!!! TH5_FILE TYPE BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine th5_file_open(this, filename, io_group_size, is_IOP)
    class(th5_file), intent(out) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: io_group_size
    logical, intent(in) :: is_IOP
    this%is_IOP = is_IOP
    call this%file%open_file(filename, io_group_size)
  end subroutine th5_file_open

  subroutine th5_file_close(this)
    class(th5_file), intent(inout) :: this
    call this%file%close_file()
  end subroutine th5_file_close

  subroutine add_sim_group(this, name, sim)
    class(th5_file), target, intent(in) :: this
    character(*), intent(in) :: name
    class(th5_sim_group), intent(out) :: sim
    integer(int64) :: gid
    sim%file => this%file
    sim%is_IOP = this%is_IOP
    sim%path = '/Simulations/' // name
    sim%groupid = this%file%create_group(sim%path)
    !! Add the standard subroups
    gid = this%file%create_group(sim%path // '/Non-series Data')
    call this%file%close_group(gid)
    gid = this%file%create_group(sim%path // '/Series Data')
    call this%file%close_group(gid)
    gid = this%file%create_group(sim%path // '/Probes')
    call this%file%close_group(gid)
  end subroutine add_sim_group

  subroutine add_unstr_mesh_group(this, name, elem_order, dim, mesh)
    class(th5_file), target, intent(in) :: this
    character(*), intent(in)  :: name
    integer, intent(in) :: elem_order, dim
    class(th5_mesh_group), intent(out) :: mesh
    integer(int64) :: gid
    mesh%file => this%file
    mesh%path = '/Meshes/' // name
    gid = this%file%create_group(mesh%path)
    call this%file%close_group(gid)
    !! Mesh group attributes
    call this%file%write_attr(mesh%path, 'Dimension', dim)
    call this%file%write_attr(mesh%path, 'Mesh Type', 'UNSTRUCTURED')
    call this%file%write_attr(mesh%path, 'Element Type', 'HEX')
    call this%file%write_attr(mesh%path, 'Element Order', 8)
  end subroutine add_unstr_mesh_group

  subroutine add_interface_mesh_group(this, name, elem_order, dim, mesh)
    class(th5_file), target, intent(in) :: this
    character(*), intent(in)  :: name
    integer, intent(in) :: elem_order, dim
    class(th5_mesh_group), intent(out) :: mesh
    integer(int64) :: gid
    mesh%file => this%file
    mesh%path = '/Meshes/' // name
    gid = this%file%create_group(mesh%path)
    call this%file%close_group(gid)
    !! Mesh group attributes
    call this%file%write_attr(mesh%path, 'Dimension', dim)
    call this%file%write_attr(mesh%path, 'Mesh Type', 'UNSTRUCTURED')
    call this%file%write_attr(mesh%path, 'Element Type', 'POLYGON')
  end subroutine add_interface_mesh_group

!!!! TH5_MESH_GROUP TYPE BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_coordinates(this, nnode, x)
    class(th5_mesh_group), intent(in) :: this
    integer, intent(in) :: nnode
    real(real64), intent(in) :: x(:,:)
    call this%file%write_dataset(this%path//'/Nodal Coordinates', x, nnode)
    call this%file%write_attr(this%path, 'Number of Nodes', nnode)
  end subroutine write_coordinates

  subroutine write_connectivity(this, ncell, connect)
    class(th5_mesh_group), intent(in) :: this
    integer, intent(in) :: ncell
    integer, intent(in) :: connect(:,:)
    character(:), allocatable :: dataset
    dataset = this%path // '/Element Connectivity'
    print*,dataset
    call this%file%write_dataset(dataset, connect, ncell)
    if(ncell > 0) then
      call this%file%write_attr(dataset, 'Offset', 1) ! who uses this?
    end if 
    call this%file%write_attr(this%path, 'Number of Elements', ncell)
  end subroutine write_connectivity

  subroutine mesh_write_attr_int32(this, name, value)
    class(th5_mesh_group), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: value
    call this%file%write_attr(this%path, name, value)
  end subroutine mesh_write_attr_int32
  
  subroutine mesh_write_attr_real64(this, name, value)
    class(th5_mesh_group), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: value
    call this%file%write_attr(this%path, name, value)
  end subroutine mesh_write_attr_real64

  
  subroutine add_mesh_link(this, mesh_name)
    class(th5_sim_group), intent(in) :: this
    character(*), intent(in) :: mesh_name
    call this%file%create_link('/Meshes/'//mesh_name, this%groupid, 'Mesh')
  end subroutine

  subroutine sim_write_repl_data_int32_r0(this, name, ldata)
    class(th5_sim_group), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: ldata
    call sim_write_repl_data_int32_r1(this, name, [ldata])
  end subroutine sim_write_repl_data_int32_r0

  subroutine sim_write_repl_data_int32_r1(this, name, ldata)
    class(th5_sim_group), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: ldata(:)
    integer :: len1, glen
    character(:), allocatable :: dataset
    dataset = this%path // '/Non-series Data/' // name
    glen = size(ldata)
    call broadcast(glen)
    len1 = merge(glen, 0, this%is_IOP)
    call this%file%write_dataset(dataset, ldata(:len1), glen)
  end subroutine sim_write_repl_data_int32_r1

  subroutine sim_write_repl_data_reall64_r0(this, name, ldata)
    class(th5_sim_group), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: ldata
    call sim_write_repl_data_reall64_r1(this, name, [ldata])
  end subroutine sim_write_repl_data_reall64_r0

  subroutine sim_write_repl_data_reall64_r1(this, name, ldata)
    class(th5_sim_group), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: ldata(:)
    integer :: len1, glen
    character(:), allocatable :: dataset
    dataset = this%path // '/Non-series Data/' // name
    glen = size(ldata)
    call broadcast(glen)
    len1 = merge(glen, 0, this%is_IOP)
    call this%file%write_dataset(dataset, ldata(:len1), glen)
  end subroutine sim_write_repl_data_reall64_r1

  subroutine sim_write_repl_data_reall64_r2(this, name, ldata)
    class(th5_sim_group), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: ldata(:,:)
    integer :: len1, len2, glen
    character(:), allocatable :: dataset
    dataset = this%path // '/Non-series Data/' // name
    len1 = size(ldata,dim=1)
    call broadcast(len1)
    glen = size(ldata,dim=2)
    call broadcast(glen)
    len2 = merge(glen, 0, this%is_IOP)
    call this%file%write_dataset(dataset, ldata(:,:len2), glen)
  end subroutine sim_write_repl_data_reall64_r2

  subroutine sim_write_dist_array_int32_r1(this, name, ldata, glen)
    class(th5_sim_group), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: ldata(:)
    integer, intent(in) :: glen
    character(:), allocatable :: path
    path = this%path // '/Non-series Data/' // name
    call this%file%write_dataset(path, ldata, glen)
  end subroutine sim_write_dist_array_int32_r1

  subroutine sim_write_dist_array_real64_r1(this, name, ldata, glen)
    class(th5_sim_group), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: ldata(:)
    integer, intent(in) :: glen
    character(:), allocatable :: path
    path = this%path // '/Non-series Data/' // name
    call this%file%write_dataset(path, ldata, glen)
  end subroutine sim_write_dist_array_real64_r1

  subroutine sim_write_dist_array_real64_r2(this, name, ldata, glen)
    class(th5_sim_group), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: ldata(:,:)
    integer, intent(in) :: glen
    character(:), allocatable :: path
    path = this%path // '/Non-series Data/' // name
    call this%file%write_dataset(path, ldata, glen)
  end subroutine sim_write_dist_array_real64_r2

  subroutine next_seq_group(this, cyc, time, seq)
    class(th5_sim_group), intent(inout) :: this
    integer, intent(in) :: cyc
    real(real64), intent(in) :: time
    class(th5_seq_group), intent(out) :: seq
    integer(int64) :: gid
    character(16) :: name
    seq%file => this%file
    this%seqno = this%seqno + 1
    write(name,'(a,1x,i0)') '/Series', this%seqno
    seq%path = this%path // '/Series Data' // trim(name)
    gid = this%file%create_group(seq%path)
    call this%file%close_group(gid)
    call this%file%write_attr(seq%path, 'cycle', cyc)
    call this%file%write_attr(seq%path, 'sequence number', this%seqno)
    call this%file%write_attr(seq%path, 'time', time)
  end subroutine next_seq_group

  subroutine create_probe(this, probe_name, probe_data, probe)
    class(th5_sim_group), intent(in) :: this
    character(*), intent(in) :: probe_name
    real(real64), intent(in) :: probe_data(:,:)
    class(th5_probe), intent(out) :: probe
    probe%file => this%file
    probe%path = this%path // '/Probes/' // probe_name
    probe%datasetid = this%file%create_probe(this%groupid, probe_name, probe_data)
  end subroutine create_probe

  subroutine sim_write_attr_int32(this, name, value)
    class(th5_sim_group), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: value
    call this%file%write_attr(this%path, name, value)
  end subroutine sim_write_attr_int32

  subroutine sim_write_attr_real64(this, name, value)
    class(th5_sim_group), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: value
    call this%file%write_attr(this%path, name, value)
  end subroutine sim_write_attr_real64

!!!! TH5_SEQ_GROUP TYPE BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine seq_write_dist_array_int8_r2(this, data_name, glen, ldata)
    class(th5_seq_group), intent(in) :: this
    character(*), intent(in) :: data_name
    integer, intent(in) :: glen  ! global length
    integer(int8), intent(in) :: ldata(:,:)
    ! FIXME: `data_name` is not always trimmed, so we need to trim it here.
    ! However the right fix is to ensure the caller passes a trimmed string.
    call this%file%write_dataset(this%path//'/'//trim(data_name), ldata, glen)
  end subroutine seq_write_dist_array_int8_r2

  subroutine seq_write_dist_array_int32_r1(this, data_name, glen, ldata)
    class(th5_seq_group), intent(in) :: this
    character(*), intent(in) :: data_name
    integer, intent(in) :: glen  ! global length
    integer(int32), intent(in) :: ldata(:)
    ! FIXME: `data_name` is not always trimmed, so we need to trim it here.
    ! However the right fix is to ensure the caller passes a trimmed string.
    call this%file%write_dataset(this%path//'/'//trim(data_name), ldata, glen)
  end subroutine seq_write_dist_array_int32_r1

  subroutine seq_write_dist_array_real64_r1(this, data_name, glen, ldata)
    class(th5_seq_group), intent(in) :: this
    character(*), intent(in) :: data_name
    integer, intent(in) :: glen  ! global length
    real(real64), intent(in) :: ldata(:) ! local data
    ! FIXME: `data_name` is not always trimmed, so we need to trim it here.
    ! However the right fix is to ensure the caller passes a trimmed string.
    call this%file%write_dataset(this%path//'/'//trim(data_name), ldata, glen)
  end subroutine seq_write_dist_array_real64_r1

  subroutine seq_write_dist_array_real64_r2(this, data_name, glen, ldata)
    class(th5_seq_group), intent(in) :: this
    character(*), intent(in) :: data_name
    integer, intent(in) :: glen  ! global length
    real(real64), intent(in) :: ldata(:,:)
    ! FIXME: `data_name` is not always trimmed, so we need to trim it here.
    ! However the right fix is to ensure the caller passes a trimmed string.
    call this%file%write_dataset(this%path//'/'//trim(data_name), ldata, glen)
  end subroutine seq_write_dist_array_real64_r2

  subroutine seq_write_attr_int32(this, attr_name, attr_value)
    class(th5_seq_group), intent(in) :: this
    character(*), intent(in) :: attr_name
    integer(int32), intent(in) :: attr_value
    call this%file%write_attr(this%path, attr_name, attr_value)
  end subroutine seq_write_attr_int32

  subroutine seq_write_attr_real64(this, attr_name, attr_value)
    class(th5_seq_group), intent(in) :: this
    character(*), intent(in) :: attr_name
    real(real64), intent(in) :: attr_value
    call this%file%write_attr(this%path, attr_name, attr_value)
  end subroutine seq_write_attr_real64

  subroutine seq_write_attr_real64_r1(this, attr_name, attr_value)
    class(th5_seq_group), intent(in) :: this
    character(*), intent(in) :: attr_name
    real(real64), intent(in) :: attr_value(:)
    call this%file%write_attr(this%path, attr_name, attr_value)
  end subroutine seq_write_attr_real64_r1

  subroutine seq_write_dataset_attr_int32(this, dataset, attr_name, attr_value)
    class(th5_seq_group), intent(in) :: this
    character(*), intent(in) :: dataset, attr_name
    integer(int32), intent(in) :: attr_value
    call this%file%write_attr(this%path // '/' // dataset, attr_name, attr_value)
  end subroutine seq_write_dataset_attr_int32

  subroutine seq_write_dataset_attr_real64(this, dataset, attr_name, attr_value)
    class(th5_seq_group), intent(in) :: this
    character(*), intent(in) :: dataset, attr_name
    real(real64), intent(in) :: attr_value
    call this%file%write_attr(this%path // '/' // dataset, attr_name, attr_value)
  end subroutine seq_write_dataset_attr_real64

  subroutine seq_write_dataset_attr_string(this, dataset, attr_name, attr_value)
    class(th5_seq_group), intent(in) :: this
    character(*), intent(in) :: dataset, attr_name
    character(*), intent(in) :: attr_value
    call this%file%write_attr(this%path // '/' // dataset, attr_name, attr_value)
  end subroutine seq_write_dataset_attr_string

!!!! TH5_PROBE TYPE BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine probe_write_data_real64_r2(this, probe_data)
    class(th5_probe), intent(in) :: this
    real(real64), intent(in), contiguous :: probe_data(:,:)
    call this%file%write_probe(this%datasetid, probe_data)
  end subroutine probe_write_data_real64_r2

  subroutine probe_write_attr_int32(this, attr_name, attr_value)
    class(th5_probe), intent(in) :: this
    character(*), intent(in) :: attr_name
    integer(int32), intent(in) :: attr_value
    call this%file%write_attr(this%path, attr_name, attr_value)
  end subroutine probe_write_attr_int32

  subroutine probe_write_attr_real64(this, attr_name, attr_value)
    class(th5_probe), intent(in) :: this
    character(*), intent(in) :: attr_name
    real(real64), intent(in) :: attr_value
    call this%file%write_attr(this%path, attr_name, attr_value)
  end subroutine probe_write_attr_real64

  subroutine probe_write_attr_string(this, attr_name, attr_value)
    class(th5_probe), intent(in) :: this
    character(*), intent(in) :: attr_name
    character(*), intent(in) :: attr_value
    call this%file%write_attr(this%path, attr_name, attr_value)
  end subroutine probe_write_attr_string

end module truchas_h5_outfile
