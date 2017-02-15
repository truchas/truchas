!! TO-DO
!!
!! Propagate return codes from scorpio to the stat argumen
!! Use stat=0 for a successful return instead of DANU_SUCCESS (=1)
!! Add standard handling of an optional STAT and ERRMSG args

#include "f90_assert.fpp"

module truchasio
  use,intrinsic :: iso_fortran_env, only: int32, real64
  use,intrinsic :: iso_c_binding, only: c_ptr, c_char, c_double, c_int, c_int8_t, &
    c_null_ptr, c_associated, c_null_char, c_f_pointer
  use pgslib_module, only: broadcast  => PGSLib_BCast
  implicit none
  private

  integer, parameter, public :: DANU_SUCCESS = 1

  type, public :: output_file
    integer :: file_id  ! h5 file id
    type(c_ptr), private :: scorp_id = c_null_ptr ! Scorpio id
    logical :: is_IOP = .false.  ! True if this PE is the IO PE.
  contains
    procedure :: open  => output_file_open
    !procedure :: close => output_file_close
    procedure :: close => output_file_final
    procedure :: mesh_add_unstructured
    procedure :: simulation_add
    procedure :: simulation_link_mesh
    !final :: output_file_final
  end type


  type, public :: simulation
    type(output_file), pointer, private :: foutput => null()
    character(:), allocatable :: path
    integer :: group_id = -1
    integer :: serno = 0  ! series counter
  contains
    procedure, private :: simulation_write_repl_rank0_int32
    procedure, private :: simulation_write_repl_rank1_int32
    procedure, private :: simulation_write_repl_rank0_real64
    procedure, private :: simulation_write_repl_rank1_real64
    procedure, private :: simulation_write_repl_rank2_real64
    generic :: write_repl_data => simulation_write_repl_rank0_int32, &
                                  simulation_write_repl_rank1_int32, &
                                  simulation_write_repl_rank0_real64, &
                                  simulation_write_repl_rank1_real64, &
                                  simulation_write_repl_rank2_real64
    procedure, private :: simulation_write_dist_rank1_int32
    procedure, private :: simulation_write_dist_rank1_real64
    procedure, private :: simulation_write_dist_rank2_real64
    generic :: write_dist_array => simulation_write_dist_rank1_int32, &
                                   simulation_write_dist_rank1_real64, &
                                   simulation_write_dist_rank2_real64
    procedure, private :: simulation_write_attr_int32
    procedure, private :: simulation_write_attr_real64
    generic :: write_attr => simulation_write_attr_int32, simulation_write_attr_real64

    procedure :: sequence_next_id
    procedure :: probe_create_data
  end type


  type, public :: output_mesh
    type(output_file), pointer, private :: foutput => null()
    character(:), allocatable :: path
  contains
    procedure :: write_coordinates
    procedure :: write_connectivity
  end type


  type, public :: sequence
    type(output_file), pointer, private :: foutput => null()
    character(:), allocatable :: path
  contains
    generic :: data_write => &
      simulation_data_write_byte_rank2, &
      simulation_data_write_integer_rank1, &
      simulation_data_write_real8_rank1, &
      simulation_data_write_real8_rank2

    procedure :: sequence_write_attr_int32
    procedure :: sequence_write_attr_real64
    procedure :: sequence_write_attr_real64_rank1
    generic :: write_attr => &
        sequence_write_attr_int32, &
        sequence_write_attr_real64, &
        sequence_write_attr_real64_rank1

    procedure :: sequence_write_dataset_attr_int32
    procedure :: sequence_write_dataset_attr_real64
    procedure :: sequence_write_dataset_attr_string
    generic :: write_dataset_attr => &
        sequence_write_dataset_attr_int32, &
        sequence_write_dataset_attr_real64, &
        sequence_write_dataset_attr_string

    ! private
    procedure, private :: simulation_data_write_byte_rank2
    procedure, private :: simulation_data_write_integer_rank1
    procedure, private :: simulation_data_write_real8_rank1
    procedure, private :: simulation_data_write_real8_rank2
  end type


  type, public :: output_probe
    type(c_ptr), private :: pid = c_null_ptr ! probe id
    type(output_file), pointer, private :: foutput => null()
    character(:), allocatable :: path
    integer :: dataset_id
  contains

    procedure, private :: probe_write_attr_int32
    procedure, private :: probe_write_attr_real64
    procedure, private :: probe_write_attr_string
    generic :: write_attr => &
        probe_write_attr_int32, &
        probe_write_attr_real64, &
        probe_write_attr_string
    procedure :: data_write => probe_data_write
  end type

  interface
    type(c_ptr) function truchas_scorpio_create_handle(filename, &
        numIOgroups) bind(c)
    import :: c_ptr, c_char, c_int
    character(kind=c_char), intent(in)  :: filename(*)
    integer(c_int), intent(in), value :: numIOgroups
    end function

    subroutine truchas_scorpio_write_dataset_1d_integer(h, &
        name, &
        vector, global_dim, local_dim) bind(c)
    import :: c_ptr, c_char, c_int
    type(c_ptr), value :: h
    character(kind=c_char), intent(in)  :: name(*)
    integer(c_int), intent(in), value :: global_dim, local_dim
    integer(c_int), intent(in) :: vector(local_dim)
    end subroutine

    subroutine truchas_scorpio_write_dataset_2d_integer(h, &
        name, &
        vector, global_dim, local_dim1, local_dim2) bind(c)
    import :: c_ptr, c_char, c_int
    type(c_ptr), value :: h
    character(kind=c_char), intent(in)  :: name(*)
    integer(c_int), intent(in), value :: global_dim, local_dim1, local_dim2
    integer(c_int), intent(in) :: vector(local_dim1, local_dim2)
    end subroutine

    subroutine truchas_scorpio_write_dataset_1d_double(h, &
        name, &
        vector, global_dim, local_dim) bind(c)
    import :: c_ptr, c_char, c_int, c_double
    type(c_ptr), value :: h
    character(kind=c_char), intent(in)  :: name(*)
    integer(c_int), intent(in), value :: global_dim, local_dim
    real(c_double), intent(in) :: vector(local_dim)
    end subroutine

    subroutine truchas_scorpio_write_dataset_2d_double(h, &
        name, &
        vector, global_dim, local_dim1, local_dim2) bind(c)
    import :: c_ptr, c_char, c_int, c_double
    type(c_ptr), value :: h
    character(kind=c_char), intent(in)  :: name(*)
    integer(c_int), intent(in), value :: global_dim, local_dim1, local_dim2
    real(c_double), intent(in) :: vector(local_dim1, local_dim2)
    end subroutine

    subroutine truchas_scorpio_write_dataset_2d_byte(h, &
        name, &
        vector, global_dim, local_dim1, local_dim2) bind(c)
    import :: c_ptr, c_char, c_int, c_double, c_int8_t
    type(c_ptr), value :: h
    character(kind=c_char), intent(in)  :: name(*)
    integer(c_int), intent(in), value :: global_dim, local_dim1, local_dim2
    integer(c_int8_t), intent(in) :: vector(local_dim1, local_dim2)
    end subroutine

    subroutine truchas_scorpio_handle_file_close(h) bind(c)
    import :: c_ptr
    type(c_ptr), value :: h
    end subroutine

    subroutine truchas_scorpio_write_attr_0d_integer(h, attr_name, attr_data, obj_name) bind(c)
    import :: c_ptr, c_char, c_int
    type(c_ptr), value :: h
    character(kind=c_char), intent(in) :: attr_name(*), obj_name(*)
    integer(c_int), intent(in) :: attr_data
    end subroutine

    subroutine truchas_scorpio_write_attr_0d_double(h, attr_name, attr_data, obj_name) bind(c)
    import :: c_ptr, c_char, c_double
    type(c_ptr), value :: h
    character(kind=c_char), intent(in) :: attr_name(*), obj_name(*)
    real(c_double), intent(in) :: attr_data
    end subroutine

    subroutine truchas_scorpio_write_attr_1d_double(h, attr_name, attr_data, ndims, adims, obj_name) bind(c)
    import :: c_ptr, c_int, c_char, c_double
    type(c_ptr), value :: h
    character(kind=c_char), intent(in) :: attr_name(*), obj_name(*)
    real(c_double), intent(in) :: attr_data(*)
    integer(c_int), value :: ndims
    integer(c_int), intent(in) :: adims(*)
    end subroutine

    subroutine truchas_scorpio_write_attr_0d_string(h, attr_name, attr_data, obj_name) bind(c)
    import :: c_ptr, c_char
    type(c_ptr), value :: h
    character(kind=c_char), intent(in) :: attr_name(*), attr_data(*), obj_name(*)
    end subroutine

    function truchas_scorpio_create_dataset_group(h, group_name) result(gid) bind(c)
    import c_ptr, c_int, c_char
    type(c_ptr), value :: h
    character(kind=c_char), intent(in) :: group_name(*)
    integer(c_int) :: gid
    end function

    subroutine truchas_scorpio_close_dataset_group(h, gid) bind(c)
    import c_ptr, c_int
    type(c_ptr), value :: h
    integer(c_int), value :: gid
    end subroutine

    subroutine truchas_scorpio_create_link(h, target, link_loc_id, link_name) bind(c)
    import c_ptr, c_char, c_int
    type(c_ptr), value :: h
    character(kind=c_char), intent(in) :: target(*), link_name(*)
    integer(c_int), value :: link_loc_id
    end subroutine

    subroutine truchas_scorpio_write_probe_data_2d_double(h, pid, dims, data) bind(c)
    import c_ptr, c_int, c_double
    type(c_ptr), value :: h
    integer(c_int), value :: pid
    integer(c_int) :: dims(*)
    real(c_double), intent(in) :: data(*)
    end subroutine

    function truchas_scorpio_create_probe_2d_double(h, sid, name, dims, data) result(pid) bind(c)
    import c_ptr, c_char, c_int, c_double
    type(c_ptr), value :: h
    character(kind=c_char), intent(in) :: name(*)
    integer(c_int), value :: sid
    integer(c_int), intent(in) :: dims(*)
    real(c_double), intent(in) :: data(*)
    integer(c_int) :: pid
    end function
  end interface

contains

! output_file

  subroutine output_file_open(this, filename, is_IOP, stat)
    class(output_file), intent(out) :: this
    character(*), intent(in) :: filename
    logical, intent(in) :: is_IOP
    integer, intent(out), optional :: stat
    this%is_IOP = is_IOP
    this%scorp_id = truchas_scorpio_create_handle(filename // c_null_char, 2)
    stat = DANU_SUCCESS
  end subroutine

  subroutine output_file_final (this)
    class(output_file), intent(inout) :: this
    if (c_associated(this%scorp_id)) then
      call truchas_scorpio_handle_file_close(this%scorp_id)
      this%scorp_id = c_null_ptr
    end if
  end subroutine

  !! The passed object is intent(out) and this causes the object to be finalized
  !! and its components default initialized on entry to the subroutine.
  subroutine output_file_close (this)
    class(output_file), intent(out) :: this
  end subroutine

  subroutine simulation_add(this, sim_name, sim, stat)

    class(output_file), target, intent(in) :: this
    character(*,kind=C_CHAR),      intent(in)  :: sim_name
    class(simulation),                 intent(out) :: sim
    integer, intent(out), optional :: stat

    integer :: gid
    character(:), allocatable :: group_name

    sim%foutput => this
    sim%path = '/Simulations/' // sim_name

    group_name = '/Simulations/' // sim_name // c_null_char
    sim%group_id = truchas_scorpio_create_dataset_group(this%scorp_id, group_name)

    !! Danu created these groups, but it seems not to be necessary, as Scorpio
    !! automatically creates them on the fly when writing datasets.  The only
    !! possible issue is if users of the file expect to find them even if empty.
    !group_name = sim%path // '/Non-series Data' // c_null_char
    !gid = truchas_scorpio_create_dataset_group(this%scorp_id, group_name)
    !call truchas_scorpio_close_dataset_group(this%scorp_id, gid)
    !
    !group_name = sim%path // '/Series Data' // c_null_char
    !gid = truchas_scorpio_create_dataset_group(this%scorp_id, group_name)
    !call truchas_scorpio_close_dataset_group(this%scorp_id, gid)
    !
    group_name = sim%path // '/Probes' // c_null_char
    gid = truchas_scorpio_create_dataset_group(this%scorp_id, group_name)
    call truchas_scorpio_close_dataset_group(this%scorp_id, gid)

    stat = DANU_SUCCESS

  end subroutine

  subroutine simulation_link_mesh(this, sim, mesh_name, stat)
    class(output_file), intent(in) :: this
    class(simulation), intent(in) :: sim
    character(*), intent(in) :: mesh_name
    character(:), allocatable :: link_target
    integer, intent(out), optional :: stat
    link_target = '/Meshes/' // mesh_name // c_null_char
    call truchas_scorpio_create_link(this%scorp_id, link_target, sim%group_id, 'Mesh'//c_null_char)
    stat = DANU_SUCCESS
  end subroutine

  subroutine mesh_add_unstructured(this,mesh_name,elem_order,mesh_dim,out_mesh,stat)
    class(output_file), target, intent(in) :: this
    character(*,kind=C_CHAR), intent(in)  :: mesh_name
    integer(C_INT), intent(in) :: elem_order, mesh_dim
    class(output_mesh), intent(out) :: out_mesh
    integer, intent(out) :: stat
    integer :: gid
    character(:), allocatable :: group_name
    out_mesh%foutput => this
    out_mesh%path = '/Meshes/' // mesh_name
    group_name = out_mesh%path // c_null_char
    gid = truchas_scorpio_create_dataset_group(this%scorp_id, group_name)
    call truchas_scorpio_close_dataset_group(this%scorp_id, gid)
    call truchas_scorpio_write_attr_0d_integer(this%scorp_id, &
        'Dimension'//c_null_char, mesh_dim, group_name)
    call truchas_scorpio_write_attr_0d_string(this%scorp_id, &
        'Mesh Type'//c_null_char, 'UNSTRUCTURED'//c_null_char, group_name)
    call truchas_scorpio_write_attr_0d_string(this%scorp_id, &
        'Element Type'//c_null_char, 'HEX'//c_null_char, group_name)
    call truchas_scorpio_write_attr_0d_integer(this%scorp_id, &
        'Element Order'//c_null_char, 8, group_name)
    stat = DANU_SUCCESS
  end subroutine

! -------------------------------------------------------------

! simulation

  subroutine simulation_write_repl_rank0_int32(this, name, ldata, stat)
    class(simulation), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: ldata
    integer, intent(out) :: stat
    call simulation_write_repl_rank1_int32(this, name, [ldata], stat)
  end subroutine

  subroutine simulation_write_repl_rank1_int32(this, name, ldata, stat)
    class(simulation), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: ldata(:)
    integer, intent(out) :: stat
    integer :: len1, glen
    character(:), allocatable :: dataset
    dataset = this%path // '/Non-series Data/' // name // c_null_char
    glen = size(ldata)
    call broadcast(glen)
    len1 = merge(glen, 0, this%foutput%is_IOP)
    call truchas_scorpio_write_dataset_1d_integer(this%foutput%scorp_id, dataset, ldata, glen, len1)
    stat = DANU_SUCCESS
  end subroutine

  subroutine simulation_write_repl_rank0_real64(this, name, ldata, stat)
    class(simulation), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: ldata
    integer, intent(out) :: stat
    call simulation_write_repl_rank1_real64(this, name, [ldata], stat)
  end subroutine

  subroutine simulation_write_repl_rank1_real64(this, name, ldata, stat)
    class(simulation), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: ldata(:)
    integer, intent(out) :: stat
    integer :: len1, glen
    character(:), allocatable :: dataset
    dataset = this%path // '/Non-series Data/' // name // c_null_char
    glen = size(ldata)
    call broadcast(glen)
    len1 = merge(glen, 0, this%foutput%is_IOP)
    call truchas_scorpio_write_dataset_1d_double(this%foutput%scorp_id, dataset, ldata, glen, len1)
    stat = DANU_SUCCESS
  end subroutine

  subroutine simulation_write_repl_rank2_real64(this, name, ldata, stat)
    class(simulation), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: ldata(:,:)
    integer, intent(out) :: stat
    integer :: len1, len2, glen
    character(:), allocatable :: dataset
    dataset = this%path // '/Non-series Data/' // name // c_null_char
    len1 = size(ldata,dim=1)
    call broadcast(len1)
    glen = size(ldata,dim=2)
    call broadcast(glen)
    len2 = merge(glen, 0, this%foutput%is_IOP)
    call truchas_scorpio_write_dataset_2d_double(this%foutput%scorp_id, dataset, ldata, glen, len1, len2)
    stat = DANU_SUCCESS
  end subroutine

  subroutine simulation_write_dist_rank1_int32(this, name, ldata, glen, stat)
    class(simulation), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: ldata(:)
    integer, intent(in) :: glen
    integer, intent(out) :: stat
    character(:), allocatable :: dataset
    dataset = this%path // '/Non-series Data/' // name // c_null_char
    call truchas_scorpio_write_dataset_1d_integer(this%foutput%scorp_id, &
        dataset, ldata, glen, size(ldata))
    stat = DANU_SUCCESS
  end subroutine

  subroutine simulation_write_dist_rank1_real64(this, name, ldata, glen, stat)
    class(simulation), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: ldata(:)
    integer, intent(in) :: glen
    integer, intent(out) :: stat
    character(:), allocatable :: dataset
    dataset = this%path // '/Non-series Data/' // name // c_null_char
    call truchas_scorpio_write_dataset_1d_double(this%foutput%scorp_id, &
        dataset, ldata, glen, size(ldata))
    stat = DANU_SUCCESS
  end subroutine

  subroutine simulation_write_dist_rank2_real64(this, name, ldata, glen, stat)
    class(simulation), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: ldata(:,:)
    integer, intent(in) :: glen
    integer, intent(out) :: stat
    character(:), allocatable :: dataset
    dataset = this%path // '/Non-series Data/' // name // c_null_char
    call truchas_scorpio_write_dataset_2d_double(this%foutput%scorp_id, &
        dataset, ldata, glen, size(ldata,1), size(ldata,2))
    stat = DANU_SUCCESS
  end subroutine

  subroutine sequence_next_id(this,cyc,time,seq,stat)
    class(simulation),                 intent(inout)  :: this
    integer(kind=C_INT),               intent(in)  :: cyc
    real(kind=C_DOUBLE),               intent(in)  :: time
    class(sequence),                   intent(out) :: seq
    integer, intent(out), optional :: stat
    integer :: gid
    character(16) :: name
    character(:), allocatable :: group_name
    seq%foutput => this%foutput
    this%serno = this%serno + 1
    write(name,'(a,1x,i0)') '/Series', this%serno
    seq%path = this%path // '/Series Data' // trim(name)

    group_name = seq%path // c_null_char
    gid = truchas_scorpio_create_dataset_group(this%foutput%scorp_id, group_name)
    call truchas_scorpio_close_dataset_group(this%foutput%scorp_id, gid)
    call truchas_scorpio_write_attr_0d_integer(this%foutput%scorp_id, &
          'cycle'//c_null_char, cyc, group_name)
    call truchas_scorpio_write_attr_0d_integer(this%foutput%scorp_id, &
          'sequence number'//c_null_char, this%serno, group_name)
    call truchas_scorpio_write_attr_0d_double(this%foutput%scorp_id, &
          'time'//c_null_char, time, group_name)
    stat = DANU_SUCCESS
  end subroutine

  subroutine probe_create_data (this, probe_name, rdata, probe, stat)
    class(simulation), intent(in) :: this
    character(*), intent(in) :: probe_name
    real(real64), intent(in) :: rdata(:,:)
    class(output_probe), intent(out) :: probe
    integer, intent(out) :: stat
    integer :: pid
    probe%foutput => this%foutput
    probe%path = this%path // '/Probes/' // probe_name
    probe%dataset_id = truchas_scorpio_create_probe_2d_double(this%foutput%scorp_id, &
                           this%group_id, probe_name//c_null_char, shape(rdata), rdata)
    stat = DANU_SUCCESS
  end subroutine

  subroutine simulation_write_attr_int32(this, name, value, stat)
    class(simulation), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: value
    integer, intent(out) :: stat
    call truchas_scorpio_write_attr_0d_integer(this%foutput%scorp_id, &
        name//c_null_char, value, this%path//c_null_char)
    stat = DANU_SUCCESS
  end subroutine

  subroutine simulation_write_attr_real64(this, name, value, stat)
    class(simulation), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: value
    integer, intent(out) :: stat
    call truchas_scorpio_write_attr_0d_double(this%foutput%scorp_id, &
        name//c_null_char, value, this%path//c_null_char)
    stat = DANU_SUCCESS
  end subroutine

! -------------------------------------------------------------

! output_mesh

  subroutine write_coordinates(this, nnode, x, stat)
    class(output_mesh), intent(in) :: this
    integer(C_INT), intent(in) :: nnode
    real(C_DOUBLE), intent(in) :: x(:,:)
    integer, intent(out) :: stat
    character(:), allocatable :: dataset
    dataset = this%path // '/Nodal Coordinates' // c_null_char
    call truchas_scorpio_write_dataset_2d_double(this%foutput%scorp_id, &
        dataset, x, nnode, size(x,1), size(x,2))
    call truchas_scorpio_write_attr_0d_integer(this%foutput%scorp_id, &
        'Number of Nodes'//c_null_char, nnode, this%path//c_null_char)
    stat = DANU_SUCCESS  ! our C layer ignores the scorpio return code
  end subroutine write_coordinates

  subroutine write_connectivity(this, ncell, idata, stat)
    class(output_mesh), intent(in) :: this
    integer(C_INT), intent(in) :: ncell
    integer(C_INT), intent(in) :: idata(:,:)
    integer, intent(out) :: stat
    character(:), allocatable :: dataset
    dataset = this%path // '/Element Connectivity' // c_null_char
    call truchas_scorpio_write_dataset_2d_integer(this%foutput%scorp_id, &
        dataset, idata, ncell, size(idata,1), size(idata,2))
    !TODO: who uses this attribute?
    call truchas_scorpio_write_attr_0d_integer(this%foutput%scorp_id, &
        'Offset'//c_null_char, 1, dataset)
    call truchas_scorpio_write_attr_0d_integer(this%foutput%scorp_id, &
        'Number of Elements'//c_null_char, ncell, this%path//c_null_char)
    stat = DANU_SUCCESS  ! our C layer ignores the scorpio return code
  end subroutine


! -------------------------------------------------------------

! sequence

  subroutine simulation_data_write_byte_rank2(this, data_name, glen, ldata, stat)
  class(sequence), intent(in) :: this
  character(*), intent(in) :: data_name
  integer,           intent(in) :: glen  ! global length
  integer(C_INT8_T), intent(in) :: ldata(:,:)
  integer, intent(out), optional :: stat
  INSIST(c_associated(this%foutput%scorp_id))
  ! FIXME: `data_name` is not always trimmed, so we need to trim it here.
  ! However the right fix is to ensure the caller passes a trimmed string.
  call truchas_scorpio_write_dataset_2d_byte(this%foutput%scorp_id, &
    this%path // "/" // trim(data_name) // c_null_char, ldata, glen, &
    size(ldata, 1), size(ldata, 2))
  stat = DANU_SUCCESS
  end subroutine

  subroutine simulation_data_write_integer_rank1(this,data_name,glen,ldata,stat)
  class(sequence), intent(in)  :: this
  character(*), intent(in)  :: data_name
  integer,                           intent(in) :: glen  ! global length
  integer(kind=C_INT), dimension(:), intent(in) :: ldata
  integer, intent(out), optional :: stat
  INSIST(c_associated(this%foutput%scorp_id))
  ! FIXME: `data_name` is not always trimmed, so we need to trim it here.
  ! However the right fix is to ensure the caller passes a trimmed string.
  call truchas_scorpio_write_dataset_1d_integer(this%foutput%scorp_id, &
    this%path // "/" // trim(data_name) // c_null_char, &
    ldata, glen, size(ldata))
  stat = DANU_SUCCESS
  end subroutine

  subroutine simulation_data_write_real8_rank1(this,data_name,glen,ldata,stat)
  class(sequence),                    intent(in) :: this
  character(*),                   intent(in) :: data_name
  integer,                            intent(in) :: glen  ! global length
  real(kind=C_DOUBLE),  dimension(:), intent(in) :: ldata ! local data
  integer, intent(out), optional :: stat
  INSIST(c_associated(this%foutput%scorp_id))
  ! FIXME: `data_name` is not always trimmed, so we need to trim it here.
  ! However the right fix is to ensure the caller passes a trimmed string.
  call truchas_scorpio_write_dataset_1d_double(this%foutput%scorp_id, &
    this%path // "/" // trim(data_name) // c_null_char, &
    ldata, glen, size(ldata))
  stat = DANU_SUCCESS
  end subroutine

  subroutine simulation_data_write_real8_rank2(this,data_name,glen,ldata,stat)
  class(sequence),                     intent(in)  :: this
  character(*),                    intent(in)  :: data_name
  integer,                             intent(in) :: glen  ! global length
  real(kind=C_DOUBLE), dimension(:,:), intent(in) :: ldata
  integer, intent(out), optional :: stat
  INSIST(c_associated(this%foutput%scorp_id))
  ! FIXME: `data_name` is not always trimmed, so we need to trim it here.
  ! However the right fix is to ensure the caller passes a trimmed string.
  call truchas_scorpio_write_dataset_2d_double(this%foutput%scorp_id, &
    this%path // "/" // trim(data_name) // c_null_char, ldata, glen, &
    size(ldata, 1), size(ldata, 2))
  stat = DANU_SUCCESS
  end subroutine

  subroutine sequence_write_attr_int32(this, attr_name, attr_value, stat)
    class(sequence), intent(in) :: this
    character(*), intent(in) :: attr_name
    integer(int32), intent(in) :: attr_value
    integer, intent(out) :: stat
    call truchas_scorpio_write_attr_0d_integer(this%foutput%scorp_id, &
        attr_name//c_null_char, attr_value, this%path//c_null_char)
    stat = DANU_SUCCESS
  end subroutine

  subroutine sequence_write_attr_real64(this, attr_name, attr_value, stat)
    class(sequence), intent(in) :: this
    character(*), intent(in) :: attr_name
    real(real64), intent(in) :: attr_value
    integer, intent(out) :: stat
    call truchas_scorpio_write_attr_0d_double(this%foutput%scorp_id, &
        attr_name//c_null_char, attr_value, this%path//c_null_char)
    stat = DANU_SUCCESS
  end subroutine

  subroutine sequence_write_attr_real64_rank1(this, attr_name, attr_value, stat)
    class(sequence), intent(in) :: this
    character(*), intent(in) :: attr_name
    real(real64), intent(in) :: attr_value(:)
    integer, intent(out) :: stat
    call truchas_scorpio_write_attr_1d_double(this%foutput%scorp_id, &
        attr_name//c_null_char, attr_value, 1, shape(attr_value), this%path//c_null_char)
    stat = DANU_SUCCESS
  end subroutine

  subroutine sequence_write_dataset_attr_int32(this, dataset, attr_name, attr_value, stat)
    class(sequence), intent(in) :: this
    character(*), intent(in) :: dataset, attr_name
    integer(int32), intent(in) :: attr_value
    integer, intent(out) :: stat
    character(:), allocatable :: path
    path = this%path // '/' // dataset // c_null_char
    call truchas_scorpio_write_attr_0d_integer(this%foutput%scorp_id, &
        attr_name//c_null_char, attr_value, path)
    stat = DANU_SUCCESS
  end subroutine

  subroutine sequence_write_dataset_attr_real64(this, dataset, attr_name, attr_value, stat)
    class(sequence), intent(in) :: this
    character(*), intent(in) :: dataset, attr_name
    real(real64), intent(in) :: attr_value
    integer, intent(out) :: stat
    character(:), allocatable :: path
    path = this%path // '/' // dataset // c_null_char
    call truchas_scorpio_write_attr_0d_double(this%foutput%scorp_id, &
        attr_name//c_null_char, attr_value, path)
    stat = DANU_SUCCESS
  end subroutine

  subroutine sequence_write_dataset_attr_string(this, dataset, attr_name, attr_value, stat)
    class(sequence), intent(in) :: this
    character(*), intent(in) :: dataset, attr_name
    character(*), intent(in) :: attr_value
    integer, intent(out) :: stat
    character(:), allocatable :: path
    path = this%path // '/' // dataset // c_null_char
    call truchas_scorpio_write_attr_0d_string(this%foutput%scorp_id, &
        attr_name//c_null_char, attr_value//c_null_char, path)
    stat = DANU_SUCCESS
  end subroutine

! -------------------------------------------------------------

! output_probe

  subroutine probe_data_write (this, rdata, stat)
    class(output_probe), intent(in) :: this
    real(real64), intent(in), contiguous :: rdata(:,:)
    integer, intent(out) :: stat
    call truchas_scorpio_write_probe_data_2d_double(this%foutput%scorp_id, &
        this%dataset_id, shape(rdata), rdata)
    stat = DANU_SUCCESS
  end subroutine

  subroutine probe_write_attr_int32(this, attr_name, attr_value, stat)
    class(output_probe), intent(in) :: this
    character(*), intent(in) :: attr_name
    integer(int32), intent(in) :: attr_value
    integer, intent(out) :: stat
    call truchas_scorpio_write_attr_0d_integer(this%foutput%scorp_id, &
        attr_name//c_null_char, attr_value, this%path//c_null_char)
    stat = DANU_SUCCESS
  end subroutine

  subroutine probe_write_attr_real64(this, attr_name, attr_value, stat)
    class(output_probe), intent(in) :: this
    character(*), intent(in) :: attr_name
    real(real64), intent(in) :: attr_value
    integer, intent(out) :: stat
    call truchas_scorpio_write_attr_0d_double(this%foutput%scorp_id, &
        attr_name//c_null_char, attr_value, this%path//c_null_char)
    stat = DANU_SUCCESS
  end subroutine

  subroutine probe_write_attr_string(this, attr_name, attr_value, stat)
    class(output_probe), intent(in) :: this
    character(*), intent(in) :: attr_name
    character(*), intent(in) :: attr_value
    integer, intent(out) :: stat
    call truchas_scorpio_write_attr_0d_string(this%foutput%scorp_id, &
        attr_name//c_null_char, attr_value//c_null_char, this%path//c_null_char)
    stat = DANU_SUCCESS
  end subroutine

end module
