#include "f90_assert.fpp"

module truchasio
  use iso_c_binding, only: c_ptr, c_char, c_double, c_int, c_int8_t, &
    c_null_ptr, c_associated, c_null_char, c_f_pointer
  use pgslib_module, only: &
    broadcast  => PGSLib_BCast, &
    collate    => PGSLib_Collate
  use danu_module, only: DANU_SUCCESS
  use danu_module, only: &
    danu_attribute_write => attribute_write, &
    danu_data_write => data_write, &
    danu_mesh_add_unstructured => mesh_add_unstructured, &
    danu_mesh_write_connectivity => mesh_write_connectivity, &
    danu_mesh_write_coordinates => mesh_write_coordinates, &
    danu_output_file_close => output_file_close, &
    danu_output_file_create => output_file_create, &
    danu_probe_create_data => probe_create_data, &
    danu_probe_data_open => probe_data_open, &
    danu_probe_data_write => probe_data_write, &
    danu_sequence_next_id => sequence_next_id, &
    danu_simulation_add => simulation_add, &
    danu_simulation_data_write => simulation_data_write, &
    danu_simulation_link_mesh => simulation_link_mesh, &
    danu_simulation_open_data => simulation_open_data
  implicit none
  private
  public DANU_SUCCESS, output_file, output_mesh, simulation, sequence, &
    output_probe, output_dataset


  type :: output_file
    type(c_ptr), private :: fid = c_null_ptr ! h5 file id
    type(c_ptr), private :: scorp_id = c_null_ptr ! Scorpio id
    logical :: is_IOP = .false.  ! True if this PE is the IO PE.
  contains
    procedure :: open  => output_file_open
    procedure :: close => output_file_close
    procedure :: mesh_add_unstructured
    procedure :: simulation_add
    procedure :: simulation_link_mesh
    final :: output_file_final
  end type


  type :: simulation
    type(c_ptr), private :: sid = c_null_ptr ! Truchas IO simulation id
    type(output_file), pointer, private :: foutput => null()
  contains
    generic :: data_write => &
      data_write_integer_rank0, &
      data_write_integer_rank1, &
      data_write_real8_rank0, &
      data_write_real8_rank1, &
      data_write_real8_rank2
    procedure :: sequence_next_id
    procedure :: probe_create_data
    procedure :: probe_data_open
    generic :: attribute_write => &
      simulation_attribute_write_integer, &
      simulation_attribute_write_real8, &
      simulation_attribute_write_character0

    ! private
    procedure, private :: data_write_integer_rank0
    procedure, private :: data_write_integer_rank1
    procedure, private :: data_write_real8_rank0
    procedure, private :: data_write_real8_rank1
    procedure, private :: data_write_real8_rank2
    procedure, private :: simulation_attribute_write_integer
    procedure, private :: simulation_attribute_write_real8
    procedure, private :: simulation_attribute_write_character0
  end type


  type :: output_mesh
    type(c_ptr), private :: mid = c_null_ptr ! mesh id
    type(output_file), pointer, private :: foutput => null()
  contains
    procedure :: write_coordinates
    procedure :: write_connectivity
  end type


  type :: sequence
    type(c_ptr), private :: nid = c_null_ptr ! sequence id
    type(output_file), pointer, private :: foutput => null()
    character(:), allocatable :: fname ! sequence path in HDF5
  contains
    procedure :: open_data => simulation_open_data
    generic :: data_write => &
      simulation_data_write_byte_rank2, &
      simulation_data_write_integer_rank1, &
      simulation_data_write_real8_rank1, &
      simulation_data_write_real8_rank2
    generic :: attribute_write => &
      sequence_attribute_write_integer, &
      sequence_attribute_write_real8, &
      sequence_attribute_write_real81, &
      sequence_attribute_write_character0

    ! private
    procedure, private :: simulation_data_write_byte_rank2
    procedure, private :: simulation_data_write_integer_rank1
    procedure, private :: simulation_data_write_real8_rank1
    procedure, private :: simulation_data_write_real8_rank2
    procedure, private :: sequence_attribute_write_integer
    procedure, private :: sequence_attribute_write_real8
    procedure, private :: sequence_attribute_write_real81
    procedure, private :: sequence_attribute_write_character0
  end type


  type :: output_probe
    type(c_ptr), private :: pid = c_null_ptr ! probe id
    type(output_file), pointer, private :: foutput => null()
  contains
    procedure :: data_write => probe_data_write
    generic :: attribute_write => &
      probe_attribute_write_integer, &
      probe_attribute_write_real8, &
      probe_attribute_write_character0

    ! private
    procedure, private :: probe_attribute_write_integer
    procedure, private :: probe_attribute_write_real8
    procedure, private :: probe_attribute_write_character0
  end type

  type :: output_dataset
    type(c_ptr), private :: did = c_null_ptr ! dataset id
    type(output_file), pointer, private :: foutput => null()
  contains
    generic :: attribute_write => &
      dataset_attribute_write_integer, &
      dataset_attribute_write_real8, &
      dataset_attribute_write_character0

    ! private
    procedure, private :: dataset_attribute_write_integer
    procedure, private :: dataset_attribute_write_real8
    procedure, private :: dataset_attribute_write_character0
  end type

  interface collate
    module procedure collate_L2, collate_I2, collate_S2, collate_D2, &
        collate_int8_2
  end interface

  interface
    type(c_ptr) function truchas_scorpio_create_handle(filename, &
        numIOgroups) bind(c)
    import :: c_ptr, c_char, c_int
    character(kind=c_char), intent(in)  :: filename(*)
    integer(c_int), intent(in), value :: numIOgroups
    end function

    subroutine truchas_scorpio_free_handle(h) bind(c)
    import :: c_ptr
    type(c_ptr), value :: h
    end subroutine

    subroutine truchas_scorpio_write_dataset_1d_integer(h, &
        name, &
        vector, global_dim, local_dim) bind(c)
    import :: c_ptr, c_char, c_int
    type(c_ptr), value :: h
    character(kind=c_char), intent(in)  :: name(*)
    integer(c_int), intent(in), value :: global_dim, local_dim
    integer(c_int), intent(in) :: vector(local_dim)
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

    integer(c_int) function truchas_scorpio_get_group_name(gid, cname) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: gid
    type(c_ptr) :: cname
    end function

    subroutine truchas_scorpio_str_free(cname) bind(c)
    import :: c_ptr
    type(c_ptr), value :: cname
    end subroutine

    type(c_ptr) function truchas_scorpio_hdf5_handle_danu_create(h) bind(c)
    import :: c_ptr
    type(c_ptr), value :: h
    end function

    subroutine truchas_scorpio_handle_initial_setup(h) bind(c)
    import :: c_ptr
    type(c_ptr), value :: h
    end subroutine

    subroutine truchas_scorpio_handle_file_close(h) bind(c)
    import :: c_ptr
    type(c_ptr), value :: h
    end subroutine

    subroutine truchas_scorpio_hdf5_handle_danu_free(fid_ptr) bind(c)
    import :: c_ptr
    type(c_ptr), value :: fid_ptr
    end subroutine

  end interface

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Specific procedures for COLLATE: extend PGSLib_Collate to rank-2 arrays
 !!

#define _LOGICAL_DATA_
#include "collate.fpp"
 
#define _INT8_DATA_
#include "collate.fpp"
 
#define _INTEGER_DATA_
#include "collate.fpp"
 
#define _SINGLE_DATA_
#include "collate.fpp"
 
#define _DOUBLE_DATA_
#include "collate.fpp"

pure integer function str_int_len(i) result(sz)
! Returns the length of the string representation of 'i'
integer, intent(in) :: i
integer, parameter :: MAX_STR = 100
character(MAX_STR) :: s
! If 's' is too short (MAX_STR too small), Fortan will abort with:
! "Fortran runtime error: End of record"
write(s, '(i0)') i
sz = len_trim(s)
end function

pure function str_int(i) result(s)
! Converts integer "i" to string
integer, intent(in) :: i
character(str_int_len(i)) :: s
write(s, '(i0)') i
end function

function stringc2f(n, cstr) result(fstr)
integer, intent(in) :: n
type(c_ptr), intent(in) :: cstr
character(:), allocatable :: fstr
character(n, kind=c_char), pointer :: fptr
call c_f_pointer(cstr, fptr)
fstr = fptr
end function

subroutine broadcast_str(s, is_IOP)
character(:), allocatable, intent(inout) :: s
logical, intent(in) :: is_IOP
integer :: n
if (is_IOP) then
  INSIST(allocated(s))
  n = len(s)
end if
call broadcast(n)
if (.not. is_IOP) then
  INSIST(.not. allocated(s))
  allocate(character(n)::s)
end if
call broadcast(s)
end subroutine

! -------------------------------------------------------------

! output_file

  subroutine output_file_open(this, filename, is_IOP, stat)
    class(output_file), intent(out) :: this
    character(*), intent(in) :: filename
    logical, intent(in) :: is_IOP
    integer, intent(out), optional :: stat
    this%is_IOP = is_IOP
    this%scorp_id = truchas_scorpio_create_handle(filename // c_null_char, 1)

    if (this%is_IOP) then
      this%fid = truchas_scorpio_hdf5_handle_danu_create(this%scorp_id)
      call truchas_scorpio_handle_initial_setup(this%scorp_id)
      !call danu_output_file_create(filename, this%fid, stat)
      stat = DANU_SUCCESS
      INSIST(c_associated(this%fid))
    end if
    call broadcast (stat)
  end subroutine

  subroutine output_file_final (this)
    type(output_file), intent(inout) :: this
    !if (c_associated(this%fid)) call danu_output_file_close(this%fid)
    if (c_associated(this%scorp_id)) then
      if (this%is_IOP) then
        call truchas_scorpio_handle_file_close(this%scorp_id)
        call truchas_scorpio_hdf5_handle_danu_free(this%fid)
      end if
      call truchas_scorpio_free_handle(this%scorp_id)
    end if
  end subroutine

  !! The passed object is intent(out) and this causes the object to be finalized
  !! and its components default initialized on entry to the subroutine.
  subroutine output_file_close (this)
    class(output_file), intent(out) :: this
    INSIST(.not. c_associated(this%fid))
  end subroutine

  subroutine simulation_add(this, sim_name, sim, stat)
  class(output_file), target, intent(in) :: this
  character(*,kind=C_CHAR),      intent(in)  :: sim_name
  class(simulation),                 intent(out) :: sim
  integer, intent(out), optional :: stat
  sim%foutput => this
  if (this%is_IOP) then
    INSIST(c_associated(this%fid))
    INSIST(.not. c_associated(sim%sid))
    call danu_simulation_add(this%fid, sim_name, sim%sid, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine simulation_link_mesh(this, sim, mesh_name, stat)
  class(output_file),                intent(in) :: this
  class(simulation),                 intent(in) :: sim
  character(*,kind=C_CHAR),      intent(in) :: mesh_name
  integer, intent(out), optional :: stat
  if (this%is_IOP) then
    INSIST(c_associated(this%fid))
    INSIST(c_associated(sim%sid))
    call danu_simulation_link_mesh(this%fid,sim%sid,mesh_name,stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine mesh_add_unstructured(this,mesh_name,elem_order,mesh_dim, &
      out_mesh,stat)
  class(output_file), target, intent(in) :: this
  character(*,kind=C_CHAR), intent(in)  :: mesh_name
  integer(kind=C_INT),          intent(in)  :: elem_order
  integer(kind=C_INT),          intent(in)  :: mesh_dim
  class(output_mesh),           intent(out) :: out_mesh
  integer, intent(out), optional :: stat
  out_mesh%foutput => this
  if (this%is_IOP) then
    INSIST(c_associated(this%fid))
    call danu_mesh_add_unstructured(this%fid,mesh_name,elem_order,mesh_dim, &
      out_mesh%mid, stat)
  end if
  call broadcast (stat)
  end subroutine

! -------------------------------------------------------------

! simulation

  subroutine data_write_integer_rank0(this,data_name,idata,stat)
  class(simulation), intent(in) :: this
  character(*),    intent(in)  :: data_name
  integer(kind=C_INT), intent(in)  :: idata
  integer, intent(out), optional :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%sid))
    call danu_data_write(this%sid, data_name, idata, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine data_write_integer_rank1(this,data_name,idata,stat)
  class(simulation), intent(in) :: this
  character(*),    intent(in)  :: data_name
  integer(kind=C_INT), dimension(:), intent(in) :: idata
  integer, intent(out), optional :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%sid))
    call danu_data_write(this%sid, data_name, idata, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine data_write_real8_rank0(this,data_name,r8data,stat)
  class(simulation), intent(in) :: this
  character(*),    intent(in)  :: data_name
  real(kind=C_DOUBLE), intent(in)  :: r8data
  integer, intent(out), optional :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%sid))
    call danu_data_write(this%sid, data_name, r8data, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine data_write_real8_rank1(this,data_name,r8data,stat)
  class(simulation), intent(in) :: this
  character(*),    intent(in)  :: data_name
  real(kind=C_DOUBLE), dimension(:), intent(in) :: r8data
  integer, intent(out), optional :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%sid))
    call danu_data_write(this%sid, data_name, r8data, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine data_write_real8_rank2(this,data_name,r8data,stat)
  class(simulation), intent(in) :: this
  character(*),    intent(in)  :: data_name
  real(kind=C_DOUBLE), dimension(:,:), intent(in) :: r8data
  integer, intent(out), optional :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%sid))
    call danu_data_write(this%sid, data_name, r8data, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine sequence_next_id(this,cyc,time,seq,stat)
  class(simulation),                 intent(in)  :: this
  integer(kind=C_INT),               intent(in)  :: cyc
  real(kind=C_DOUBLE),               intent(in)  :: time
  class(sequence),                   intent(out) :: seq
  integer, intent(out), optional :: stat
  integer :: cname_size
  type(c_ptr) :: cname
  seq%foutput => this%foutput
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%sid))
    call danu_sequence_next_id(this%sid, cyc, time, seq%nid, stat)
    cname_size = truchas_scorpio_get_group_name(seq%nid, cname)
    seq%fname = stringc2f(cname_size, cname)
    call truchas_scorpio_str_free(cname)
  end if
  call broadcast_str(seq%fname, this%foutput%is_IOP)
  call broadcast (stat)
  end subroutine

  subroutine probe_data_open(this,probe_name,probe,stat)
  class(simulation),                 intent(in)  :: this
  character(*,kind=C_CHAR),      intent(in)  :: probe_name
  class(output_probe),               intent(out) :: probe
  integer, intent(out), optional :: stat
  probe%foutput => this%foutput
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%sid))
    call danu_probe_data_open(this%sid, probe_name, probe%pid, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine probe_create_data(this,probe_name,rdata,probe,stat)
  class(simulation),                   intent(in)  :: this
  character(*,kind=C_CHAR),        intent(in)  :: probe_name
  real(kind=C_DOUBLE), dimension(:,:), intent(in)  :: rdata
  class(output_probe),                 intent(out) :: probe
  integer, intent(out), optional :: stat
  probe%foutput => this%foutput
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%sid))
    call danu_probe_create_data(this%sid, probe_name, rdata, probe%pid, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine simulation_attribute_write_integer(this,attr_name,value,stat)
  class(simulation), intent(in)  :: this
  character(*),  intent(in)  :: attr_name
  integer,           intent(in)  :: value
  integer, intent(out), optional :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%sid))
    call danu_attribute_write(this%sid, attr_name, value, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine simulation_attribute_write_real8(this,attr_name,value,stat)
  class(simulation),   intent(in)  :: this
  character(*),    intent(in)  :: attr_name
  real(kind=C_DOUBLE), intent(in)  :: value
  integer, intent(out), optional   :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%sid))
    call danu_attribute_write(this%sid, attr_name, value, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine simulation_attribute_write_character0(this,attr_name,char_data,stat)
  class(simulation),            intent(in)  :: this
  character(*),             intent(in)  :: attr_name
  character(*,kind=C_CHAR), intent(in)  :: char_data
  integer, intent(out), optional   :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%sid))
    call danu_attribute_write(this%sid, attr_name, char_data, stat)
  end if
  call broadcast (stat)
  end subroutine


! -------------------------------------------------------------

! output_mesh

  subroutine write_coordinates(this,nnodes,x,y,z,stat)
  class(output_mesh), intent(in) :: this
  integer(kind=C_INT),               intent(in)  :: nnodes
  real(kind=C_DOUBLE), dimension(:), intent(in)  :: x
  real(kind=C_DOUBLE), dimension(:), intent(in)  :: y
  real(kind=C_DOUBLE), dimension(:), intent(in)  :: z
  integer, intent(out), optional :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%mid))
    call danu_mesh_write_coordinates(this%mid, nnodes, x, y, z, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine write_connectivity(this,num,idata,stat)
  class(output_mesh), intent(in) :: this
  integer(kind=C_INT),                 intent(in) :: num
  integer(kind=C_INT), dimension(:,:), intent(in) :: idata
  integer, intent(out), optional :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%mid))
    call danu_mesh_write_connectivity(this%mid, num, idata, stat)
  end if
  call broadcast (stat)
  end subroutine


! -------------------------------------------------------------

! sequence

  subroutine simulation_open_data(this,dataname,dataset,stat)
  class(sequence),                   intent(in)  :: this
  character(*,kind=C_CHAR),      intent(in)  :: dataname
  class(output_dataset),             intent(out) :: dataset
  integer, intent(out), optional :: stat
  dataset%foutput => this%foutput
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%nid))
    call danu_simulation_open_data(this%nid, dataname, dataset%did, stat)
  end if
  call broadcast (stat)
  end subroutine

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
    this%fname // "/" // trim(data_name) // c_null_char, ldata, glen, &
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
    this%fname // "/" // trim(data_name) // c_null_char, &
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
    this%fname // "/" // trim(data_name) // c_null_char, &
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
    this%fname // "/" // trim(data_name) // c_null_char, ldata, glen, &
    size(ldata, 1), size(ldata, 2))
  stat = DANU_SUCCESS
  end subroutine

  subroutine sequence_attribute_write_integer(this,attr_name,value,stat)
  class(sequence), intent(in)  :: this
  character(*),  intent(in)  :: attr_name
  integer,           intent(in)  :: value
  integer, intent(out), optional :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%nid))
    call danu_attribute_write(this%nid, attr_name, value, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine sequence_attribute_write_real8(this,attr_name,value,stat)
  class(sequence), intent(in)  :: this
  character(*),    intent(in)  :: attr_name
  real(kind=C_DOUBLE), intent(in)  :: value
  integer, intent(out), optional   :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%nid))
    call danu_attribute_write(this%nid, attr_name, value, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine sequence_attribute_write_real81(this,attr_name,value,stat)
  class(sequence), intent(in)  :: this
  character(*),    intent(in)  :: attr_name
  real(kind=C_DOUBLE), intent(in)  :: value(:)
  integer, intent(out), optional   :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%nid))
    call danu_attribute_write(this%nid, attr_name, value, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine sequence_attribute_write_character0(this,attr_name,char_data,stat)
  class(sequence),          intent(in)  :: this
  character(*),             intent(in)  :: attr_name
  character(*,kind=C_CHAR), intent(in)  :: char_data
  integer, intent(out), optional   :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%nid))
    call danu_attribute_write(this%nid, attr_name, char_data, stat)
  end if
  call broadcast (stat)
  end subroutine

! -------------------------------------------------------------

! output_probe

  subroutine probe_data_write(this,rdata,stat)
  class(output_probe),                  intent(in)  :: this
  real(kind=C_DOUBLE), dimension(:,:),  intent(in)  :: rdata
  integer, intent(out), optional :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%pid))
    call danu_probe_data_write(this%pid, rdata, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine probe_attribute_write_integer(this,attr_name,value,stat)
  class(output_probe), intent(in)  :: this
  character(*),  intent(in)  :: attr_name
  integer,           intent(in)  :: value
  integer, intent(out), optional :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%pid))
    call danu_attribute_write(this%pid, attr_name, value, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine probe_attribute_write_real8(this,attr_name,value,stat)
  class(output_probe), intent(in)  :: this
  character(*),    intent(in)  :: attr_name
  real(kind=C_DOUBLE), intent(in)  :: value
  integer, intent(out), optional   :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%pid))
    call danu_attribute_write(this%pid, attr_name, value, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine probe_attribute_write_character0(this,attr_name,char_data,stat)
  class(output_probe),          intent(in)  :: this
  character(*),             intent(in)  :: attr_name
  character(*,kind=C_CHAR), intent(in)  :: char_data
  integer, intent(out), optional   :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%pid))
    call danu_attribute_write(this%pid, attr_name, char_data, stat)
  end if
  call broadcast (stat)
  end subroutine

! -------------------------------------------------------------

! dataset

  subroutine dataset_attribute_write_integer(this,attr_name,value,stat)
  class(output_dataset),    intent(in)  :: this
  character(*),  intent(in)  :: attr_name
  integer,           intent(in)  :: value
  integer, intent(out), optional :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%did))
    call danu_attribute_write(this%did, attr_name, value, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine dataset_attribute_write_real8(this,attr_name,value,stat)
  class(output_dataset),      intent(in)  :: this
  character(*),    intent(in)  :: attr_name
  real(kind=C_DOUBLE), intent(in)  :: value
  integer, intent(out), optional   :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%did))
    call danu_attribute_write(this%did, attr_name, value, stat)
  end if
  call broadcast (stat)
  end subroutine

  subroutine dataset_attribute_write_character0(this,attr_name,char_data,stat)
  class(output_dataset),               intent(in)  :: this
  character(*),             intent(in)  :: attr_name
  character(*,kind=C_CHAR), intent(in)  :: char_data
  integer, intent(out), optional   :: stat
  if (this%foutput%is_IOP) then
    INSIST(c_associated(this%did))
    call danu_attribute_write(this%did, attr_name, char_data, stat)
  end if
  call broadcast (stat)
  end subroutine

end module
