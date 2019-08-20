!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module grid_mapping_wrapper

  use,intrinsic :: iso_c_binding
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: iso_fortran_env, only: error_unit
  use exodus_mesh_type
  implicit none
  private

  public :: mesh_map, map_cell_field_c

  type :: base_mesh
    integer :: ncell, nnode
    real(r8), allocatable :: coord(:,:)   ! node coordinates
    integer,  allocatable :: connect(:,:) ! cell connectivity
    integer,  allocatable :: blockid(:)   ! cell block IDs
    logical :: is_tet_mesh = .false.
  contains
    procedure, private :: base_mesh_init_exo
    procedure, private :: base_mesh_init_raw
    generic :: init => base_mesh_init_exo, base_mesh_init_raw
  end type base_mesh

  type, bind(c) :: map_data
    integer(c_int) :: ncell, nnode
    type(c_ptr) :: coord, connect, blockid, mesh_map
  end type map_data

contains

  !! Compute a mesh map from base mesh (ultimately an h5 file coming from python)
  !! to an exodus mesh. Returns both the exodus mesh data and the mesh map.
  function mesh_map(nnode, ncell, connect, blockid, coord, exodus_filename, scale_factor) bind(c)

    use exodus_mesh_io, only: read_exodus_mesh
    use grid_mapping_module, only: gm_mesh, compute_int_volumes, grid_int_vols

    integer(c_int), intent(in), value :: nnode, ncell
    integer(c_int), intent(in) :: connect(8,ncell), blockid(ncell)
    real(c_double), intent(in) :: coord(3,nnode)
    character(kind=c_char), intent(in) :: exodus_filename(*)
    real(c_double), intent(in), value :: scale_factor
    type(map_data) :: mesh_map

    integer :: stat
    character(:), allocatable :: fpath, errmsg
    type(exodus_mesh) :: exomesh
    type(base_mesh) :: base_src
    type(gm_mesh) :: gm_src, gm_dest
    type(base_mesh), pointer :: base_dest => null()
    type(grid_int_vols), pointer :: map => null()

    !! Allocate pointers to be returned through the C bindings.
    allocate(base_dest, map)

    !! Read and scale the exodus mesh.
    fpath = trim(c_to_f_str(exodus_filename))
    call read_exodus_mesh(fpath, exomesh, stat, errmsg)
    if (stat /= 0) then
      write(error_unit,'(3a)') 'Error reading ', trim(fpath), ': ' // errmsg
      stop 1
    end if

    if (scale_factor /= 1) &
        exomesh%coord = scale_factor * exomesh%coord

    !! Calculate the mesh mapping data.
    call base_src%init(nnode, ncell, connect, blockid, coord)
    call base_dest%init(exomesh)
    call gm_mesh_init(base_src, gm_src)
    call gm_mesh_init(base_dest, gm_dest)
    call compute_int_volumes(gm_src, gm_dest, map)

    !! Store data in output variable
    mesh_map%ncell = base_dest%ncell
    mesh_map%nnode = base_dest%nnode
    mesh_map%coord = c_loc(base_dest%coord)
    mesh_map%connect = c_loc(base_dest%connect)
    mesh_map%blockid = c_loc(base_dest%blockid)
    mesh_map%mesh_map = c_loc(map)

  end function mesh_map


  subroutine map_cell_field_c(mesh_map_c, form, in_size, out_size, defval, src, dest) bind(c)

    use grid_mapping_module, only: map_cell_field, grid_int_vols

    type(c_ptr), intent(in), value :: mesh_map_c
    integer(c_int), intent(in), value :: form, in_size, out_size
    real(c_double), intent(in), value :: defval
    real(c_double), intent(in) :: src(in_size)
    real(c_double), intent(out) :: dest(out_size)

    integer, parameter :: MAP_FORM_DEFAULT  = 0
    integer, parameter :: MAP_FORM_CONSTANT = 1
    integer, parameter :: MAP_FORM_INTEGRAL = 2

    type(grid_int_vols), pointer :: mesh_map => null()

    call c_f_pointer(mesh_map_c, mesh_map)

    select case (form)
    case (MAP_FORM_CONSTANT)  ! THIS PRESERVES A CONSTANT-VALUED FIELD
      call map_cell_field(src, dest, mesh_map, strict=.true., defval=defval, &
          preserve_constants=.true.)
    case (MAP_FORM_INTEGRAL)  ! THIS PRESERVES THE INTEGRAL OF THE FIELD
      call map_cell_field(src, dest, mesh_map, strict=.true., defval=defval, &
          exactly_conservative=.true.)
    case (MAP_FORM_DEFAULT)   ! I DON'T ACTUALLY KNOW WHAT THIS DOES
      call map_cell_field(src, dest, mesh_map, strict=.true., defval=defval)
    case default
      print *, 'Error: Map cell form not recognized.'
      stop 1
    end select

  end subroutine map_cell_field_c


  !! =========== PRIVATE ROUTINES ===========

  !! Convert a null-terminated C char* to a Fortran character scalar.
  function c_to_f_str(cstr) result(fstr)

    character(kind=c_char), intent(in) :: cstr(*)
    character(:), allocatable :: fstr

    integer :: i

    i = 0
    do while (cstr(i+1) /= c_null_char)
      i = i + 1
    end do

    allocate(character(len=i) :: fstr)
    fstr = transfer(cstr(:i), fstr)

  end function c_to_f_str


  !! Wire a 'base mesh' to a 'gm mesh'.
  subroutine gm_mesh_init(src, dest)

    use grid_mapping_module, only: gm_mesh

    type(base_mesh), intent(in) :: src
    type(gm_mesh), intent(out) :: dest

    dest%nnod = src%nnode
    dest%nelt = src%ncell
    dest%pos_node  = src%coord
    if (src%is_tet_mesh) then
      dest%node_elt  = src%connect(2:5,:)
    else
      dest%node_elt  = src%connect
    end if
    dest%block_elt = src%blockid

  end subroutine gm_mesh_init

  !! Initialize a BASE_MESH from an ExodusII mesh
  subroutine base_mesh_init_exo(this, mesh)

    class(base_mesh),   intent(out) :: this
    type(exodus_mesh), intent(in)  :: mesh

    !! Truchas degenerate hex node numbering to ExodusII element node numbering.
    integer, parameter ::   TET_NODE_MAP(8) = [1,1,2,3,4,4,4,4]
    integer, parameter :: WEDGE_NODE_MAP(8) = [1,4,5,2,3,6,6,3]

    !! ExodusII element types identified by their number of nodes.
    integer, parameter :: TET=4, WEDGE=6, HEX=8

    integer :: n, i, j, nodes_per_elem
    logical :: non_hex_mesh

    this%nnode = mesh%num_node
    this%ncell = mesh%num_elem

    if (mesh%num_dim /= 3) then
      write(error_unit,'(a,i0)') 'Error: target mesh is not 3D: num_dim=', mesh%num_dim
      stop 1
    end if
    allocate(this%coord(3,this%nnode))
    this%coord = mesh%coord

    allocate(this%connect(8,this%ncell), this%blockid(this%ncell))

    !! Translate Exodus element connectivity to Truchas convention.
    n = 0
    non_hex_mesh = .false.
    this%is_tet_mesh = .true.
    do i = 1, mesh%num_eblk
      nodes_per_elem = size(mesh%eblk(i)%connect,dim=1)
      !! Translate Exodus element connectivity to Truchas convention.  See the ExodusII
      !! section in MESH_READ from MESH_INPUT_MODULE and the later section that converts
      !! non-hex elements into degenerate hexes as the reference for what must be done here.
      select case (nodes_per_elem)
      case (HEX)
        this%is_tet_mesh = .false.
        do j = 1, mesh%eblk(i)%num_elem
          n = n + 1
          this%blockid(n) = mesh%eblk(i)%ID
          this%connect(:,n) = mesh%eblk(i)%connect(:,j)
        end do
      case (TET)
        non_hex_mesh = .true.
        do j = 1, mesh%eblk(i)%num_elem
          n = n + 1
          this%blockid(n) = mesh%eblk(i)%ID
          this%connect(:,n) = mesh%eblk(i)%connect(TET_NODE_MAP,j)
        end do
      case (WEDGE)
        this%is_tet_mesh = .false.
        non_hex_mesh = .true.
        do j = 1, mesh%eblk(i)%num_elem
          n = n + 1
          this%blockid(n) = mesh%eblk(i)%ID
          this%connect(:,n) = mesh%eblk(i)%connect(WEDGE_NODE_MAP,j)
        end do
      case default
        write(error_unit,'(a,i0,a)') 'Error: unknown element type in target mesh: ', &
            nodes_per_elem, '-node element'
        stop 1
      end select
    end do
    if (n /= this%ncell) then
      print *, 'Internal Error: n /= this%ncell.'
      stop 1
    end if

    if (non_hex_mesh .and. .not.this%is_tet_mesh) &
        write(error_unit,'(a)') 'Warning: target mesh contains degenerate hex elements'

  end subroutine base_mesh_init_exo


  !! Initialize a BASE_MESH from raw arrays (ultimately from h5 file)
  subroutine base_mesh_init_raw(this, nnode, ncell, connect, blockid, coord)
    class(base_mesh), intent(out) :: this
    integer, intent(in) :: nnode, ncell, connect(:,:), blockid(:)
    real(r8), intent(in) :: coord(:,:)
    this%nnode = nnode
    this%ncell = ncell
    this%coord = coord
    this%connect = connect
    this%blockid = blockid
    this%is_tet_mesh = all(this%connect(1,:) == this%connect(2,:))
  end subroutine base_mesh_init_raw

end module grid_mapping_wrapper
