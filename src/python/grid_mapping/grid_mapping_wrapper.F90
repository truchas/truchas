!!
!! GRID_MAPPING_WRAPPER
!!
!! This provides an interface to the Kuprat and Portage mesh mappers for the
!! Python mapped restart utility.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! December 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module grid_mapping_wrapper

  use,intrinsic :: iso_c_binding
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: iso_fortran_env, only: error_unit
  use unstr_mesh_type
  use ext_exodus_mesh_type
  use data_mapper_class
  implicit none
  private

  public :: mapper_init, map_field

  type, public, bind(c) :: map_data
    integer(c_int) :: ncell, nnode
    type(c_ptr) :: coord, connect, blockid, mapper
  end type map_data

  type :: mapper_switch
    type(unstr_mesh), pointer :: src => null()
    type(unstr_mesh), pointer :: dest => null()
    integer, allocatable :: connect(:,:), blockid(:)
    class(data_mapper), allocatable :: mapper
  contains
    final :: delete
  end type mapper_switch

  logical :: mpi_initialized = .false.

contains

  !! Final subroutine for MAPPER_SWITCH objects
  subroutine delete(this)
    type(mapper_switch), intent(inout) :: this
    if (associated(this%src)) deallocate(this%src)
    if (associated(this%dest)) deallocate(this%dest)
  end subroutine delete


  !! Compute a mesh map from base mesh (ultimately an h5 file coming from python)
  !! to an exodus mesh. Returns both the exodus mesh data and the mapper.
  function mapper_init(nnode, ncell, connect, blockid, coord, exodus_filename, scale_factor, &
      use_portage_backend) &
    result(mesh_map) bind(c)

    use kuprat_mapper_type
#ifdef USE_PORTAGE
    use portage_mapper_type
#endif

    integer(c_int), intent(in), value :: nnode, ncell
    integer(c_int), intent(in) :: connect(8,ncell), blockid(ncell)
    real(c_double), intent(in) :: coord(3,nnode)
    character(kind=c_char), intent(in) :: exodus_filename(*)
    real(c_double), intent(in), value :: scale_factor
    logical(c_bool), intent(in), value :: use_portage_backend
    type(map_data) :: mesh_map

    type(mapper_switch), pointer :: this => null()
    character(:), allocatable :: filename

    if (.not.mpi_initialized) call mpi_init

    allocate(this)

#ifdef USE_PORTAGE
    if (use_portage_backend) then
      print "(a)", "Using Portage mapper..."
      allocate(portage_mapper :: this%mapper)
    else
      print "(a)", "Using Kuprat mapper..."
      allocate(kuprat_mapper :: this%mapper)
    end if
#else
    if (use_portage_backend) then
      write(error_unit,"(a)") "ERROR: Portage is not supported in the build of Truchas. Aborting."
      error stop
    end if
    print "(a)", "Using Kuprat mapper..."
    allocate(kuprat_mapper :: this%mapper)
#endif

    call mesh_from_raw(coord, connect, blockid, this%src)
    call mesh_from_file(c_to_f_str(exodus_filename), scale_factor, this%dest)
    call compute_connect_blockid(this)
    call this%mapper%init(this%src, this%dest)

    !! Store data in output variable
    mesh_map%ncell = this%dest%ncell
    mesh_map%nnode = this%dest%nnode
    mesh_map%coord = c_loc(this%dest%x)
    mesh_map%connect = c_loc(this%connect)
    mesh_map%blockid = c_loc(this%blockid)
    mesh_map%mapper = c_loc(this)

  end function mapper_init


  subroutine map_field(mapper, map_type, src_size, dest_size, defval, src, dest) bind(c)
    type(c_ptr), intent(in), value :: mapper
    integer(c_int), intent(in), value :: map_type, src_size, dest_size
    real(c_double), intent(in), value :: defval
    real(c_double), intent(in) :: src(src_size)
    real(c_double), intent(out) :: dest(dest_size)
    type(mapper_switch), pointer :: this => null()
    call c_f_pointer(mapper, this)
    call this%mapper%map_field(src, dest, defval, map_type)
  end subroutine map_field


  subroutine mpi_init()

    use parallel_util_module, only: parallel_init
    use parallel_communication
    use pgslib_module
    use truchas_env, only: prefix
    use truchas_logging_services

    character(PGSLib_CL_MAX_TOKEN_LENGTH), pointer :: argv(:) => null()

    call parallel_init(argv)
    call init_parallel_communication
    prefix = 'grid-mapper' ! TLS will write to 'grid-mapper.log'
    call TLS_initialize

    mpi_initialized = .true.

  end subroutine mpi_init


  subroutine compute_connect_blockid(this)

    type(mapper_switch), intent(inout) :: this

    integer, parameter :: TET_NODE_MAP(8) = [1,1,2,3,4,4,4,4]
    integer, parameter :: PYR_NODE_MAP(8) = [1,2,3,4,5,5,5,5]
    integer, parameter :: PRI_NODE_MAP(8) = [1,4,5,2,3,6,6,3]

    integer :: j

    allocate(this%connect(8,this%dest%ncell), this%blockid(this%dest%ncell))
    do j = 1, this%dest%ncell
      associate (cnode => this%dest%cnode(this%dest%xcnode(j):this%dest%xcnode(j+1)-1))
        select case (size(cnode))
        case (4)  ! tet
          this%connect(:,j) = cnode(TET_NODE_MAP)
        case (5)  ! pyramid
          this%connect(:,j) = cnode(PYR_NODE_MAP)
        case (6)  ! prism
          this%connect(:,j) = cnode(PRI_NODE_MAP)
        case (8)  ! hex
          this%connect(:,j) = cnode
        case default
          this%connect(:,j) = -1
          write(error_unit,"(a)") "ERROR: Invalid mesh element."
          error stop
          !INSIST(.false.)
        end select
      end associate

      associate (bitmask => this%dest%cell_set_mask(j))
        !INSIST(popcnt(bitmask) == 1)
        this%blockid(j) = this%dest%cell_set_id(trailz(bitmask))
      end associate
    end do

  end subroutine compute_connect_blockid


  subroutine mesh_from_raw(coord, connect, blockid, mesh)
    use parameter_list_type
    use unstr_mesh_factory
    use ext_exodus_mesh_type
    real(r8), intent(in) :: coord(:,:)
    integer, intent(in) :: connect(:,:), blockid(:)
    type(unstr_mesh), intent(out), pointer :: mesh
    integer :: stat
    character(:), allocatable :: errmsg
    type(parameter_list) :: params
    type(ext_exodus_mesh) :: exomesh
    call init_exodus_mesh(coord, connect, blockid, exomesh)
    mesh => new_unstr_mesh_aux(exomesh, params, stat, errmsg)
    if (stat /= 0) then
      write(error_unit,'(a)') errmsg
      error stop
    end if
  end subroutine mesh_from_raw


  subroutine mesh_from_file(exodus_filename, scale_factor, mesh)

    use parameter_list_type
    use unstr_mesh_factory

    character(*), intent(in) :: exodus_filename
    real(r8), intent(in) :: scale_factor
    type(unstr_mesh), intent(out), pointer :: mesh

    integer :: stat
    character(:), allocatable :: errmsg
    type(parameter_list) :: params

    call params%set("mesh-file", exodus_filename)
    call params%set("coord-scale-factor", scale_factor)
    call params%set("exodus-block-modulus", 10000)
    mesh => new_unstr_mesh(params, stat, errmsg)
    if (stat /= 0) then
      write(error_unit,'(a)') errmsg
      error stop
    end if

  end subroutine mesh_from_file


  subroutine init_exodus_mesh(coord, connect, blockid, mesh)

    real(r8), intent(in) :: coord(:,:)
    integer, intent(in) :: connect(:,:), blockid(:)
    type(ext_exodus_mesh), intent(out) :: mesh

    integer, parameter :: exodus_block_modulus = 10000

    mesh%num_dim = 3
    mesh%num_node = size(coord, dim=2)
    mesh%num_elem = size(connect, dim=2)
    mesh%coord = coord
    mesh%num_nset = 0
    mesh%num_sset = 0
    allocate(mesh%nset(mesh%num_nset), mesh%sset(mesh%num_sset))
    call identify_blocks
    call mesh%set_no_links

  contains

    subroutine identify_blocks()

      integer, parameter :: TET_NODE_MAP(4) = [1,3,4,5]
      integer, parameter :: PYR_NODE_MAP(5) = [1,2,3,4,5]
      integer, parameter :: PRI_NODE_MAP(6) = [1,4,5,2,3,6]

      integer :: j, n, k, eblk_id(100), eblk_n(100), eblkid

      k = 1
      eblk_n(1) = nnode(connect(:,1))
      eblk_id(1) = exodus_blockid(blockid(1), eblk_n(1))
      do j = 2, mesh%num_elem
        n = nnode(connect(:,j))
        eblkid = exodus_blockid(blockid(j), n)
        if (all(eblk_id(:k) /= eblkid)) then
          k = k + 1
          eblk_id(k) = eblkid
          eblk_n(k) = n
        end if
      end do

      mesh%num_eblk = k
      allocate(mesh%eblk(mesh%num_eblk))
      do k = 1, mesh%num_eblk
        associate (eblk => mesh%eblk(k))
          ! Merge exodus blocks, as done in new_unstr_mesh_file
          eblk%id = modulo(eblk_id(k), exodus_block_modulus)
          eblk%num_nodes_per_elem = eblk_n(k)

          n = 0
          do j = 1, mesh%num_elem
            if (blockid(j) == eblk%id .and. nnode(connect(:,j)) == eblk%num_nodes_per_elem) &
                n = n + 1
          end do
          eblk%num_elem = n
          allocate(eblk%connect(eblk%num_nodes_per_elem,eblk%num_elem))

          n = 0
          do j = 1, mesh%num_elem
            if (blockid(j) == eblk%id .and. nnode(connect(:,j)) == eblk%num_nodes_per_elem) then
              n = n + 1
              select case (eblk%num_nodes_per_elem)
              case (4)
                eblk%connect(:,n) = connect(TET_NODE_MAP,j)
              case (5)
                eblk%connect(:,n) = connect(PYR_NODE_MAP,j)
              case (6)
                eblk%connect(:,n) = connect(PRI_NODE_MAP,j)
              case (8)
                eblk%connect(:,n) = connect(:,j)
              case default
                eblk%connect(:,n) = -1
                write(error_unit,"(a)") "ERROR: Invalid mesh element."
                error stop
              end select
            end if
          end do
        end associate
      end do

      ! Ensure all cells were identified
      if (mesh%num_elem /= sum(mesh%eblk(:)%num_elem)) then
        write(error_unit,"(a)") "ERROR: Mismatch during block binning."
        error stop
      end if

    end subroutine identify_blocks

    pure integer function nnode(cnode)
      integer, intent(in) :: cnode(:)
      if (cnode(1) == cnode(2)) then
        nnode = 4 ! tet
      else if (cnode(7) == cnode(8)) then
        nnode = 5 ! pyramid
      else if (cnode(6) == cnode(7)) then
        nnode = 6 ! wedge
      else
        nnode = 8 ! hex
      end if
    end function nnode

    integer function exodus_blockid(blockid, nnode)
      integer, intent(in) :: blockid, nnode
      select case (nnode)
      case (4) ! TET4
        exodus_blockid = 1*exodus_block_modulus + blockid
      case (5) ! PYR5
        exodus_blockid = 2*exodus_block_modulus + blockid
      case (6) ! WED6
        exodus_blockid = 3*exodus_block_modulus + blockid
      case (8) ! HEX8
        exodus_blockid = 0*exodus_block_modulus + blockid
      case default
        write(error_unit,"(a)") "ERROR: Invalid mesh element."
        error stop
      end select
    end function exodus_blockid

  end subroutine init_exodus_mesh


  !! Convert a null-terminated C char* to a Fortran character scalar.
  pure function c_to_f_str(cstr) result(fstr)
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

end module grid_mapping_wrapper
