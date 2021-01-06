!!
!! MESH_MANAGER
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! July 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"
!#define DEBUG_MESH

module mesh_manager

  use parameter_list_type
  use truchas_logging_services
  use simpl_mesh_type
  implicit none
  private

  public :: enable_mesh, named_mesh_ptr, unstr_mesh_ptr, simpl_mesh_ptr
  public :: read_truchas_mesh_namelists, init_mesh_manager
  public :: get_main_mesh_size  ! HACK

  public :: simpl_mesh ! re-export; do not want this -- FIXME

  interface init_mesh_manager
    module procedure init_mesh_manager, init_mesh_manager_params
  end interface

  type(parameter_list), target, save :: meshes

  !! The value of the 'mesh' parameter will be of this private type.
  !! A PARAMETER_LIST stores (shallow) copies of parameter values, and thus
  !! cannot properly store a mesh object.  Instead we have it store a pointer
  !! to a mesh object by boxing the pointer in a TYPE(ANY_MESH) variable and
  !! and storing that variable.
  type :: any_mesh
    class(*), pointer :: mesh => null()
  end type any_mesh

contains

  !! NNC, Jun 2014.  A new initialization routine that takes a parameter list
  !! to specify the mesh files and parameters.  This is an alternative to
  !! using READ_MESH_NAMELISTS to read the data from an input file, and was
  !! added to simplify the development of unit tests requiring a mesh.  The
  !! expected syntax of the parameter list is a list of sublists.  The name
  !! of a sublist is taken as the name of the mesh, and the sublist has the
  !! following parameters:
  !!
  !!                 'mesh-file' : string (required)
  !!        'coord-scale-factor' : real-scalar (optional; default 1.0)
  !!    'interface-side-set-ids' : integer-array (optional)

  subroutine init_mesh_manager_params (params)

    use kinds, only: r8
    use string_utilities, only: raise_case

    type(parameter_list) :: params

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist1, plist2
    integer  :: exodus_block_modulus
    real(r8) :: coord_scale_factor
    character(:), allocatable :: mesh_file
#ifdef INTEL_COMPILER_WORKAROUND
    character(:), allocatable :: name
#endif
    integer, allocatable :: side_sets(:)

    !! Copy the input parameter list to the modules private parameter list.
    !! We want case-insensitive mesh names in the module, so we convert them
    !! to upper case when copying; later searching will use the upper-cased
    !! version of the specified name.
    piter = parameter_list_iterator(params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist1 => piter%sublist()
#ifdef INTEL_COMPILER_WORKAROUND
      !Intel tracking ID: DPD200357444
      name = piter%name()
      name = raise_case(name)
      plist2 => meshes%sublist(name)
#else
      plist2 => meshes%sublist(raise_case(piter%name()))
#endif
      call plist1%get ('mesh-file', mesh_file)
      call plist2%set ('mesh-file', mesh_file)
      call plist2%set ('enabled', .true.)

      if (plist1%is_parameter('coord-scale-factor')) then
        call plist1%get ('coord-scale-factor', coord_scale_factor)
        call plist2%set ('coord-scale-factor', coord_scale_factor)
      end if

      if (plist1%is_parameter('interface-side-set-ids')) then
        call plist1%get ('interface-side-set-ids', side_sets)
        call plist2%set ('interface-side-set-ids', side_sets)
      end if

      if (plist1%is_parameter('exodus-block-modulus')) then
        call plist1%get ('exodus-block-modulus', exodus_block_modulus)
        call plist2%set ('exodus-block-modulus', exodus_block_modulus)
      end if

      call piter%next
    end do

    !! Now call the original initialization procedure.
    call init_mesh_manager

  end subroutine init_mesh_manager_params

  subroutine init_mesh_manager

    use unstr_mesh_factory
    use simpl_mesh_factory
#ifdef DEBUG_MESH
    use unstr_mesh_gmv
#endif
    use truchas_timers

    integer :: stat
    logical :: enabled, em_mesh
    character(:), allocatable :: errmsg
    type(parameter_list), pointer :: plist
    type(parameter_list_iterator) :: piter
    type(unstr_mesh), pointer :: umesh
    type(simpl_mesh), pointer :: smesh
#if defined(GNU_PR49213)
    type(any_mesh) :: box
#endif

    piter = parameter_list_iterator(meshes)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      call plist%get ('enabled', enabled, default=.false.)
      if (enabled) then
        call start_timer('mesh-'//piter%name())
        call TLS_info ('')
        call TLS_info ('Initializing mesh "' // piter%name() // '" ...')
        call plist%get ('em-mesh', em_mesh, default=.false.)
        if (em_mesh) then
          smesh => new_simpl_mesh(plist, stat, errmsg)
          if (stat /= 0) call TLS_fatal (errmsg)
#if defined(GNU_PR49213)
          box%mesh => smesh
          call plist%set ('mesh', box)
#else
          call plist%set ('mesh', any_mesh(smesh))
#endif
          call smesh%write_profile
        else
          umesh => new_unstr_mesh(plist, stat, errmsg)
          if (stat /= 0) call TLS_fatal (errmsg)
#if defined(GNU_PR49213)
          box%mesh => umesh
          call plist%set ('mesh', box)
#else
          call plist%set ('mesh', any_mesh(umesh))
#endif
          call umesh%write_profile
          call umesh%check_bndry_face_set
#ifdef DEBUG_MESH
          call gmv_open (piter%name()//'.gmv')
          call gmv_write_unstr_mesh (umesh)
          call gmv_close
#endif
        end if
        call TLS_info ('  Mesh "' // trim(piter%name()) // '" initialized')
        call stop_timer('mesh-'//piter%name())
      end if
      call piter%next
    end do

  end subroutine init_mesh_manager

  !! This procedure reads the MESH and ALTMESH namelists and uses the data to
  !! initialize the module variable MESHES parameter list.

  subroutine read_truchas_mesh_namelists (lun)

    use altmesh_namelist
    use input_utilities, only: NULL_I

    integer, intent(in) :: lun

    type(parameter_list), pointer :: plist

    !! The MESH namelist: the "MAIN" mesh used by most physics
    plist => meshes%sublist('MAIN')
    call plist%set ('mesh', any_mesh())
    call plist%set ('enabled', .true.)  ! always enabled
    call read_mesh_namelist (lun, plist)

    !! The ALTMESH namelist: the "ALT" mesh used by EM
    call read_altmesh_namelist (lun)
    if (altmesh_exists) then
      plist => meshes%sublist('ALT')
      call plist%set ('mesh', any_mesh())
      call plist%set ('mesh-file', trim(altmesh_file))
      call plist%set ('coord-scale-factor', altmesh_coordinate_scale_factor)
      call plist%set ('rotation-angles', rotation_angles)
      call plist%set ('em-mesh', .true.)
      call plist%set ('partitioner', trim(partitioner))
      select case (partitioner)
      case ('file')
        call plist%set ('partition-file', trim(partition_file))
        call plist%set ('first-partition', first_partition)
#ifdef USE_METIS
      case ('metis')
        plist => plist%sublist('metis-options')
        if (metis_ctype   /= NULL_I) call plist%set('ctype',   metis_ctype)
        if (metis_ncuts   /= NULL_I) call plist%set('ncuts',   metis_ncuts)
        if (metis_niter   /= NULL_I) call plist%set('niter',   metis_niter)
        if (metis_ufactor /= NULL_I) call plist%set('ufactor', metis_ufactor)
        if (metis_minconn /= NULL_I) call plist%set('minconn', metis_minconn)
        if (metis_contig  /= NULL_I) call plist%set('contig',  metis_contig)
        if (metis_seed    /= NULL_I) call plist%set('seed',    metis_seed)
        if (metis_dbglvl  /= NULL_I) call plist%set('dbglvl',  metis_dbglvl)
#endif
      end select
    end if

  end subroutine read_truchas_mesh_namelists

  !! This auxiliary procedure reads the MESH namelist and stuffs the results
  !! into the passed parameter list PARAMS.

  subroutine read_mesh_namelist(lun, params)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
    use string_utilities, only: i_to_c, lower_case
    use truchas_env, only: input_dir
    use parallel_communication, only: is_IOP, broadcast

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer :: ios
    logical :: found
    character(80) :: iom
    type(parameter_list), pointer :: plist

    !! Namelist variables
    character(16)  :: partitioner
    character(511) :: mesh_file, partition_file
    real(r8) :: coordinate_scale_factor, rotation_angles(3)
    integer :: exodus_block_modulus, gap_element_blocks(50), interface_side_sets(127), first_partition
    namelist /mesh/ mesh_file, coordinate_scale_factor, rotation_angles, &
                    exodus_block_modulus, gap_element_blocks, interface_side_sets, &
                    partitioner, partition_file, first_partition

    !! Metis parameters
    integer :: metis_ctype, metis_ncuts, metis_niter, metis_ufactor, &
               metis_minconn, metis_contig, metis_seed, metis_dbglvl
    namelist /mesh/ metis_ctype, metis_ncuts, metis_niter, metis_ufactor, &
                    metis_minconn, metis_contig, metis_seed, metis_dbglvl

    !! Namelist variables for the internal mesh
    type :: coord_grid
      real(r8) :: coarse_grid(0:100)
      integer  :: intervals(100)
      real(r8) :: ratio(100)
    end type
    type(coord_grid) :: x_axis, y_axis, z_axis
    real(r8) :: noise_factor
    namelist /mesh/ x_axis, y_axis, z_axis, noise_factor

    call TLS_info('')
    call TLS_info('Reading MESH namelist ...')

    !! Locate the MESH namelist (required)
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'mesh', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0)  call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('MESH namelist not found')

    !! Default values
    mesh_file = NULL_C
    coordinate_scale_factor = 1.0_r8
    rotation_angles = 0
    exodus_block_modulus = 10000
    gap_element_blocks = NULL_I
    interface_side_sets = NULL_I
    partitioner = NULL_C
    partition_file = NULL_C
    first_partition = NULL_I
    call coord_grid_default(x_axis)
    call coord_grid_default(y_axis)
    call coord_grid_default(z_axis)
    noise_factor = NULL_R

    metis_ctype   = NULL_I
    metis_ncuts   = NULL_I
    metis_niter   = NULL_I
    metis_ufactor = NULL_I
    metis_minconn = NULL_I
    metis_contig  = NULL_I
    metis_seed    = -314159
    metis_dbglvl  = NULL_I

    !! Read the MESH namelist
    if (is_IOP) read(lun,nml=mesh,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading MESH namelist: ' // trim(iom))

    !! Broadcast the namelist variables
    call broadcast(mesh_file)
    call broadcast(coordinate_scale_factor)
    call broadcast(rotation_angles)
    call broadcast(exodus_block_modulus)
    call broadcast(gap_element_blocks)
    call broadcast(interface_side_sets)
    call broadcast(partitioner)
    call broadcast(partition_file)
    call broadcast(first_partition)
    call coord_grid_broadcast(x_axis)
    call coord_grid_broadcast(y_axis)
    call coord_grid_broadcast(z_axis)
    call broadcast(noise_factor)

    call broadcast(metis_ctype)
    call broadcast(metis_ncuts)
    call broadcast(metis_niter)
    call broadcast(metis_ufactor)
    call broadcast(metis_minconn)
    call broadcast(metis_contig)
    call broadcast(metis_seed)
    call broadcast(metis_dbglvl)

    !! Check and process the namelist variables, and stuff them into the return PARAMS.

    if (mesh_file /= NULL_C) then
      call external_mesh_input
    else
      call TLS_info('  MESH_FILE not specified; using inputs for an internally generated mesh')
      call internal_mesh_input
    end if

    if (coordinate_scale_factor <= 0.0_r8) call TLS_fatal('COORDINATE_SCALE_FACTOR must be > 0')
    call params%set('coord-scale-factor', coordinate_scale_factor)

    call params%set('rotation-angles', rotation_angles)

    if (partitioner == NULL_C) partitioner = 'chaco'
    select case (lower_case(partitioner))
    case ('chaco')
    case ('metis')
#ifdef USE_METIS
      plist => params%sublist('metis-options')
      if (metis_ctype   /= NULL_I) call plist%set('ctype',   metis_ctype)
      if (metis_ncuts   /= NULL_I) call plist%set('ncuts',   metis_ncuts)
      if (metis_niter   /= NULL_I) call plist%set('niter',   metis_niter)
      if (metis_ufactor /= NULL_I) call plist%set('ufactor', metis_ufactor)
      if (metis_minconn /= NULL_I) call plist%set('minconn', metis_minconn)
      if (metis_contig  /= NULL_I) call plist%set('contig',  metis_contig)
      if (metis_seed    /= NULL_I) call plist%set('seed',    metis_seed)
      if (metis_dbglvl  /= NULL_I) call plist%set('dbglvl',  metis_dbglvl)
#else
      call TLS_fatal('PARTITIONER = "metis" is not supported by this Truchas build')
#endif

    case ('block')
    case ('file')
      if (partition_file == NULL_C) call TLS_fatal('PARTITION_FILE not specified')
      if (partition_file(1:1) /= '/') partition_file = trim(input_dir) // trim(partition_file)
      if (is_IOP) inquire(file=partition_file,exist=found)
      call broadcast(found)
      if (.not.found) call TLS_fatal('PARTITION_FILE not found: ' // trim(partition_file))
      call params%set('partition-file', trim(partition_file))
      if (first_partition == NULL_I) first_partition = 0
      if (.not.any(first_partition == [0,1])) call TLS_fatal('FIRST_PARTITION must be 0 or 1')
      call params%set('first-partition', first_partition)
    case default
      call TLS_fatal('unknown value for PARTITIONER: ' // trim(partitioner))
    end select
    call params%set('partitioner', lower_case(partitioner))

  contains

    subroutine coord_grid_default(this)
      class(coord_grid), intent(inout) :: this
      this%coarse_grid = NULL_R
      this%intervals = NULL_I
      this%ratio = NULL_R
    end subroutine coord_grid_default

    subroutine coord_grid_broadcast(this)
      class(coord_grid), intent(inout) :: this
      call broadcast(this%coarse_grid)
      call broadcast(this%intervals)
      call broadcast(this%ratio)
    end subroutine coord_grid_broadcast

    subroutine external_mesh_input

      integer, allocatable :: iarray(:)

      if (mesh_file(1:1) /= '/') then ! not an absolute path
        mesh_file = trim(input_dir) // trim(mesh_file)
      end if
      if (is_IOP) inquire(file=mesh_file,exist=found)
      call broadcast(found)
      if (.not.found) call TLS_fatal('MESH_FILE not found: ' // trim(mesh_file))
      call params%set('mesh-file', trim(mesh_file))

      if (exodus_block_modulus < 0) call TLS_fatal('EXODUS_BLOCK_MODULUS must be >= 0')
      call params%set('exodus-block-modulus', exodus_block_modulus)

      iarray = pack(gap_element_blocks, mask=(gap_element_blocks /= NULL_I))
      if (size(iarray) > 0) call params%set('gap-element-block-ids', iarray)

      iarray = pack(interface_side_sets, mask=(interface_side_sets /= NULL_I))
      if (size(iarray) > 0) call params%set('interface-side-set-ids', iarray)

      !Temporary: provides data for the temporary get_main_mesh_size procedure
      block
        use exodus_truchas_hack, only: read_exodus_mesh_size
        integer :: ncells, nnodes, stat
        if (is_IOP) call read_exodus_mesh_size(trim(mesh_file), nnodes, ncells, stat)
        call broadcast(stat)
        if (stat /= 0) call TLS_fatal('error reading MESH_FILE "' // trim(mesh_file) // '"')
        call broadcast(nnodes)
        call broadcast(ncells)
        call params%set('nnodes', nnodes)
        call params%set('ncells', ncells)
      end block

    end subroutine external_mesh_input

    subroutine internal_mesh_input

      type(parameter_list), pointer :: plist

      plist => params%sublist('x-axis')
      call set_axis_params(plist, x_axis, 'X_AXIS')

      plist => params%sublist('y-axis')
      call set_axis_params(plist, y_axis, 'Y_AXIS')

      plist => params%sublist('z-axis')
      call set_axis_params(plist, z_axis, 'Z_AXIS')

      !Temporary: provides data for the temporary get_main_mesh_size procedure
      block
        integer :: array(3), ncells, nnodes
        array(1) = sum(x_axis%intervals, mask=(x_axis%intervals /= NULL_I))
        array(2) = sum(y_axis%intervals, mask=(y_axis%intervals /= NULL_I))
        array(3) = sum(z_axis%intervals, mask=(z_axis%intervals /= NULL_I))
        ncells = product(array)
        nnodes = product(array+1)
        call params%set('ncells', ncells)
        call params%set('nnodes', nnodes)
      end block

    end subroutine internal_mesh_input

    subroutine set_axis_params(params, axis, axis_name)

      type(parameter_list), intent(inout) :: params
      type(coord_grid), intent(in) :: axis
      character(*), intent(in) :: axis_name

      character(:), allocatable :: varname
      real(r8), allocatable :: rarray(:)
      integer,  allocatable :: iarray(:)
      integer :: n

      varname = axis_name // '%COARSE_GRID'
      if (all(axis%coarse_grid == NULL_R)) call TLS_fatal(varname // ' not specified')
      rarray = pack(axis%coarse_grid, mask=(axis%coarse_grid /= NULL_R))
      n = size(rarray)
      if (n < 2) then
        call TLS_fatal(varname // ' requires at least 2 values')
      else if (any(rarray(2:n) <= rarray(1:n-1))) then
        call TLS_fatal(varname // ' values must be strictly increasing')
      else
        call params%set('coarse-grid', rarray)
      end if

      n = n - 1 ! expected size for the following arrays

      varname = axis_name // '%INTERVALS'
      if (all(axis%intervals == NULL_I)) call TLS_fatal(varname // ' not specified')
      iarray = pack(axis%intervals, mask=(axis%intervals /= NULL_I))
      if (size(iarray) /= n) then
        call TLS_fatal(i_to_c(n) // ' values required for ' // varname)
      else if (any(iarray < 1)) then
        call TLS_fatal(varname // ' values must be > 0')
      else
        call params%set('intervals', iarray)
      end if

      varname = axis_name // '%RATIO'
      if (any(axis%ratio /= NULL_R)) then
        rarray = pack(axis%ratio, mask=(axis%ratio /= NULL_R))
        if (size(rarray) /= n) then
          call TLS_fatal(i_to_c(n) // ' values required for ' // varname)
        else if (any(rarray < 0)) then
          call TLS_fatal(varname // ' values must be > 0')
        else
          call params%set('ratio', rarray)
        end if
      end if

    end subroutine set_axis_params

  end subroutine read_mesh_namelist

  subroutine enable_mesh (name, exists)
    use string_utilities, only: raise_case
    character(*), intent(in) :: name
    logical, intent(out) :: exists
    type(parameter_list), pointer :: plist
    exists = meshes%is_sublist(raise_case(name))
    if (exists) then
      plist => meshes%sublist(raise_case(name))
      call plist%set ('enabled', .true.)
    end if
  end subroutine enable_mesh

  !! Returns a pointer to the CLASS(BASE_MESH) mesh object that corresponds to
  !! NAME. A null pointer is returned if NAME is not recognized or if the mesh
  !! that corresponds to NAME is not of class BASE_MESH.
  function named_mesh_ptr (name)
    use base_mesh_class
    character(*), intent(in) :: name
    class(base_mesh), pointer :: named_mesh_ptr
    class(*), pointer :: csptr
    csptr => mesh_ptr(name)
    named_mesh_ptr => null()
    if (.not.associated(csptr)) return
    select type (csptr)
    class is (base_mesh)
      named_mesh_ptr => csptr
    end select
  end function named_mesh_ptr

  !! Returns a pointer to the TYPE(UNSTR_MESH) mesh object that corresponds to
  !! NAME. A null pointer is returned if NAME is not recognized or if the mesh
  !! that corresponds to NAME is not of type UNSTR_MESH.
  function unstr_mesh_ptr (name)
    use unstr_mesh_type
    character(*), intent(in) :: name
    type(unstr_mesh), pointer :: unstr_mesh_ptr
    class(*), pointer :: csptr
    csptr => mesh_ptr(name)
    unstr_mesh_ptr => null()
    if (.not.associated(csptr)) return
    select type (csptr)
    class is (unstr_mesh)
      unstr_mesh_ptr => csptr
    end select
  end function unstr_mesh_ptr

  !! Returns a pointer to the TYPE(SIMPL_MESH) mesh object that corresponds to
  !! NAME. A null pointer is returned if NAME is not recognized or if the mesh
  !! that corresponds to NAME is not of type SIMPL_MESH.
  function simpl_mesh_ptr (name)
    use simpl_mesh_type
    character(*), intent(in) :: name
    type(simpl_mesh), pointer :: simpl_mesh_ptr
    class(*), pointer :: csptr
    csptr => mesh_ptr(name)
    simpl_mesh_ptr => null()
    if (.not.associated(csptr)) return
    select type (csptr)
    class is (simpl_mesh)
      simpl_mesh_ptr => csptr
    end select
  end function simpl_mesh_ptr

  !! Auxiliary function returns a CLASS(*) pointer to the arbitrary mesh object
  !! that corresponds to NAME.  NAME is case-insensitive.  A null pointer is
  !! returned if NAME is not recognized.
  function mesh_ptr (name)
    use string_utilities, only: raise_case
    character(*), intent(in) :: name
    class(*), pointer :: mesh_ptr
    type(parameter_list), pointer :: plist
    class(*), allocatable :: cs
    mesh_ptr => null()
    if (meshes%is_sublist(raise_case(name))) then
      plist => meshes%sublist(raise_case(name))
      call plist%get_any ('mesh', cs)
      select type (cs)  ! unbox the mesh pointer
      type is (any_mesh)
        mesh_ptr => cs%mesh
      class default
        INSIST(.false.)
      end select
    end if
  end function mesh_ptr

  !! NNC, Feb 2016.  This awful hack is necessitated by LIN_SOLVER_INPUT which
  !! initializes default solver parameters in the opaque Ubik structure array.
  !! Several of these parameters depend on the number of cells and nodes in the
  !! mesh, and it is doing this *before* the mesh has actually been read; i.e.
  !! it is doing it at the wrong time.  The legacy mesh accomodated this design
  !! flaw by pre-reading the sizes at the time the MESH namelist was read.  So
  !! we have to do the same thing for the new mesh.  A difference here is that
  !! we are storing this metadata with the mesh_manager and not the new mesh.

  subroutine get_main_mesh_size (ncells_tot, nnodes_tot)
    use parameter_list_type
    type(parameter_list), pointer :: plist
    integer, intent(out) :: ncells_tot, nnodes_tot
    if (meshes%is_sublist('MAIN')) then
      plist => meshes%sublist('MAIN')
      call plist%get ('ncells', ncells_tot)
      call plist%get ('nnodes', nnodes_tot)
    else
      ncells_tot = -1
      nnodes_tot = -1
    end if
  end subroutine get_main_mesh_size

end module mesh_manager
