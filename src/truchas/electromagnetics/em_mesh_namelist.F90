!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module em_mesh_namelist

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  implicit none
  private

  public :: em_mesh_namelist_exists, read_em_mesh_namelist

contains

  logical function em_mesh_namelist_exists(lun) result(found)
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist
    use string_utilities, only: i_to_c
    integer, intent(in) :: lun
    integer :: ios
    if (is_IOP) rewind(lun)
    if (is_IOP) call seek_to_namelist(lun, 'em_mesh', found, iostat=ios)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
  end function
  
  subroutine read_em_mesh_namelist(lun, params)

    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
    use string_utilities, only: i_to_c, lower_case
    use truchas_env, only: input_dir
    use parallel_communication, only: is_IOP, broadcast
    use parameter_list_type

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer :: ios
    logical :: found
    character(128) :: iom
    type(parameter_list), pointer :: plist

    !! Namelist variables
    character(16) :: partitioner
    character(511) :: mesh_file, partition_file
    real(r8) :: coord_scale_factor, rotation_angles(3)
    integer :: first_partition
    namelist /em_mesh/ mesh_file, coord_scale_factor, rotation_angles, &
        partitioner, partition_file, first_partition

    !! Metis parameters
    integer :: metis_ptype, metis_iptype, metis_ctype, metis_ncuts, metis_niter, &
               metis_ufactor, metis_minconn, metis_contig, metis_seed, metis_dbglvl
    namelist /em_mesh/ metis_ptype, metis_iptype, metis_ctype, metis_ncuts, metis_niter, &
                       metis_ufactor, metis_minconn, metis_contig, metis_seed, metis_dbglvl

    !! Namelist variables for the internal mesh
    type :: coord_grid
      real(r8) :: coarse_grid(0:100)
      integer  :: intervals(100)
      real(r8) :: ratio(100)
    end type
    type(coord_grid) :: x_axis, y_axis, z_axis
    real(r8) :: noise_factor
    namelist /em_mesh/ x_axis, y_axis, z_axis, noise_factor

    call TLS_info ('Reading EM_MESH namelist ...')

    if (is_IOP) rewind(lun)

    !! Locate the EM_MESH namelist (required)
    if (is_IOP) call seek_to_namelist(lun, 'em_mesh', found, iostat=ios)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('EM_MESH namelist not found')

    !! Default values
    mesh_file = NULL_C
    coord_scale_factor = NULL_R
    rotation_angles = NULL_R
    partitioner = NULL_C
    partition_file = NULL_C
    first_partition = NULL_I
    call coord_grid_default(x_axis)
    call coord_grid_default(y_axis)
    call coord_grid_default(z_axis)
    noise_factor = NULL_R

    metis_ptype   = NULL_I
    metis_iptype  = NULL_I
    metis_ctype   = NULL_I
    metis_ncuts   = NULL_I
    metis_niter   = NULL_I
    metis_ufactor = NULL_I
    metis_minconn = NULL_I
    metis_contig  = NULL_I
    metis_seed    = -314159
    metis_dbglvl  = NULL_I
    
    !! Read the EM_MESH namelist
    if (is_IOP) read(lun,nml=em_mesh,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading EM_MESH namelist: ' // trim(iom))

    !! Broadcast the namelist variables.
    call broadcast(mesh_file)
    call broadcast(coord_scale_factor)
    call broadcast(rotation_angles)
    call broadcast(partitioner)
    call broadcast(partition_file)
    call broadcast(first_partition)
    call coord_grid_broadcast(x_axis)
    call coord_grid_broadcast(y_axis)
    call coord_grid_broadcast(z_axis)
    call broadcast(noise_factor)

    call broadcast(metis_ptype)
    call broadcast(metis_iptype)
    call broadcast(metis_ctype)
    call broadcast(metis_ncuts)
    call broadcast(metis_niter)
    call broadcast(metis_ufactor)
    call broadcast(metis_minconn)
    call broadcast(metis_contig)
    call broadcast(metis_seed)
    call broadcast(metis_dbglvl)

    !! Check the input variables for errors.
    if (mesh_file /= NULL_C) then
      call external_mesh_input
    else
      call TLS_info('  MESH_FILE not specified; using inputs for an internally generated mesh')
      call internal_mesh_input
    end if

    if (coord_scale_factor /= NULL_R) then
      if (coord_scale_factor <= 0.0_r8) call TLS_fatal('COORD_SCALE_FACTOR must be > 0.0')
      call params%set('coord-scale-factor', coord_scale_factor)
    end if

    if (any(rotation_angles /= NULL_R)) then
      if (any(rotation_angles == NULL_R)) call TLS_fatal('ROTATION_ANGLES requires 3 values')
      call params%set('rotation-angles', rotation_angles)
    end if

    if (partitioner == NULL_C) partitioner = 'metis'
    select case (lower_case(partitioner))
    case ('metis')
      plist => params%sublist('metis-options')
      if (metis_ptype   /= NULL_I) call plist%set('ptype',   metis_ptype)
      if (metis_iptype  /= NULL_I) call plist%set('iptype',  metis_iptype)
      if (metis_ctype   /= NULL_I) call plist%set('ctype',   metis_ctype)
      if (metis_ncuts   /= NULL_I) call plist%set('ncuts',   metis_ncuts)
      if (metis_niter   /= NULL_I) call plist%set('niter',   metis_niter)
      if (metis_ufactor /= NULL_I) call plist%set('ufactor', metis_ufactor)
      if (metis_minconn /= NULL_I) call plist%set('minconn', metis_minconn)
      if (metis_contig  /= NULL_I) call plist%set('contig',  metis_contig)
      if (metis_seed    /= NULL_I) call plist%set('seed',    metis_seed)
      if (metis_dbglvl  /= NULL_I) call plist%set('dbglvl',  metis_dbglvl)
    case ('block')
    case ('file')
      if (partition_file == NULL_C) call TLS_fatal('PARTITION_FILE not specified')
      if (partition_file(1:1) /= '/') partition_file = trim(input_dir) // trim(partition_file)
      if (is_IOP) inquire(file=partition_file,exist=found)
      call broadcast(found)
      if (.not.found) call TLS_fatal('PARTITION_FILE not found: ' // trim(partition_file))
      call params%set('partition-file', trim(partition_file))
      if (first_partition /= NULL_I) then
        if (all(first_partition /= [0,1])) call TLS_fatal('FIRST_PARTITION must be 0 or 1')
        call params%set('first-partition', first_partition)
      end if
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

      if (mesh_file(1:1) /= '/') then ! not an absolute path
        mesh_file = trim(input_dir) // trim(mesh_file)
      end if
      if (is_IOP) inquire(file=mesh_file,exist=found)
      call broadcast(found)
      if (.not.found) call TLS_fatal('MESH_FILE not found: ' // trim(mesh_file))
      call params%set('mesh-file', trim(mesh_file))

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

  end subroutine read_em_mesh_namelist

end module em_mesh_namelist

