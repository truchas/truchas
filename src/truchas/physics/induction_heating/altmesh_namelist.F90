!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module altmesh_namelist

  use kinds, only: r8
  implicit none
  private

  public :: read_altmesh_namelist

  logical, public :: altmesh_exists = .false.
  character(16), public :: partitioner
  character(511), public :: altmesh_file, grid_transfer_file, partition_file
  real(r8), public :: altmesh_coordinate_scale_factor, rotation_angles(3)
  integer, public :: first_partition
  character(16), public :: data_mapper_kind

  namelist /altmesh/ altmesh_file, altmesh_coordinate_scale_factor, rotation_angles, &
      grid_transfer_file, partitioner, partition_file, first_partition, data_mapper_kind

  !! Metis parameters
  integer, public :: metis_ptype, metis_iptype, metis_ctype, metis_ncuts, metis_niter, &
      metis_ufactor, metis_minconn, metis_contig, metis_seed, metis_dbglvl
  namelist /altmesh/ metis_ptype, metis_iptype, metis_ctype, metis_ncuts, metis_niter, &
                  metis_ufactor, metis_minconn, metis_contig, metis_seed, metis_dbglvl

contains

  subroutine read_altmesh_namelist(lun)

    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I
    use string_utilities, only: i_to_c, lower_case
    use truchas_env, only: input_dir
    use parallel_communication, only: is_IOP, broadcast
    use truchas_logging_services

    integer, intent(in) :: lun

    logical :: found
    integer :: ios
    character(128) :: iom

    !! Locate the ALTMESH namelist (optional)
    altmesh_exists = .false.
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'ALTMESH', found=altmesh_exists, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(altmesh_exists)
    if (.not.altmesh_exists) return

    !! Read the ALTMESH namelist, assigning default values first.
    call TLS_info ('Reading ALTMESH namelist ...')
    if (is_IOP) then
      altmesh_file = NULL_C
      altmesh_coordinate_scale_factor = 1.0_r8
      rotation_angles = 0
      grid_transfer_file = NULL_C
      partitioner = NULL_C
      partition_file = NULL_C
      first_partition = NULL_I
      data_mapper_kind = 'default'
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
      read(lun,nml=altmesh,iostat=ios,iomsg=iom)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading ALTMESH namelist: ' // trim(iom))

    !! Broadcast the namelist variables.
    call broadcast(altmesh_exists)
    call broadcast(altmesh_file)
    call broadcast(altmesh_coordinate_scale_factor)
    call broadcast(rotation_angles)
    call broadcast(grid_transfer_file)
    call broadcast(partitioner)
    call broadcast(partition_file)
    call broadcast(first_partition)
    call broadcast(data_mapper_kind)
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
    if (altmesh_file == NULL_C) call TLS_fatal ('ALTMESH_FILE not specified')
    if (altmesh_file(1:1) /= '/') altmesh_file = trim(input_dir) // trim(altmesh_file)
    if (is_IOP) inquire(file=altmesh_file,exist=found)
    call broadcast (found)
    if (.not.found) call TLS_fatal ('ALTMESH_FILE not found: ' // trim(altmesh_file))

    if (grid_transfer_file /= NULL_C) then
      grid_transfer_file = adjustl(grid_transfer_file)
      if (grid_transfer_file(1:1) /= '/') then
        grid_transfer_file = trim(input_dir) // trim(grid_transfer_file)
      end if
    end if

    if (altmesh_coordinate_scale_factor <= 0.0_r8) &
        call TLS_fatal('ALTMESH_COORDINATE_SCALE_FACTOR must be > 0')

    if (partitioner == NULL_C) partitioner = 'chaco'
    select case (lower_case(partitioner))
    case ('chaco')
    case ('metis')
#ifndef USE_METIS
      call TLS_fatal('PARTITIONER = "metis" is not supported by this Truchas build')
#endif
    case ('block')
    case ('file')
      if (partition_file == NULL_C) call TLS_fatal('PARTITION_FILE not specified')
      if (partition_file(1:1) /= '/') partition_file = trim(input_dir) // trim(partition_file)
      if (is_IOP) inquire(file=partition_file,exist=found)
      call broadcast(found)
      if (.not.found) call TLS_fatal('PARTITION_FILE not found: ' // trim(partition_file))
      if (first_partition == NULL_I) first_partition = 0
      if (.not.any(first_partition == [0,1])) call TLS_fatal ('FIRST_PARTITION must be 0 or 1')
    case default
      call TLS_fatal ('unknown value for PARTITIONER: ' // trim(partitioner))
    end select

    select case (data_mapper_kind)
    case ('default')
    case ('portage')
#ifndef USE_PORTAGE
      call TLS_fatal('DATA_MAPPER_KIND = "portage" is not supported by this Truchas build')
#endif
    case default
      call TLS_fatal('invalid value for DATA_MAPPER_KIND: ' // trim(data_mapper_kind))
    end select

  end subroutine read_altmesh_namelist

end module altmesh_namelist

