MODULE INTERFACE_OUTPUT_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define variables procedures for the output of interface
  !   information.
  !
  !   Public Interface:
  !
  !     * call INTERFACE_OUTPUT (Ro)
  !
  !        Write interface information out to the .int file: interface
  !        unit normals, the interface plane constant (intercept Ro),
  !        and the interface cell vertices.
  !
  ! Contains: INTERFACE_OUTPUT
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module, only: mops
  implicit none
  private

  public :: INTERFACE_OUTPUT

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! OUTPUTS namelist interface variables.
  integer, dimension(mops), save, public :: Int_Output_Dt_Multiplier
 
  ! Logical flags.
  logical, save, public :: interface_dump, time_for_int_dump
  
  integer, save :: int_lun = -1
  logical, save :: int_file_opened = .false.

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
CONTAINS

  SUBROUTINE INTERFACE_OUTPUT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Write interface information: Unit normals, the plane constant
    !   (intercept Rho), and the cell vertices containing the interface.
    !
    !=======================================================================
    use parallel_info_module, only: p_info
    use interface_module, only: Int_Geom
    use mesh_module,      only: orthogonal_mesh
    use parameter_module, only: ndim, nicells, nvc
    use pgslib_module,    only: PGSLib_GLOBAL_SUM, pgslib_collate
    use time_step_module, only: cycle_number, t
    use output_utilities, only: ANNOUNCE_FILE_WRITE
#ifdef SUPPORTS_NEWUNIT
    use truchas_env, only: output_file_name
#else
    use truchas_env, only: output_file_name, new_unit
#endif
 
    ! Local Variables
    integer  :: v, n, nicells_tot
    real(r8) :: dump_time
    real(r8), allocatable :: array(:)

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
    if (p_info%IOP .and. .not.int_file_opened) then
#ifdef SUPPORTS_NEWUNIT
      open(newunit=int_lun,file=output_file_name('int'),form='unformatted',action='write',status='replace')
#else
      call new_unit (int_lun)
      open(unit=int_lun,file=output_file_name('int'),form='unformatted',action='write',status='replace')
#endif
      int_file_opened = .true.
    end if

    nicells_tot = PGSLIB_Global_Sum(nicells)

    if (p_info%IOP) then
      allocate(array(nicells_tot))
    else
      allocate(array(0))
    end if

    ! Write time, cycle, and number of interface cells
    dump_time = t

    if (p_info%IOP) then
      write(int_lun) dump_time
      write(int_lun) cycle_number
      write(int_lun) nicells_tot
      write(int_lun) orthogonal_mesh
    end if

    ! Write interface cell vertices
    do v = 1,nvc
       do n = 1,ndim
          call pgslib_collate (array, INT_GEOM%Cell_Coord(n,v))
          if (p_info%IOP) write(int_lun) array
       end do
    end do

    ! Write interface intercepts and normals
    call pgslib_collate (array, INT_GEOM%Rho)
    if (p_info%IOP) write(int_lun) array


    do n = 1,ndim
       call pgslib_collate (array, INT_GEOM%Normal(n))
       if (p_info%IOP) write(int_lun) array
    end do

    call ANNOUNCE_FILE_WRITE ('Interface dump', output_file_name('int'))

    ! Set interface dump flag back to .false.
    !-mf now set in volume_track_module
    !-mf-jim Oct 12 2005 
    !time_for_int_dump = .false.

  END SUBROUTINE INTERFACE_OUTPUT

END MODULE INTERFACE_OUTPUT_MODULE
