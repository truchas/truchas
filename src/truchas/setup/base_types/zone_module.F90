!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ZONE_MODULE
  !=======================================================================
  ! Purpose:
  !    Define the CELL_AVG derived type
  !
  ! Contains:
  !    STREAM_ZONE_OUT (Stream, Zone)
  !    ZONE_COLLATE (Zone)
  !    STREAM_ZONE_IN (Stream, Zone)
  !    ZONE_DISTRIBUTE (Zone_Tot)
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !=======================================================================
  use kinds, only: r8
  use legacy_mesh_api, only: ndim
  implicit none
  private

  ! public variables and types
  public :: CELL_AVG, Zone, COLLATE, PERMUTE_ZONE, read_zone_data

  INTERFACE COLLATE
     MODULE PROCEDURE COLLATE_ZONE
  END INTERFACE

  !-----------------------------------------------------------------------------

  ! CELL_AVG structure - hold physical state variables
  type CELL_AVG
     real(r8)                      :: Rho       ! current cell average density
     real(r8)                      :: Rho_Old   ! past cell average density
     real(r8)                      :: Temp      ! current temperature
     real(r8)                      :: Temp_Old  ! past temperature
     real(r8)                      :: Enthalpy      ! current enthalpy
     real(r8)                      :: Enthalpy_Old  ! past enthalpy
     real(r8)                      :: P         ! current pressure
     real(r8), dimension(ndim)     :: Vc        ! current cell-centered velocity
     real(r8), dimension(ndim)     :: Vc_old    ! past cell-centered velocity
  end type CELL_AVG

  type(CELL_AVG), dimension(:), pointer :: Zone

CONTAINS

  FUNCTION ZONE_COLLATE (Zone)
    !==================================================================
    ! Purpose(s):
    !   Collate a distributed zone into a single large zone on IO PE
    !==================================================================
    use parallel_communication, only: is_IOP
    use legacy_mesh_api, only: ncells_tot, ncells

    ! Arguments
    type(CELL_AVG),          dimension(ncells), intent(IN) :: Zone
    type(CELL_AVG), pointer, dimension(:)                  :: ZONE_COLLATE

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (is_IOP) then
       ALLOCATE(Zone_Collate(ncells_tot))
    else
       ALLOCATE(Zone_Collate(0))
    end if

    call COLLATE(Zone_Collate, Zone)

  END FUNCTION ZONE_COLLATE

  SUBROUTINE COLLATE_ZONE (Collated_Zone, Local_Zone)
    !==================================================================
    ! Purpose(s):
    !   Collate a distributed zone into a single large zone on IO PE
    !==================================================================
    use parallel_communication,        only: collate

    ! Arguments
    type(CELL_AVG), dimension(:), intent(IN    ) :: Local_Zone
    type(CELL_AVG), dimension(:), intent(   OUT) :: Collated_Zone

    ! Local variables
    integer :: n

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    call collate (Collated_Zone%Rho,          Local_Zone%Rho)
    call collate (Collated_Zone%Rho_Old,      Local_Zone%Rho_Old)
    call collate (Collated_Zone%Temp,         Local_Zone%Temp)
    call collate (Collated_Zone%Temp_Old,     Local_Zone%Temp_Old)
    call collate (Collated_Zone%Enthalpy,     Local_Zone%Enthalpy)
    call collate (Collated_Zone%Enthalpy_Old, Local_Zone%Enthalpy_Old)
    call collate (Collated_Zone%P,            Local_Zone%P)

    do n = 1,ndim
       call collate (Collated_Zone%Vc(n),      Local_Zone%Vc(n))
       call collate (Collated_Zone%Vc_Old(n),  Local_Zone%Vc_Old(n))
    end do

  END SUBROUTINE COLLATE_ZONE

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_ZONE_DATA
 !!
 !! Neil N. Carlson <nnc@lanl.gov>
 !! 6 Jan 2005, revised 21 Apr 2005
 !!
 !! This subroutine reads the zone data from a restart file opened (and pre-
 !! positioned) on UNIT, and initializes the module structure ZONE with this
 !! data (properly distributed and permuted).  VERSION is the version number
 !! of the restart file format.
 !!
 !! NB: We have no data to initialize the *_OLD components of the zone
 !! structure, so we simply set those to be the current values.  We wouldn't
 !! normally need to do this (they are overwritten at the start of a time
 !! step) except the output routine wants to write them at the start :-(
 !!

  subroutine read_zone_data (unit, version)

    use legacy_mesh_api, only: pcell => unpermute_mesh_vector
    use restart_utilities, only: read_dist_array

    integer, intent(in) :: unit, version

    integer :: n

    call read_dist_array (unit, zone%rho,      pcell, 'READ_ZONE_DATA: error reading RHO record')
    call read_dist_array (unit, zone%temp,     pcell, 'READ_ZONE_DATA: error reading TEMP record')
    call read_dist_array (unit, zone%enthalpy, pcell, 'READ_ZONE_DATA: error reading ENTHALPY record')
    call read_dist_array (unit, zone%p,        pcell, 'READ_ZONE_DATA: error reading PRESSURE record')

    zone%rho_old      = zone%rho
    zone%temp_old     = zone%temp
    zone%enthalpy_old = zone%enthalpy

    do n = 1, ndim
      call read_dist_array (unit, zone%vc(n), pcell, 'READ_ZONE_DATA: error reading VC records')
      zone%vc_old(n) = zone%vc(n)
    end do

  end subroutine read_zone_data

  SUBROUTINE PERMUTE_ZONE (Permuted_Zone, Orig_Zone, Permuter, SCOPE)
    !==================================================================
    ! Purpose(s):
    !   Permute zone according the Permuter vector
    !==================================================================
    use parallel_scope

    ! Arguments
    type(CELL_AVG), dimension(:), intent(IN   ) :: Orig_Zone
    type(CELL_AVG), dimension(:), intent(  OUT) :: Permuted_Zone
    integer, dimension(:), intent(IN) :: Permuter
    type (PL_SCOPE), OPTIONAL,    intent(IN   ) :: SCOPE

    ! Local variables
    type (PL_SCOPE) :: Desired_Scope

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Default scope is global
    if (PRESENT(SCOPE)) then
       Desired_Scope = SCOPE
    else
       Desired_Scope = GLOBAL_SCOPE
    end if

    if (DESIRED_SCOPE == GLOBAL_SCOPE) then
       call PERMUTE_ZONE_GLOBAL(Permuted_Zone, Orig_Zone, Permuter)
    end if

    if (DESIRED_SCOPE == LOCAL_SCOPE) then
       call PERMUTE_ZONE_LOCAL (Permuted_Zone, Orig_Zone, Permuter)
    end if

  end SUBROUTINE PERMUTE_ZONE

  SUBROUTINE PERMUTE_ZONE_GLOBAL (Permuted_Zone, Orig_Zone, Permuter)
    !==================================================================
    ! Purpose(s):
    !   Permute zone according the Permuter vector, global version
    !==================================================================
    use pgslib_module,    only: PGSLib_Permute,    &
                                PGSLIB_Deallocate_Trace, &
                                PGSLib_GS_Trace

    ! Arguments
    type(CELL_AVG), dimension(:), intent(IN   ) :: Orig_Zone
    type(CELL_AVG), dimension(:), intent(  OUT) :: Permuted_Zone
    integer, dimension(:), intent(IN) :: Permuter

    ! Local variables
    integer :: n
    type (PGSLib_GS_Trace), POINTER :: Zone_Trace

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    NULLIFY(Zone_Trace)

    call PGSLib_PERMUTE (DEST   = Permuted_Zone%Rho,          &
                         SOURCE = Orig_Zone%Rho,              &
                         INDEX  = Permuter,       &
                         TRACE  = Zone_Trace)
    call PGSLib_PERMUTE (DEST   = Permuted_Zone%Rho_Old,      &
                         SOURCE = Orig_Zone%Rho_Old,              &
                         INDEX  = Permuter,       &
                         TRACE  = Zone_Trace)
    call PGSLib_PERMUTE (DEST   = Permuted_Zone%Temp,         &
                         SOURCE = Orig_Zone%Temp,                 &
                         INDEX  = Permuter,       &
                         TRACE  = Zone_Trace)
    call PGSLib_PERMUTE (DEST   = Permuted_Zone%Temp_Old,     &
                         SOURCE = Orig_Zone%Temp_Old,             &
                         INDEX  = Permuter,       &
                         TRACE  = Zone_Trace)
    call PGSLib_PERMUTE (DEST   = Permuted_Zone%Enthalpy,     &
                         SOURCE = Orig_Zone%Enthalpy,             &
                         INDEX  = Permuter,       &
                         TRACE  = Zone_Trace)
    call PGSLib_PERMUTE (DEST   = Permuted_Zone%Enthalpy_Old, &
                         SOURCE = Orig_Zone%Enthalpy_Old,         &
                         INDEX  = Permuter,       &
                         TRACE  = Zone_Trace)
    call PGSLib_PERMUTE (DEST   = Permuted_Zone%P,            &
                         SOURCE = Orig_Zone%P,                    &
                         INDEX  = Permuter,       &
                         TRACE  = Zone_Trace)

    do n = 1,ndim
       call PGSLib_PERMUTE (DEST   = Permuted_Zone%Vc(n),   &
                            SOURCE = Orig_Zone%Vc(n),           &
                            INDEX  = Permuter,  &
                            TRACE  = Zone_Trace)
       call PGSLib_PERMUTE (DEST   = Permuted_Zone%Vc_Old(n),   &
                            SOURCE = Orig_Zone%Vc_Old(n),       &
                            INDEX  = Permuter,  &
                            TRACE  = Zone_Trace)
    end do

    ! Done with the trace
    call PGSLib_DEALLOCATE_TRACE (Zone_Trace)

  END SUBROUTINE PERMUTE_ZONE_GLOBAL

  SUBROUTINE PERMUTE_ZONE_LOCAL (Permuted_Zone, Orig_Zone, Permuter)
    !==================================================================
    ! Purpose(s):
    !   Permute zone according the Permuter vector, local version
    !   The Permuter vector refers to local indices.
    !   The input and output vectors must have the same size.
    !==================================================================

    ! Arguments
    type(CELL_AVG), dimension(:), intent(IN   ) :: Orig_Zone
    type(CELL_AVG), dimension(:), intent(  OUT) :: Permuted_Zone
    integer, dimension(:), intent(IN) :: Permuter

    ! Local variables
    integer :: cell

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do cell = 1, SIZE(Permuter)
       Permuted_Zone(Permuter(cell)) = Orig_Zone(Cell)
    end do

  END SUBROUTINE PERMUTE_ZONE_LOCAL

END MODULE ZONE_MODULE
