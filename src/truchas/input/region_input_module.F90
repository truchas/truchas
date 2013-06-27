MODULE REGION_INPUT_MODULE

  use kinds, only: r8
  implicit none
  private

  ! Public Subroutines
  public :: REGION_CHECK, REGION_DEFAULT, REGION_INPUT, REGION_READ

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE REGION_CHECK () 

    return

  END SUBROUTINE REGION_CHECK

  SUBROUTINE REGION_READ (lun)

    use region_data, only: nregion
    use parameter_module, only: mregion
    use parallel_info_module, only: p_info
    
    integer, intent(in) :: lun

    Integer :: ir
    logical :: region_namelist

    ! Initialize the number of bodies
    nregion = 0

    if (p_info%IOP) then
       ! Rewind to prepare for reading regions
       rewind lun
    endif

    ! read REGIONs
    do ir=1,mregion
       region_namelist = .false.
       call REGION_INPUT(lun, region_namelist)
       if (.not. region_namelist) exit
       call REGION_CHECK()
    end do

    return

  END SUBROUTINE REGION_READ

  SUBROUTINE REGION_DEFAULT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Default REGION namelist.
    !
    !=======================================================================
   
    return

  END SUBROUTINE REGION_DEFAULT

  SUBROUTINE REGION_INPUT (lun, region_namelist)

    use region_data,            only: x1,y1,z1,x2,y2,z2,flow_off,Regions, &
                                      nregion
    use input_utilities,        only: seek_to_namelist
    use parallel_info_module,   only: p_info

    ! Argument List
    integer, intent(in) :: lun
    logical :: region_namelist

    logical :: fatal
    integer :: ioerror

    namelist /REGION/ x1,y1,z1,x2,y2,z2,flow_off

    ! some default values for REGION
    x1 = 0
    y1 = 0
    z1 = 0

    x2 = 0
    y2 = 0
    z2 = 0

    flow_off = .false.

    IO_PE_ONLY: if (p_info%IOP) then
       call seek_to_namelist (lun, 'REGION', found=region_namelist)
       if (region_namelist) then
          read (lun, NML = REGION, IOSTAT = ioerror)
          if (ioerror /= 0) then ! If read error, then didn't read namelist
             region_namelist = .false.
             fatal         = .true.
          end if
       end if
    end if IO_PE_ONLY

    ! Broadcast data just read in to all PE's.
    call REGION_INPUT_PARALLEL (region_namelist)

    if (region_namelist) then

      nregion = nregion + 1

      Regions(nregion)%x1(1) = x1
      Regions(nregion)%x1(2) = y1
      Regions(nregion)%x1(3) = z1

      Regions(nregion)%x2(1) = x2
      Regions(nregion)%x2(2) = y2
      Regions(nregion)%x2(3) = z2

      Regions(nregion)%flow_off = flow_off

    endif

  END SUBROUTINE REGION_INPUT

  SUBROUTINE REGION_INPUT_PARALLEL (region_namelist)
    !======================================================================
    ! Purpose(s):
    !
    !   Broadcast items in BODY namelist from IO PE to all other PE's.
    !   Also broadcast body_namelist flag.
    !
    !======================================================================
    use region_data, only: x1,y1,z1,x2,y2,z2,flow_off
    use parallel_info_module, only: p_info
    use pgslib_module, only: PGSLIB_BCAST

    ! Argument List

    logical :: region_namelist

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Broadcast Data
    if (.NOT. p_info%UseGlobalServices) then
       call PGSLIB_BCAST (x1)     
       call PGSLIB_BCAST (y1)
       call PGSLIB_BCAST (z1)
       call PGSLIB_BCAST (x2)
       call PGSLIB_BCAST (y2)
       call PGSLIB_BCAST (z2)
       call PGSLIB_BCAST (flow_off)
       call PGSLIB_BCAST (region_namelist)
    end if

  END SUBROUTINE REGION_INPUT_PARALLEL

END MODULE REGION_INPUT_MODULE
