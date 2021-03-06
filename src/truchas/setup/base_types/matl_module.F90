!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE MATL_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define the MATL and MATL_SLOT derived types and the procedures
  !   necessary to allocate/deallocate them.
  !
  !   Public Interface(s):
  !
  !     * call GATHER_VOF (m, Vof), GATHER_VOF_OLD (m, Vof_Old)
  !         Gather material m volume fractions from Matl; place them in Vof.
  !
  !     * call SCATTER_VOF (m, Vof)
  !         Scatter material m volume fractions in Vof into Matl.
  !
  !     * call SLOT_DECREASE (ABC, slot, slot_new)
  !         Decrease the ABC size from slot to slot_new.
  !
  !     * call SLOT_INCREASE (ABC, slot, slot_new)
  !         Increase the ABC size from slot to slot_new.
  !
  !     * call SLOT_COMPRESS (ABC, slot)
  !         Compress the slot structure of ABC(slot) to remove zeroes
  !
  !     * call SLOT_SET (ABC, s)
  !         Initialize slot s of ABC.
  !
  ! Contains: GATHER_VOF, GATHER_VOF_OLD
  !           SCATTER_VOF
  !           SLOT_DECREASE
  !           SLOT_INCREASE
  !           SLOT_COMPRESS
  !           SLOT_SET
  !
  ! Author(s): Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module, only: max_relation_forms, max_slots
  use truchas_logging_services
  implicit none
  private

  ! Public Variables
  public :: MATERIAL, MATL_SLOT, MATL, COLLATE
  ! Public Procedures
  public :: GATHER_VOF, SCATTER_VOF, &
            SLOT_DECREASE, SLOT_INCREASE, SLOT_COMPRESS,      &
            SLOT_SET, GATHER_VOF_OLD, SLOT_RESIZE

  INTERFACE COLLATE
     MODULE PROCEDURE COLLATE_MATL
  END INTERFACE

  interface ASSIGNMENT(=)
     module procedure MATL_ASSIGN
  end interface

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Variables needed in processing the namelist input
  character(LEN = 80), dimension(max_relation_forms), public, save :: &
                       Relation_Forms, P_Change_Forms

  ! Define MATERIAL Structure
  type MATERIAL
     ! Material Identifier
     integer :: Id = 0

     ! Volume Fraction
     real(r8)   :: Vof = 0.0_r8

     ! Old time step Volume Fraction
     real(r8)   :: Vof_Old = 0.0_r8

  end type MATERIAL

  ! Define MATL_SLOT Structure
  type MATL_SLOT
     type(MATERIAL), dimension(:), pointer :: Cell => NULL()
  end type MATL_SLOT

  type(MATL_SLOT), dimension(max_slots), save, Target  :: Matl

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE GATHER_VOF (m, Vof)
    !=======================================================================
    ! Purpose(s):
    !   Gather material m volume fractions from the Matl derived
    !   type and place them into array Vof.
    !=======================================================================
    use parameter_module, only: mat_slot

    ! Argument List
    integer,  intent(IN)  :: m
    real(r8), intent(OUT) :: Vof(:)

    ! Local Variables
    integer :: s

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize
    Vof = 0.0_r8

    ! Loop over slots, gathering material m volume fractions
    do s = 1,mat_slot
       where (Matl(s)%Cell%Id == m) Vof = Matl(s)%Cell%Vof
    end do

  END SUBROUTINE GATHER_VOF

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE GATHER_VOF_OLD (m, Vof_Old)
    !=======================================================================
    ! Purpose(s):
    !   Gather material m volume fractions from the Matl derived
    !   type and place them into array Vof.
    !=======================================================================
    use parameter_module, only: mat_slot
    use legacy_mesh_api,  only: ncells

    ! Argument List
    integer, intent(IN) :: m
    real(r8), dimension(ncells), intent(OUT) :: Vof_Old

    ! Local Variables
    integer :: s

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize
    Vof_Old = 0.0_r8

    ! Loop over slots, gathering material m volume fractions
    do s = 1,mat_slot
       where (Matl(s)%Cell%Id == m) Vof_Old = Matl(s)%Cell%Vof_Old
    end do

  END SUBROUTINE GATHER_VOF_OLD

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SCATTER_VOF (m, Vof)
    !=======================================================================
    ! Purpose(s):
    !   Scatter material m volume fractions in array Vof into the
    !   appropriate slots of the Matl derived type.
    !=======================================================================
    use parameter_module, only: mat_slot
    use legacy_mesh_api,  only: ncells

    ! Argument List
    integer, intent(IN) :: m
    real(r8), dimension(ncells), intent(IN) :: Vof

    ! Local Variables
    integer :: s

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Loop over slots, scattering material m volume fractions
    do s = 1,mat_slot
       where (Matl(s)%Cell%Id == m) Matl(s)%Cell%Vof = Vof
    end do

  END SUBROUTINE SCATTER_VOF

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SLOT_DECREASE (ABC, slot, slot_new)
    !=======================================================================
    ! Purpose(s):
    !   Decrease the ABC structure size from slot to slot_new.
    !=======================================================================
    use parameter_module, only: max_slots

    ! Argument List
    type(MATL_SLOT), dimension(max_slots), intent(INOUT) :: ABC
    integer, intent(INOUT) :: slot, slot_new

    ! Local Variables
    integer :: s
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Check to make sure slot_new < slot and >= 0; if not, terminate
    if (slot_new >= slot .or. slot_new < 0) then
       write (message, 10) slot_new, slot
10     format ('SLOT_DECREASE: unacceptable value of ',i2,' for slot_new; must be < ',i2,' and >= 0')
       call TLS_panic (message)
    end if

    ! Deallocate a mesh size array for each appropriate slot
    do s = slot,slot_new+1,-1

       if (ASSOCIATED(ABC(s)%Cell)) DEALLOCATE (ABC(s)%Cell)

    end do

    slot = slot_new

  END SUBROUTINE SLOT_DECREASE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SLOT_INCREASE (ABC, slot, slot_new)
    !=======================================================================
    ! Purpose(s):
    !   Increase the ABC structure size from slot to slot_new.
    !=======================================================================
    use parameter_module, only: max_slots
    use legacy_mesh_api,  only: ncells

    ! Argument List
    type(MATL_SLOT), dimension(max_slots), intent(INOUT) :: ABC
    integer, intent(INOUT) :: slot, slot_new

    ! Local Variables
    integer :: memstat, s
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Check to make slot < slot_new <= max_slots; if not, terminate
    if (slot >= slot_new .or. slot_new > max_slots) then
       write (message, 10) slot_new, slot, max_slots
10     format ('SLOT_INCREASE: unacceptable value of ',i2,' for slot_new; must be > ',i2,' and <= ',i2)
       call TLS_panic (message)
    end if

    ! Bump the starting slot number by 1
    slot = slot + 1

    ! Allocate a mesh size array for each appropriate slot s
    do s = slot,slot_new

       ALLOCATE (ABC(s)%Cell(ncells), STAT = memstat)
       if (memstat /= 0) call TLS_panic ('SLOT_INCREASE: memory allocation error')

       ! Initialize the mesh size array
       call SLOT_SET (ABC, s)

    end do

    slot = slot_new

  END SUBROUTINE SLOT_INCREASE
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  subroutine slot_resize (abc, slot, slot_new)

    use legacy_mesh_api, only: ncells

    type(matl_slot), intent(inout) :: abc(:)
    integer,         intent(inout) :: slot
    integer,         intent(in)    :: slot_new

    integer :: s
    character(128) :: message

    !! Check for invalid values for SLOT_NEW.
    if (slot_new < 0 .or. slot_new >= size(abc)) then
       write (message, 10) slot_new, size(abc)
10     format ('SLOT_RESIZE: unacceptable value of ',i2,' for slot_new; must be >= 0 and <= ',i2)
       call TLS_panic (message)
    end if

    if (slot_new < slot) then
      do s = slot_new+1, slot
        deallocate(abc(s)%cell)
      end do
    else if (slot_new > slot) then
      do s = slot+1, slot_new
        allocate(abc(s)%cell(ncells))
        call slot_set (abc, s)
      end do
    end if

    slot = slot_new

  end subroutine slot_resize

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SLOT_COMPRESS (ABC, slot)
    !=======================================================================
    ! Purpose(s):
    !
    !   This routine compresses the slot structure to remove unused holes.
    !
    !=======================================================================
    use parameter_module, only: max_slots
    use legacy_mesh_api,  only: ncells

    ! Argument List
    type(MATL_SLOT), dimension(max_slots), intent(INOUT) :: ABC
    integer, intent(INOUT) :: slot

    ! Local Variables
    integer :: s_old, s_new
    logical, dimension(ncells) :: Move_slot
    integer, dimension(ncells) :: Num_mat 

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Zero the starting number of materials in each cell.
    Num_mat = 0

    ! Examine each slot to determine if it contains a material that should
    ! have its storage compressed.
    do s_old = 1,slot

       Move_slot = ABC(s_old)%Cell%Id /= 0

       ! If this slot has a valid material, increment the number of
       ! materials in this cell.
       where (Move_slot) Num_mat = Num_mat + 1

       ! If it does have a material, put its properties in the slot
       ! corresponding to the current number of materials in this cell.
       do s_new = 1,slot
          if (s_new == s_old) cycle
          where (Move_slot .and. s_new == Num_mat)
             ABC(s_new)%Cell%Id       = ABC(s_old)%Cell%Id
             ABC(s_new)%Cell%Vof      = ABC(s_old)%Cell%Vof
             ABC(s_new)%Cell%Vof_Old  = ABC(s_old)%Cell%Vof_Old
             ABC(s_old)%Cell%Id       = 0
             ABC(s_old)%Cell%Vof      = 0.0_r8
             ABC(s_old)%Cell%Vof_Old  = 0.0_r8
          end where
       end do

    end do

  END SUBROUTINE SLOT_COMPRESS

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SLOT_SET (ABC, s)
    !=======================================================================
    ! Purpose(s):
    !   Initialize the ABC structure for slot s.
    !=======================================================================
    use parameter_module, only: max_slots

    ! Argument List
    type(MATL_SLOT), dimension(max_slots), intent(INOUT) :: ABC
    integer, intent(IN) :: s

    ! Local Variables
    integer :: Id_def  = 0
    real(r8)   :: Vof_def = 0.0_r8

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Assign defaults
    ABC(s)%Cell%Id      = Id_def
    ABC(s)%Cell%Vof     = Vof_def
    ABC(s)%Cell%Vof_Old = Vof_def

  END SUBROUTINE SLOT_SET

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE COLLATE_MATL (Collated_Matl, Local_Matl)
    !==================================================================
    ! Purpose(s):
    !   Collate a distributed matl into a single large matl on IO PE
    !==================================================================
    use parameter_module,     only: mat_slot, max_slots
    use pgslib_module,        only: PGSLib_COLLATE

    ! Arguments
    type(MATL_SLOT), dimension(max_slots), intent(INOUT) :: Collated_Matl
    type(MATL_SLOT), dimension(max_slots), intent(IN   ) :: Local_Matl
    
    ! Local variables
    integer :: s

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do s = 1, mat_slot
       call PGSLib_COLLATE(Collated_Matl(s)%Cell%Id,      Local_Matl(s)%Cell%Id )
       call PGSLib_COLLATE(Collated_Matl(s)%Cell%Vof,     Local_Matl(s)%Cell%Vof )
       call PGSLib_COLLATE(Collated_Matl(s)%Cell%Vof_Old, Local_Matl(s)%Cell%Vof_Old )
    end do

  END SUBROUTINE COLLATE_MATL
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE MATL_ASSIGN (Matl_Dest, Matl_Src)
    !====================================================================
    ! Purposes(s):
    !   Assign matl data structure across =.  This involves
    !   copy of appropriate Cell pointers.
    !
    !   WARNING: This routine assumes that matl data structures have
    !   max_slots slots, of which mat_slot are used.
    !
    !====================================================================
    use parameter_module, only: mat_slot, max_slots

    ! Arguments
    type(MATL_SLOT), dimension(max_slots), intent(IN)  :: Matl_Src
    type(MATL_SLOT), dimension(max_slots), intent(INOUT) :: Matl_Dest

    ! Local variables
    integer :: s

    do s = 1, mat_slot
       Matl_Dest(s)%Cell = Matl_Src(s)%Cell
    end do

  END SUBROUTINE MATL_ASSIGN

END MODULE MATL_MODULE
