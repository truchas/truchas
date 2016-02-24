!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module BC_SET_MODULE
  !-----------------------------------------------------------------------------
  ! Purpose:
  !   define procedures to control the application of boundary conditions
  !
  ! Contains: InitializeBoundaryConditions
  !           FinalizeBoundaryConditions
  !           SetBoundaryCondition
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !            Bryan R. Lally, LANL ESA-EPE (lally@lanl.gov)
  !-----------------------------------------------------------------------------
  Implicit None

  Private

  ! public procedures
  Public :: InitializeBoundaryConditions, &
     FinalizeBoundaryConditions,   &
     SetBoundaryCondition

  !-----------------------------------------------------------------------------

Contains

  subroutine InitializeBC_Struct(BC_Struct)

    ! Initialize a BC_Structure
    use bc_type_module,       only: BC_STRUCTURE

    type(BC_STRUCTURE), intent(OUT) :: BC_Struct
    BC_Struct%ID =''
    BC_Struct%Exist = -1
    BC_Struct%IFlag = -1
    BC_Struct%Shift = -1
    NULLIFY(BC_Struct%BCID)
    NULLIFY(BC_Struct%BCOperator)
    NULLIFY(BC_Struct%Flag)
    NULLIFY(BC_Struct%UseFunc)
    NULLIFY(BC_Struct%Value1)
    NULLIFY(BC_Struct%Value2)
    NULLIFY(BC_Struct%Value3)
    NULLIFY(BC_Struct%Value_Array)
  END subroutine InitializeBC_Struct

  subroutine DeallocateBC_Struct(BC_Struct)

    ! Deallocate a BC_Structure
    use bc_type_module,       only: BC_STRUCTURE

    type(BC_STRUCTURE), intent(INOUT) :: BC_Struct

    IF (ASSOCIATED(BC_Struct%BCID))          DEALLOCATE(BC_Struct%BCID)
    IF (ASSOCIATED(BC_Struct%BCOperator))    DEALLOCATE(BC_Struct%BCOperator)
    IF (ASSOCIATED(BC_Struct%Flag))          DEALLOCATE(BC_Struct%Flag)
    IF (ASSOCIATED(BC_Struct%UseFunc))       DEALLOCATE(BC_Struct%UseFunc)
    IF (ASSOCIATED(BC_Struct%Value1))        DEALLOCATE(BC_Struct%Value1)
    IF (ASSOCIATED(BC_Struct%Value2))        DEALLOCATE(BC_Struct%Value2)
    IF (ASSOCIATED(BC_Struct%Value3))        DEALLOCATE(BC_Struct%Value3)
    IF (ASSOCIATED(BC_Struct%Value_Array))   DEALLOCATE(BC_Struct%Value_Array)

    call InitializeBC_Struct(BC_Struct)

  END subroutine DeallocateBC_Struct

  !-----------------------------------------------------------------------------

  Subroutine InitializeBoundaryConditions (tflag, cflag, fflag)
    !---------------------------------------------------------------------------
    ! Purpose:
    !   initialize boundary conditions
    !---------------------------------------------------------------------------
    use bc_data_module,       only: BC_C_EXIST, BC_C_IFLAG, BC_C_SHIFT, &
                                    BC_P_EXIST, BC_P_IFLAG, BC_P_SHIFT, &
                                    BC_T_EXIST, BC_T_IFLAG, BC_T_SHIFT, &
                                    BC_V_EXIST, BC_V_IFLAG, BC_V_SHIFT
    use bc_type_module,       only: BC_C, BC_P, BC_T, BC_V
    use legacy_mesh_api,      only: ncells, nfc

    ! Argument List
    logical, intent(IN) :: cflag ! Concentration flag
    logical, intent(IN) :: fflag ! Fluid flow flag
    logical, intent(IN) :: tflag ! Temperature flag

    !---------------------------------------------------------------------------

    ! boundary conditions
    ! Concentration
    call InitializeBC_Struct(BC_C)
    BC_C%ID    = 'C'
    BC_C%Exist = BC_C_EXIST
    BC_C%IFlag = BC_C_IFLAG
    BC_C%Shift = BC_C_SHIFT

!!$    NULLIFY(BC_C%Flag)
!!$    NULLIFY(BC_C%Value1)
!!$    NULLIFY(BC_C%Value2)
!!$    NULLIFY(BC_C%Value3)
!!$    NULLIFY(BC_C%UseFunc)

    ! Pressure
    call InitializeBC_Struct(BC_P)
    BC_P%ID    = 'P'
    BC_P%Exist = BC_P_EXIST
    BC_P%IFlag = BC_P_IFLAG
    BC_P%Shift = BC_P_SHIFT

!!$    NULLIFY(BC_P%Flag)
!!$    NULLIFY(BC_P%Value1)
!!$    NULLIFY(BC_P%Value2)
!!$    NULLIFY(BC_P%Value3)
!!$    NULLIFY(BC_P%UseFunc)

    ! Temperature
    call InitializeBC_Struct(BC_T)
    BC_T%ID    = 'T'
    BC_T%Exist = BC_T_EXIST
    BC_T%IFlag = BC_T_IFLAG
    BC_T%Shift = BC_T_SHIFT

!!$    NULLIFY(BC_T%Flag)
!!$    NULLIFY(BC_T%Value1)
!!$    NULLIFY(BC_T%Value2)
!!$    NULLIFY(BC_T%Value3)
!!$    NULLIFY(BC_T%Value_Array)
!!$    NULLIFY(BC_T%UseFunc)

    ! Velocity
    call InitializeBC_Struct(BC_V)
    BC_V%ID    = 'V'
    BC_V%Exist = BC_V_EXIST
    BC_V%IFlag = BC_V_IFLAG
    BC_V%Shift = BC_V_SHIFT

!!$    NULLIFY(BC_V%Flag)
!!$    NULLIFY(BC_V%Value1)
!!$    NULLIFY(BC_V%Value2)
!!$    NULLIFY(BC_V%Value3)
!!$    NULLIFY(BC_V%UseFunc)

    ! allocate arrays
    ! Concentration
    if (cflag) then
       allocate(BC_C%BCID(nfc,ncells))
       BC_C%BCID = -1

       allocate(BC_C%BCOperator(nfc,ncells))
       BC_C%BCOperator = 0

       allocate(BC_C%Flag(ncells))
       BC_C%Flag = 0
    end if

    ! Fluid Flow - Pressure and Velocity
    if (fflag) then
       allocate(BC_P%BCID(nfc,ncells))
       BC_P%BCID = -1

       allocate(BC_P%BCOperator(nfc,ncells))
       BC_P%BCOperator = 0

       allocate(BC_P%Flag(ncells))
       BC_P%Flag = 0
       allocate(BC_V%Flag(ncells))
       BC_V%Flag = 0
    end if

    ! Temperature
    if (tflag) then
       allocate(BC_T%BCID(nfc,ncells))
       BC_T%BCID = -1

       allocate(BC_T%BCOperator(nfc,ncells))
       BC_T%BCOperator = 0

       allocate(BC_T%Flag(ncells))
       BC_T%Flag = 0
    end if

  End Subroutine InitializeBoundaryConditions

  !-----------------------------------------------------------------------------

  Subroutine FinalizeBoundaryConditions ()
    !---------------------------------------------------------------------------
    ! Purpose:
    !   finalize boundary conditions
    !---------------------------------------------------------------------------
    Use bc_type_module, Only: BC_C, BC_P, BC_T, BC_V

    Implicit None

    !---------------------------------------------------------------------------

    ! Concentration
    call DeallocateBC_Struct(BC_C)
!!$    if (ASSOCIATED(BC_C%Flag))   DEALLOCATE(BC_C%Flag)
!!$    if (ASSOCIATED(BC_C%Value1)) DEALLOCATE(BC_C%Value1)
!!$    if (ASSOCIATED(BC_C%Value2)) DEALLOCATE(BC_C%Value2)
!!$    if (ASSOCIATED(BC_C%Value3)) DEALLOCATE(BC_C%Value3)
!!$    if (ASSOCIATED(BC_C%UseFunc)) DEALLOCATE(BC_C%UseFunc)

    ! Pressure
    call DeallocateBC_Struct(BC_P)
!!$    if (ASSOCIATED(BC_P%Flag))   DEALLOCATE(BC_P%Flag)
!!$    if (ASSOCIATED(BC_P%Value1)) DEALLOCATE(BC_P%Value1)
!!$    if (ASSOCIATED(BC_P%Value2)) DEALLOCATE(BC_P%Value2)
!!$    if (ASSOCIATED(BC_P%Value3)) DEALLOCATE(BC_P%Value3)
!!$    if (ASSOCIATED(BC_P%UseFunc)) DEALLOCATE(BC_P%UseFunc)

    ! Temperature
    call DeallocateBC_Struct(BC_T)
!!$    if (ASSOCIATED(BC_T%Flag))   DEALLOCATE(BC_T%Flag)
!!$    if (ASSOCIATED(BC_T%Value1)) DEALLOCATE(BC_T%Value1)
!!$    if (ASSOCIATED(BC_T%Value2)) DEALLOCATE(BC_T%Value2)
!!$    if (ASSOCIATED(BC_T%Value3)) DEALLOCATE(BC_T%Value3)
!!$    if (ASSOCIATED(BC_T%Value_Array)) DEALLOCATE(BC_T%Value_Array)
!!$    if (ASSOCIATED(BC_T%UseFunc)) DEALLOCATE(BC_T%UseFunc)

    ! Velocity
    call DeallocateBC_Struct(BC_V)
!!$    if (ASSOCIATED(BC_V%Flag))   DEALLOCATE(BC_V%Flag)
!!$    if (ASSOCIATED(BC_V%Value1)) DEALLOCATE(BC_V%Value1)
!!$    if (ASSOCIATED(BC_V%Value2)) DEALLOCATE(BC_V%Value2)
!!$    if (ASSOCIATED(BC_V%Value3)) DEALLOCATE(BC_V%Value3)
!!$    if (ASSOCIATED(BC_V%UseFunc)) DEALLOCATE(BC_V%UseFunc)

    return

  End Subroutine FinalizeBoundaryConditions

  !-----------------------------------------------------------------------------

  Subroutine SetBoundaryCondition (BCID, BC, mask, face, bctype, UseFunc, value1, value2, value3, value_array)
    !---------------------------------------------------------------------------
    ! Purpose:
    !   set boundary conditions
    !---------------------------------------------------------------------------
    use bc_data_module,       only: BC_C_DIRICHLET, BC_C_HNEUMANN, &
                                    BC_C_NEUMANN,                  &

                                    BC_P_DIRICHLET, BC_P_HNEUMANN, &
                                    BC_P_NEUMANN,                  &

                                    BC_T_DIRICHLET, BC_T_HNEUMANN, &
                                    BC_T_HTC,  &
                                    BC_T_NEUMANN, BC_T_RADIATION,  &
                                    BC_T_VFRADIATION, BC_T_HTC_GAP, &
                                    BC_V_DIRICHLET, BC_V_HNEUMANN, &
                                    BC_V_NEUMANN
    use bc_enum_types,        only: BC_DIRICHLET_OP,    &
                                    BC_NEUMANN_OP,      &
                                    BC_HNEUMANN_OP,     &
                                    BC_RADIATION_OP,    &
                                    BC_VFRADIATION_OP,  &
                                    BC_NO_OP,           &
                                    BC_HTC_EXTERNAL_OP, &
                                    BC_HTC_INTERNAL_OP, &
                                    BC_HTC_GAP_OP

    use bc_type_module,       only: BC_STRUCTURE
    use legacy_mesh_api,      only: ncells, nfc, Mesh
    use kinds, only: r8
    use pgslib_module,        only: PGSLIB_GLOBAL_COUNT

    ! Argument List
    integer,                                     intent(in)           :: BCID
    type   (BC_STRUCTURE),                       intent(INOUT)        :: BC
    logical, dimension(ncells), intent(IN)           :: mask
    integer,                    intent(IN)           :: face
    integer,                    intent(IN)           :: bctype
    logical,                    intent(IN), optional :: UseFunc
    real(r8),                   intent(IN), optional :: value1
    real(r8),                   intent(IN), optional :: value2
    real(r8),                   intent(IN), optional :: value3
    real(r8), dimension(:),     intent(IN), optional :: value_array

    ! Local Variables
    integer                    :: bit_mask
    integer                    :: bits
    integer                    :: bitsi
    integer                    :: k, length_value_array
    integer, dimension(ncells) :: bit_mask_array
    integer, dimension(ncells) :: bits_array
    integer, dimension(ncells) :: bitsi_array
    integer                                     :: tmpOp
    integer                                     :: n_interior_faces 
    
    !---------------------------------------------------------------------------

    ! select case and type
    tmpOp = BC_No_Op
    select case (BC%ID)
       ! Concentration
       case ('C')
         select case (bctype)
            case (BC_C_DIRICHLET)
              bits = BC_C_DIRICHLET
              tmpOp = BC_Dirichlet_Op
            case (BC_C_HNEUMANN)
              bits = BC_C_HNEUMANN
              tmpOp = BC_HNEUMANN_OP
            case (BC_C_NEUMANN)
              bits = BC_C_NEUMANN
              tmpOp = BC_NEUMANN_OP
         end select

       ! Pressure
       case ('P')
         select case (bctype)
            case (BC_P_DIRICHLET)
              bits = BC_P_DIRICHLET
              tmpOp = BC_Dirichlet_Op
            case (BC_P_HNEUMANN)
              bits = BC_P_HNEUMANN
              tmpOp = BC_HNEUMANN_OP
            case (BC_P_NEUMANN)
              bits = BC_P_NEUMANN
              tmpOp = BC_NEUMANN_OP
         end select

      ! Temperature
      case ('T')
         select case (bctype)
            case (BC_T_DIRICHLET)
              bits = BC_T_DIRICHLET
              tmpOp = BC_Dirichlet_Op
            case (BC_T_HNEUMANN)
              bits = BC_T_HNEUMANN
              tmpOp = BC_HNEUMANN_OP
            case (BC_T_NEUMANN)
              bits = BC_T_NEUMANN
              tmpOp = BC_NEUMANN_OP
            case (BC_T_HTC)
              bits = BC_T_HTC
              n_interior_faces = PGSLib_Global_COUNT(Mask .and. Mesh%Ngbr_cell(face) /= 0)
              if (n_interior_faces == 0) then
                 tmpOp = BC_HTC_EXTERNAL_OP 
              else
                 tmpOp = BC_HTC_INTERNAL_OP 
              endif
            case (BC_T_RADIATION)
              bits = BC_T_RADIATION
              tmpOp = BC_RADIATION_OP
            case (BC_T_VFRADIATION)
              bits = BC_T_VFRADIATION
              tmpOp = BC_VFRADIATION_OP
            case (BC_T_HTC_GAP)
              bits = BC_T_HTC_GAP
              tmpOp = BC_HTC_GAP_OP
         end select

       ! Velocity
       case ('V')
         select case (bctype)
            case (BC_V_DIRICHLET)
              bits = BC_V_DIRICHLET
              tmpOp = BC_Dirichlet_Op
            case (BC_V_HNEUMANN)
              bits = BC_V_HNEUMANN
              tmpOp = BC_HNEUMANN_OP
            case (BC_V_NEUMANN)
              bits = BC_V_NEUMANN
              tmpOp = BC_NEUMANN_OP
         end select
    end select


    where(mask)
       BC%BCOperator(face,:) = tmpOp
    end where

!    if (internal) then
!       bits = IOR(bits,BC%IFlag)
!    end if

    ! shift flag bits
    bit_mask = NOT(ISHFT(IOR(BC%EXIST,BC%IFLAG),(face - 1) * BC%Shift))
    bits  = ISHFT(bits,    (face - 1) * BC%Shift)
    bitsi = ISHFT(BC%IFlag,(face - 1) * BC%Shift)
    bit_mask_array = bit_mask
    bits_array     = bits
    bitsi_array    = 0
 
    ! search for internal boundaries
    Where (Mesh%ngbr_cell(face) /= 0)
       bitsi_array = bitsi
    End Where

    ! Tag the masked faces with the ID of the original BC condition
    where (mask)
       BC%BCID(face,:) = BCID
    end where

    where (mask)
       BC%Flag = Iand(BC%Flag,bit_mask_array)              ! clear what was there
       BC%Flag = IOR (BC%Flag,IOR(bits_array,bitsi_array)) ! or in new bits
    end where

    ! allocate and set boundary conditions

    ! Use function to determine value for boundary condition variable
    if (PRESENT(UseFunc)) then
       if (.not. ASSOCIATED(BC%UseFunc)) then
          allocate(BC%UseFunc(nfc,ncells))
       end if
       where (mask)
          BC%UseFunc(face,:) = UseFunc
       end where
    end if

    ! Value Number 1
    if (PRESENT(value1)) then
       if (.not. ASSOCIATED(BC%Value1)) then
          allocate(BC%Value1(nfc,ncells))
       end if
       where (mask)
          BC%Value1(face,:) = value1
       end where
    end if

    ! Value Number 2
    if (PRESENT(value2)) then
       if (.not. ASSOCIATED(BC%Value2)) then
          allocate(BC%Value2(nfc,ncells))
       end if
       Where (mask)
          BC%Value2(face,:) = value2
       End Where
    end if

    ! Value Number 3
    if (PRESENT(value3)) then
       if (.not. ASSOCIATED(BC%Value3)) then
          allocate(BC%Value3(nfc,ncells))
       end if
       Where (mask)
          BC%Value3(face,:) = value3
       End Where
    end if

    ! Value Array
    if (PRESENT(value_array)) then
       length_value_array = SIZE (value_array)
       if (.not. ASSOCIATED(BC%Value_Array)) then
          allocate(BC%Value_Array(length_value_array,nfc,ncells))
          BC%Value_Array = 0.0_r8
       end if
       do k = 1,length_value_array
           Where (mask)
              BC%Value_Array(k,face,:) = value_array(k)
           End Where
       end do
    end if

  End Subroutine SetBoundaryCondition

  !-----------------------------------------------------------------------------

End Module BC_SET_MODULE
