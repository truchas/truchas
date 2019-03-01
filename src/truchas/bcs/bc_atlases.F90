!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module BC_ATLASES
  !-----------------------------------------------------------------------------
  ! Purpose:
  !   Provide the data structures to support boundary condition atlas
  !   operations in Telluride.
  !
  ! Provides:
  !
  ! Documentation mostly in the documentation directory.
  !
  ! Author: Robert Ferrell (ferrell@cpca.com)
  !-----------------------------------------------------------------------------

  use bc_atlases_util,       only: CANONICALIZE,  &
                                   COLLATE,       &
                                   ORDER,         &
                                   PERMUTE,       &
                                   RENUMBER_CELLS

  use bc_atlases_data_types, only: BC_Atlas_Data, &
                                   BC_Atlas_Spec, &
                                   BC_Atlas,      &
                                   BC_Get_Face,   &
                                   BC_Get_Cell,   &
                                   BC_Get_Offset, &
                                   BC_Get_Length, &
                                   BC_Get_Values, &
                                   BC_Get_ValueIndex,&
                                   BC_Get_Usefunction,&
                                   BC_Get_Positions, &
                                   BC_Get_Data,   &
                                   BC_Set_Data_Size, &
                                   BC_NumberOfCharts,&
                                   BC_Get_Name,   &
                                   BC_Set_Name,   &
                                   BC_Get_DOF,    &
                                   BC_Set_DOF,    &
                                   Get_Scope,     &
                                   Set_Scope,     &
                                   Data_Size,     &
                                   SIZE,          &
                                   PERMUTE, CLONE, REDISTRIBUTE, &
                                   Dimensionality,   &
                                   INITIALIZE, FREE, ALLOC, REALLOC, APPEND,&
                                   ASSIGNMENT(=)
                                   
  

  PUBLIC :: CLONE
  PUBLIC :: Dimensionality
  PUBLIC :: INITIALIZE, FREE, ALLOC, REALLOC, APPEND
  PUBLIC :: ASSIGNMENT(=)
  PUBLIC :: GET_SCOPE, SET_SCOPE

END Module BC_ATLASES

