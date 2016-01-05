!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module BC_ENUM_TYPES
  !-----------------------------------------------------------------------------
  ! Purpose:
  !   Provide the constant parameters to support boundary condition
  !   operations
  !
  ! Provides:
  !
  ! Documentation mostly in the documentation directory.
  !
  ! Author: Robert Ferrell (ferrell@cpca.com)
  !-----------------------------------------------------------------------------
  Implicit None
  PRIVATE
  
  PUBLIC :: BC_Position_Data, BC_Values
  PUBLIC :: BC_REFLECTIVE_OP, BC_DIRICHLET_OP, BC_NEUMANN_OP, BC_HNEUMANN_OP, BC_NO_OP
  PUBLIC :: BC_RADIATION_OP, BC_VFRADIATION_OP, BC_HTC_EXTERNAL_OP, BC_HTC_INTERNAL_OP, &
            BC_HTC_GAP_OP, BC_EXTERIOR_Op, BC_HTC_RADIATION_EXT_OP, BC_HTC_RADIATION_INT_OP
  PUBLIC :: BC_X_TRACTION_OP, BC_Y_TRACTION_OP, BC_Z_TRACTION_OP
  PUBLIC :: BC_X_DISPLACEMENT_OP, BC_Y_DISPLACEMENT_OP, BC_Z_DISPLACEMENT_OP
  PUBLIC :: BC_NORM_DISPLACEMENT_OP, BC_NORM_TRACTION_OP, BC_FREE_INTERFACE_OP
  PUBLIC :: BC_NORMAL_CONSTRAINT_OP, BC_CONTACT_OP
  PUBLIC :: BC_MAX_OPERATORS
  PUBLIC :: BC_INVALID_FACE, BC_INVALID_CELL, BC_INVALID_LENGTH, BC_INVALID_OFFSET
  PUBLIC :: BC_INVALID_DIMENSIONALITY, BC_INVALID_SIZE, BC_INVALID_ATLASINDEX
  PUBLIC :: BC_INVALID_COUNT, BC_INVALID_ID
  PUBLIC :: BC_PRESSURE_ID, BC_TEMPERATURE_ID, BC_DISPLACEMENT_ID
  PUBLIC :: BC_Cells_Scope_Local, BC_Cells_Scope_Global, BC_Cells_No_Scope
  PUBLIC :: BC_OP_INVALID_STATE, BC_OP_ACTIVE, BC_OP_DEACTIVE
  PUBLIC :: BC_STRING_LEN, BC_NO_NAME

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Type definitions
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Data items provided by this module
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! For now, a fixed number of arrays, with "enums" to identify
  ! contents of each array.

  integer, parameter :: BC_Position_Data        = 1
  integer, parameter :: BC_Values               = 2
  integer, parameter :: BC_INVALID_ID           = -1
  integer, parameter :: BC_INVALID_COUNT        = -1
  integer, parameter :: BC_STRING_LEN           = 1024
  character(BC_STRING_LEN), parameter:: BC_NO_NAME = 'NO NAME'

  ! Supported BC Operators

  ! For now, fixed maximum number of opererators in each BC_Specifier
  integer, parameter :: BC_MAX_OPERATORS = 30

  integer, parameter :: BC_DIRICHLET_OP         = 1
  integer, parameter :: BC_HNEUMANN_OP          = 2
  integer, parameter :: BC_NEUMANN_OP           = 3
  integer, parameter :: BC_REFLECTIVE_OP        = 4
  integer, parameter :: BC_RADIATION_OP         = 5
  integer, parameter :: BC_HTC_EXTERNAL_OP      = 6
  integer, parameter :: BC_HTC_INTERNAL_OP      = 7
  integer, parameter :: BC_EXTERIOR_OP          = 8
  integer, parameter :: BC_NO_OP                = -1
  ! Solid Mechanics BCs
  integer, parameter :: BC_X_TRACTION_OP        = 9
  integer, parameter :: BC_Y_TRACTION_OP        = 10
  integer, parameter :: BC_Z_TRACTION_OP        = 11
  integer, parameter :: BC_X_DISPLACEMENT_OP    = 12
  integer, parameter :: BC_Y_DISPLACEMENT_OP    = 13
  integer, parameter :: BC_Z_DISPLACEMENT_OP    = 14

  integer, parameter :: BC_NORM_DISPLACEMENT_OP = 15
  integer, parameter :: BC_NORM_TRACTION_OP     = 16
  integer, parameter :: BC_FREE_INTERFACE_OP    = 17
  integer, parameter :: BC_NORMAL_CONSTRAINT_OP = 18
  integer, parameter :: BC_CONTACT_OP           = 19

  ! Enclosure radiation mode
  integer, parameter :: BC_VFRADIATION_OP       = 20

  ! HTC at at gap interface
  integer, parameter :: BC_HTC_GAP_OP           = 21
  ! Radiation plus Internal HTC
  integer, parameter :: BC_HTC_RADIATION_EXT_OP = 22
  integer, parameter :: BC_HTC_RADIATION_INT_OP = 23


  ! Parameters for identifying variables
  integer, parameter :: BC_PRESSURE_ID          = 1
  integer, parameter :: BC_TEMPERATURE_ID       = 2
  integer, parameter :: BC_DISPLACEMENT_ID      = 3

  ! Parameters for identifying scope
  integer, parameter :: BC_Cells_No_Scope       = -1
  integer, parameter :: BC_Cells_Scope_Global   = 1
  integer, parameter :: BC_Cells_Scope_Local    = 2

  ! Parameters for defining state of an operator
  integer, parameter :: BC_OP_INVALID_STATE = -1
  integer, parameter :: BC_OP_ACTIVE        =  1
  integer, parameter :: BC_OP_DEACTIVE      =  2
  

  ! Defaults used in various routines

  integer, parameter :: BC_INVALID_FACE   = -5
  integer, parameter :: BC_INVALID_CELL   = -6
  integer, parameter :: BC_INVALID_OFFSET = -1
  integer, parameter :: BC_INVALID_LENGTH = -7
  integer, parameter :: BC_INVALID_SIZE   = -8
  integer, parameter :: BC_INVALID_DIMENSIONALITY = -9
  integer, parameter :: BC_INVALID_ATLASINDEX = -10

END Module BC_ENUM_TYPES

