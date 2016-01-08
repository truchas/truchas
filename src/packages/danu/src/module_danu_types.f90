!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  DANU
!     Fortran module that defines the types
!     needed by Danu
!
module danu_types

    use iso_c_binding

    implicit none

    public


    integer, parameter :: Danu_INT32  = SELECTED_INT_KIND(6)
    integer, parameter :: Danu_INT64  = SELECTED_INT_KIND(16)
    integer, parameter :: Danu_REAL32 = SELECTED_REAL_KIND(6)
    integer, parameter :: Danu_REAL64 = SELECTED_INT_KIND(12)

    integer, parameter :: Danu_INT    = C_INT
    integer, parameter :: Danu_LONG   = C_LONG
    integer, parameter :: Danu_FLOAT  = C_FLOAT 
    integer, parameter :: Danu_DOUBLE = C_DOUBLE 
    integer, parameter :: Danu_CHAR   = C_CHAR
    integer, parameter :: Danu_SIZE_T = C_SIZE_T

    ! HDF5 types
    ! In the C code these are pointers to structs that hold
    ! the values. This eliminates the need to determine the correct 
    ! integer kind type to hold 
    type(C_PTR)        :: Danu_HID_PTR_T
    type(C_PTR)        :: Danu_HSIZE_T

end module danu_types    

    
