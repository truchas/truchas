!*******************************************************************************
! Copyright Notice
!  + 2010-2012 North Carolina State University
!  + 2010-2012 Pacific Northwest National Laboratory
! 
! This file is part of SCORPIO.
! 
! SCORPIO is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! any later version.
! 
! SCORPIO is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with SCORPIO.  If not, see <http://www.gnu.org/licenses/>.
! 
!*******************************************************************************
!     Fortran Parallel IO interface - Constants declaration

  INTEGER SCORPIO_UNIFORM_CONTIGUOUS_READ, SCORPIO_NONUNIFORM_CONTIGUOUS_READ
  INTEGER SCORPIO_NONCONTIGUOUS_READ, SCORPIO_EVERYONE_ENTIRE_DATASET_READ
  INTEGER SCORPIO_UNIFORM_CONTIGUOUS_WRITE, SCORPIO_NONUNIFORM_CONTIGUOUS_WRITE
  INTEGER SCORPIO_NONCONTIGUOUS_WRITE

  PARAMETER (SCORPIO_UNIFORM_CONTIGUOUS_READ=0, SCORPIO_NONUNIFORM_CONTIGUOUS_READ=1)
  PARAMETER (SCORPIO_NONCONTIGUOUS_READ=2, SCORPIO_EVERYONE_ENTIRE_DATASET_READ=3)
  PARAMETER (SCORPIO_UNIFORM_CONTIGUOUS_WRITE=4, SCORPIO_NONUNIFORM_CONTIGUOUS_WRITE=5)
  PARAMETER (SCORPIO_NONCONTIGUOUS_WRITE=6)

  INTEGER SCORPIO_FILE_CREATE, SCORPIO_FILE_READONLY, SCORPIO_FILE_READWRITE
  PARAMETER (SCORPIO_FILE_CREATE=0, SCORPIO_FILE_READONLY=1, SCORPIO_FILE_READWRITE=2)

  INTEGER SCORPIO_INTEGER, SCORPIO_DOUBLE, SCORPIO_FLOAT, SCORPIO_LONG, SCORPIO_CHAR, SCORPIO_BYTE
  PARAMETER (SCORPIO_INTEGER=0, SCORPIO_DOUBLE=1, SCORPIO_FLOAT=2, SCORPIO_LONG=3, SCORPIO_CHAR=4, SCORPIO_BYTE=5)
