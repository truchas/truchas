Program MAIN
  !-----------------------------------------------------------------------------
  ! TRUCHAS: A three-dimensional, high-resolution Eulerian algorithm
  ! designed to model incompressible flows of multiple immiscible
  ! fluids, coupled with heat transfer, and possibly phase change. The
  ! fluids are delineated with interfaces that can have surface
  ! tension.
  !
  ! Copyright 2007-2013. Los Alamos National Security, LLC.
  !                                                                  
  ! This material was produced under U.S. Government contract DE-AC52-06NA25396
  ! for Los Alamos National Laboratory (LANL), which is operated by Los Alamos
  ! National Security, LLC for the U.S. Department of Energy. The U.S. Government
  ! has rights to use, reproduce, and distribute this software.  NEITHER THE
  ! GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS
  ! OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software
  ! is modified to produce derivative works, such modified software should be
  ! clearly marked, so as not to confuse it with the version available from LANL.
  ! 
  ! Additionally, this program is free software; you can redistribute it and/or
  ! modify it under the terms of the GNU General Public License as published by
  ! the Free Software Foundation; either version 2 of the License, or (at your
  ! option) any later version. Accordingly, this program is distributed in the
  ! hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
  ! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  ! See the GNU General Public License for more details.
  !
  ! For information regarding Truchas, visit the Truchas home page at:
  !
  !         http://telluride.lanl.gov
  !
  ! or contact the development team at:
  !
  !         telluride-support@lanl.gov
  !-----------------------------------------------------------------------------

!  use,intrinsic :: ieee_exceptions
  Use drivers, Only: CODE
  Implicit None

  !-----------------------------------------------------------------------------

!   call ieee_set_halting_mode ((/ieee_invalid, ieee_overflow, ieee_divide_by_zero/), .true.)
   Call CODE ()

End Program MAIN
