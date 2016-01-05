!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module parallel_communication
  use kinds
  public :: is_IOP, global_dot_product
  logical :: is_IOP = .true.
contains
  function global_dot_product (a, b) result (dp)
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: dp
    dp = dot_product(a,b)
  end function
end module
