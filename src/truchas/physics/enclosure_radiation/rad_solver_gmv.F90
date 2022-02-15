!!
!! RAD_SOLVER_GMV
!!
!! This is a quick-n-dirty high level layer over GMV's gmvwrite C library
!! that provides some procedures for writing a radiation enclosure surface
!! mesh and face-based fields over that mesh to a GMV-format graphics file.
!! It is intended for ad hoc use in standalone test programs and debugging
!! situations.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Revised May 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! USAGE
!!
!!  This extends RAD_ENCL_GMV; refer to that module for usage info.  This
!!  module extends the following generic subroutines for RAD_SOLVER type
!!  objects.
!!
!!  CALL GMV_WRITE_ENCLOSURE (THIS) writes the distributed radiation enclosure
!!    surface mesh contained within the RAD_SOLVER type object THIS to the
!!    graphics file.  When executing in parallel, the partitioning of faces
!!    is written as flag data.
!!
!!  CALL GMV_WRITE_VARIABLE (THIS, VAR, NAME) writes the distributed variable
!!    data VAR to the graphics file.  THIS is a RAD_SOLVER type object, and
!!    VAR a face field defined on the enclosure mesh contained within THIS.
!!    NAME is an arbitrary string used to label the variable; only the first
!!    8 characters are significant.  This may be called multiple times.
!!

#include "f90_assert.fpp"

module rad_solver_gmv

  use kinds, only: r8
  use rad_encl_gmv
  use rad_solver_type
  implicit none
  private

  !! Re-export everything from RAD_ENCL_GMV with added specific subroutines
  !! for the generic GMV_WRITE_ENCLOSURE and GMV_WRITE_VARIABLE
  public :: gmv_open, gmv_close, gmv_write_enclosure
  public :: gmv_begin_variables, gmv_end_variables, gmv_write_variable

  interface gmv_write_enclosure
    procedure gmv_write_encl
  end interface

  interface gmv_write_variable
    procedure gmv_write_var
  end interface

contains

  subroutine gmv_write_encl (this)
    use rad_solver_type
    type(rad_solver), intent(in) :: this
    call gmv_write_enclosure (this%encl)
  end subroutine gmv_write_encl

  subroutine gmv_write_var (this, u, name)
    type(rad_solver), intent(in) :: this
    real(r8), intent(in) :: u(:)
    character(*), intent(in) :: name
    call gmv_write_variable (this%encl, u, name)
  end subroutine gmv_write_var

end module rad_solver_gmv
