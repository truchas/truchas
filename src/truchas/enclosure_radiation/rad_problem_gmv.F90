!!
!! RAD_PROBLEM_GMV
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
!!  This extends RAD_ENCL_GMV and RAD_SOLVER_GMV; refer to those modules for
!!  usage info, especially the former.  This module extends the following
!!  generic subroutines for RAD_PROBLEM type objects.
!!
!!  CALL GMV_WRITE_ENCLOSURE (THIS) writes the distributed radiation enclosure
!!    surface mesh contained within the RAD_PROBLEM type object THIS to the
!!    graphics file.  When executing in parallel, the partitioning of faces
!!    is written as flag data.
!!
!!  CALL GMV_WRITE_VARIABLE (THIS, VAR, NAME) writes the distributed variable
!!    data VAR to the graphics file.  THIS is a RAD_PROBLEM type object.
!!    VAR a face-based field defined on the associated heat conduction mesh
!!    and its value on the radiation enclosure surface mesh contained within
!!    THIS is extracted.  NAME is an arbitrary string used to label the
!!    variable; only the first 8 characters are significant.  This may be
!!    called multiple times.
!!

#include "f90_assert.fpp"

module rad_problem_gmv

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use rad_encl_gmv
  use rad_problem_type
  use parallel_permutations
  implicit none
  private

  !! Re-export everything from RAD_SOLVER_GMV with added specific subroutines
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
    type(rad_problem), intent(in) :: this
    call gmv_write_enclosure (this%encl)
  end subroutine gmv_write_encl

  subroutine gmv_write_var (this, var, name)
    type(rad_problem), intent(in) :: this
    real(r8), intent(in) :: var(:)
    character(*), intent(in) :: name
    real(r8) :: var_er(this%nface_er)
    ASSERT(size(var) == this%nface_hc)
    call reorder (this%perm_er_to_hc, var_er, var)
    call gmv_write_variable (this%encl, var_er, name)
  end subroutine gmv_write_var

end module rad_problem_gmv
