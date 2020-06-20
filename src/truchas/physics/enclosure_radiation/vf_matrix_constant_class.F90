!!
!! VF_MATRIX_CONSTANT_CLASS
!!
!! A common interface for time-independent view factor implementations that
!! extends the VF_MATRIX class. This extension specifies an initialization
!! procedure which is expected to acquire the view factors from a radiation
!! enclosure file.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! June 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vf_matrix_constant_class

  use vf_matrix_class
  use rad_encl_file_type
  implicit none
  private

  type, abstract, extends(vf_matrix), public :: vf_matrix_constant
  contains
    procedure(init), deferred :: init
  end type

  abstract interface
    subroutine init(this, file)
      import :: vf_matrix_constant, rad_encl_file
      class(vf_matrix_constant), intent(out) :: this
      type(rad_encl_file), intent(in) :: file
    end subroutine
  end interface

end module vf_matrix_constant_class
