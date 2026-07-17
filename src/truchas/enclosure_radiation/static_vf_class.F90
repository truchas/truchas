!!
!! STATIC_VF_CLASS
!!
!! A common interface for time-independent view factor implementations that
!! extends the ENCL_VF class. This extension specifies an initialization
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

module static_vf_class

  use encl_vf_class
  use rad_encl_file_type
  implicit none
  private

  type, abstract, extends(encl_vf), public :: static_vf
  contains
    procedure(init), deferred :: init
  end type

  abstract interface
    subroutine init(this, file)
      import :: static_vf, rad_encl_file
      class(static_vf), intent(out) :: this
      type(rad_encl_file), intent(in) :: file
    end subroutine
  end interface

end module static_vf_class
