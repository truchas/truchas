!!
!! DATA_MAPPER_CLASS
!!
!! This defines the interface used by Truchas to map cell-centered fields
!! between the main heat transfer mesh and the tetrahedral EM mesh.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module data_mapper_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: data_mapper
  contains
    procedure(init), deferred :: init
    procedure(map_field), deferred :: map_field
  end type

  abstract interface
    subroutine init(this, mesh1, mesh2)
      use base_mesh_class
      import data_mapper
      class(data_mapper), intent(out) :: this
      class(base_mesh), intent(in) :: mesh1, mesh2
    end subroutine
    subroutine map_field(this, src, dest, defval, map_type, pullback)
      import data_mapper, r8
      class(data_mapper), intent(in) :: this
      real(r8), intent(in) :: src(:)
      real(r8), intent(out) :: dest(:)
      real(r8), intent(in) :: defval
      integer, intent(in) :: map_type
      logical, intent(in), optional :: pullback
    end subroutine
  end interface

  !! Mapping type options
  integer, parameter, public :: LOCALLY_CONSERVATIVE  = 1
  integer, parameter, public :: LOCALLY_BOUNDED       = 2
  integer, parameter, public :: GLOBALLY_CONSERVATIVE = 3

end module data_mapper_class
