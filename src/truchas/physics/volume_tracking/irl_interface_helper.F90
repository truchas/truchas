!!
!! Provides helper functions to faciliate Truchas interace with the
!! Interface Reconstruction Library (IRL)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module irl_interface_helper

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use irl_fortran_interface
  implicit none
  
  
  !procedure truchas_tet_to_irl
  !procedure truchas_pyramid_to_irl
  !procedure truchas_wedge_to_irl
  !procedure truchas_hex_to_irl
  
  integer, parameter, private :: &
    truchas_irl_tet_mapping(4) = [1,2,3,4]

  integer, parameter, private :: &
    truchas_irl_hex_mapping(8) = [6,7,8,5,2,3,4,1]
    
  integer, parameter, private :: &
    truchas_irl_dod_mapping(8) = [6,7,8,5,2,3,4,1]

contains

  subroutine truchas_tet_to_irl(a_truchas_tet, a_irl_tet)
  
    real(r8), intent(in) :: a_truchas_tet(:,:)
    type(Tet_type), intent(inout) :: a_irl_tet

    call construct(a_irl_tet, a_truchas_tet)

  end subroutine truchas_tet_to_irl
  
  subroutine truchas_hex_to_irl(a_truchas_hex, a_irl_hex)
  
    real(r8), intent(in) :: a_truchas_hex(:,:)
    type(Hex_type), intent(inout) :: a_irl_hex
       
    call construct(a_irl_hex, a_truchas_hex(:,truchas_irl_hex_mapping(:)))
    
  end subroutine truchas_hex_to_irl
  
  subroutine truchas_dod_to_irl(a_truchas_dod, a_irl_dod)
  
    real(r8), intent(in) :: a_truchas_dod(:,:)
    type(Dod_type), intent(inout) :: a_irl_dod
       
    call construct(a_irl_dod, a_truchas_dod(:,truchas_irl_dod_mapping(:)))
    
  end subroutine truchas_dod_to_irl

end module irl_interface_helper
