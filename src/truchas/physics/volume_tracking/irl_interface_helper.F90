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
  

  interface truchas_poly_to_irl
    module procedure truchas_tet_to_irl
    module procedure truchas_pyramid_to_irl
    module procedure truchas_wedge_to_irl
    module procedure truchas_hex_to_irl
    module procedure truchas_sym_tet_to_irl
    module procedure truchas_sym_pyramid_to_irl
    module procedure truchas_sym_wedge_to_irl
    module procedure truchas_sym_hex_to_irl    
    module procedure truchas_octa_to_irl
    module procedure truchas_dod_to_irl
    module procedure truchas_capdod_LLLL_to_irl
    module procedure truchas_capdod_LLLT_to_irl
    module procedure truchas_capdod_LTLT_to_irl
    module procedure truchas_capdod_LLTT_to_irl
    module procedure truchas_capdod_LTTT_to_irl
    module procedure truchas_capdod_TTTT_to_irl
    module procedure truchas_capocta_LLL_to_irl
    module procedure truchas_capocta_LLT_to_irl
    module procedure truchas_capocta_LTT_to_irl
    module procedure truchas_capocta_TTT_to_irl    
  end interface truchas_poly_to_irl

  

  
  integer, parameter, private :: &
       truchas_irl_tet_mapping(4) = [2,3,4,1]
    
  integer, parameter, private :: &
       truchas_irl_pyramid_mapping(5) = [4,3,2,1,5] 
    
  integer, parameter, private :: &
       truchas_irl_wedge_mapping(6) = [5,6,4,2,3,1]

  integer, parameter, private :: &
       truchas_irl_hex_mapping(8) = [6,7,8,5,2,3,4,1]

  integer, parameter, private :: &
       truchas_irl_sym_tet_mapping(8) = [2,3,4,1,6,8,7,5]
    
  integer, parameter, private :: &
       truchas_irl_sym_pyramid_mapping(10) = [4,3,2,1,5,10,8,7,6,9] 
    
  integer, parameter, private :: &
       truchas_irl_sym_wedge_mapping(11) = [5,6,4,2,3,1,11,8,9,7,10]

  integer, parameter, private :: &
       truchas_irl_sym_hex_mapping(14) = [6,7,8,5,2,3,4,1,14,10,11,12,9,13]
  
  integer, parameter, private :: &
       truchas_irl_octa_mapping(6) = [5,6,4,2,3,1]
    
  integer, parameter, private :: &
       truchas_irl_dod_mapping(8) = [6,7,8,5,2,3,4,1]

  integer, parameter, private :: &
       truchas_irl_capdod_mapping(9) = [6,7,8,5,2,3,4,1,9]

  integer, parameter, public :: &
       truchas_irl_capocta_mapping(7) = [5,6,4,2,3,1,7]  

contains

  subroutine truchas_tet_to_irl(a_truchas_tet, a_irl_tet)
  
    real(r8), intent(in) :: a_truchas_tet(:,:)
    type(Tet_type), intent(inout) :: a_irl_tet

    call construct(a_irl_tet, a_truchas_tet(:,truchas_irl_tet_mapping))

  end subroutine truchas_tet_to_irl
  
  subroutine truchas_pyramid_to_irl(a_truchas_pyramid, a_irl_pyramid)
  
    real(r8), intent(in) :: a_truchas_pyramid(:,:)
    type(Pyrmd_type), intent(inout) :: a_irl_pyramid

    call construct(a_irl_pyramid, a_truchas_pyramid(:,truchas_irl_pyramid_mapping))

  end subroutine truchas_pyramid_to_irl
  
  subroutine truchas_wedge_to_irl(a_truchas_wedge, a_irl_wedge)
  
    real(r8), intent(in) :: a_truchas_wedge(:,:)
    type(TriPrism_type), intent(inout) :: a_irl_wedge
       
    call construct(a_irl_wedge, a_truchas_wedge(:,truchas_irl_wedge_mapping))
    
  end subroutine truchas_wedge_to_irl

  subroutine truchas_hex_to_irl(a_truchas_hex, a_irl_hex)
  
    real(r8), intent(in) :: a_truchas_hex(:,:)
    type(Hex_type), intent(inout) :: a_irl_hex
       
    call construct(a_irl_hex, a_truchas_hex(:,truchas_irl_hex_mapping))
    
  end subroutine truchas_hex_to_irl  

  subroutine truchas_sym_tet_to_irl(a_truchas_tet, a_irl_tet)
  
    real(r8), intent(in) :: a_truchas_tet(:,:)
    type(SymTet_type), intent(inout) :: a_irl_tet

    call construct(a_irl_tet, a_truchas_tet(:,truchas_irl_sym_tet_mapping))

  end subroutine truchas_sym_tet_to_irl
  
  subroutine truchas_sym_pyramid_to_irl(a_truchas_pyramid, a_irl_pyramid)
  
    real(r8), intent(in) :: a_truchas_pyramid(:,:)
    type(SymPyrmd_type), intent(inout) :: a_irl_pyramid

    call construct(a_irl_pyramid, a_truchas_pyramid(:,truchas_irl_sym_pyramid_mapping))

  end subroutine truchas_sym_pyramid_to_irl
  
  subroutine truchas_sym_wedge_to_irl(a_truchas_wedge, a_irl_wedge)
  
    real(r8), intent(in) :: a_truchas_wedge(:,:)
    type(SymTriPrism_type), intent(inout) :: a_irl_wedge
       
    call construct(a_irl_wedge, a_truchas_wedge(:,truchas_irl_sym_wedge_mapping))
    
  end subroutine truchas_sym_wedge_to_irl

  subroutine truchas_sym_hex_to_irl(a_truchas_hex, a_irl_hex)
  
    real(r8), intent(in) :: a_truchas_hex(:,:)
    type(SymHex_type), intent(inout) :: a_irl_hex
       
    call construct(a_irl_hex, a_truchas_hex(:,truchas_irl_sym_hex_mapping))
    
  end subroutine truchas_sym_hex_to_irl
  
  subroutine truchas_octa_to_irl(a_truchas_octa, a_irl_octa)
  
    real(r8), intent(in) :: a_truchas_octa(:,:)
    type(Octa_type), intent(inout) :: a_irl_octa
       
    call construct(a_irl_octa, a_truchas_octa(:,truchas_irl_octa_mapping))
    
  end subroutine truchas_octa_to_irl
  
  subroutine truchas_dod_to_irl(a_truchas_dod, a_irl_dod)
  
    real(r8), intent(in) :: a_truchas_dod(:,:)
    type(Dod_type), intent(inout) :: a_irl_dod
       
    call construct(a_irl_dod, a_truchas_dod(:,truchas_irl_dod_mapping))
    
  end subroutine truchas_dod_to_irl

  subroutine truchas_capdod_LLLL_to_irl(a_truchas_capdod, a_irl_capdod)
  
    real(r8), intent(in) :: a_truchas_capdod(:,:)
    type(CapDod_LLLL_type), intent(inout) :: a_irl_capdod
       
    call construct(a_irl_capdod, a_truchas_capdod(:,truchas_irl_capdod_mapping))
    
  end subroutine truchas_capdod_LLLL_to_irl

  subroutine truchas_capdod_LLLT_to_irl(a_truchas_capdod, a_irl_capdod)
  
    real(r8), intent(in) :: a_truchas_capdod(:,:)
    type(CapDod_LLLT_type), intent(inout) :: a_irl_capdod
       
    call construct(a_irl_capdod, a_truchas_capdod(:,truchas_irl_capdod_mapping))
    
  end subroutine truchas_capdod_LLLT_to_irl

  subroutine truchas_capdod_LTLT_to_irl(a_truchas_capdod, a_irl_capdod)
  
    real(r8), intent(in) :: a_truchas_capdod(:,:)
    type(CapDod_LTLT_type), intent(inout) :: a_irl_capdod
       
    call construct(a_irl_capdod, a_truchas_capdod(:,truchas_irl_capdod_mapping))
    
  end subroutine truchas_capdod_LTLT_to_irl


  subroutine truchas_capdod_LLTT_to_irl(a_truchas_capdod, a_irl_capdod)
  
    real(r8), intent(in) :: a_truchas_capdod(:,:)
    type(CapDod_LLTT_type), intent(inout) :: a_irl_capdod
       
    call construct(a_irl_capdod, a_truchas_capdod(:,truchas_irl_capdod_mapping))
    
  end subroutine truchas_capdod_LLTT_to_irl


  subroutine truchas_capdod_LTTT_to_irl(a_truchas_capdod, a_irl_capdod)

    real(r8), intent(in) :: a_truchas_capdod(:,:)
    type(CapDod_LTTT_type), intent(inout) :: a_irl_capdod
       
    call construct(a_irl_capdod, a_truchas_capdod(:,truchas_irl_capdod_mapping))
    
  end subroutine truchas_capdod_LTTT_to_irl


  subroutine truchas_capdod_TTTT_to_irl(a_truchas_capdod, a_irl_capdod)
  
    real(r8), intent(in) :: a_truchas_capdod(:,:)
    type(CapDod_TTTT_type), intent(inout) :: a_irl_capdod
       
    call construct(a_irl_capdod, a_truchas_capdod(:,truchas_irl_capdod_mapping))
    
  end subroutine truchas_capdod_TTTT_to_irl

  
  subroutine truchas_capocta_LLL_to_irl(a_truchas_capocta, a_irl_capocta)
  
    real(r8), intent(in) :: a_truchas_capocta(:,:)
    type(CapOcta_LLL_type), intent(inout) :: a_irl_capocta
       
    call construct(a_irl_capocta, a_truchas_capocta(:,truchas_irl_capocta_mapping))
    
  end subroutine truchas_capocta_LLL_to_irl

  
  subroutine truchas_capocta_LLT_to_irl(a_truchas_capocta, a_irl_capocta)
  
    real(r8), intent(in) :: a_truchas_capocta(:,:)
    type(CapOcta_LLT_type), intent(inout) :: a_irl_capocta
       
    call construct(a_irl_capocta, a_truchas_capocta(:,truchas_irl_capocta_mapping))
    
  end subroutine truchas_capocta_LLT_to_irl

  subroutine truchas_capocta_LTT_to_irl(a_truchas_capocta, a_irl_capocta)
  
    real(r8), intent(in) :: a_truchas_capocta(:,:)
    type(CapOcta_LTT_type), intent(inout) :: a_irl_capocta
       
    call construct(a_irl_capocta, a_truchas_capocta(:,truchas_irl_capocta_mapping))
    
  end subroutine truchas_capocta_LTT_to_irl

  
  subroutine truchas_capocta_TTT_to_irl(a_truchas_capocta, a_irl_capocta)
  
    real(r8), intent(in) :: a_truchas_capocta(:,:)
    type(CapOcta_TTT_type), intent(inout) :: a_irl_capocta
       
    call construct(a_irl_capocta, a_truchas_capocta(:,truchas_irl_capocta_mapping))
    
  end subroutine truchas_capocta_TTT_to_irl    
  

end module irl_interface_helper
