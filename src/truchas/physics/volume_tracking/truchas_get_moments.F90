!!
!! Traverse the mesh to obtain volume moments for a flux volume
!! in an unsplit VOF framework.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module truchas_getMoments_mod

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use irl_fortran_interface
  implicit none

  interface truchas_getMoments
     module procedure truchas_getMoments_CapDod_LLLL
     module procedure truchas_getMoments_CapDod_LLLT
     module procedure truchas_getMoments_CapDod_LTLT
     module procedure truchas_getMoments_CapDod_LLTT
     module procedure truchas_getMoments_CapDod_LTTT
     module procedure truchas_getMoments_CapDod_TTTT
     module procedure truchas_getMoments_CapOcta_LLL
     module procedure truchas_getMoments_CapOcta_LLT
     module procedure truchas_getMoments_CapOcta_LTT
     module procedure truchas_getMoments_CapOcta_TTT     
  end interface truchas_getMoments

contains

  subroutine truchas_getMoments_CapDod_LLLL(a_mesh, a_nmat, a_flux_volume, a_starting_index, a_reconstruction, a_flux)

    use traversal_tracker_type, only : traversal_tracker
    use unstr_mesh_type
    use tagged_volumes_type
    
    type(unstr_mesh), pointer, intent(in) :: a_mesh
    integer, intent(in) :: a_nmat
    type(CapDod_LLLL_type), intent(in) :: a_flux_volume
    integer, intent(in) :: a_starting_index(2)
    type(LocSepLink_type), intent(in) :: a_reconstruction(:)
    type(tagged_volumes), intent(inout) :: a_flux

    type(traversal_tracker) :: tracker
    integer :: current_cell
    type(SepVol_type) :: irl_cell_flux
    real(r8) :: cell_flux(a_nmat)
    integer :: n

    if(a_starting_index(2) == 0) then
      call tracker%init([a_starting_index(1)])
    else
      call tracker%init(a_starting_index)      
    end if    
    call new(irl_cell_flux)
    call a_flux%clear()
    do while(tracker%still_cells_to_visit())

       current_cell = tracker%get_next_cell()
       
       call setMinimumVolToTrack(a_mesh%volume(current_cell)*1.0e-15_r8)
       call getNormMoments(a_flux_volume, a_reconstruction(current_cell), irl_cell_flux)

       cell_flux(1) = getVolume(irl_cell_flux, 0)
       cell_flux(2) = getVolume(irl_cell_flux, 1)       
       
       if(sum(abs(cell_flux)) > 1.0e-15*a_mesh%volume(current_cell)) then
          call a_flux%add_flux(current_cell, cell_flux)         
          associate( cn => a_mesh%cnhbr(a_mesh%xcnhbr(current_cell):a_mesh%xcnhbr(current_cell+1)-1))
            do n = 1, size(cn)
               if(cn(n) /= 0 .and. tracker%cell_not_encountered(cn(n))) then
                  call tracker%add_cell(cn(n))
               end if
            end do
          end associate   
       end if               
    end do
    
  end subroutine truchas_getMoments_CapDod_LLLL

  subroutine truchas_getMoments_CapDod_LLLT(a_mesh, a_nmat, a_flux_volume, a_starting_index, a_reconstruction, a_flux)

    use traversal_tracker_type, only : traversal_tracker    
    use unstr_mesh_type
    use tagged_volumes_type
    
    type(unstr_mesh), pointer, intent(in) :: a_mesh
    integer, intent(in) :: a_nmat
    type(CapDod_LLLT_type), intent(in) :: a_flux_volume
    integer, intent(in) :: a_starting_index(2)
    type(LocSepLink_type), intent(in) :: a_reconstruction(:)
    type(tagged_volumes), intent(inout) :: a_flux

    type(traversal_tracker) :: tracker
    integer :: current_cell
    type(SepVol_type) :: irl_cell_flux
    real(r8) :: cell_flux(a_nmat)
    integer :: n

    if(a_starting_index(2) == 0) then
      call tracker%init([a_starting_index(1)])
    else
      call tracker%init(a_starting_index)      
    end if
    call new(irl_cell_flux)    
    call a_flux%clear()
    do while(tracker%still_cells_to_visit())

       current_cell = tracker%get_next_cell()
       
       call setMinimumVolToTrack(a_mesh%volume(current_cell)*1.0e-15_r8)
       call getNormMoments(a_flux_volume, a_reconstruction(current_cell), irl_cell_flux)

       cell_flux(1) = getVolume(irl_cell_flux, 0)
       cell_flux(2) = getVolume(irl_cell_flux, 1)       
   
       if(sum(abs(cell_flux)) > 1.0e-15*a_mesh%volume(current_cell)) then
          call a_flux%add_flux(current_cell, cell_flux)         
          associate( cn => a_mesh%cnhbr(a_mesh%xcnhbr(current_cell):a_mesh%xcnhbr(current_cell+1)-1))
            do n = 1, size(cn)
               if(cn(n) /= 0 .and. tracker%cell_not_encountered(cn(n))) then
                  call tracker%add_cell(cn(n))
               end if
            end do
          end associate   
       end if               
    end do
    
  end subroutine truchas_getMoments_CapDod_LLLT

  subroutine truchas_getMoments_CapDod_LTLT(a_mesh, a_nmat, a_flux_volume, a_starting_index, a_reconstruction, a_flux)

    use traversal_tracker_type, only : traversal_tracker    
    use unstr_mesh_type
    use tagged_volumes_type
    
    type(unstr_mesh), pointer, intent(in) :: a_mesh
    integer, intent(in) :: a_nmat
    type(CapDod_LTLT_type), intent(in) :: a_flux_volume
    integer, intent(in) :: a_starting_index(2)
    type(LocSepLink_type), intent(in) :: a_reconstruction(:)
    type(tagged_volumes), intent(inout) :: a_flux

    type(traversal_tracker) :: tracker
    integer :: current_cell
    type(SepVol_type) :: irl_cell_flux
    real(r8) :: cell_flux(a_nmat)
    integer :: n

    
    if(a_starting_index(2) == 0) then
      call tracker%init([a_starting_index(1)])
    else
      call tracker%init(a_starting_index)      
    end if    
    call new(irl_cell_flux)    
    call a_flux%clear()
    do while(tracker%still_cells_to_visit())

       current_cell = tracker%get_next_cell()
       
       call setMinimumVolToTrack(a_mesh%volume(current_cell)*1.0e-15_r8)
       call getNormMoments(a_flux_volume, a_reconstruction(current_cell), irl_cell_flux)

       cell_flux(1) = getVolume(irl_cell_flux, 0)
       cell_flux(2) = getVolume(irl_cell_flux, 1)       
   
       if(sum(abs(cell_flux)) > 1.0e-15*a_mesh%volume(current_cell)) then
          call a_flux%add_flux(current_cell, cell_flux)         
          associate( cn => a_mesh%cnhbr(a_mesh%xcnhbr(current_cell):a_mesh%xcnhbr(current_cell+1)-1))
            do n = 1, size(cn)
               if(cn(n) /= 0 .and. tracker%cell_not_encountered(cn(n))) then
                  call tracker%add_cell(cn(n))
               end if
            end do
          end associate   
       end if               
    end do
    
  end subroutine truchas_getMoments_CapDod_LTLT

  subroutine truchas_getMoments_CapDod_LLTT(a_mesh, a_nmat, a_flux_volume, a_starting_index, a_reconstruction, a_flux)

    use traversal_tracker_type, only : traversal_tracker    
    use unstr_mesh_type
    use tagged_volumes_type
    
    type(unstr_mesh), pointer, intent(in) :: a_mesh
    integer, intent(in) :: a_nmat
    type(CapDod_LLTT_type), intent(in) :: a_flux_volume
    integer, intent(in) :: a_starting_index(2)
    type(LocSepLink_type), intent(in) :: a_reconstruction(:)
    type(tagged_volumes), intent(inout) :: a_flux

    type(traversal_tracker) :: tracker
    integer :: current_cell
    type(SepVol_type) :: irl_cell_flux
    real(r8) :: cell_flux(a_nmat)
    integer :: n

    
    if(a_starting_index(2) == 0) then
      call tracker%init([a_starting_index(1)])
    else
      call tracker%init(a_starting_index)      
    end if    
    call new(irl_cell_flux)    
    call a_flux%clear()
    do while(tracker%still_cells_to_visit())

       current_cell = tracker%get_next_cell()
       
       call setMinimumVolToTrack(a_mesh%volume(current_cell)*1.0e-15_r8)
       call getNormMoments(a_flux_volume, a_reconstruction(current_cell), irl_cell_flux)

       cell_flux(1) = getVolume(irl_cell_flux, 0)
       cell_flux(2) = getVolume(irl_cell_flux, 1)       
   
       if(sum(abs(cell_flux)) > 1.0e-15*a_mesh%volume(current_cell)) then
          call a_flux%add_flux(current_cell, cell_flux)         
          associate( cn => a_mesh%cnhbr(a_mesh%xcnhbr(current_cell):a_mesh%xcnhbr(current_cell+1)-1))
            do n = 1, size(cn)
               if(cn(n) /= 0 .and. tracker%cell_not_encountered(cn(n))) then
                  call tracker%add_cell(cn(n))
               end if
            end do
          end associate   
       end if               
    end do
    
  end subroutine truchas_getMoments_CapDod_LLTT

  subroutine truchas_getMoments_CapDod_LTTT(a_mesh, a_nmat, a_flux_volume, a_starting_index, a_reconstruction, a_flux)

    use traversal_tracker_type, only : traversal_tracker    
    use unstr_mesh_type
    use tagged_volumes_type
    
    type(unstr_mesh), pointer, intent(in) :: a_mesh
    integer, intent(in) :: a_nmat
    type(CapDod_LTTT_type), intent(in) :: a_flux_volume
    integer, intent(in) :: a_starting_index(2)
    type(LocSepLink_type), intent(in) :: a_reconstruction(:)
    type(tagged_volumes), intent(inout) :: a_flux

    type(traversal_tracker) :: tracker
    integer :: current_cell
    type(SepVol_type) :: irl_cell_flux
    real(r8) :: cell_flux(a_nmat)
    integer :: n

    
    if(a_starting_index(2) == 0) then
      call tracker%init([a_starting_index(1)])
    else
      call tracker%init(a_starting_index)      
    end if    
    call new(irl_cell_flux)    
    call a_flux%clear()
    do while(tracker%still_cells_to_visit())

       current_cell = tracker%get_next_cell()
       
       call setMinimumVolToTrack(a_mesh%volume(current_cell)*1.0e-15_r8)
       call getNormMoments(a_flux_volume, a_reconstruction(current_cell), irl_cell_flux)

       cell_flux(1) = getVolume(irl_cell_flux, 0)
       cell_flux(2) = getVolume(irl_cell_flux, 1)       
   
       if(sum(abs(cell_flux)) > 1.0e-15*a_mesh%volume(current_cell)) then
          call a_flux%add_flux(current_cell, cell_flux)         
          associate( cn => a_mesh%cnhbr(a_mesh%xcnhbr(current_cell):a_mesh%xcnhbr(current_cell+1)-1))
            do n = 1, size(cn)
               if(cn(n) /= 0 .and. tracker%cell_not_encountered(cn(n))) then
                  call tracker%add_cell(cn(n))
               end if
            end do
          end associate   
       end if               
    end do
    
  end subroutine truchas_getMoments_CapDod_LTTT  

  subroutine truchas_getMoments_CapDod_TTTT(a_mesh, a_nmat, a_flux_volume, a_starting_index, a_reconstruction, a_flux)

    use traversal_tracker_type, only : traversal_tracker    
    use unstr_mesh_type
    use tagged_volumes_type
    
    type(unstr_mesh), pointer, intent(in) :: a_mesh
    integer, intent(in) :: a_nmat
    type(CapDod_TTTT_type), intent(in) :: a_flux_volume
    integer, intent(in) :: a_starting_index(2)
    type(LocSepLink_type), intent(in) :: a_reconstruction(:)
    type(tagged_volumes), intent(inout) :: a_flux

    type(traversal_tracker) :: tracker
    integer :: current_cell
    type(SepVol_type) :: irl_cell_flux
    real(r8) :: cell_flux(a_nmat)
    integer :: n

    
    if(a_starting_index(2) == 0) then
      call tracker%init([a_starting_index(1)])
    else
      call tracker%init(a_starting_index)      
    end if    
    call new(irl_cell_flux)    
    call a_flux%clear()
    do while(tracker%still_cells_to_visit())

       current_cell = tracker%get_next_cell()
       
       call setMinimumVolToTrack(a_mesh%volume(current_cell)*1.0e-15_r8)
       call getNormMoments(a_flux_volume, a_reconstruction(current_cell), irl_cell_flux)

       cell_flux(1) = getVolume(irl_cell_flux, 0)
       cell_flux(2) = getVolume(irl_cell_flux, 1)       
   
       if(sum(abs(cell_flux)) > 1.0e-15*a_mesh%volume(current_cell)) then
          call a_flux%add_flux(current_cell, cell_flux)         
          associate( cn => a_mesh%cnhbr(a_mesh%xcnhbr(current_cell):a_mesh%xcnhbr(current_cell+1)-1))
            do n = 1, size(cn)
               if(cn(n) /= 0 .and. tracker%cell_not_encountered(cn(n))) then
                  call tracker%add_cell(cn(n))
               end if
            end do
          end associate   
       end if               
    end do
    
  end subroutine truchas_getMoments_CapDod_TTTT

  subroutine truchas_getMoments_CapOcta_LLL(a_mesh, a_nmat, a_flux_volume, a_starting_index, a_reconstruction, a_flux)

    use traversal_tracker_type, only : traversal_tracker    
    use unstr_mesh_type
    use tagged_volumes_type
    
    type(unstr_mesh), pointer, intent(in) :: a_mesh
    integer, intent(in) :: a_nmat
    type(CapOcta_LLL_type), intent(in) :: a_flux_volume
    integer, intent(in) :: a_starting_index(2)
    type(LocSepLink_type), intent(in) :: a_reconstruction(:)
    type(tagged_volumes), intent(inout) :: a_flux

    type(traversal_tracker) :: tracker
    integer :: current_cell
    type(SepVol_type) :: irl_cell_flux
    real(r8) :: cell_flux(a_nmat)
    integer :: n

    
    if(a_starting_index(2) == 0) then
      call tracker%init([a_starting_index(1)])
    else
      call tracker%init(a_starting_index)      
    end if    
    call new(irl_cell_flux)    
    call a_flux%clear()
    do while(tracker%still_cells_to_visit())

       current_cell = tracker%get_next_cell()
       
       call setMinimumVolToTrack(a_mesh%volume(current_cell)*1.0e-15_r8)
       call getNormMoments(a_flux_volume, a_reconstruction(current_cell), irl_cell_flux)

       cell_flux(1) = getVolume(irl_cell_flux, 0)
       cell_flux(2) = getVolume(irl_cell_flux, 1)       
   
       if(sum(abs(cell_flux)) > 1.0e-15*a_mesh%volume(current_cell)) then
          call a_flux%add_flux(current_cell, cell_flux)         
          associate( cn => a_mesh%cnhbr(a_mesh%xcnhbr(current_cell):a_mesh%xcnhbr(current_cell+1)-1))
            do n = 1, size(cn)
               if(cn(n) /= 0 .and. tracker%cell_not_encountered(cn(n))) then
                  call tracker%add_cell(cn(n))
               end if
            end do
          end associate   
       end if               
    end do
    
  end subroutine truchas_getMoments_CapOcta_LLL

  subroutine truchas_getMoments_CapOcta_LLT(a_mesh, a_nmat, a_flux_volume, a_starting_index, a_reconstruction, a_flux)

    use traversal_tracker_type, only : traversal_tracker    
    use unstr_mesh_type
    use tagged_volumes_type
    
    type(unstr_mesh), pointer, intent(in) :: a_mesh
    integer, intent(in) :: a_nmat
    type(CapOcta_LLT_type), intent(in) :: a_flux_volume
    integer, intent(in) :: a_starting_index(2)
    type(LocSepLink_type), intent(in) :: a_reconstruction(:)
    type(tagged_volumes), intent(inout) :: a_flux

    type(traversal_tracker) :: tracker
    integer :: current_cell
    type(SepVol_type) :: irl_cell_flux
    real(r8) :: cell_flux(a_nmat)
    integer :: n

    
    if(a_starting_index(2) == 0) then
      call tracker%init([a_starting_index(1)])
    else
      call tracker%init(a_starting_index)      
    end if    
    call new(irl_cell_flux)    
    call a_flux%clear()
    do while(tracker%still_cells_to_visit())

       current_cell = tracker%get_next_cell()
       
       call setMinimumVolToTrack(a_mesh%volume(current_cell)*1.0e-15_r8)
       call getNormMoments(a_flux_volume, a_reconstruction(current_cell), irl_cell_flux)

       cell_flux(1) = getVolume(irl_cell_flux, 0)
       cell_flux(2) = getVolume(irl_cell_flux, 1)       
   
       if(sum(abs(cell_flux)) > 1.0e-15*a_mesh%volume(current_cell)) then
          call a_flux%add_flux(current_cell, cell_flux)         
          associate( cn => a_mesh%cnhbr(a_mesh%xcnhbr(current_cell):a_mesh%xcnhbr(current_cell+1)-1))
            do n = 1, size(cn)
               if(cn(n) /= 0 .and. tracker%cell_not_encountered(cn(n))) then
                  call tracker%add_cell(cn(n))
               end if
            end do
          end associate   
       end if               
    end do
    
  end subroutine truchas_getMoments_CapOcta_LLT

  subroutine truchas_getMoments_CapOcta_LTT(a_mesh, a_nmat, a_flux_volume, a_starting_index, a_reconstruction, a_flux)

    use traversal_tracker_type, only : traversal_tracker    
    use unstr_mesh_type
    use tagged_volumes_type
    
    type(unstr_mesh), pointer, intent(in) :: a_mesh
    integer, intent(in) :: a_nmat
    type(CapOcta_LTT_type), intent(inout) :: a_flux_volume
    integer, intent(in) :: a_starting_index(2)
    type(LocSepLink_type), intent(in) :: a_reconstruction(:)
    type(tagged_volumes), intent(inout) :: a_flux

    type(traversal_tracker) :: tracker
    integer :: current_cell
    type(SepVol_type) :: irl_cell_flux
    real(r8) :: cell_flux(a_nmat)
    integer :: n

    
    if(a_starting_index(2) == 0) then
      call tracker%init([a_starting_index(1)])
    else
      call tracker%init(a_starting_index)      
    end if    
    call new(irl_cell_flux)    
    call a_flux%clear()
    
    do while(tracker%still_cells_to_visit())

       current_cell = tracker%get_next_cell()
       
       call setMinimumVolToTrack(a_mesh%volume(current_cell)*1.0e-15_r8)
       call getNormMoments(a_flux_volume, a_reconstruction(current_cell), irl_cell_flux)

       cell_flux(1) = getVolume(irl_cell_flux, 0)
       cell_flux(2) = getVolume(irl_cell_flux, 1)   
   
       if(sum(abs(cell_flux)) > 1.0e-15*a_mesh%volume(current_cell)) then
          call a_flux%add_flux(current_cell, cell_flux)         
          associate( cn => a_mesh%cnhbr(a_mesh%xcnhbr(current_cell):a_mesh%xcnhbr(current_cell+1)-1))
            do n = 1, size(cn)
               if(cn(n) /= 0 .and. tracker%cell_not_encountered(cn(n))) then
                  call tracker%add_cell(cn(n))
               end if
            end do
          end associate   
       end if               
    end do
    
  end subroutine truchas_getMoments_CapOcta_LTT

  subroutine truchas_getMoments_CapOcta_TTT(a_mesh, a_nmat, a_flux_volume, a_starting_index, a_reconstruction, a_flux)

    use traversal_tracker_type, only : traversal_tracker    
    use unstr_mesh_type
    use tagged_volumes_type
    
    type(unstr_mesh), pointer, intent(in) :: a_mesh
    integer, intent(in) :: a_nmat
    type(CapOcta_TTT_type), intent(in) :: a_flux_volume
    integer, intent(in) :: a_starting_index(2)
    type(LocSepLink_type), intent(in) :: a_reconstruction(:)
    type(tagged_volumes), intent(inout) :: a_flux

    type(traversal_tracker) :: tracker
    integer :: current_cell
    type(SepVol_type) :: irl_cell_flux
    real(r8) :: cell_flux(a_nmat)
    integer :: n

    
    if(a_starting_index(2) == 0) then
      call tracker%init([a_starting_index(1)])
    else
      call tracker%init(a_starting_index)      
    end if    
    call new(irl_cell_flux)    
    call a_flux%clear()
    do while(tracker%still_cells_to_visit())

       current_cell = tracker%get_next_cell()
       
       call setMinimumVolToTrack(a_mesh%volume(current_cell)*1.0e-15_r8)
       call getNormMoments(a_flux_volume, a_reconstruction(current_cell), irl_cell_flux)

       cell_flux(1) = getVolume(irl_cell_flux, 0)
       cell_flux(2) = getVolume(irl_cell_flux, 1)       
   
       if(sum(abs(cell_flux)) > 1.0e-15*a_mesh%volume(current_cell)) then
          call a_flux%add_flux(current_cell, cell_flux)         
          associate( cn => a_mesh%cnhbr(a_mesh%xcnhbr(current_cell):a_mesh%xcnhbr(current_cell+1)-1))
            do n = 1, size(cn)
               if(cn(n) /= 0 .and. tracker%cell_not_encountered(cn(n))) then
                  call tracker%add_cell(cn(n))
               end if
            end do
          end associate   
       end if               
    end do
    
  end subroutine truchas_getMoments_CapOcta_TTT
  
end module truchas_getMoments_mod
