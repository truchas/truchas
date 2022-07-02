!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

MODULE MATL_UTILITIES

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  IMPLICIT NONE
  private

  public :: MATL_GET_VOF  ! These two enable removing Matl from the volume tracker.
  public :: MATL_SET_VOF  !    - Markus
  public :: define_matl, read_matl_data, matl_get_cell_vof

  public :: UPDATE_MATL   ! Moved here from heat_transfer.F90 by Sriram

CONTAINS

  SUBROUTINE MATL_GET_VOF (VOF)
    ! ==========================================================================
    ! Sriram originally wrote the skeleton for a routine that would allow access
    ! to any part of Matl.  But I don't think that that'll every be necessary,
    ! as the only part that's really used much are the current Vof values.  And
    ! so I've taken the liberty of simplifying this to a routine that fills an
    ! array (nmat,ncells) with current VOF values.
    !
    ! Written by:
    ! Markus Bussmann (University of Toronto)
    !===========================================================================
    use matl_module, only: Matl, ncells, nmat, mat_slot

    ! Arguments
    real(r8), dimension(nmat,ncells), intent(INOUT) :: VOF

    ! Local Variables
    integer :: m, n, s

    VOF = 0.0_r8

    do n = 1,ncells
       do s = 1,mat_slot
          m = Matl(s)%Cell(n)%Id
          if (m > 0) then
             VOF(m,n) = Matl(s)%Cell(n)%Vof
          end if
       end do
    end do

  END SUBROUTINE MATL_GET_VOF
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! MATL_GET_CELL_VOF
 !!
 !! Neil Carlson <nnc@lanl.gov> 28 Jul 2009
 !!
 !! Returns the material volume fractions for cell N in the array VOF.
 !!

  subroutine matl_get_cell_vof (n, vof)
  
    use matl_module, only: matl, ncells, nmat, mat_slot
    
    integer, intent(in) :: n
    real(r8), intent(out) :: vof(:)
    
    integer :: s, m
    
    ASSERT(size(vof) == nmat)
    ASSERT(n >= 1 .and. n <= ncells)
    
    vof = 0.0
    do s = 1, mat_slot
      m = matl(s)%cell(n)%id
      if (m > 0) vof(m) = matl(s)%cell(n)%vof
    end do
    
  end subroutine matl_get_cell_vof

  SUBROUTINE MATL_SET_VOF (VOF)
    ! ==========================================================================
    ! Equivalent to MATL_GET_VOF, but this routine writes VOF data back into Matl.
    ! A little more complicated than the GET routine, because the number of slots
    ! in Matl may not be the same as the maximum number of materials in a cell.
    !
    ! Written by:
    ! Markus Bussmann (University of Toronto)
    !===========================================================================
    use matl_module, only: Matl, SLOT_COMPRESS, SLOT_DECREASE, SLOT_INCREASE, ncells, &
        nmat, mat_slot
    use parallel_communication, only: global_maxval

    ! Arguments
    real(r8), dimension(nmat,ncells), intent(INOUT) :: VOF

    ! Local Variables
    integer, dimension(ncells) :: mat
    integer :: m, n, s, max_mat

    ! Begin by compressing the current slot structure.
    call SLOT_COMPRESS (Matl, mat_slot)

    ! Determine the maximum number of materials in any one cell.
    mat = 0
    do m = 1,nmat
       where (Vof(m,:) > 0.0_r8) mat = mat + 1
    end do
    max_mat = global_maxval(mat)

    if (max_mat < mat_slot) then
       CALL SLOT_DECREASE (Matl, mat_slot, max_mat)
       mat_slot = max_mat
    else if (max_mat > mat_slot) then
       CALL SLOT_INCREASE (Matl, mat_slot, max_mat)
       mat_slot = max_mat
    end if

    ! Now insert VOF values that are zero, to open slots for new materials.
    do n = 1,ncells
       MATERIALS_1: do m = 1,nmat
          if (Vof(m,n) == 0.0_r8) then
             do s = 1,mat_slot
                if (Matl(s)%Cell(n)%Id == m) then
                   Matl(s)%Cell(n)%Id = 0
                   Matl(s)%Cell(n)%Vof = 0.0_r8
                   CYCLE MATERIALS_1
                end if
             end do
          end if
       end do MATERIALS_1
    end do
 
    ! And now update existing values, and insert new values.
    do n = 1,ncells
       MATERIALS_2: do m = 1, nmat

          if (Vof(m,n) == 0.0_r8) cycle MATERIALS_2

          ! Look for an existing Matl slot for material number m.
          do s = 1,mat_slot
             if (Matl(s)%Cell(n)%Id == m) then
                Matl(s)%Cell(n)%Vof = VOF(m,n)
                cycle MATERIALS_2
             end if
          end do

          ! Otherwise look for an open slot that needs to be filled.
          do s = 1,mat_slot
             if (Matl(s)%Cell(n)%Id == 0) then
                Matl(s)%Cell(n)%Id = m
                Matl(s)%Cell(n)%Vof = VOF(m,n)
                cycle MATERIALS_2
             end if
          end do

       end do MATERIALS_2
    end do

  END SUBROUTINE MATL_SET_VOF

  SUBROUTINE UPDATE_MATL (VF_NEW)
    !===========================================================================
    ! Purpose:
    !
    !  Update the matl pointer with the new  volume fractions and  densitiy
    !  values
    !===========================================================================
    use matl_module, only: Matl, SLOT_INCREASE, SLOT_DECREASE, ncells, &
        mat_slot, mat_slot_new, nmat
    use parallel_communication, only: global_maxval

    ! Arguments
    real(r8), dimension(0:nmat, 1:ncells), intent(IN) :: VF_New

    ! Local Variables
    integer :: i, m, s, slots_needed
    real(r8), dimension(ncells) :: nslots_cell

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Check the number of slots and increase if needed
    nslots_cell = 0
    do m = 1,nmat
       where (VF_New(m,:) > 0.0_r8) nslots_cell = nslots_cell + 1
    end do
    slots_needed = global_maxval(nslots_cell)
    if (slots_needed > mat_slot) then
       mat_slot_new = slots_needed
       call SLOT_INCREASE (Matl, mat_slot, mat_slot_new)
    end if

    ! Updating the matl database with the new volume fractions
    do i=1,ncells
       nslots_cell(i)=1

       ! Resetting the information for the materials by slots
       do s = 1,mat_slot
          matl(s)%cell(i)%id = 0
          matl(s)%cell(i)%vof = 0.0
       end do
       ! Updating the information in slots, by materials
       s = nslots_cell(i)
       do m = 1,nmat
          if (VF_New(m,i) > 0) then
             Matl(s)%Cell(i)%Id = m
             Matl(s)%Cell(i)%Vof  = VF_New(m,i)
             s = s+1
          end if
       end do
       nslots_cell(i) = s-1
    end do

    ! Check the number of slots and DECREASE if needed
    slots_needed = global_maxval(nslots_cell)
    if (slots_needed < mat_slot) then
       mat_slot_new = slots_needed
       call SLOT_DECREASE(Matl, mat_slot, mat_slot_new)
    end if

  END SUBROUTINE UPDATE_MATL

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DEFINE_MATL
 !!
 !! Neil Carlson <nnc@lanl.gov>
 !! 14 Apr 2005
 !!
 !! Given the volume fraction array VF, covering all materials and cells, this
 !! routine initializes the data in the MATL structure.  The first index of
 !! the VF array is the material index, and the second is the cell index.
 !!
 !! NB: MATL is intent(inout) because it has pointer components that may be
 !! associated.  None of its values are used, however, and all are overwritten.
 !!

  subroutine define_matl (vf)

    use matl_module, only: matl, material, slot_resize, ncells, nmat, mat_slot
    use parallel_communication, only: global_maxval

    real(r8), intent(in)    :: vf(:,:)

    integer :: j, m, s
    type(material) :: null_mat  ! This is default initialized.
    type(material), allocatable :: mlist(:)

    ASSERT( size(vf,1) <= nmat )
    ASSERT( size(vf,2) == ncells )

    !! Find the max number of materials in any one cell and resize MATL accordingly.
    m = global_maxval(count(vf > 0.0_r8, dim=1))
    call slot_resize (matl, mat_slot, m)

    !! Set up the material list; only the VOF values change from cell to cell.
    allocate(mlist(size(vf,1)))  ! This is default initialized.
    do m = 1, size(mlist)
      mlist(m)%ID  = m
    end do

    !! Define MATL.
    do j = 1, ncells
      mlist%vof = vf(:,j) ! stuff in the VF values for this cell.
      !! Pack the valid part of the material list into MATL; can't use the pack intrinsic :(
      s = 1
      do m = 1, size(mlist)
        if (vf(m,j) <= 0.0_r8) cycle
        matl(s)%cell(j) = mlist(m)
        s = s + 1
      end do
      !! Fill the rest of MATL with null material data.
      do s = s, mat_slot
        matl(s)%cell(j) = null_mat
      end do
    end do

    deallocate(mlist)

  end subroutine define_matl

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_MATL_DATA
 !!
 !! Neil N. Carlson <nnc@lanl.gov>
 !! 18 Apr 2005
 !!
 !! This subroutine reads the volume fraction data from a restart file opened
 !! (and pre-positioned) on UNIT, and uses this data (properly distributed and
 !! permuted) to define the module structure MATL.  VERSION is the version
 !! number of the restart file format.
 !!
 !! NB: It is implicitly assumed that the material indices in the restart file
 !! directly correspond to the material indices generated by the input file.
 !! We require that the number of materials agree between the input and restart
 !! files; this constraint could probably be relaxed to allow the input to
 !! specify more (new additional) materials.
 !!

  subroutine read_matl_data (unit, version)

    use matl_module, only: nmat
    use restart_utilities, only: read_var, read_dist_array, halt
    use string_utilities, only: i_to_c
    use base_mesh_class
    use mesh_manager, only: named_mesh_ptr

    integer, intent(in) :: unit, version

    integer :: n
    real(r8), allocatable :: vf(:,:)
    class(base_mesh), pointer :: mesh

    mesh => named_mesh_ptr('MAIN')
    INSIST(associated(mesh))

    !! Read the number of materials defined in the restart file.
    call read_var (unit, n, 'READ_MATL_DATA: error reading NMAT record')
    if (n /= nmat) call halt ('READ_MATL_DATA: incompatible NMAT value: ' // i_to_c(n))

    !! Read the volume fraction array.
    allocate(vf(n,mesh%ncell_onP))
    call read_dist_array (unit, vf, mesh%xcell(:mesh%ncell_onP), 'READ_MATL_DATA: error reading VF records')

    !! Derive the MATL structure from the volume fraction array.
    call define_matl (vf)
    deallocate(vf)

  end subroutine read_matl_data

END MODULE MATL_UTILITIES
