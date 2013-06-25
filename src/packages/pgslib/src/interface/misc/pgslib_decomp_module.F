MODULE PGSLib_Decomp_Module
  !======================================================================
  !  PURPOSE
  !    Misc routines for determining decompostions
  !======================================================================

  ! $Id: pgslib_decomp_module.F,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $

  USE PGSLib_Type_MODULE

  implicit none
  PRIVATE
  PUBLIC :: PGSLib_Block_Decompose, PGSLib_Block_SIZE

  INTERFACE PGSLib_Block_Decompose
     MODULE PROCEDURE block_parameters
  END INTERFACE

  INTERFACE PGSLib_Block_Size
     MODULE PROCEDURE block_size
  END INTERFACE

CONTAINS

  function block_size(N_Tot, PE)
    !======================================================================
    !  PURPOSE
    !    For N_tot total elements return number of elements on this PE,
    !    assuming block decomposition.
    !======================================================================
    USE pgslib_utility_module

    ! arguments
    integer (PGSLib_Int_Type), intent(IN) :: N_Tot
    integer (PGSLib_Int_Type), intent(IN), optional :: PE

    ! function return    
    integer (PGSLib_Int_Type) :: block_size

    ! local variables
    integer (PGSLib_Int_Type) :: max_block, min_block, critPE, thisPE

    if (present(PE)) then
       thisPE = PE
    else
       thisPE = PGSLib_Inquire_thisPE()
    end if
    call PGSLib_Block_Decompose(max_block, min_block, critPE, N_Tot)

    if (thisPE <= critPE) then
       block_size = max_block
    else
       block_size = min_block
    end if

    return
  end function block_size

  SUBROUTINE block_parameters(max_block, min_block, critPE, N_Tot)
    !=======================================================================
    ! PURPOSE -
    !   For a given N_Tot, determine block distribution for N_tot.
    !   PE's <= critPE have block size max_block.  PE's > critPT
    !   have blcok size = min_block.
    !   The algorithm used here has (max_block - min_block) = 1
    !=======================================================================
    USE pgslib_utility_module

    ! Arguments
    integer (PGSLib_Int_Type), intent(OUT) :: max_block, min_block, critPE
    integer (PGSLib_Int_Type), intent(IN ) :: N_Tot

    ! Local variables
    integer (PGSLib_Int_Type) :: nPE


    nPE = PGSLib_Inquire_nPE()
    ! First find the average number of data items in a block, then truncate
    ! to an integer, to get min_block
    min_block = FLOOR((REAL(N_tot)/REAL(nPE)))

    ! Max_block is one larger than min_block
    max_block = min_block + 1

    ! critPE is the cutoff PE
    ! (Counting of PEs is one based at F90 level)
    critPE = (N_Tot - min_block*nPE) 

    return
  end SUBROUTINE block_parameters


end MODULE PGSLib_Decomp_Module
