!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE FISCHER_MODULE
 !=======================================================================
  ! Purpose(s):
  !
  !   Implements the Fischer intial guess technique from
  !   
  !   P. F. Fischer, "Projection Techniques for iterative solution of Ax=b with
  !    successive right-hand sides.  Comput. Methods Appl. Mech. Engrg. v. 163
  !    pp. 193--204, 1998
  !
  !   Public Interface:
  !
  !     * call FISCHER_INITIAL_GUESS( x_guess, b )
  !       call FISCHER_UPDATE_SPACE( x_soln, b )
  !
  !
  ! Contains: 
  !           FISCHER_INITIAL_GUESS
  !           FISCHER_UPDATE_SPACE
  !
  !=======================================================================

  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private 

  real(r8), pointer, public, save, dimension(:,:) :: b_tilde
  real(r8), pointer, public, save, dimension(:,:) :: x_tilde
  real(r8), pointer, public, save, dimension(:)   :: x_guess

  integer, parameter :: max_num_vecs = 6
  integer :: cur_num_vecs

  public:: FISCHER_INITIALIZE, FISCHER_INITIAL_GUESS, FISCHER_UPDATE_SPACE

CONTAINS

  SUBROUTINE FISCHER_INITIALIZE ()
    use legacy_mesh_api, only: ncells

    integer :: status

    ALLOCATE( b_tilde(max_num_vecs,ncells), &
         x_tilde(max_num_vecs,ncells), x_guess(ncells), &
         STAT=status)

    if (status /= 0) call TLS_panic ('FISCHER: allocation failed')

    cur_num_vecs = 0

    x_guess(:) = 0.0_r8

  END SUBROUTINE FISCHER_INITIALIZE

  SUBROUTINE FISCHER_INITIAL_GUESS( x, b )
    !=======================================================================
    ! Purpose(s):
    !
    !   Given the new right hand side b, this computes a new guess for
    !   the solution x based on previous stored solution values
    !
    !   CAUTION: The space of vectors is updated after the true solution
    !    to Ax=b is found.  That updated requires the guess we calculate
    !    here.  As such, the guess is stored in the module (x_guess).  
    !    The usual code flow will be
    !
    !      call Fischer_Initial_Guess( x', b)
    !      solve Ax=b using x'
    !      call Fischer_update_space( x, b )
    !
    !    If Fischer_Initial_Guess is called repeatedly before the update
    !    then the update will be wrong.
    !
    !=======================================================================
    use legacy_mesh_api, only: ncells
    use pgslib_module, only: PGSLib_Global_DOT_PRODUCT

    real(r8), dimension(ncells) :: x
    real(r8), dimension(ncells) :: b

    integer :: i
    real(r8), dimension(max_num_vecs) :: alpha

    do i=1, cur_num_vecs
       alpha(i) = PGSLib_Global_DOT_PRODUCT(b,b_tilde(i,:))
    enddo

    x_guess(:) = 0

    do i=1, cur_num_vecs
       x_guess(:) = x_guess(:) + alpha(i)*x_tilde(i,:)
    enddo

    x(:) = x_guess(:)

  END SUBROUTINE FISCHER_INITIAL_GUESS


  SUBROUTINE FISCHER_UPDATE_SPACE( x_soln, MatVec, b )
    !=======================================================================
    ! Purpose(s):
    !
    !   After a PPE solution, the set of projection vectors (i.e. previous
    !   solutions and right hand sides) is updated.
    !
    !=======================================================================
    use legacy_mesh_api, only: ncells
    use cutoffs_module, only: alittle
    use UbikSolve_module
    use pgslib_module, only: PGSLib_Global_DOT_PRODUCT

    ! For MatVec definition
#include "solver_function_prototypes.fpp"

    real(r8), dimension(ncells) :: x_soln
    real(r8), dimension(ncells) :: b
    real(r8) :: b_norm

    real(r8), dimension(max_num_vecs) :: alpha

    type(Ubik_vector_type), target :: x_vec

    integer :: i, status

    ! MatVec takes a Ubik_vector_type as the vector to be multiplied,
    ! so set one up to contain x
    call Ubik_create (x_vec, overlap_only=.true.)

    if (cur_num_vecs .eq. max_num_vecs) then

       ! \tilde b = M x
       x_tilde(1,:) = x_soln(:)
       call Ubik_set_values_ptr(x_vec, x_tilde(1,:))
       call MatVec( x_vec, b_tilde(1,:), status )

       ! \tilde b <- \tilde b/ \| \tilde b \|
       ! \tilde x <- \tilde x/ \| \tilde b \|
       b_norm = SQRT(PGSLib_Global_DOT_PRODUCT(b_tilde(1,:),b_tilde(1,:)))

       if (b_norm .gt. alittle) then
          b_tilde(1,:) = b_tilde(1,:)/b_norm
          x_tilde(1,:) = x_tilde(1,:)/b_norm

          cur_num_vecs = 1
       endif
    else
       cur_num_vecs = cur_num_vecs + 1

       x_tilde(cur_num_vecs,:) = x_soln(:) - x_guess(:)

       ! \tilde b = M \tilde x
       call Ubik_set_values_ptr(x_vec, x_tilde(cur_num_vecs,:))
       call MatVec( x_vec, b_tilde(cur_num_vecs,:), status )

       do i=1, cur_num_vecs-1
          alpha(i) = PGSLib_Global_DOT_PRODUCT(b_tilde(cur_num_vecs,:),b_tilde(i,:))
       enddo

       do i=1, cur_num_vecs-1
          b_tilde(cur_num_vecs,:) = b_tilde(cur_num_vecs,:) - alpha(i)*b_tilde(i,:)
          x_tilde(cur_num_vecs,:) = x_tilde(cur_num_vecs,:) - alpha(i)*x_tilde(i,:)
       enddo

       b_norm = SQRT(PGSLib_Global_DOT_PRODUCT(b_tilde(cur_num_vecs,:),b_tilde(cur_num_vecs,:)))

       if (b_norm .gt. alittle) then
          b_tilde(cur_num_vecs,:) = b_tilde(cur_num_vecs,:)/b_norm
          x_tilde(cur_num_vecs,:) = x_tilde(cur_num_vecs,:)/b_norm
       else
          ! Disregard the vector
          cur_num_vecs = cur_num_vecs - 1
       endif
    endif

  END SUBROUTINE FISCHER_UPDATE_SPACE

END MODULE FISCHER_MODULE
