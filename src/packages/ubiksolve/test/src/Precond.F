! -*- Mode: f90 -*-

!*********************************************************************
!
! $Id: Precond.F,v 1.12 2003/10/19 22:27:12 turner Exp $
! 
!*********************************************************************
module Precond_module
  use Ubik_module
  implicit none

  interface GeneratePrecond
     module procedure GeneratePrecond_Full
     module procedure GeneratePrecond_ELL
  end interface

  private
  public :: &
       precond_control_type, &
       precond_control, &
       SetPrecondDefaults, &
       ApplyPrecond, &
       GeneratePrecond

  type precond_control_type

     ! precond - controls preconditioning
     !    1  ==>  m-step Jacobi iteration
     !              (also uses istep and omega)
     !    2  ==>  m-step SSOR iteration
     !              (also uses istep and omega)
     !    3  ==>  Incomplete Cholesky
     !              (also uses omega)
     !    4  ==>  Incomplete LU
     !              (also uses omega)
     ! 
     ! Recommendations:
     !   - precond = 1 and precond = 2 are equivalent to precond = 3 with
     !     steps = 1 and omega = 1.0.  All three are also equivalent to simply
     !     diagonally scaling the original system prior to solution.  Hence
     !     neither of these are recommended (precond = 1 and precond = 2 are
     !     really only useful as templates while developing other
     !     preconditioners).
     ! 
     !   - M-step Jacobi (precond = 3), with steps > 1, vectorizes well, but
     !     is typically not very effective.
     ! 
     !   - M-step SSOR is often very effective, but is difficult to
     !     vectorize (the implementation in UbikSolve is not vectorized).
     ! 
     !   - For the incomplete factorization preconditioners, omega
     !     controls the degree to which fill-in elements are added back
     !     into the diagonal:
     !       omega = 0  ==>  standard incomplete factorization
     !       omega = 1  ==>  modified incomplete factorization
     !       0<omega<1  ==>  relaxed incomplete factorization
     ! 
     !   - The incomplete factorization preconditioners can be quite
     !     effective in reducing the number of iterations, but are
     !     often so expensive to compute as to negate any benefit in
     !     using them.  In addition both the generation and application
     !     of incomplete factorization preconditioners is difficult to
     !     vectorize, and the implementations in UbikSolve are not 
     !     vectorized.  One way to make them more useful when solving a
     !     number of linear systems with similar coefficients (e.g. a 
     !     time-dependent problem) is to amortize the cost of computing
     !     them over several systems.  That is, call the solver routine
     !     initially with precond_in = 0, and have the solver generate and
     !     use one of the incomplete factorization preconditioners.
     !     Then on subsequent calls to the solver, set precond_in = 1, so
     !     that the already-computed preconditioner will be reused.
     ! 
     !   - While incomplete Cholesky (IC) factorization preconditioners
     !     are appropriate only for symmetric coefficients, incomplete
     !     LU (ILU) factorization preconditioners can be used for
     !     either symmetric or nonsymmetric coefficients.  Typically,
     !     for a symmetric coefficient, one of the IC factorization
     !     should be preferrable to one of the ILU factorization, since
     !     the number of operations required to generate an IC
     !     factorization would be about 1/2 that required for an ILU
     !     factorization.  However, if the coefficient is stored in ELL
     !     format, the nature of the data structure causes generation
     !     of an IC factorization to be highly inefficient.  Hence, for
     !     coefficients stored in ELL format, ILU factorization
     !     preconditioners are preferrable to IC factorization
     !     preconditioners, even if the coefficient is symmetric.
     integer(Ubik_int_type) :: precond

     ! steps - number of Jacobi or SSOR steps to take when
     !         precond = 3 or 4 [DEFAULT = 1]
     integer(Ubik_int_type) :: steps

     ! omega - relaxation parameter [DEFAULT = 1.0]
     real(Ubik_real_type) :: omega

     real(Ubik_real_type) :: tiny

  end type precond_control_type

  type(precond_control_type) :: precond_control

contains

  !*********************************************************************
  ! Sets reasonable default values for control parameters.
  ! 
  ! <PARAMETER LIST>
  ! 
  !  In/Out:
  !   precond_control - structure containing preconditioner control parameters
  !*********************************************************************
  subroutine SetPrecondDefaults (precond_control)
    type(precond_control_type), intent(out) :: precond_control

    precond_control%precond = Ubik_Precond_none ! no preconditioning
    precond_control%steps = 1                 ! take 1 step in m-step preconditioners
    precond_control%omega = 1.0               ! relaxation parameter
    precond_control%tiny = TINY(0.0)         ! smallest number to worry about

    return
  end subroutine SetPrecondDefaults

  subroutine Pre_Full (coeff, precond, b, x_vec, status)
    real(Ubik_real_type), intent(in), dimension(:,:) :: coeff, precond

#include "Precond-guts.F"

  end subroutine Pre_Full

  !*********************************************************************

  subroutine Pre_ELL (coeff, precond, b, x_vec, status)
    type(Ubik_ELL_matrix), intent(in) :: coeff, precond

#include "Precond-guts.F"

  end subroutine Pre_ELL

  !*********************************************************************
  ! 
  ! Applies preconditioner as appropriate (solve Ap*x = b).
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   n - number of unknowns
  !   b - b-vector (source)
  ! 
  !  Output:
  !   x - x-vector (unknown)
  !   status - return status:
  !     -3  ==>  internal error
  !     -1  ==>  invalid argument(s)
  !      0  ==>  success
  !
  !*********************************************************************
  subroutine ApplyPrecond (b, x_vec, status)
    use Full_data_module
    use ELL_data_module
    use test_data_module

    real(Ubik_real_type), intent(in), dimension(:), target :: b
    type(Ubik_vector_type), intent(inout) :: x_vec
    integer(Ubik_int_type), intent(out) :: status

    select case (prob)
    case (Saad_Full)
       call Pre_Full (coeff_full, precond_full, b, x_vec, status)
    case (Saad_ELL)
       call Pre_ELL (coeff_ell, precond_ell, b, x_vec, status)
    case default
       x_vec = b
    end select

  end subroutine ApplyPrecond

  !*********************************************************************
  ! 
  ! Computes an approximation to a or an approximation to inv(a)
  ! as appropriate, to be used as a preconditioner.
  ! 
  ! Note that for both implementations of Jacobi preconditioning,
  ! (control%precond=1 and control%precond=2), the main diagonal is
  ! assumed to be in the first column.
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   control - structure containing control parameters
  !   coeff - coeff. matrix
  ! 
  !  Out:
  !   precond - preconditioner
  ! 
  ! <SUBROUTINES REQUIRED>
  ! 
  !  Ubik_IC
  !  Ubik_ILU
  !  Ubik_WriteMatrix
  !
  ! <UNDOCUMENTED FEATURE>
  !
  !  Two simple implementations of Jacobi preconditioning have been
  !  implemented, primarily as templates for development of other
  !  preconditioners.
  !
  !  These preconditioners can be accessed by setting ipre to -2 or -3:
  !      -2 ==> Jacobi approximation of the coefficient
  !      -3 ==> Jacobi approximation of the inverse of coefficient
  !
  !*********************************************************************
  subroutine GeneratePrecond_Full (control, coeff, precond)
    real(Ubik_real_type), intent(in), dimension(:,:) :: coeff
    real(Ubik_real_type), intent(out), dimension(:,:) :: precond

#include "GeneratePrecond-guts.F"

  end subroutine GeneratePrecond_Full

  !*********************************************************************

  subroutine GeneratePrecond_ELL (control, coeff, precond)
    type(Ubik_ELL_matrix), intent(in) :: coeff
    type(Ubik_ELL_matrix), intent(out) :: precond

#include "GeneratePrecond-guts.F"

  end subroutine GeneratePrecond_ELL

end module Precond_module
