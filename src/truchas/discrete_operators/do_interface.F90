MODULE DO_INTERFACE
  !=======================================================================
  ! Purpose(s):
  !   Provides model developer interface into discrete_operators
  !   public data structures, subroutines and functions.
  !   Developers should use discrete_operators interfaces available through 
  !   this module ONLY.
  !
  !   Public Data Structures: DO_Specifier
  !                           DO_Diag_Specifier
  !                           DO_DiagListSpec
  !                           cMat_row_type
  !
  !   Public Parameters:      DO_SOLVE_DEFAULT
  !                           DO_SOLVE_ORTHO
  !                           DO_SOLVE_LU_LSLR
  !                           DO_SOLVE_SVD_LSLR
  !                           DO_NUM_ST
  !
  !   Public Interfaces:      DO_INIT_SS
  !                           DO_DESTROY_SS
  !                           DO_DESTROY_cM_compressed
  !                           DO_DESTROY_cM_full
  !                           DO_GRADIENT_FACE
  !                           DO_FACE_SOLVE
  !                           DO_UPDATE_WEIGHTS
  !                           DO_GoodSolution
  !                           DO_GoodPhiSolution
  !
  ! DO_ use documentation:
  !
  ! The interface has been broken into "multiple steps" which (hopefully) will make it a little more 
  ! transparent to use.
  !
  ! The user is now free to create as many "solve specifiers" as they wish in a model. A solve specifier 
  ! is a data structure which points to all of the detail that has to be kept track of for a given type 
  ! of solve (i.e. ORTHO or one of the LSLR techniques). Prior to the creation step (i.e. prior to the 
  ! call to do_init_ss), the solve specifier must be a NULL pointer of type DO_Specifier. There is a 
  ! choice of 3 solve techniques: DO_SOLVE_ORTHO (our "good" old ORTHO solve), DO_SOLVE_LU_LSLR and 
  ! DO_SOLVE_SVD_LSLR (LSLR solution via "old" LU with full pivot and LSLR via "new" SVD).
  ! 
  ! Given the untested nature of the SVD sovle, it is desirable to start work as much as possible with 
  ! both LU and SVD techniques to help verify that (when they work), they are equivalent. The SVD solution
  ! has a number of properties easier to understand and manipulate than the LU and so it would also make
  ! sense to eventually make SVD the default.
  ! 
  ! The solve specifier is initialized by a call to do_init_ss, performed only once per solve specifier.
  ! An existing solve specifier may be "destroyed" using DO_DESTROY_SS and then re-initialized if desired.
  ! 
  ! The call to do_update_weights updates the solve technique associated with a give solve specifier based
  !  on a set of updated weights. The update algorithm is designed to update only those components for 
  ! which the weights have changed (since initialization or the last call to do_update_weights). Note that
  ! this capability has also had limited testing owing to the lack of problems in the test suite actually
  ! using LSLR.
  ! 
  ! As before, DO_GRADIENT_FACE produces the gradient and phi face values. One advantage of the SVD 
  ! implementation is that it is now possible to write routines which return only the gradient or 
  ! face value, never having to have calculated the unwanted component. But there hasn't been any real 
  ! call to implement this yet and it would be SVD only.
  ! 
  ! I have also implemented functions: GoodPhiSolution (which face Phi values are "good") and GoodSolution
  ! (which of the Grad & Phi solutions is "good"). Note that "good" currently means "were any components 
  ! pivoted out of the solution" (LU) or "were any SVD weights less than 10.e-6 orders of magnitude less 
  ! than the maximum weight (SVD).
  !
  ! There's also code written, by not currently active to set the "goodness" based on a calculation of 
  ! "standard uncertainty". Since there's been no real research into how masking based on standard 
  ! uncertainty would affect solution quality, the code is currently commented out.
  ! 
  ! I should note that there's a lot of recent work as of the writing of this documentation indicating 
  ! that solution quality is also dependent on the actual field being solved for. In such cases, all the 
  ! operator indicators may say a solution is "good", but when compared to a known, analytic test 
  ! function, the solutions are in fact "bad". It is an open question whether these pathologic cases have 
  ! any relevance to real problems.
  !
  ! J Durachta, 01/27/04
  !
  ! Author(s): Doug Kothe (dbk@lanl.gov)
  !            Jeff Durachta (durachta@verizon.net)
  !            Robert Ferrell (ferrell@cpca.com)
  !
  !=======================================================================
  use do_base_types,         only: DO_Specifier,DO_Diag_Specifier, DO_DiagListSpec, &
                                   DO_SOLVE_DEFAULT,DO_SOLVE_ORTHO, &
                                   DO_SOLVE_LU_LSLR,DO_SOLVE_SVD_LSLR,DO_NUM_ST, &
                                   cMat_row_type
  use do_discrete_operators, only: DO_GRADIENT_FACE,DO_FACE_SOLVE,DO_UPDATE_WEIGHTS, &
                                   DO_GoodSolution,DO_GoodPhiSolution
  use do_solve_specifier,    only: DO_INIT_SS, DO_DESTROY_SS, &
                                   DO_GET_cM_compressed, DO_DESTROY_cM_compressed, &
                                   DO_GET_cM_full, DO_DESTROY_cM_full

  implicit none
  private

  ! Public Data Structures
  public :: DO_Specifier
  public :: DO_Diag_Specifier
  public :: DO_DiagListSpec
  public :: cMat_row_type

  ! Public Parameters
  public :: DO_SOLVE_DEFAULT
  public :: DO_SOLVE_ORTHO
  public :: DO_SOLVE_LU_LSLR
  public :: DO_SOLVE_SVD_LSLR
  public :: DO_NUM_ST

  ! Public Subroutines
  public :: DO_INIT_SS, DO_DESTROY_SS, DO_DESTROY_cM_compressed, DO_DESTROY_cM_full
  public :: DO_GRADIENT_FACE, DO_FACE_SOLVE, DO_UPDATE_WEIGHTS

  ! Public Functions
  public :: DO_GoodSolution, DO_GoodPhiSolution
  public :: DO_GET_cM_compressed, DO_GET_cM_full

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

END MODULE DO_INTERFACE
