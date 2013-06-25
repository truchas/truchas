! -*- Mode: f90 -*-

!*********************************************************************
!
! $Id: UbikMatrix.F,v 1.5 2003/10/19 22:12:48 turner Exp $
!
! Primary interface to UbikSolve.  Host codes should access UbikSolve data types,
! routines, etc. only through this module (or one like it).
!
! The prefix Ubik_ is prepended to all publicly exposed entities in order to
! help guard against namespace collision.  If it should be a problem in
! spite of this, alternative interfaces to UbikSolve can easily be created by
! copying this module and replacing Ubik_ with some other Fortran-legal string.
!
! The "only" form of "use" is used, so all entities are listed, even if not
! renamed (some have Ubik_ prepended already).  This makes 
!
! Note that we need not worry about public/private here, since that is
! controlled within the modules comprising UbikSolve - everything they expose
! externally can be exposed to host codes.
!
!*********************************************************************
module UbikMatrix_module
  use UbikVector_module

  use ubik_constants_module, only: &
       Ubik_Full => Full, &
       Ubik_ELL => ELL

  use Ubik_Full_module, only: &
       Ubik_MatVec => MatVec, &
       Ubik_Diag => Diag, &
       Ubik_Norm => Norm, &
       Ubik_Dump => Dump, &
       Ubik_EnsureFormat => EnsureFormat, &
       Ubik_IC => IC, &
       Ubik_ILU => ILU, &
       Ubik_Jacobi => Jacobi, &
       Ubik_Load => Load, &
       Ubik_LUdecomp => LUdecomp, &
       Ubik_ScaleMatrixByDiagonal => ScaleMatrixByDiagonal, &
       Ubik_ScaleMatrixByRow => ScaleMatrixByRow, &     
       Ubik_ScaleSystemByDiagonal => ScaleSystemByDiagonal, &
       Ubik_ScaleSystemByRow => ScaleSystemByRow, &     
       Ubik_SolveLTriang => SolveLTriang, &
       Ubik_SolveUTriang => SolveUTriang, &
       Ubik_SOR => SOR, &
       Ubik_Spew => Spew, &
       Ubik_Sum

  use Ubik_CSC_module, only: &
       Ubik_CSC_matrix => CSC_matrix, &
       assignment(=), &
       Ubik_create => create, &
       Ubik_destroy => destroy, &
       Ubik_alloc_failure => alloc_failure, &
       Ubik_set_status_alloc_failure => set_status_alloc_failure, &
       Ubik_Load => Load

  use Ubik_ELL_module, only: &
       Ubik_ELL_matrix => ELL_matrix, &
       assignment(=), &
       Ubik_create => create, &
       Ubik_destroy => destroy, &
       Ubik_alloc_failure => alloc_failure, &
       Ubik_set_status_alloc_failure => set_status_alloc_failure, &
       Ubik_MatVec => MatVec, &
       SIZE, &
       Ubik_Diag => Diag, &
       Ubik_DiagonalDominance => DiagonalDominance, &
       Ubik_Norm => Norm, &
       Ubik_NumberOfNonpositiveRows => NumberOfNonpositiveRows, &
       Ubik_NumberOfNonsymmElements => NumberOfNonsymmElements, &
       Ubik_Dump => Dump, &
       Ubik_EnsureFormat => EnsureFormat, &
       Ubik_IC => IC, &
       Ubik_ILU => ILU, &
       Ubik_Jacobi => Jacobi, &
       Ubik_Load => Load, &
       Ubik_LUdecomp => LUdecomp, &
       Ubik_ScaleMatrixByDiagonal => ScaleMatrixByDiagonal, &
       Ubik_ScaleMatrixByRow => ScaleMatrixByRow, &     
       Ubik_ScaleSystemByDiagonal => ScaleSystemByDiagonal, &
       Ubik_ScaleSystemByRow => ScaleSystemByRow, &     
       Ubik_SolveLTriang => SolveLTriang, &
       Ubik_SolveUTriang => SolveUTriang, &
       Ubik_SOR => SOR, &
       Ubik_Spew => Spew

  use Ubik_Matrix_module, only: &
       assignment(=)

end module UbikMatrix_module
