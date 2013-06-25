!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Global Reduction Routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! This file contains the DUMMY, serial versions.  This file
!     is replaced by a .c file containing MPI calls in the
!     actual library.
!
!


! $Id: redux-f.F,v 1.2 2001/03/22 00:26:13 ferrell Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! MINVAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine PGSLib_global_min_int_c(MinA)
        integer, INTENT(INOUT):: MinA

        RETURN
        END subroutine PGSLib_global_min_int_c

        subroutine PGSLib_global_min_float_c(MinA)
        REAL, INTENT(INOUT):: MinA

        RETURN
        END subroutine PGSLib_global_min_float_c

        subroutine PGSLib_global_min_double_c(MinA)
        real (KIND=KIND(1.0D0)), INTENT(INOUT):: MinA

        RETURN
        end subroutine PGSLib_global_min_double_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! MAXVAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine PGSLib_global_max_int_c(MaxA)
        integer, INTENT(INOUT):: MaxA

        RETURN
        END subroutine PGSLib_global_max_int_c

        subroutine PGSLib_global_max_float_c(MaxA)
        REAL, INTENT(INOUT):: MaxA

        RETURN
        END subroutine PGSLib_global_max_float_c

        subroutine PGSLib_global_max_double_c(MaxA)
        real (KIND=KIND(1.0D0)), INTENT(INOUT):: MaxA

        RETURN
        end subroutine PGSLib_global_max_double_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! SUM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine PGSLib_global_sum_int_C(SumA)
        integer, INTENT(INOUT):: SumA

        RETURN
        END subroutine PGSLib_global_sum_int_C

        subroutine PGSLib_global_sum_float_C(SumA)
        REAL, INTENT(INOUT):: SumA

        RETURN
        END subroutine PGSLib_global_sum_float_C

        subroutine PGSLib_global_sum_double_C(SumA)
        real (KIND=KIND(1.0D0)), INTENT(INOUT):: SumA

        RETURN
        end subroutine PGSLib_global_sum_double_C

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! ALL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine PGSLib_global_all_log_C(AllA)
        INTEGER, INTENT(INOUT):: AllA

        RETURN
        END subroutine PGSLib_global_all_log_C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! ANY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine PGSLib_global_any_log_C(AnyA)
        INTEGER, INTENT(INOUT):: AnyA

        RETURN
        END subroutine PGSLib_global_any_log_C

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! MINLOC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine pgslib_global_minloc_int_c(minv, glbl_index)
          integer, intent (IN OUT) :: glbl_index
          integer, intent (IN OUT) :: minv

          RETURN
        END subroutine pgslib_global_minloc_int_c

        subroutine pgslib_global_minloc_float_c(minv, glbl_index)
          integer, intent (IN OUT) :: glbl_index
          real,    intent (IN OUT) :: minv

          RETURN
        END subroutine pgslib_global_minloc_float_c

        subroutine pgslib_global_minloc_double_c(minv, glbl_index)
          integer,            intent (IN OUT) :: glbl_index
          real (kind(1.0d0)), intent (IN OUT) :: minv

          RETURN
        END subroutine pgslib_global_minloc_double_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! MAXLOC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine pgslib_global_maxloc_int_c(minv, glbl_index)
          integer, intent (IN OUT) :: glbl_index
          integer, intent (IN OUT) :: minv

          RETURN
        END subroutine pgslib_global_maxloc_int_c

        subroutine pgslib_global_maxloc_float_c(minv, glbl_index)
          integer, intent (IN OUT) :: glbl_index
          real,    intent (IN OUT) :: minv

          RETURN
        END subroutine pgslib_global_maxloc_float_c

        subroutine pgslib_global_maxloc_double_c(minv, glbl_index)
          integer,            intent (IN OUT) :: glbl_index
          real (kind(1.0d0)), intent (IN OUT) :: minv

          RETURN
        END subroutine pgslib_global_maxloc_double_c

