!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  DUMMY IO ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! THESE WILL BE REPLACE WITH C ROUTINES !!!!!!!!!!!!!!!!!!!!!

! $Id: input-f.F,v 1.1.1.1 2000/10/11 22:44:23 ferrell Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Broadcast ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_bcast_int_scalar_c(scalar)
           IMPLICIT NONE
           INTEGER :: scalar
           RETURN
         END SUBROUTINE PGSLib_bcast_int_scalar_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_bcast_float_scalar_c(scalar)
           IMPLICIT NONE
           REAL ::scalar
           RETURN
         END SUBROUTINE PGSLib_bcast_float_scalar_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_bcast_double_scalar_c(scalar)
           IMPLICIT NONE
           REAL (KIND=KIND(1.0D0)) :: scalar
           RETURN
         END SUBROUTINE PGSLib_bcast_double_scalar_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_bcast_log_scalar_c(scalar)
           IMPLICIT NONE
           LOGICAL :: scalar
           RETURN
         END SUBROUTINE PGSLib_bcast_log_scalar_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Broadcast of Vectors
! These are dummy routines which are C routines in the actual library.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_bcast_int_vector_c(vector, vec_len)
           IMPLICIT NONE
           INTEGER, dimension(*):: vector
           integer vec_len

           return
         END SUBROUTINE PGSLib_bcast_int_vector_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_bcast_float_vector_c(vector, vec_len)
           IMPLICIT NONE
           REAL, dimension(*) ::vector
           integer vec_len

           return
         END SUBROUTINE PGSLib_bcast_float_vector_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_bcast_double_vector_c(vector, vec_len)
           IMPLICIT NONE
           REAL (KIND=KIND(1.0D0)), dimension(*):: vector
           integer vec_len

           return
         END SUBROUTINE PGSLib_bcast_double_vector_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_bcast_log_vector_c(vector, vec_len)
           IMPLICIT NONE
           LOGICAL, dimension(*):: vector
           integer vec_len

           return
         END SUBROUTINE PGSLib_bcast_log_vector_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_bcast_char_vector_c(vector, vec_len)
           IMPLICIT NONE
           CHARACTER (*) :: vector
           integer vec_len
           
           return
         END SUBROUTINE PGSLib_bcast_char_vector_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         DISTRIBUTE ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Dist_Scalar Routines

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_Dist_int_scalar_c(scalar_out, scalarv_in)
           IMPLICIT NONE
           INTEGER, INTENT(OUT) :: scalar_out
           INTEGER, INTENT(IN), DIMENSION(*):: scalarv_in

           scalar_out = scalarv_in( 1)
           RETURN
         END SUBROUTINE PGSLib_Dist_int_scalar_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_Dist_float_scalar_c(scalar_out, scalarv_in)
           IMPLICIT NONE
           REAL, INTENT(OUT) ::scalar_out
           REAL, INTENT(IN), DIMENSION(*) ::scalarv_in

           scalar_out = scalarv_in( 1)
           RETURN
         END SUBROUTINE PGSLib_Dist_float_scalar_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_Dist_double_scalar_c(scalar_out, scalarv_in)
           IMPLICIT NONE
           REAL (KIND=KIND(1.0D0)), INTENT(OUT) :: scalar_out
           REAL (KIND=KIND(1.0D0)), INTENT(IN), DIMENSION(*) :: scalarv_in

           scalar_out = scalarv_in( 1)
           RETURN
         END SUBROUTINE PGSLib_Dist_double_scalar_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_Dist_log_scalar_c(scalar_out, scalarv_in)
           IMPLICIT NONE
           LOGICAL, INTENT(OUT) :: scalar_out
           LOGICAL, INTENT(IN), DIMENSION(*) :: scalarv_in

           scalar_out = scalarv_in( 1)
           RETURN
         END SUBROUTINE PGSLib_Dist_log_scalar_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Dist Vector Routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_dist_int_vector_c(vector_out, out_len, vector_in, lengths)
           IMPLICIT NONE
           INTEGER, INTENT(IN):: out_len
           INTEGER, INTENT(OUT), DIMENSION(out_len):: vector_out
           INTEGER, INTENT(IN), DIMENSION(*):: vector_in
           INTEGER, INTENT(IN), DIMENSION(*):: lengths

           ! Local variables

           ! This only works on serial (or shared mem) version
           vector_out = vector_in(1:Lengths(1))

           RETURN
         END SUBROUTINE PGSLib_dist_int_vector_c

         SUBROUTINE PGSLib_dist_float_vector_c(vector_out, out_len, vector_in, lengths)
           IMPLICIT NONE
           INTEGER, INTENT(IN):: out_len
           REAL,    INTENT(OUT), DIMENSION(out_len):: vector_out
           REAL,    INTENT(IN), DIMENSION(*):: vector_in
           INTEGER, INTENT(IN), DIMENSION(*):: lengths

           ! Local variables
           vector_out = vector_in(1:Lengths(1))

           RETURN
         END SUBROUTINE PGSLib_dist_float_vector_c

         SUBROUTINE PGSLib_dist_double_vector_c(vector_out, out_len, vector_in, lengths)
           IMPLICIT NONE
           INTEGER, INTENT(IN):: out_len
           REAL (KIND=KIND(1.0D0)), INTENT(OUT), DIMENSION(out_len):: vector_out
           REAL (KIND=KIND(1.0D0)), INTENT(IN), DIMENSION(*):: vector_in
           INTEGER (KIND(1)), INTENT(IN), DIMENSION(*):: lengths

           ! Local variables

           vector_out = vector_in(1:Lengths(1))

           RETURN
         END SUBROUTINE PGSLib_dist_double_vector_c

         SUBROUTINE PGSLib_dist_log_vector_c(vector_out, out_len, vector_in, lengths)
           IMPLICIT NONE
           INTEGER, INTENT(IN):: out_len
           LOGICAL,  INTENT(OUT), DIMENSION(out_len):: vector_out
           LOGICAL, INTENT(IN), DIMENSION(*):: vector_in
           INTEGER, INTENT(IN), DIMENSION(*):: lengths

           ! Local variables

           vector_out = vector_in(1:Lengths(1))

           RETURN
         END SUBROUTINE PGSLib_dist_log_vector_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         COLLATE ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Collate_Scalar Routines

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_Collate_int_scalar_c(scalarv_out, scalar_in)
           IMPLICIT NONE
           INTEGER, INTENT(OUT), DIMENSION(*) :: scalarv_out
           INTEGER, INTENT(IN)                :: scalar_in

           ! This only works for serial version
           scalarv_out( 1) = scalar_in

           RETURN
         END SUBROUTINE PGSLib_Collate_int_scalar_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_Collate_float_scalar_c(scalarv_out, scalar_in)
           IMPLICIT NONE
           REAL, INTENT(OUT), DIMENSION(*) ::scalarv_out
           REAL, INTENT(IN)                ::scalar_in

           ! This only works for serial version
           scalarv_out( 1) = scalar_in

           RETURN
         END SUBROUTINE PGSLib_Collate_float_scalar_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_Collate_double_scalar_c(scalarv_out, scalar_in)
           IMPLICIT NONE
           REAL (KIND=KIND(1.0D0)), INTENT(OUT), DIMENSION(*)  :: scalarv_out
           REAL (KIND=KIND(1.0D0)), INTENT(IN)                 :: scalar_in

           ! This only works for serial version
           scalarv_out( 1) = scalar_in

           RETURN
         END SUBROUTINE PGSLib_Collate_double_scalar_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_Collate_log_scalar_c(scalarv_out, scalar_in)
           IMPLICIT NONE
           LOGICAL, INTENT(OUT) , DIMENSION(*) :: scalarv_out
           LOGICAL, INTENT(IN)                 :: scalar_in

           ! This only works for serial version
           scalarv_out( 1) = scalar_in

           RETURN
         END SUBROUTINE PGSLib_Collate_log_scalar_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Collate Vector Routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE PGSLib_Collate_int_vector_c(vector_out, Lengths, vector_in, In_Length)
           IMPLICIT NONE
           INTEGER, INTENT(OUT), DIMENSION(*):: vector_out
           INTEGER, INTENT(IN),  DIMENSION(*):: lengths
           INTEGER, INTENT(IN),  DIMENSION(*):: vector_in
           INTEGER, INTENT(IN):: In_Length

           ! Local variables

           vector_out(1:Lengths(1)) = vector_in(1:Lengths(1))
           
           RETURN
         END SUBROUTINE PGSLib_Collate_int_vector_c

         SUBROUTINE PGSLib_Collate_float_vector_c(vector_out, Lengths, vector_in, In_Length)
           IMPLICIT NONE
           REAL,   INTENT(OUT), DIMENSION(*):: vector_out
           INTEGER, INTENT(IN),  DIMENSION(*):: lengths
           REAL   ,   INTENT(IN),  DIMENSION(*):: vector_in
           INTEGER, INTENT(IN):: In_Length


           ! Local variables

           vector_out(1:Lengths(1)) = vector_in(1:Lengths(1))
           
           RETURN
         END SUBROUTINE PGSLib_Collate_float_vector_c

         SUBROUTINE PGSLib_Collate_double_vector_c(vector_out, Lengths, vector_in, In_Length)
           IMPLICIT NONE
           REAL (KIND=KIND(1.0D0)), INTENT(OUT), DIMENSION(*):: vector_out
           INTEGER, INTENT(IN),  DIMENSION(*):: lengths
           REAL (KIND=KIND(1.0D0)), INTENT(IN),  DIMENSION(*):: vector_in
           INTEGER, INTENT(IN):: In_Length


           ! Local variables
           INTEGER :: LowIndex, HighIndex

           ! This only works on serial (or shared mem) version
           ! Plucks the correct segment from vector_in.
           LowIndex  = SUM(Lengths(1:0 ) ) + 1
           HighIndex = SUM(Lengths(1: 1) )
!           vector_out(Lengths(LowIndex:HighIndex)) = vector_in
           vector_out(1:Lengths(1)) = vector_in(1:Lengths(1))
           
           RETURN
         END SUBROUTINE PGSLib_Collate_double_vector_c

         SUBROUTINE PGSLib_Collate_log_vector_c(vector_out, Lengths, vector_in, In_Length)
           IMPLICIT NONE
           LOGICAL, INTENT(OUT), DIMENSION(*):: vector_out
           INTEGER, INTENT(IN),  DIMENSION(*):: lengths
           LOGICAL, INTENT(IN),  DIMENSION(*):: vector_in
           INTEGER, INTENT(IN):: In_Length


           ! Local variables

           vector_out(1:Lengths(1)) = vector_in(1:Lengths(1))
           
           RETURN
         END SUBROUTINE PGSLib_Collate_log_vector_c

         SUBROUTINE PGSLib_Collate_char_vector_c(vector_out, Lengths, vector_in, In_Length)
           IMPLICIT NONE
           CHARACTER (*), INTENT(  OUT), DIMENSION(*):: vector_out
           INTEGER       , INTENT(IN   ), DIMENSION(*):: lengths
           CHARACTER (*), INTENT(IN   ), DIMENSION(*):: vector_in
           INTEGER       , INTENT(IN   )              :: In_Length

           
           vector_out(1:Lengths(1)) = vector_in(1:Lengths(1))
           
           RETURN
         END SUBROUTINE PGSLib_Collate_char_vector_c

