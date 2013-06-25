MODULE VAR_VECTOR_TYPES
  !=======================================================================
  ! Purpose(s):
  !
  !   Define structures for supporting "Ragged Arrays".  Here we define
  !   the "varying vectors" which will be the elements of a Ragged Array.
  !
  !   Public Interface(s):
  !
  !   Public Type(s):
  ! 
  !     * INT_VAR_VECTOR
  !
  !     * REAL_VAR_VECTOR
  !
  !     * LOG_VAR_VECTOR
  !
  !
  ! Contains: INT_VAR_VECTOR
  !           REAL_VAR_VECTOR
  !           LOG_VAR_VECTOR
  !
  ! Author(s): Robert C. Ferrell, CPCA (ferrell@cpca.com)
  !=======================================================================

  use kind_module,      only: int_kind, log_kind, real_kind

  implicit none

  private

  ! Public variables and types
  public :: INT_VAR_VECTOR,    &
            REAL_VAR_VECTOR,   &
            LOG_VAR_VECTOR

  ! File Version
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Define varying vectors
  type INT_VAR_VECTOR
     integer (int_kind)                        :: L
     integer (int_kind), POINTER, DIMENSION(:) :: V => NULL()
     integer (int_kind), POINTER, DIMENSION(:) :: CONTAINER => NULL()
  end type INT_VAR_VECTOR
  
  type REAL_VAR_VECTOR
     integer (int_kind)                        :: L
     REAL (REAL_kind),   POINTER, DIMENSION(:) :: V => NULL()
     REAL (REAL_kind),   POINTER, DIMENSION(:) :: CONTAINER => NULL()
  end type REAL_VAR_VECTOR
  
  type LOG_VAR_VECTOR
     integer (int_kind)                        :: L
     logical (log_kind), POINTER, DIMENSION(:) :: V => NULL()
     logical (log_kind), POINTER, DIMENSION(:) :: CONTAINER => NULL()
  end type LOG_VAR_VECTOR
  
end MODULE VAR_VECTOR_TYPES
