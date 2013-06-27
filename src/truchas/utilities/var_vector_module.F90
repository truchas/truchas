MODULE Var_Vector_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define routines to support operations on varying vectors, especially
  !   "Ragged Arrays" (arrays of varying vectors)
  !
  ! Public Subroutines
  !   * CREATE 
  !       Takes a scalar or array of unallocated varying vectors and a conformal
  !        SIZE arguments.  Allocates the vectors according to the SIZE argument
  !
  !   * DESTROY
  !       Takes a scalar or array of allocated varying vectors and deallocates 
  !       each vector.
  !
  ! Public Functions
  !   * SIZE
  !       Takes a scalar or array of allocated or unallocated varying vectors
  !       and returns a conformal scalar or array of sizes.  Unallocated 
  !       varying vectors have negative size.
  !
  !   * VECTOR
  !       1. Takes a ragged array (array of varying vectors) as an argument.
  !       Returns a pointer to a single large array which contains
  !       all the values of the ragged array.  Modifying these values may
  !       or may not effect the values in the ragged array.
  !       2. Takes a varying vector SCALAR as an argument.  Returns a 
  !       pointer to a vector of data which contains the values.
  !
  !   * FLATTEN
  !       Synonmous with VECTOR
  !
  !   * STUFF
  !       1. Takes a ragged array (array of varying vectors) and a
  !       flat array as an argument.
  !       STUFFs the values of the source array into the ragged array.
  !       2. Takes a scalar varying vector and a vector as arguments.
  !       STUFFs the vector into the varying vector data field.  The size
  !       of the source vector argument and the destination varying vector
  !       must conform.
  !
  ! Public Operators
  !   * = 
  !       Copies a scalar or array of varying vectors.  LHS must be
  !       conformal with RHS.
  !
  !   * +
  !       Adds conformal scalars or arrays of varying vectors
  !
  ! Author(s): Robert C. Ferrell, CPCA (ferrell@cpca.com)
  !
  !=======================================================================
  use kinds, only: r8
  use ArrayAllocate_Module
  use truchas_logging_services
  use var_vector_types
  implicit none
  private

  ! Public variables and types
  public :: CREATE,        &
            DESTROY,       &
            SIZES,         &
            VECTOR,        &
            FLATTEN,       &
            STUFF,         &
            ASSIGNMENT(=) 
!            OPERATOR(+)

  ! Pulbic types provided by var_vector_types
  public :: int_var_vector,  &
            real_var_vector, &
            log_var_vector

  INTERFACE CREATE
     module procedure V_V_CREATE_V_INT
     module procedure V_V_CREATE_V_REAL
     module procedure V_V_CREATE_V_LOG
  end INTERFACE

  INTERFACE DESTROY
     module procedure V_V_DESTROY_V_INT
     module procedure V_V_DESTROY_V_REAL
     module procedure V_V_DESTROY_V_LOG
  end INTERFACE

  INTERFACE SIZES
     module procedure V_V_SIZES_V_INT
     module procedure V_V_SIZES_V_REAL
     module procedure V_V_SIZES_V_LOG
     module procedure V_V_SIZES_S_INT
     module procedure V_V_SIZES_S_REAL
     module procedure V_V_SIZES_S_LOG
  end INTERFACE
  
  INTERFACE FLATTEN
     module procedure V_V_FLATTEN_V_INT
     module procedure V_V_FLATTEN_V_REAL
     module procedure V_V_FLATTEN_V_LOG
     module procedure V_V_FLATTEN_S_INT
     module procedure V_V_FLATTEN_S_REAL
     module procedure V_V_FLATTEN_S_LOG
  end INTERFACE

  INTERFACE VECTOR
     module procedure V_V_FLATTEN_V_INT
     module procedure V_V_FLATTEN_V_REAL
     module procedure V_V_FLATTEN_V_LOG
     module procedure V_V_FLATTEN_S_INT
     module procedure V_V_FLATTEN_S_REAL
     module procedure V_V_FLATTEN_S_LOG
  end INTERFACE
  
  INTERFACE STUFF
     module procedure V_V_STUFF_V_INT
     module procedure V_V_STUFF_V_REAL
     module procedure V_V_STUFF_V_LOG
     module procedure V_V_STUFF_S_INT
     module procedure V_V_STUFF_S_REAL
     module procedure V_V_STUFF_S_LOG
  end INTERFACE


  INTERFACE ASSIGNMENT ( = )
     module procedure V_V_STUFF_V_INT
     module procedure V_V_STUFF_V_REAL
     module procedure V_V_STUFF_V_LOG
     module procedure V_V_STUFF_S_INT
     module procedure V_V_STUFF_S_REAL
     Module procedure V_V_STUFF_S_LOG
  end INTERFACE
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  !***** V_V_CREATE Routines *****

  !====================================================================
  ! Purpose(s):
  !   Allocate ARRAY of varying vectors according to SIZES
  !   NOTE: This allocates each of the vectors.  The array 
  !   is already allocated.
  !====================================================================

  ! Integer

#define _ROUTINE_NAME_   V_V_CREATE_V_INT
#define _DATA_TYPE_      integer
#define _VAR_DATA_TYPE_  INT_VAR_VECTOR
#include "v_v_create.fpp"
  
  ! Real

#define _ROUTINE_NAME_   V_V_CREATE_V_REAL
#define _DATA_TYPE_      real(r8)
#define _VAR_DATA_TYPE_  REAL_VAR_VECTOR
#include "v_v_create.fpp"
  
  ! Logical

#define _ROUTINE_NAME_   V_V_CREATE_V_LOG
#define _DATA_TYPE_      logical
#define _VAR_DATA_TYPE_  LOG_VAR_VECTOR
#include "v_v_create.fpp"
  

  ! ***** V_V_DESTROY routines

  !====================================================================
  ! Purpose(s):
  !   Deallocate ARRAY of varying vectors.  
  !   NOTE: This deallocates each of the varying vectors.  It does not 
  !   deallocate the ARRAY.
  !====================================================================
    
  ! Integer
#define _ROUTINE_NAME_  V_V_DESTROY_V_INT
#define _DATA_TYPE_     integer
#define _VAR_DATA_TYPE_ INT_VAR_VECTOR
#include "v_v_destroy.fpp"

  ! Real
#define _ROUTINE_NAME_  V_V_DESTROY_V_REAL
#define _DATA_TYPE_     real(r8)
#define _VAR_DATA_TYPE_ REAL_VAR_VECTOR
#include "v_v_destroy.fpp"

  ! Logical
#define _ROUTINE_NAME_  V_V_DESTROY_V_LOG
#define _DATA_TYPE_     logical
#define _VAR_DATA_TYPE_ LOG_VAR_VECTOR
#include "v_v_destroy.fpp"

  ! ***** V_V_SIZES_V Routines
  !====================================================================
  ! Purpose(s):
  !   Return an array of the sizes of var_vectors
  !====================================================================
  
  ! Integer
#define _ROUTINE_NAME_  V_V_SIZES_V_INT
#define _DATA_TYPE_     integer
#define _VAR_DATA_TYPE_ INT_VAR_VECTOR
#include "v_v_sizes_v.fpp"

  ! REAL
#define _ROUTINE_NAME_  V_V_SIZES_V_REAL
#define _DATA_TYPE_     real(r8)
#define _VAR_DATA_TYPE_ REAL_VAR_VECTOR
#include "v_v_sizes_v.fpp"

  ! Logical
#define _ROUTINE_NAME_  V_V_SIZES_V_LOG
#define _DATA_TYPE_     logical
#define _VAR_DATA_TYPE_ LOG_VAR_VECTOR
#include "v_v_sizes_v.fpp"

  ! ***** V_V_SIZES_S Routines
  !====================================================================
  ! Purpose(s):
  !   Return the size of the var_vector
  !====================================================================
  
  ! Integer
#define _ROUTINE_NAME_  V_V_SIZES_S_INT
#define _DATA_TYPE_     integer
#define _VAR_DATA_TYPE_ INT_VAR_VECTOR
#include "v_v_sizes_s.fpp"

  ! REAL
#define _ROUTINE_NAME_  V_V_SIZES_S_REAL
#define _DATA_TYPE_     real(r8)
#define _VAR_DATA_TYPE_ REAL_VAR_VECTOR
#include "v_v_sizes_s.fpp"

  ! Logical
#define _ROUTINE_NAME_  V_V_SIZES_S_LOG
#define _DATA_TYPE_     logical
#define _VAR_DATA_TYPE_ LOG_VAR_VECTOR
#include "v_v_sizes_s.fpp"

  !***** FLATTEN routines
  !====================================================================
  ! Purpose(s):
  !   Return a pointer to a single large array which contains
  !   all the varying vectors in the ragged array.
  !====================================================================
  
  ! Integer
#define _ROUTINE_NAME_  V_V_FLATTEN_V_INT
#define _DATA_TYPE_     integer
#define _VAR_DATA_TYPE_ INT_VAR_VECTOR
#include "v_v_flatten_v.fpp"

  ! REAL
#define _ROUTINE_NAME_  V_V_FLATTEN_V_REAL
#define _DATA_TYPE_     real(r8)
#define _VAR_DATA_TYPE_ REAL_VAR_VECTOR
#include "v_v_flatten_v.fpp"

  ! Logical
#define _ROUTINE_NAME_  V_V_FLATTEN_V_LOG
#define _DATA_TYPE_     logical
#define _VAR_DATA_TYPE_ LOG_VAR_VECTOR
#include "v_v_flatten_v.fpp"

  !====================================================================
  ! Purpose(s):
  !   Return a pointer to vector of values in the varying vector
  !====================================================================

  ! Integer
#define _ROUTINE_NAME_  V_V_FLATTEN_S_INT
#define _DATA_TYPE_     integer
#define _VAR_DATA_TYPE_ INT_VAR_VECTOR
#include "v_v_flatten_s.fpp"

  ! Real
#define _ROUTINE_NAME_  V_V_FLATTEN_S_REAL
#define _DATA_TYPE_     real(r8)
#define _VAR_DATA_TYPE_ REAL_VAR_VECTOR
#include "v_v_flatten_s.fpp"

  ! Logical
#define _ROUTINE_NAME_  V_V_FLATTEN_S_LOG
#define _DATA_TYPE_     logical
#define _VAR_DATA_TYPE_ LOG_VAR_VECTOR
#include "v_v_flatten_s.fpp"

  !***** STUFF routines
  !====================================================================
  ! Purpose(s):
  !   STUFF the values of SOURCE into RAGGED_ARRAY.  In this imlementation
  !   this is efficient because RAGGED_ARRAY has a large container.
  !   If RAGGED_ARRAY is not allocated, then it is an error to call this routine.
  !   It is required that SIZE(SOURCE) == SIZES(RAGGED_ARRAY)
  !====================================================================

  ! Integer
#define _ROUTINE_NAME_  V_V_STUFF_V_INT
#define _DATA_TYPE_     integer
#define _VAR_DATA_TYPE_ INT_VAR_VECTOR
#include "v_v_stuff_v.fpp"  

  ! Real
#define _ROUTINE_NAME_  V_V_STUFF_V_REAL
#define _DATA_TYPE_     real(r8)
#define _VAR_DATA_TYPE_ REAL_VAR_VECTOR
#include "v_v_stuff_v.fpp"  

  ! Logical
#define _ROUTINE_NAME_  V_V_STUFF_V_LOG
#define _DATA_TYPE_     logical
#define _VAR_DATA_TYPE_ LOG_VAR_VECTOR
#include "v_v_stuff_v.fpp"  
  

  ! ***** STUFF SCALAR
  !====================================================================
  ! Purpose(s):
  !   STUFF the values of SOURCE into V_V_SCALAR.  This hides the
  !   implementatioh of varying vectors.
  !   If V_V_SCALAR is not allocated, then it is an error to call this routine.
  !   It is required that SIZE(SOURCE) == SIZES(V_V_SCALAR)
  !====================================================================
  
  ! Integer
#define _ROUTINE_NAME_  V_V_STUFF_S_INT
#define _DATA_TYPE_     integer
#define _VAR_DATA_TYPE_ INT_VAR_VECTOR
#include "v_v_stuff_s.fpp"

  ! Real
#define _ROUTINE_NAME_  V_V_STUFF_S_REAL
#define _DATA_TYPE_     real(r8)
#define _VAR_DATA_TYPE_ REAL_VAR_VECTOR
#include "v_v_stuff_s.fpp"

  ! Logical
#define _ROUTINE_NAME_  V_V_STUFF_S_LOG
#define _DATA_TYPE_     logical
#define _VAR_DATA_TYPE_ LOG_VAR_VECTOR
#include "v_v_stuff_s.fpp"

END MODULE Var_Vector_MODULE
