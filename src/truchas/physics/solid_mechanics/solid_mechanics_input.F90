!!
!!  SOLID_MECHANICS_INPUT
!!
!! This module owns the logical flags solid_mechanics and
!! stress_reduced_integration. Both flags are used outside
!! of solid_mechanics, but should be elevated to a multiprocessor
!! component at some point in the future.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  SOLID_MECHANICS
!!    Logical flag that indicates if the solid mechanics package
!!    is active. Was originally defined in the SOLID_MECHANICS_DATA
!!    module.
!!
!!  SOLID_MECHANICS_BODY_FORCE
!!    Logical flag that controls the addition of body forces in 
!!    SOLID_MECHANICS_MODULE. Set in PHYSICS_INPUT and is set to
!!    false by default. Was originally defined in SOLID_MECHANICS_DATA
!!    module.
!!
!!  STRESS_REDUCED_INTEGRATION
!!    Logical flag that controls the stress integration. If .true. only
!!    one stress per cell, otherwise calculate the elastic stress at the
!!    centroid of each control volume. It is NOT supported at this time
!!    Should always be defined .false.. The NUMERICS_INPUT_MODULE sets this
!!    parameter. Originally defined in NODE_OPERATOR_MODULE. 
!!
!!  These input parameters were originally defined in SOLID_MECHANICS_DATA
!!  Comments here are based on reading code and comment in the 
!!  NUMERICS_INPUT_MODULE and the usage in the SM package.
!!
!!  DISPLACEMENT_LINEAR_SOLUTION
!!    String variable defined in NUMERICS_INPUT_MODULE that controls which
!!    Ubik linear solver is selected for displacement. Appears to be hardcoded 
!!    to 'default' and error is thrown for any other choice.
!!
!!  DISPLACEMENT_NONLINEAR_SOLUTION
!!    String variable defined in NUMERICS_INPUT_MODULE that controls which
!!    Ubik nonlinear solver is selected for displacement. Appears to be !!
!!    hardcoded to 'default' and error is thrown for any other choice.
!!      
!!  CONTACT_DISTANCE
!!  CONTACT_NORM_TRAC
!!    Parameters used to compute lambda in GET_LAMBDA. Not in contact if
!!    n (surface normal)  \dot u_diff is larger than CONTACT_DISTANCE and 
!!    normal traction (?) is greater than CONTACT_NORM_TRAC. Value of 
!!    lambda depends on both if either condition fails.
!!  
!!  CONTACT_PENALTY
!!    Parameter used for pre-conditioning in MECH_PRECOND_DISP_CONSTRAINTS
!!    
!!  STRAIN_LIMIT
!!    From the comments in NUMERICS_INPUT_MODULE: Minimum plastic strain limit.
!!   
!!  UBIK_DISPLACEMENT
!!    Ubik linear solver control parameter. Used in NUMERICS_INPUT_MODULE.
!!
!!  NK_DIAPLCEMENT
!!    NK (nonlinear) solver control parameter. Used in NUMERICS_INPUT_MODULE.
!!    
module solid_mechanics_input

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  public

  logical :: solid_mechanics_body_force
  logical :: stress_reduced_integration

  !! Numeric input parameters
      
  !! Contact parameters
  real(r8) :: contact_distance
  real(r8) :: contact_norm_trac
  real(r8) :: contact_penalty
  
  !! Maximum plsatic strain increment
  real(r8) :: strain_limit

end module solid_mechanics_input

