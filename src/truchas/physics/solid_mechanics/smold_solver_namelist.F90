!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module smold_solver_namelist

  use parameter_list_type
  implicit none
  private

  public :: read_smold_solver_namelist

  type(parameter_list), public :: linear_params
  type(parameter_list), public :: nonlinear_params

contains

  subroutine read_smold_solver_namelist(lun)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist, NULL_R, NULL_I, NULL_C
    use string_utilities, only: i_to_c
    use truchas_logging_services
    use solid_mechanics_input, only: displacement_nonlinear_solution

    integer, intent(in) :: lun

    logical :: found
    integer :: ios
    character(64) :: iom

    !! deprecated namelist parameters
    character(64) :: method, preconditioning_scope, preconditioning_method

    character(64) :: name, linear_solver_name
    real(r8) :: convergence_criterion, nlk_vector_tolerance
    integer :: maximum_iterations, nlk_max_vectors
    namelist /nonlinear_solver/ name, linear_solver_name, &
        convergence_criterion, maximum_iterations, nlk_vector_tolerance, nlk_max_vectors, &
        method

    integer :: preconditioning_steps
    real(r8) :: relaxation_parameter
    namelist /linear_solver/ name, preconditioning_steps, relaxation_parameter, &
        method, preconditioning_scope, preconditioning_method
    
    call TLS_info('')
    call TLS_info('Reading NONLINEAR_SOLVER namelist ...')

    !! Locate the NONLINEAR_SOLVER namelist associated with displacement_nonlinear_solution
    found = .false.
    if (is_IOP) then
      rewind(lun)
      do while (.not.found)
        !! Default values
        name = NULL_C
        linear_solver_name = NULL_C
        convergence_criterion = NULL_R
        maximum_iterations = NULL_I
        nlk_vector_tolerance = NULL_R
        nlk_max_vectors = NULL_I

        call seek_to_namelist(lun, 'nonlinear_solver', found, iostat=ios)
        if (ios /= 0 .or. .not.found) exit

        read(lun, nml=nonlinear_solver, iostat=ios, iomsg=iom)
        if (ios /= 0) exit

        found = trim(name) == trim(displacement_nonlinear_solution)
      end do
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios) // ':  ' // trim(iom))
    call broadcast(found)
    if (.not.found) call TLS_fatal('NONLINEAR_SOLVER namelist not found')

    call broadcast(maximum_iterations)
    call broadcast(convergence_criterion)
    call broadcast(nlk_max_vectors)
    call broadcast(nlk_vector_tolerance)
    if (maximum_iterations /= NULL_I) call nonlinear_params%set('nlk-max-iter', maximum_iterations)
    if (convergence_criterion /= NULL_R) call nonlinear_params%set('nlk-tol', convergence_criterion)
    if (nlk_max_vectors /= NULL_I) call nonlinear_params%set('nlk-max-vec', nlk_max_vectors)
    if (nlk_vector_tolerance /= NULL_R) call nonlinear_params%set('nlk-vec-tol', nlk_vector_tolerance)


    call TLS_info('')
    call TLS_info('Reading LINEAR_SOLVER namelist ...')

    !! Locate the LINEAR_SOLVER namelist associated with linear_solver_name
    found = .false.
    if (is_IOP) then
      rewind(lun)
      do while (.not.found)
        !! Default values
        name = NULL_C
        preconditioning_steps = NULL_I
        relaxation_parameter = NULL_R

        call seek_to_namelist(lun, 'linear_solver', found, iostat=ios)
        if (ios /= 0 .or. .not.found) exit

        read(lun, nml=linear_solver, iostat=ios, iomsg=iom)
        if (ios /= 0) exit

        found = trim(name) == trim(linear_solver_name)
      end do
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios) // ':  ' // trim(iom))
    call broadcast(found)
    if (.not.found) call TLS_fatal('LINEAR_SOLVER namelist not found')

    call broadcast(preconditioning_steps)
    call broadcast(relaxation_parameter)
    if (preconditioning_steps /= NULL_I) call linear_params%set('precon-num-iter', preconditioning_steps)
    if (relaxation_parameter /= NULL_R) call linear_params%set('precon-relaxation-parameter', relaxation_parameter)

  end subroutine read_smold_solver_namelist

end module smold_solver_namelist
