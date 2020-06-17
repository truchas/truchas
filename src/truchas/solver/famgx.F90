!!
!! FAMGX
!!
!!  A Fortran interface to a subset of the AmgX linear solver library.
!!
!!  Zach Jibben <zjibben@lanl.gov
!!  June 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module famgx

  use amgx_c_binding
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
  use,intrinsic :: iso_c_binding, only: amgx_obj => c_ptr
  use,intrinsic :: iso_c_binding, only: amgx_null_obj => c_null_ptr
#ifndef INTEL_SRN04330341
  use,intrinsic :: iso_c_binding, only: amgx_associated => c_associated
#endif
  implicit none
  private

  public :: amgx_obj, amgx_null_obj, amgx_associated

  public :: famgx_initialize, famgx_initialize_plugins
  public :: famgx_finalize, famgx_finalize_plugins

  public :: famgx_config_create
  public :: famgx_resources_create
  public :: famgx_matrix_create
  public :: famgx_vector_create
  public :: famgx_solver_create
  public :: famgx_config_destroy
  public :: famgx_resources_destroy
  public :: famgx_vector_destroy
  public :: famgx_matrix_destroy
  public :: famgx_solver_destroy

  public :: famgx_matrix_upload_all_global
  public :: famgx_solver_setup
  public :: famgx_pin_memory
  public :: famgx_vector_upload
  public :: famgx_vector_download
  public :: famgx_vector_set_zero
  public :: famgx_solver_solve_with_0_initial_guess
  public :: famgx_solver_get_status

#ifdef INTEL_SRN04330341
  interface amgx_associated
    procedure amgx_associated1, amgx_associated2
  end interface
#endif

contains

#ifdef INTEL_SRN04330341
  ! this only covers the present use case of c_associated in the amgx context
  logical function amgx_associated1(c_ptr1)
    use,intrinsic :: iso_c_binding, only: c_associated
    type(c_ptr), intent(in) :: c_ptr1
    amgx_associated1 = c_associated(c_ptr1)
  end function
  logical function amgx_associated2(c_ptr1, c_ptr2)
    use,intrinsic :: iso_c_binding, only: c_associated
    type(c_ptr), intent(in) :: c_ptr1, c_ptr2
    amgx_associated2 = c_associated(c_ptr1, c_ptr2)
  end function
#endif

  subroutine famgx_initialize(ierr)
    integer, intent(out) :: ierr
    ierr = AMGX_initialize()
  end subroutine famgx_initialize

  subroutine famgx_initialize_plugins(ierr)
    integer, intent(out) :: ierr
    ierr = AMGX_initialize_plugins()
  end subroutine famgx_initialize_plugins

  subroutine famgx_finalize(ierr)
    integer, intent(out) :: ierr
    ierr = AMGX_finalize()
  end subroutine famgx_finalize

  subroutine famgx_finalize_plugins(ierr)
    integer, intent(out) :: ierr
    ierr = AMGX_finalize_plugins()
  end subroutine famgx_finalize_plugins

  subroutine famgx_config_create(config, options, ierr)
    use,intrinsic :: iso_c_binding, only: c_char
    type(amgx_obj), intent(out) :: config
    character(*), intent(in) :: options
    integer, intent(out) :: ierr
    ierr = AMGX_config_create(config, f_c_string(options))
  end subroutine famgx_config_create

  subroutine famgx_resources_create(resources, config, ndevices, gpu_ids, ierr)
    type(amgx_obj), intent(out) :: resources
    type(amgx_obj), intent(in) :: config
    integer, intent(in) :: ndevices, gpu_ids(:)
    integer, intent(out) :: ierr
    ierr = AMGX_resources_create_ext(resources, config, ndevices, gpu_ids)
  end subroutine famgx_resources_create

  subroutine famgx_matrix_create(matrix, resources, ierr)
    type(amgx_obj), intent(out) :: matrix
    type(amgx_obj), intent(in) :: resources
    integer, intent(out) :: ierr
    ierr = AMGX_matrix_create_ext(matrix, resources)
  end subroutine famgx_matrix_create

  subroutine famgx_vector_create(vector, resources, ierr)
    type(amgx_obj), intent(out) :: vector
    type(amgx_obj), intent(in) :: resources
    integer, intent(out) :: ierr
    ierr = AMGX_vector_create_ext(vector, resources)
  end subroutine famgx_vector_create

  subroutine famgx_solver_create(solver, resources, config, ierr)
    type(amgx_obj), intent(out) :: solver
    type(amgx_obj), intent(in) :: resources, config
    integer, intent(out) :: ierr
    ierr = AMGX_solver_create_ext(solver, config, resources)
  end subroutine famgx_solver_create

  subroutine famgx_config_destroy(config, ierr)
    type(amgx_obj), intent(inout) :: config
    integer, intent(out) :: ierr
    ierr = AMGX_config_destroy(config)
    config = amgx_null_obj
  end subroutine famgx_config_destroy

  subroutine famgx_resources_destroy(resources, ierr)
    type(amgx_obj), intent(inout) :: resources
    integer, intent(out) :: ierr
    ierr = AMGX_resources_destroy(resources)
    resources = amgx_null_obj
  end subroutine famgx_resources_destroy

  subroutine famgx_vector_destroy(vector, ierr)
    type(amgx_obj), intent(inout) :: vector
    integer, intent(out) :: ierr
    ierr = AMGX_vector_destroy(vector)
    vector = amgx_null_obj
  end subroutine famgx_vector_destroy

  subroutine famgx_matrix_destroy(matrix, ierr)
    type(amgx_obj), intent(inout) :: matrix
    integer, intent(out) :: ierr
    ierr = AMGX_matrix_destroy(matrix)
    matrix = amgx_null_obj
  end subroutine famgx_matrix_destroy

  subroutine famgx_solver_destroy(solver, ierr)
    type(amgx_obj), intent(inout) :: solver
    integer, intent(out) :: ierr
    ierr = AMGX_solver_destroy(solver)
    solver = amgx_null_obj
  end subroutine famgx_solver_destroy

  subroutine famgx_solver_setup(solver, matrix, ierr)
    type(amgx_obj), intent(out) :: solver
    type(amgx_obj), intent(in) :: matrix
    integer, intent(out) :: ierr
    ierr = AMGX_solver_setup(solver, matrix)
  end subroutine famgx_solver_setup

  subroutine famgx_pin_memory(x, ierr)
    use,intrinsic :: iso_c_binding, only: c_loc
    real(r8), intent(inout) :: x(:)
    integer, intent(out) :: ierr
    ierr = AMGX_pin_memory(c_loc(x))
  end subroutine famgx_pin_memory

  subroutine famgx_vector_upload(dev, len, block_dim, host, ierr)
    type(amgx_obj), intent(in) :: dev
    integer, intent(in) :: len, block_dim
    real(r8), intent(in) :: host(:)
    integer, intent(out) :: ierr
    ierr = AMGX_vector_upload(dev, len, block_dim, host)
  end subroutine famgx_vector_upload

  subroutine famgx_vector_download(dev, host, ierr)
    type(amgx_obj), intent(in) :: dev
    real(r8), intent(out) :: host(:)
    integer, intent(out) :: ierr
    ierr = AMGX_vector_download(dev, host)
  end subroutine famgx_vector_download

  subroutine famgx_vector_set_zero(dev, len, block_dim, ierr)
    type(amgx_obj), intent(in) :: dev
    integer, intent(in) :: len, block_dim
    integer, intent(out) :: ierr
    ierr = AMGX_vector_set_zero(dev, len, block_dim)
  end subroutine famgx_vector_set_zero

  subroutine famgx_solver_solve_with_0_initial_guess(solver, rhs, soln, ierr)
    type(amgx_obj), intent(in) :: solver, rhs, soln
    integer, intent(out) :: ierr
    ierr = AMGX_solver_solve_with_0_initial_guess(solver, rhs, soln)
  end subroutine famgx_solver_solve_with_0_initial_guess

  subroutine famgx_solver_get_status(solver, status, ierr)
    type(amgx_obj), intent(in) :: solver
    integer, intent(out) :: status, ierr
    ierr = AMGX_solver_get_status(solver, status)
  end subroutine famgx_solver_get_status

  !! Note diag_data is hardcoded to NULL
  subroutine famgx_matrix_upload_all_global(matrix, nrow_global, nrow_onP, len_onP, &
      block_dimx, block_dimy, row_offsets, col_indices_global, data, &
      allocated_halo_depth, num_import_rings, partition_vector, ierr)
    type(amgx_obj), intent(in) :: matrix
    integer, intent(in) :: nrow_global, nrow_onP, len_onP, block_dimx, block_dimy, allocated_halo_depth, num_import_rings
    integer, intent(in) :: row_offsets(:), col_indices_global(:), partition_vector(:)
    real(r8), intent(in) :: data(:)
    integer, intent(out) :: ierr
    ierr = AMGX_matrix_upload_all_global(matrix, nrow_global, nrow_onP, len_onP, &
        block_dimx, block_dimy, row_offsets, col_indices_global, data, amgx_null_obj, &
        allocated_halo_depth, num_import_rings, partition_vector)
  end subroutine famgx_matrix_upload_all_global

  !! Create a null-terminated C string from an input Fortran character scalar.
  !! This function imitates the Fortran 202x F_C_STRING intrinsic, and can
  !! eventually be removed. Cf. J3 papers 18-258r2 and 19-197r3.
  function f_c_string(fstr, trim_str) result(cstr)

    use,intrinsic :: iso_c_binding, only: c_char, c_null_char

    character(*), intent(in) :: fstr
    logical, intent(in), optional :: trim_str
    character(kind=c_char, len=:), allocatable :: cstr

    logical :: trimh

    trimh = .true.
    if (present(trim_str)) trimh = trim_str

    if (trimh) then
      cstr = trim(fstr) // c_null_char
    else
      cstr = fstr // c_null_char
    end if

  end function f_c_string

end module famgx
