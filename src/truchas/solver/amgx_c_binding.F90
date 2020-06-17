!!
!! AMGX_C_BINDING
!!
!! Raw bindings to a subset of the AmgX library C interface.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! June 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module amgx_c_binding

  use,intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_size_t
  implicit none
  private


  !! Init && finalize
  public :: AMGX_initialize, AMGX_initialize_plugins
  public :: AMGX_finalize, AMGX_finalize_plugins
  interface
    function AMGX_initialize() &
        result(ierr) bind(c, name="AMGX_initialize")
      import c_int
      integer(c_int) :: ierr
    end function
    function AMGX_initialize_plugins() &
        result(ierr) bind(c, name="AMGX_initialize_plugins")
      import c_int
      integer(c_int) :: ierr
    end function
    function AMGX_finalize() &
        result(ierr) bind(c, name="AMGX_finalize")
      import c_int
      integer(c_int) :: ierr
    end function
    function AMGX_finalize_plugins() &
        result(ierr) bind(c, name="AMGX_finalize_plugins")
      import c_int
      integer(c_int) :: ierr
    end function
  end interface

  !! Resource create & destroy
  public :: AMGX_config_create
  public :: AMGX_resources_create_ext
  public :: AMGX_matrix_create_ext
  public :: AMGX_vector_create_ext
  public :: AMGX_solver_create_ext
  public :: AMGX_config_destroy
  public :: AMGX_resources_destroy
  public :: AMGX_vector_destroy
  public :: AMGX_matrix_destroy
  public :: AMGX_solver_destroy
  interface
    function AMGX_config_create(config, options) &
        result(ierr) bind(c, name="AMGX_config_create")
      import c_ptr, c_char, c_int
      type(c_ptr) :: config
      character(kind=c_char, len=*) :: options
      integer(c_int) :: ierr
    end function
    function AMGX_resources_create_ext(resources, config, ndevices, gpu_ids) &
        result(ierr) bind(c, name="AMGX_resources_create_ext")
      import c_ptr, c_int
      type(c_ptr) :: resources
      type(c_ptr), value :: config
      integer(c_int), value :: ndevices
      integer(c_int) :: gpu_ids(*)
      integer(c_int) :: ierr
    end function
    function AMGX_matrix_create_ext(matrix, resources) &
        result(ierr) bind(c, name="AMGX_matrix_create_ext")
      import c_ptr, c_int
      type(c_ptr) :: matrix
      type(c_ptr), value :: resources
      integer(c_int) :: ierr
    end function
    function AMGX_vector_create_ext(vector, resources) &
        result(ierr) bind(c, name="AMGX_vector_create_ext")
      import c_ptr, c_int
      type(c_ptr) :: vector
      type(c_ptr), value :: resources
      integer(c_int) :: ierr
    end function
    function AMGX_solver_create_ext(solver, resources, config) &
        result(ierr) bind(c, name="AMGX_solver_create_ext")
      import c_ptr, c_int
      type(c_ptr) :: solver
      type(c_ptr), value :: resources, config
      integer(c_int) :: ierr
    end function

    function AMGX_config_destroy(config) &
        result(ierr) bind(c, name="AMGX_config_destroy")
      import c_ptr, c_int
      type(c_ptr), value :: config
      integer(c_int) :: ierr
    end function
    function AMGX_resources_destroy(resources) &
        result(ierr) bind(c, name="AMGX_resources_destroy")
      import c_ptr, c_int
      type(c_ptr), value :: resources
      integer(c_int) :: ierr
    end function
    function AMGX_vector_destroy(vector) &
        result(ierr) bind(c, name="AMGX_vector_destroy")
      import c_ptr, c_int
      type(c_ptr), value :: vector
      integer(c_int) :: ierr
    end function
    function AMGX_matrix_destroy(matrix) &
        result(ierr) bind(c, name="AMGX_matrix_destroy")
      import c_ptr, c_int
      type(c_ptr), value :: matrix
      integer(c_int) :: ierr
    end function
    function AMGX_solver_destroy(solver) &
        result(ierr) bind(c, name="AMGX_solver_destroy")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int) :: ierr
    end function
  end interface


  public :: AMGX_matrix_upload_all_global
  public :: AMGX_solver_setup
  public :: AMGX_pin_memory
  public :: AMGX_vector_upload
  public :: AMGX_vector_download
  public :: AMGX_vector_set_zero
  public :: AMGX_solver_solve_with_0_initial_guess
  public :: AMGX_solver_get_status
  interface
    function AMGX_matrix_upload_all_global(matrix, nrow_global, nrow_onP, len_onP, &
        block_dimx, block_dimy, row_offsets, col_indices_global, data, diag_data, &
        allocated_halo_depth, num_import_rings, partition_vector) &
        result(ierr) bind(c, name="AMGX_matrix_upload_all_global")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: matrix, diag_data
      integer(c_int), value :: nrow_global, nrow_onP, len_onP, block_dimx, block_dimy, allocated_halo_depth, num_import_rings
      integer(c_int) :: row_offsets(*), col_indices_global(*), partition_vector(*)
      real(c_double) :: data(*) !, diag_data(*)
      integer(c_int) :: ierr
    end function
    function AMGX_solver_setup(solver, matrix) &
        result(ierr) bind(c, name="AMGX_solver_setup")
      import c_ptr, c_int
      type(c_ptr), value :: solver, matrix
      integer(c_int) :: ierr
    end function
    function AMGX_pin_memory(amgx_obj) &
        result(ierr) bind(c, name="AMGX_pin_memory")
      import c_ptr, c_int
      type(c_ptr), value :: amgx_obj
      integer(c_int) :: ierr
    end function
    function AMGX_vector_upload(dev, len, block_dim, host) &
        result(ierr) bind(c, name="AMGX_vector_upload")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: dev
      real(c_double) :: host(*)
      integer(c_int), value :: len, block_dim
      integer(c_int) :: ierr
    end function
    function AMGX_vector_download(dev, host) &
        result(ierr) bind(c, name="AMGX_vector_download")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: dev
      real(c_double) :: host(*)
      integer(c_int) :: ierr
    end function
    function AMGX_vector_set_zero(dev, len, block_dim) &
        result(ierr) bind(c, name="AMGX_vector_set_zero")
      import c_ptr, c_int
      type(c_ptr), value :: dev
      integer(c_int), value :: len, block_dim
      integer(c_int) :: ierr
    end function
    function AMGX_solver_solve_with_0_initial_guess(solver, rhs, soln) &
        result(ierr) bind(c, name="AMGX_solver_solve_with_0_initial_guess")
      import c_ptr, c_int
      type(c_ptr), value :: solver, rhs, soln
      integer(c_int) :: ierr
    end function
    function AMGX_solver_get_status(solver, status) &
        result(ierr) bind(c, name="AMGX_solver_get_status")
      import c_ptr, c_int
      type(c_ptr), value :: solver
      integer(c_int) :: status
      integer(c_int) :: ierr
    end function
  end interface

end module amgx_c_binding
