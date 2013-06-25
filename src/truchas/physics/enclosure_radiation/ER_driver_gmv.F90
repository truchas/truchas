#include "f90_assert.fpp"

module ER_driver_gmv

  use ER_solver_gmv, only: ERD_gmv_open => ERS_gmv_open, ERD_gmv_close => ERS_gmv_close, &
                           ERD_gmv_begin_variables => ERS_gmv_begin_variables, &
                           ERD_gmv_end_variables => ERS_gmv_end_variables, &
                           ERS_gmv_write_enclosure, ERS_gmv_write_var
  implicit none
  public

contains

  subroutine ERD_gmv_write_enclosure (this)

    use ER_driver, only: ERD_problem

    type(ERD_problem), intent(in) :: this

    call ERS_gmv_write_enclosure (this%sol)

  end subroutine ERD_gmv_write_enclosure


  subroutine ERD_gmv_write_var (this, var, name)

    use kinds, only: r8
    use ER_driver, only: ERD_problem
    use parallel_permutations

    type(ERD_problem), intent(in) :: this
    real(r8), intent(in) :: var(:)
    character(len=*), intent(in) :: name

    real(r8) :: var_er(this%nface_er)

    ASSERT(size(var) == this%nface_hc)

    call reorder (this%perm_er_to_hc, var_er, var)
    call ERS_gmv_write_var (this%sol, var_er, name)

  end subroutine ERD_gmv_write_var

end module ER_driver_gmv
