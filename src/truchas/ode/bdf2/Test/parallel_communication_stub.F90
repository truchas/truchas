module parallel_communication
  use kinds
  public :: is_IOP, global_dot_product
  logical :: is_IOP = .true.
contains
  function global_dot_product (a, b) result (dp)
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: dp
    dp = dot_product(a,b)
  end function
end module
