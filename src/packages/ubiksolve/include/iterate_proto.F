! -*- Mode: f90 -*-

interface
   subroutine Iterate (b, x, status)
     use ubik_types_module
     use ubik_vector_type_module
     real(real_type), intent(in), dimension(:) :: b
     real(real_type), intent(inout), dimension(:) :: x
     integer(int_type), intent(out) :: status
   end subroutine Iterate
end interface
