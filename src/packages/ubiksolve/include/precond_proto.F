! -*- Mode: f90 -*-

interface
   subroutine ApplyPrecond (b, x, status)
     use ubik_types_module
     use ubik_vector_type_module
     real(real_type), intent(in), dimension(:), target :: b
     type(vector_type), intent(inout) :: x
     integer(int_type), intent(out) :: status
   end subroutine ApplyPrecond
end interface
